import copy

import numpy as np
from torch_geometric.loader import DataLoader
from tqdm import tqdm

from confidence.dataset import ListDataset
from utils import so3, torus
from utils.sampling import randomize_position, sampling
import torch
from utils.diffusion_utils import get_t_schedule


def loss_function(tr_pred, rot_pred, tor_pred, data, t_to_sigma, device, tr_weight=1, rot_weight=1,
                  tor_weight=1, apply_mean=True, no_torsion=False):
    tr_sigma, rot_sigma, tor_sigma = t_to_sigma(
        *[torch.cat([d.complex_t[noise_type] for d in data]) if device.type == 'cuda' else data.complex_t[noise_type]
          for noise_type in ['tr', 'rot', 'tor']])
    mean_dims = (0, 1) if apply_mean else 1

    # translation component
    tr_score = torch.cat([d.tr_score for d in data], dim=0) if device.type == 'cuda' else data.tr_score
    tr_sigma = tr_sigma.unsqueeze(-1)
    tr_loss = ((tr_pred.cpu() - tr_score) ** 2 * tr_sigma ** 2).mean(dim=mean_dims)
    tr_base_loss = (tr_score ** 2 * tr_sigma ** 2).mean(dim=mean_dims).detach()

    # rotation component
    rot_score = torch.cat([d.rot_score for d in data], dim=0) if device.type == 'cuda' else data.rot_score
    rot_score_norm = so3.score_norm(rot_sigma.cpu()).unsqueeze(-1)
    rot_loss = (((rot_pred.cpu() - rot_score) / rot_score_norm) ** 2).mean(dim=mean_dims)
    rot_base_loss = ((rot_score / rot_score_norm) ** 2).mean(dim=mean_dims).detach()

    # torsion component
    if not no_torsion:
        edge_tor_sigma = torch.from_numpy(
            np.concatenate([d.tor_sigma_edge for d in data] if device.type == 'cuda' else data.tor_sigma_edge))
        tor_score = torch.cat([d.tor_score for d in data], dim=0) if device.type == 'cuda' else data.tor_score
        tor_score_norm2 = torch.tensor(torus.score_norm(edge_tor_sigma.cpu().numpy())).float()
        tor_loss = ((tor_pred.cpu() - tor_score) ** 2 / tor_score_norm2)
        tor_base_loss = ((tor_score ** 2 / tor_score_norm2)).detach()
        if apply_mean:
            tor_loss, tor_base_loss = tor_loss.mean() * torch.ones(1, dtype=torch.float), tor_base_loss.mean() * torch.ones(1, dtype=torch.float)
        else:
            index = torch.cat([torch.ones(d['ligand'].edge_mask.sum()) * i for i, d in
                               enumerate(data)]).long() if device.type == 'cuda' else data['ligand'].batch[
                data['ligand', 'ligand'].edge_index[0][data['ligand'].edge_mask]]
            num_graphs = len(data) if device.type == 'cuda' else data.num_graphs
            t_l, t_b_l, c = torch.zeros(num_graphs), torch.zeros(num_graphs), torch.zeros(num_graphs)
            c.index_add_(0, index, torch.ones(tor_loss.shape))
            c = c + 0.0001
            t_l.index_add_(0, index, tor_loss)
            t_b_l.index_add_(0, index, tor_base_loss)
            tor_loss, tor_base_loss = t_l / c, t_b_l / c
    else:
        if apply_mean:
            tor_loss, tor_base_loss = torch.zeros(1, dtype=torch.float), torch.zeros(1, dtype=torch.float)
        else:
            tor_loss, tor_base_loss = torch.zeros(len(rot_loss), dtype=torch.float), torch.zeros(len(rot_loss), dtype=torch.float)

    loss = tr_loss * tr_weight + rot_loss * rot_weight + tor_loss * tor_weight
    return loss, tr_loss.detach(), rot_loss.detach(), tor_loss.detach(), tr_base_loss, rot_base_loss, tor_base_loss


class AverageMeter():
    def __init__(self, types, unpooled_metrics=False, intervals=1):
        self.types = types
        self.intervals = intervals
        self.count = 0 if intervals == 1 else torch.zeros(len(types), intervals)
        self.acc = {t: torch.zeros(intervals) for t in types}
        self.unpooled_metrics = unpooled_metrics

    def add(self, vals, interval_idx=None):
        if self.intervals == 1:
            self.count += 1 if vals[0].dim() == 0 else len(vals[0])
            for type_idx, v in enumerate(vals):
                self.acc[self.types[type_idx]] += v.sum() if self.unpooled_metrics else v
        else:
            for type_idx, v in enumerate(vals):
                self.count[type_idx].index_add_(0, interval_idx[type_idx], torch.ones(len(v)))
                if not torch.allclose(v, torch.tensor(0.0)):
                    self.acc[self.types[type_idx]].index_add_(0, interval_idx[type_idx], v)

    def summary(self):
        if self.intervals == 1:
            out = {k: v.item() / self.count for k, v in self.acc.items()}
            return out
        else:
            out = {}
            for i in range(self.intervals):
                for type_idx, k in enumerate(self.types):
                    out['int' + str(i) + '_' + k] = (
                            list(self.acc.values())[type_idx][i] / self.count[type_idx][i]).item()
            return out


def train_epoch(model, loader, optimizer, device, t_to_sigma, loss_fn, ema_weigths):
    model.train()
    meter = AverageMeter(['loss', 'tr_loss', 'rot_loss', 'tor_loss', 'tr_base_loss', 'rot_base_loss', 'tor_base_loss'])

    for data in tqdm(loader, total=len(loader)):
        if device.type == 'cuda' and len(data) == 1 or device.type == 'cpu' and data.num_graphs == 1:
            print("Skipping batch of size 1 since otherwise batchnorm would not work.")
        optimizer.zero_grad()
        try:
            tr_pred, rot_pred, tor_pred = model(data)
            loss, tr_loss, rot_loss, tor_loss, tr_base_loss, rot_base_loss, tor_base_loss = \
                loss_fn(tr_pred, rot_pred, tor_pred, data=data, t_to_sigma=t_to_sigma, device=device)
            loss.backward()
            optimizer.step()
            ema_weigths.update(model.parameters())
            meter.add([loss.cpu().detach(), tr_loss, rot_loss, tor_loss, tr_base_loss, rot_base_loss, tor_base_loss])
        except RuntimeError as e:
            if 'out of memory' in str(e):
                print('| WARNING: ran out of memory, skipping batch')
                for p in model.parameters():
                    if p.grad is not None:
                        del p.grad  # free some memory
                torch.cuda.empty_cache()
                continue
            elif 'Input mismatch' in str(e):
                print('| WARNING: weird torch_cluster error, skipping batch')
                for p in model.parameters():
                    if p.grad is not None:
                        del p.grad  # free some memory
                torch.cuda.empty_cache()
                continue
            else:
                raise e

    return meter.summary()


def test_epoch(model, loader, device, t_to_sigma, loss_fn, test_sigma_intervals=False):
    model.eval()
    meter = AverageMeter(['loss', 'tr_loss', 'rot_loss', 'tor_loss', 'tr_base_loss', 'rot_base_loss', 'tor_base_loss'],
                         unpooled_metrics=True)

    if test_sigma_intervals:
        meter_all = AverageMeter(
            ['loss', 'tr_loss', 'rot_loss', 'tor_loss', 'tr_base_loss', 'rot_base_loss', 'tor_base_loss'],
            unpooled_metrics=True, intervals=10)

    for data in tqdm(loader, total=len(loader)):
        try:
            with torch.no_grad():
                tr_pred, rot_pred, tor_pred = model(data)

            loss, tr_loss, rot_loss, tor_loss, tr_base_loss, rot_base_loss, tor_base_loss = \
                loss_fn(tr_pred, rot_pred, tor_pred, data=data, t_to_sigma=t_to_sigma, apply_mean=False, device=device)
            meter.add([loss.cpu().detach(), tr_loss, rot_loss, tor_loss, tr_base_loss, rot_base_loss, tor_base_loss])

            if test_sigma_intervals > 0:
                complex_t_tr, complex_t_rot, complex_t_tor = [torch.cat([d.complex_t[noise_type] for d in data]) for
                                                              noise_type in ['tr', 'rot', 'tor']]
                sigma_index_tr = torch.round(complex_t_tr.cpu() * (10 - 1)).long()
                sigma_index_rot = torch.round(complex_t_rot.cpu() * (10 - 1)).long()
                sigma_index_tor = torch.round(complex_t_tor.cpu() * (10 - 1)).long()
                meter_all.add(
                    [loss.cpu().detach(), tr_loss, rot_loss, tor_loss, tr_base_loss, rot_base_loss, tor_base_loss],
                    [sigma_index_tr, sigma_index_tr, sigma_index_rot, sigma_index_tor, sigma_index_tr, sigma_index_rot,
                     sigma_index_tor, sigma_index_tr])

        except RuntimeError as e:
            if 'out of memory' in str(e):
                print('| WARNING: ran out of memory, skipping batch')
                for p in model.parameters():
                    if p.grad is not None:
                        del p.grad  # free some memory
                torch.cuda.empty_cache()
                continue
            elif 'Input mismatch' in str(e):
                print('| WARNING: weird torch_cluster error, skipping batch')
                for p in model.parameters():
                    if p.grad is not None:
                        del p.grad  # free some memory
                torch.cuda.empty_cache()
                continue
            else:
                raise e

    out = meter.summary()
    if test_sigma_intervals > 0: out.update(meter_all.summary())
    return out


def inference_epoch(model, complex_graphs, device, t_to_sigma, args):
    t_schedule = get_t_schedule(inference_steps=args.inference_steps)
    tr_schedule, rot_schedule, tor_schedule = t_schedule, t_schedule, t_schedule

    dataset = ListDataset(complex_graphs)
    loader = DataLoader(dataset=dataset, batch_size=1, shuffle=False)
    rmsds = []

    for orig_complex_graph in tqdm(loader):
        data_list = [copy.deepcopy(orig_complex_graph)]
        randomize_position(data_list, args.no_torsion, False, args.tr_sigma_max)

        predictions_list = None
        failed_convergence_counter = 0
        while predictions_list == None:
            try:
                predictions_list, confidences = sampling(data_list=data_list, model=model.module if device.type=='cuda' else model,
                                                         inference_steps=args.inference_steps,
                                                         tr_schedule=tr_schedule, rot_schedule=rot_schedule,
                                                         tor_schedule=tor_schedule,
                                                         device=device, t_to_sigma=t_to_sigma, model_args=args)
            except Exception as e:
                if 'failed to converge' in str(e):
                    failed_convergence_counter += 1
                    if failed_convergence_counter > 5:
                        print('| WARNING: SVD failed to converge 5 times - skipping the complex')
                        break
                    print('| WARNING: SVD failed to converge - trying again with a new sample')
                else:
                    raise e
        if failed_convergence_counter > 5: continue
        if args.no_torsion:
            orig_complex_graph['ligand'].orig_pos = (orig_complex_graph['ligand'].pos.cpu().numpy() +
                                                     orig_complex_graph.original_center.cpu().numpy())

        filterHs = torch.not_equal(predictions_list[0]['ligand'].x[:, 0], 0).cpu().numpy()

        if isinstance(orig_complex_graph['ligand'].orig_pos, list):
            orig_complex_graph['ligand'].orig_pos = orig_complex_graph['ligand'].orig_pos[0]

        ligand_pos = np.asarray(
            [complex_graph['ligand'].pos.cpu().numpy()[filterHs] for complex_graph in predictions_list])
        orig_ligand_pos = np.expand_dims(
            orig_complex_graph['ligand'].orig_pos[filterHs] - orig_complex_graph.original_center.cpu().numpy(), axis=0)
        rmsd = np.sqrt(((ligand_pos - orig_ligand_pos) ** 2).sum(axis=2).mean(axis=1))
        rmsds.append(rmsd)

    rmsds = np.array(rmsds)
    losses = {'rmsds_lt2': (100 * (rmsds < 2).sum() / len(rmsds)),
              'rmsds_lt5': (100 * (rmsds < 5).sum() / len(rmsds))}
    return losses
