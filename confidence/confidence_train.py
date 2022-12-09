import gc
import math
import os

import shutil

from argparse import Namespace, ArgumentParser, FileType
import torch.nn.functional as F

import wandb
import torch
from sklearn.metrics import roc_auc_score
from torch_geometric.loader import DataListLoader, DataLoader
from tqdm import tqdm

from confidence.dataset import ConfidenceDataset
from utils.training import AverageMeter

torch.multiprocessing.set_sharing_strategy('file_system')

import yaml
from utils.utils import save_yaml_file, get_optimizer_and_scheduler, get_model


parser = ArgumentParser()
parser.add_argument('--config', type=FileType(mode='r'), default=None)
parser.add_argument('--original_model_dir', type=str, default='workdir', help='Path to folder with trained model and hyperparameters')
parser.add_argument('--restart_dir', type=str, default=None, help='')
parser.add_argument('--use_original_model_cache', action='store_true', default=False, help='If this is true, the same dataset as in the original model will be used. Otherwise, the dataset parameters are used.')
parser.add_argument('--data_dir', type=str, default='data/PDBBind_processed/', help='Folder containing original structures')
parser.add_argument('--ckpt', type=str, default='best_model.pt', help='Checkpoint to use inside the folder')
parser.add_argument('--model_save_frequency', type=int, default=0, help='Frequency with which to save the last model. If 0, then only the early stopping criterion best model is saved and overwritten.')
parser.add_argument('--best_model_save_frequency', type=int, default=0, help='Frequency with which to save the best model. If 0, then only the early stopping criterion best model is saved and overwritten.')
parser.add_argument('--run_name', type=str, default='test_confidence', help='')
parser.add_argument('--project', type=str, default='diffdock_confidence', help='')
parser.add_argument('--split_train', type=str, default='data/splits/timesplit_no_lig_overlap_train', help='Path of file defining the split')
parser.add_argument('--split_val', type=str, default='data/splits/timesplit_no_lig_overlap_val', help='Path of file defining the split')
parser.add_argument('--split_test', type=str, default='data/splits/timesplit_test', help='Path of file defining the split')

# Inference parameters for creating the positions and rmsds that the confidence predictor will be trained on.
parser.add_argument('--cache_path', type=str, default='data/cacheNew', help='Folder from where to load/restore cached dataset')
parser.add_argument('--cache_ids_to_combine', nargs='+', type=str, default=None, help='RMSD value below which a prediction is considered a postitive. This can also be multiple cutoffs.')
parser.add_argument('--cache_creation_id', type=int, default=None, help='number of times that inference is run on the full dataset before concatenating it and coming up with the full confidence dataset')
parser.add_argument('--wandb', action='store_true', default=False, help='')
parser.add_argument('--inference_steps', type=int, default=2, help='Number of denoising steps')
parser.add_argument('--samples_per_complex', type=int, default=3, help='')
parser.add_argument('--balance', action='store_true', default=False, help='If this is true than we do not force the samples seen during training to be the same amount of negatives as positives')
parser.add_argument('--rmsd_prediction', action='store_true', default=False, help='')
parser.add_argument('--rmsd_classification_cutoff', nargs='+', type=float, default=2, help='RMSD value below which a prediction is considered a postitive. This can also be multiple cutoffs.')

parser.add_argument('--log_dir', type=str, default='workdir', help='')
parser.add_argument('--main_metric', type=str, default='accuracy', help='Metric to track for early stopping. Mostly [loss, accuracy, ROC AUC]')
parser.add_argument('--main_metric_goal', type=str, default='max', help='Can be [min, max]')
parser.add_argument('--transfer_weights', action='store_true', default=False, help='')
parser.add_argument('--batch_size', type=int, default=5, help='')
parser.add_argument('--lr', type=float, default=1e-3, help='')
parser.add_argument('--w_decay', type=float, default=0.0, help='')
parser.add_argument('--scheduler', type=str, default='plateau', help='')
parser.add_argument('--scheduler_patience', type=int, default=20, help='')
parser.add_argument('--n_epochs', type=int, default=5, help='')

# Dataset
parser.add_argument('--limit_complexes', type=int, default=0, help='')
parser.add_argument('--all_atoms', action='store_true', default=True, help='')
parser.add_argument('--multiplicity', type=int, default=1, help='')
parser.add_argument('--chain_cutoff', type=float, default=10, help='')
parser.add_argument('--receptor_radius', type=float, default=30, help='')
parser.add_argument('--c_alpha_max_neighbors', type=int, default=10, help='')
parser.add_argument('--atom_radius', type=float, default=5, help='')
parser.add_argument('--atom_max_neighbors', type=int, default=8, help='')
parser.add_argument('--matching_popsize', type=int, default=20, help='')
parser.add_argument('--matching_maxiter', type=int, default=20, help='')
parser.add_argument('--max_lig_size', type=int, default=None, help='Maximum number of heavy atoms')
parser.add_argument('--remove_hs', action='store_true', default=False, help='remove Hs')
parser.add_argument('--num_conformers', type=int, default=1, help='')
parser.add_argument('--esm_embeddings_path', type=str, default=None,help='If this is set then the LM embeddings at that path will be used for the receptor features')
parser.add_argument('--no_torsion', action='store_true', default=False, help='')

# Model
parser.add_argument('--num_conv_layers', type=int, default=2, help='Number of interaction layers')
parser.add_argument('--max_radius', type=float, default=5.0, help='Radius cutoff for geometric graph')
parser.add_argument('--scale_by_sigma', action='store_true', default=True, help='Whether to normalise the score')
parser.add_argument('--ns', type=int, default=16, help='Number of hidden features per node of order 0')
parser.add_argument('--nv', type=int, default=4, help='Number of hidden features per node of order >0')
parser.add_argument('--distance_embed_dim', type=int, default=32, help='')
parser.add_argument('--cross_distance_embed_dim', type=int, default=32, help='')
parser.add_argument('--no_batch_norm', action='store_true', default=False, help='If set, it removes the batch norm')
parser.add_argument('--use_second_order_repr', action='store_true', default=False, help='Whether to use only up to first order representations or also second')
parser.add_argument('--cross_max_distance', type=float, default=80, help='')
parser.add_argument('--dynamic_max_cross', action='store_true', default=False, help='')
parser.add_argument('--dropout', type=float, default=0.0, help='MLP dropout')
parser.add_argument('--embedding_type', type=str, default="sinusoidal", help='')
parser.add_argument('--sigma_embed_dim', type=int, default=32, help='')
parser.add_argument('--embedding_scale', type=int, default=10000, help='')
parser.add_argument('--confidence_no_batchnorm', action='store_true', default=False, help='')
parser.add_argument('--confidence_dropout', type=float, default=0.0, help='MLP dropout in confidence readout')

args = parser.parse_args()
if args.config:
    config_dict = yaml.load(args.config, Loader=yaml.FullLoader)
    arg_dict = args.__dict__
    for key, value in config_dict.items():
        if isinstance(value, list):
            for v in value:
                arg_dict[key].append(v)
        else:
            arg_dict[key] = value
    args.config = args.config.name
assert(args.main_metric_goal == 'max' or args.main_metric_goal == 'min')

def train_epoch(model, loader, optimizer, rmsd_prediction):
    model.train()
    meter = AverageMeter(['confidence_loss'])

    for data in tqdm(loader, total=len(loader)):
        if device.type == 'cuda' and len(data) % torch.cuda.device_count() == 1 or device.type == 'cpu' and data.num_graphs == 1:
            print("Skipping batch of size 1 since otherwise batchnorm would not work.")
        optimizer.zero_grad()
        try:
            pred = model(data)
            if rmsd_prediction:
                labels = torch.cat([graph.rmsd for graph in data]).to(device) if isinstance(data, list) else data.rmsd
                confidence_loss = F.mse_loss(pred, labels)
            else:
                if isinstance(args.rmsd_classification_cutoff, list):
                    labels = torch.cat([graph.y_binned for graph in data]).to(device) if isinstance(data, list) else data.y_binned
                    confidence_loss = F.cross_entropy(pred, labels)
                else:
                    labels = torch.cat([graph.y for graph in data]).to(device) if isinstance(data, list) else data.y
                    confidence_loss = F.binary_cross_entropy_with_logits(pred, labels)
            confidence_loss.backward()
            optimizer.step()
            meter.add([confidence_loss.cpu().detach()])
        except RuntimeError as e:
            if 'out of memory' in str(e):
                print('| WARNING: ran out of memory, skipping batch')
                for p in model.parameters():
                    if p.grad is not None:
                        del p.grad  # free some memory
                torch.cuda.empty_cache()
                gc.collect()
                continue
            else:
                raise e

    return meter.summary()

def test_epoch(model, loader, rmsd_prediction):
    model.eval()
    meter = AverageMeter(['loss'], unpooled_metrics=True) if rmsd_prediction else AverageMeter(['confidence_loss', 'accuracy', 'ROC AUC'], unpooled_metrics=True)
    all_labels = []
    for data in tqdm(loader, total=len(loader)):
        try:
            with torch.no_grad():
                pred = model(data)
            affinity_loss = torch.tensor(0.0, dtype=torch.float, device=pred[0].device)
            accuracy = torch.tensor(0.0, dtype=torch.float, device=pred[0].device)
            if rmsd_prediction:
                labels = torch.cat([graph.rmsd for graph in data]).to(device) if isinstance(data, list) else data.rmsd
                confidence_loss = F.mse_loss(pred, labels)
                meter.add([confidence_loss.cpu().detach()])
            else:
                if isinstance(args.rmsd_classification_cutoff, list):
                    labels = torch.cat([graph.y_binned for graph in data]).to(device) if isinstance(data,list) else data.y_binned
                    confidence_loss = F.cross_entropy(pred, labels)
                else:
                    labels = torch.cat([graph.y for graph in data]).to(device) if isinstance(data, list) else data.y
                    confidence_loss = F.binary_cross_entropy_with_logits(pred, labels)
                    accuracy = torch.mean((labels == (pred > 0).float()).float())
                try:
                    roc_auc = roc_auc_score(labels.detach().cpu().numpy(), pred.detach().cpu().numpy())
                except ValueError as e:
                    if 'Only one class present in y_true. ROC AUC score is not defined in that case.' in str(e):
                        roc_auc = 0
                    else:
                        raise e
            meter.add([confidence_loss.cpu().detach(), accuracy.cpu().detach(), torch.tensor(roc_auc)])
            all_labels.append(labels)

        except RuntimeError as e:
            if 'out of memory' in str(e):
                print('| WARNING: ran out of memory, skipping batch')
                for p in model.parameters():
                    if p.grad is not None:
                        del p.grad  # free some memory
                torch.cuda.empty_cache()
                continue
            else:
                raise e

    all_labels = torch.cat(all_labels)

    if rmsd_prediction:
        baseline_metric = ((all_labels - all_labels.mean()).abs()).mean()
    else:
        baseline_metric = all_labels.sum() / len(all_labels)
    results = meter.summary()
    results.update({'baseline_metric': baseline_metric})
    return meter.summary(), baseline_metric


def train(args, model, optimizer, scheduler, train_loader, val_loader, run_dir):
    best_val_metric = math.inf if args.main_metric_goal == 'min' else 0
    best_epoch = 0

    print("Starting training...")
    for epoch in range(args.n_epochs):
        logs = {}
        train_metrics = train_epoch(model, train_loader, optimizer, args.rmsd_prediction)
        print("Epoch {}: Training loss {:.4f}".format(epoch, train_metrics['confidence_loss']))

        val_metrics, baseline_metric = test_epoch(model, val_loader, args.rmsd_prediction)
        if args.rmsd_prediction:
            print("Epoch {}: Validation loss {:.4f}".format(epoch, val_metrics['confidence_loss']))
        else:
            print("Epoch {}: Validation loss {:.4f}  accuracy {:.4f}".format(epoch, val_metrics['confidence_loss'], val_metrics['accuracy']))

        if args.wandb:
            logs.update({'valinf_' + k: v for k, v in val_metrics.items()}, step=epoch + 1)
            logs.update({'train_' + k: v for k, v in train_metrics.items()}, step=epoch + 1)
            logs.update({'mean_rmsd' if args.rmsd_prediction else 'fraction_positives': baseline_metric,
                         'current_lr': optimizer.param_groups[0]['lr']})
            wandb.log(logs, step=epoch + 1)

        if scheduler:
            scheduler.step(val_metrics[args.main_metric])

        state_dict = model.module.state_dict() if device.type == 'cuda' else model.state_dict()

        if args.main_metric_goal == 'min' and val_metrics[args.main_metric] < best_val_metric or \
                args.main_metric_goal == 'max' and val_metrics[args.main_metric] > best_val_metric:
            best_val_metric = val_metrics[args.main_metric]
            best_epoch = epoch
            torch.save(state_dict, os.path.join(run_dir, 'best_model.pt'))
        if args.model_save_frequency > 0 and (epoch + 1) % args.model_save_frequency == 0:
            torch.save(state_dict, os.path.join(run_dir, f'model_epoch{epoch+1}.pt'))
        if args.best_model_save_frequency > 0 and (epoch + 1) % args.best_model_save_frequency == 0:
            shutil.copyfile(os.path.join(run_dir, 'best_model.pt'), os.path.join(run_dir, f'best_model_epoch{epoch+1}.pt'))

        torch.save({
            'epoch': epoch,
            'model': state_dict,
            'optimizer': optimizer.state_dict(),
        }, os.path.join(run_dir, 'last_model.pt'))

    print("Best Validation accuracy {} on Epoch {}".format(best_val_metric, best_epoch))


def construct_loader_confidence(args, device):
    common_args = {'cache_path': args.cache_path, 'original_model_dir': args.original_model_dir, 'device': device,
                   'inference_steps': args.inference_steps, 'samples_per_complex': args.samples_per_complex,
                   'limit_complexes': args.limit_complexes, 'all_atoms': args.all_atoms, 'balance': args.balance,
                   'rmsd_classification_cutoff': args.rmsd_classification_cutoff, 'use_original_model_cache': args.use_original_model_cache,
                   'cache_creation_id': args.cache_creation_id, "cache_ids_to_combine": args.cache_ids_to_combine,
                   "model_ckpt": args.ckpt}
    loader_class = DataListLoader if torch.cuda.is_available() else DataLoader

    exception_flag = False
    try:
        train_dataset = ConfidenceDataset(split="train", args=args, **common_args)
        train_loader = loader_class(dataset=train_dataset, batch_size=args.batch_size, shuffle=True)
    except Exception as e:
        if 'The generated ligand positions with cache_id do not exist:' in str(e):
            print("HAPPENING | Encountered the following exception when loading the confidence train dataset:")
            print(str(e))
            print("HAPPENING | We are still continuing because we want to try to generate the validation dataset if it has not been created yet:")
            exception_flag = True
        else: raise e

    val_dataset = ConfidenceDataset(split="val", args=args, **common_args)
    val_loader = loader_class(dataset=val_dataset, batch_size=args.batch_size, shuffle=True)

    if exception_flag: raise Exception('We encountered the exception during train dataset loading: ', e)
    return train_loader, val_loader


if __name__ == '__main__':
    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    with open(f'{args.original_model_dir}/model_parameters.yml') as f:
        score_model_args = Namespace(**yaml.full_load(f))

    # construct loader
    train_loader, val_loader = construct_loader_confidence(args, device)
    model = get_model(score_model_args if args.transfer_weights else args, device, t_to_sigma=None, confidence_mode=True)
    optimizer, scheduler = get_optimizer_and_scheduler(args, model, scheduler_mode=args.main_metric_goal)

    if args.transfer_weights:
        print("HAPPENING | Transferring weights from original_model_dir to the new model after using original_model_dir's arguments to construct the new model.")
        checkpoint = torch.load(os.path.join(args.original_model_dir,args.ckpt), map_location=device)
        model_state_dict = model.state_dict()
        transfer_weights_dict = {k: v for k, v in checkpoint.items() if k in list(model_state_dict.keys())}
        model_state_dict.update(transfer_weights_dict)  # update the layers with the pretrained weights
        model.load_state_dict(model_state_dict)

    elif args.restart_dir:
        dict = torch.load(f'{args.restart_dir}/last_model.pt', map_location=torch.device('cpu'))
        model.module.load_state_dict(dict['model'], strict=True)
        optimizer.load_state_dict(dict['optimizer'])
        print("Restarting from epoch", dict['epoch'])

    numel = sum([p.numel() for p in model.parameters()])
    print('Model with', numel, 'parameters')

    if args.wandb:
        wandb.init(
            entity='entity',
            settings=wandb.Settings(start_method="fork"),
            project=args.project,
            name=args.run_name,
            config=args
        )
        wandb.log({'numel': numel})

    # record parameters
    run_dir = os.path.join(args.log_dir, args.run_name)
    yaml_file_name = os.path.join(run_dir, 'model_parameters.yml')
    save_yaml_file(yaml_file_name, args.__dict__)
    args.device = device

    train(args, model, optimizer, scheduler, train_loader, val_loader, run_dir)
