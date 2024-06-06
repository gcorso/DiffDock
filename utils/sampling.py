import copy
import random

import numpy as np
import torch
from torch_geometric.data import Batch
from torch_geometric.loader import DataLoader

from utils.diffusion_utils import modify_conformer, set_time, modify_conformer_batch
from utils.torsion import modify_conformer_torsion_angles
from scipy.spatial.transform import Rotation as R
from utils.utils import crop_beyond
from utils.logging_utils import get_logger


def randomize_position(data_list, no_torsion, no_random, tr_sigma_max, pocket_knowledge=False, pocket_cutoff=7,
                       initial_noise_std_proportion=-1.0, choose_residue=False):
    # in place modification of the list
    center_pocket = data_list[0]['receptor'].pos.mean(dim=0)
    if pocket_knowledge:
        complex = data_list[0]
        d = torch.cdist(complex['receptor'].pos, torch.from_numpy(complex['ligand'].orig_pos[0]).float() - complex.original_center)
        label = torch.any(d < pocket_cutoff, dim=1)

        if torch.any(label):
            center_pocket = complex['receptor'].pos[label].mean(dim=0)
        else:
            print("No pocket residue below minimum distance ", pocket_cutoff, "taking closest at", torch.min(d))
            center_pocket = complex['receptor'].pos[torch.argmin(torch.min(d, dim=1)[0])]

    if not no_torsion:
        # randomize torsion angles
        for complex_graph in data_list:
            torsion_updates = np.random.uniform(low=-np.pi, high=np.pi, size=complex_graph['ligand'].edge_mask.sum())
            complex_graph['ligand'].pos = \
                modify_conformer_torsion_angles(complex_graph['ligand'].pos,
                                                complex_graph['ligand', 'ligand'].edge_index.T[
                                                    complex_graph['ligand'].edge_mask],
                                                complex_graph['ligand'].mask_rotate[0], torsion_updates)

    for complex_graph in data_list:
        # randomize position
        molecule_center = torch.mean(complex_graph['ligand'].pos, dim=0, keepdim=True)
        random_rotation = torch.from_numpy(R.random().as_matrix()).float()
        complex_graph['ligand'].pos = (complex_graph['ligand'].pos - molecule_center) @ random_rotation.T + center_pocket
        # base_rmsd = np.sqrt(np.sum((complex_graph['ligand'].pos.cpu().numpy() - orig_complex_graph['ligand'].pos.numpy()) ** 2, axis=1).mean())

        if not no_random:  # note for now the torsion angles are still randomised
            if choose_residue:
                idx = random.randint(0, len(complex_graph['receptor'].pos)-1)
                tr_update = torch.normal(mean=complex_graph['receptor'].pos[idx:idx+1], std=0.01)
            elif initial_noise_std_proportion >= 0.0:
                std_rec = torch.sqrt(torch.mean(torch.sum(complex_graph['receptor'].pos ** 2, dim=1)))
                tr_update = torch.normal(mean=0, std=std_rec * initial_noise_std_proportion / 1.73, size=(1, 3))
            else:
                # if initial_noise_std_proportion < 0.0, we use the tr_sigma_max multiplied by -initial_noise_std_proportion
                tr_update = torch.normal(mean=0, std=-initial_noise_std_proportion * tr_sigma_max, size=(1, 3))
            complex_graph['ligand'].pos += tr_update


def is_iterable(arr):
    try:
        some_object_iterator = iter(arr)
        return True
    except TypeError as te:
        return False


def sampling(data_list, model, inference_steps, tr_schedule, rot_schedule, tor_schedule, device, t_to_sigma, model_args,
             no_random=False, ode=False, visualization_list=None, confidence_model=None, confidence_data_list=None, confidence_model_args=None,
             t_schedule=None, batch_size=32, no_final_step_noise=False, pivot=None, return_full_trajectory=False,
             temp_sampling=1.0, temp_psi=0.0, temp_sigma_data=0.5, return_features=False):
    N = len(data_list)
    trajectory = []
    logger = get_logger()
    if return_features:
        lig_features, rec_features = [], []
        assert batch_size >= N, "Not implemented yet"

    loader = DataLoader(data_list, batch_size=batch_size)
    assert not (return_full_trajectory or return_features or pivot), "Not implemented yet in new inference version"

    mask_rotate = torch.from_numpy(data_list[0]['ligand'].mask_rotate[0]).to(device)

    confidence = None
    if confidence_model is not None:
        confidence_loader = iter(DataLoader(confidence_data_list, batch_size=batch_size))
        confidence = []

    with torch.no_grad():
        for batch_id, complex_graph_batch in enumerate(loader):
            b = complex_graph_batch.num_graphs
            n = len(complex_graph_batch['ligand'].pos) // b
            complex_graph_batch = complex_graph_batch.to(device)

            for t_idx in range(inference_steps):
                t_tr, t_rot, t_tor = tr_schedule[t_idx], rot_schedule[t_idx], tor_schedule[t_idx]
                dt_tr = tr_schedule[t_idx] - tr_schedule[t_idx + 1] if t_idx < inference_steps - 1 else tr_schedule[t_idx]
                dt_rot = rot_schedule[t_idx] - rot_schedule[t_idx + 1] if t_idx < inference_steps - 1 else rot_schedule[t_idx]
                dt_tor = tor_schedule[t_idx] - tor_schedule[t_idx + 1] if t_idx < inference_steps - 1 else tor_schedule[t_idx]

                tr_sigma, rot_sigma, tor_sigma = t_to_sigma(t_tr, t_rot, t_tor)

                if hasattr(model_args, 'crop_beyond') and model_args.crop_beyond is not None:
                    #print('Cropping beyond', tr_sigma * 3 + model_args.crop_beyond, 'for score model')
                    mod_complex_graph_batch = copy.deepcopy(complex_graph_batch).to_data_list()
                    for batch in mod_complex_graph_batch:
                        crop_beyond(batch, tr_sigma * 3 + model_args.crop_beyond, model_args.all_atoms)
                    mod_complex_graph_batch = Batch.from_data_list(mod_complex_graph_batch)
                else:
                    mod_complex_graph_batch = complex_graph_batch

                set_time(mod_complex_graph_batch, t_schedule[t_idx] if t_schedule is not None else None, t_tr, t_rot, t_tor, b,
                         'all_atoms' in model_args and model_args.all_atoms, device)

                tr_score, rot_score, tor_score = model(mod_complex_graph_batch)[:3]
                mean_scores = torch.mean(tr_score, dim=-1)
                num_nans = torch.sum(torch.isnan(mean_scores))
                if num_nans > 0:
                    name = complex_graph_batch['name']
                    if isinstance(name, list):
                        name = name[0]
                    logger.warning(f"Complex {name} Batch {batch_id+1} Inference Iteration {t_idx}: "
                                   f"{num_nans} / {mean_scores.numel()} samples failed")

                    # Set the nan values to a small value, just want to disturb slightly
                    # Hopefully won't get nan the next iteration
                    tr_score.nan_to_num_(nan=(eps := 0.01*torch.nanmean(tr_score.abs())), posinf=eps, neginf=-eps)
                    rot_score.nan_to_num_(nan=(eps := 0.01*torch.nanmean(rot_score.abs())), posinf=eps, neginf=-eps)
                    tor_score.nan_to_num_(nan=(eps := 0.01*torch.nanmean(tor_score.abs())), posinf=eps, neginf=-eps)
                    del eps

                tr_g = tr_sigma * torch.sqrt(torch.tensor(2 * np.log(model_args.tr_sigma_max / model_args.tr_sigma_min)))
                rot_g = rot_sigma * torch.sqrt(torch.tensor(2 * np.log(model_args.rot_sigma_max / model_args.rot_sigma_min)))

                if ode:
                    tr_perturb = (0.5 * tr_g ** 2 * dt_tr * tr_score)
                    rot_perturb = (0.5 * rot_score * dt_rot * rot_g ** 2)
                else:
                    tr_z = torch.zeros((min(batch_size, N), 3), device=device) if no_random or (no_final_step_noise and t_idx == inference_steps - 1) \
                        else torch.normal(mean=0, std=1, size=(min(batch_size, N), 3), device=device)
                    tr_perturb = (tr_g ** 2 * dt_tr * tr_score + tr_g * np.sqrt(dt_tr) * tr_z)

                    rot_z = torch.zeros((min(batch_size, N), 3), device=device) if no_random or (no_final_step_noise and t_idx == inference_steps - 1) \
                        else torch.normal(mean=0, std=1, size=(min(batch_size, N), 3), device=device)
                    rot_perturb = (rot_score * dt_rot * rot_g ** 2 + rot_g * np.sqrt(dt_rot) * rot_z)

                if not model_args.no_torsion:
                    tor_g = tor_sigma * torch.sqrt(torch.tensor(2 * np.log(model_args.tor_sigma_max / model_args.tor_sigma_min)))
                    if ode:
                        tor_perturb = (0.5 * tor_g ** 2 * dt_tor * tor_score)
                    else:
                        tor_z = torch.zeros(tor_score.shape, device=device) if no_random or (no_final_step_noise and t_idx == inference_steps - 1) \
                            else torch.normal(mean=0, std=1, size=tor_score.shape, device=device)
                        tor_perturb = (tor_g ** 2 * dt_tor * tor_score + tor_g * np.sqrt(dt_tor) * tor_z)
                    torsions_per_molecule = tor_perturb.shape[0] // b
                else:
                    tor_perturb = None

                if not is_iterable(temp_sampling):
                    temp_sampling = [temp_sampling] * 3
                if not is_iterable(temp_psi):
                    temp_psi = [temp_psi] * 3

                if not is_iterable(temp_sampling): temp_sampling = [temp_sampling] * 3
                if not is_iterable(temp_psi): temp_psi = [temp_psi] * 3
                if not is_iterable(temp_sigma_data): temp_sigma_data = [temp_sigma_data] * 3

                assert len(temp_sampling) == 3
                assert len(temp_psi) == 3
                assert len(temp_sigma_data) == 3

                if temp_sampling[0] != 1.0:
                    tr_sigma_data = np.exp(temp_sigma_data[0] * np.log(model_args.tr_sigma_max) + (1 - temp_sigma_data[0]) * np.log(model_args.tr_sigma_min))
                    lambda_tr = (tr_sigma_data + tr_sigma) / (tr_sigma_data + tr_sigma / temp_sampling[0])
                    tr_perturb = (tr_g ** 2 * dt_tr * (lambda_tr + temp_sampling[0] * temp_psi[0] / 2) * tr_score + tr_g * np.sqrt(dt_tr * (1 + temp_psi[0])) * tr_z)

                if temp_sampling[1] != 1.0:
                    rot_sigma_data = np.exp(temp_sigma_data[1] * np.log(model_args.rot_sigma_max) + (1 - temp_sigma_data[1]) * np.log(model_args.rot_sigma_min))
                    lambda_rot = (rot_sigma_data + rot_sigma) / (rot_sigma_data + rot_sigma / temp_sampling[1])
                    rot_perturb = (rot_g ** 2 * dt_rot * (lambda_rot + temp_sampling[1] * temp_psi[1] / 2) * rot_score + rot_g * np.sqrt(dt_rot * (1 + temp_psi[1])) * rot_z)

                if temp_sampling[2] != 1.0:
                    tor_sigma_data = np.exp(temp_sigma_data[2] * np.log(model_args.tor_sigma_max) + (1 - temp_sigma_data[2]) * np.log(model_args.tor_sigma_min))
                    lambda_tor = (tor_sigma_data + tor_sigma) / (tor_sigma_data + tor_sigma / temp_sampling[2])
                    tor_perturb = (tor_g ** 2 * dt_tor * (lambda_tor + temp_sampling[2] * temp_psi[2] / 2) * tor_score + tor_g * np.sqrt(dt_tor * (1 + temp_psi[2])) * tor_z)

                # Apply noise
                complex_graph_batch['ligand'].pos = \
                    modify_conformer_batch(complex_graph_batch['ligand'].pos, complex_graph_batch, tr_perturb, rot_perturb,
                                           tor_perturb if not model_args.no_torsion else None, mask_rotate)

                if visualization_list is not None:
                    for idx_b in range(b):
                        visualization_list[batch_id * batch_size + idx_b].add((
                                complex_graph_batch['ligand'].pos[idx_b*n:n*(idx_b+1)].detach().cpu() +
                                data_list[batch_id * batch_size + idx_b].original_center.detach().cpu()),
                                part=1, order=t_idx + 2)

            for i in range(b):
               data_list[batch_id * batch_size + i]['ligand'].pos = complex_graph_batch['ligand'].pos[i*n:n*(i+1)]

            if visualization_list is not None:
                for idx, visualization in enumerate(visualization_list):
                    visualization.add((data_list[idx]['ligand'].pos.detach().cpu() + data_list[idx].original_center.detach().cpu()),
                                      part=1, order=2)

            if confidence_model is not None:
                if confidence_data_list is not None:
                    confidence_complex_graph_batch = next(confidence_loader)
                    confidence_complex_graph_batch['ligand'].pos = complex_graph_batch['ligand'].pos.cpu()

                    if hasattr(confidence_model_args, 'crop_beyond') and confidence_model_args.crop_beyond is not None:
                        confidence_complex_graph_batch = confidence_complex_graph_batch.to_data_list()
                        for batch in confidence_complex_graph_batch:
                            crop_beyond(batch, confidence_model_args.crop_beyond, confidence_model_args.all_atoms)
                        confidence_complex_graph_batch = Batch.from_data_list(confidence_complex_graph_batch)

                    confidence_complex_graph_batch = confidence_complex_graph_batch.to(device)
                    set_time(confidence_complex_graph_batch, 0, 0, 0, 0, b, confidence_model_args.all_atoms, device)
                    out = confidence_model(confidence_complex_graph_batch)
                else:
                    out = confidence_model(complex_graph_batch)

                if type(out) is tuple:
                    out = out[0]
                confidence.append(out)

    if confidence_model is not None:
        confidence = torch.cat(confidence, dim=0)
        confidence = torch.nan_to_num(confidence, nan=-1000)

    if return_full_trajectory:
        return data_list, confidence, trajectory
    elif return_features:
        lig_features = torch.cat(lig_features, dim=0)
        rec_features = torch.cat(rec_features, dim=0)
        return data_list, confidence, lig_features, rec_features

    return data_list, confidence


def compute_affinity(data_list, affinity_model, affinity_data_list, device, parallel, all_atoms, include_miscellaneous_atoms):

    with torch.no_grad():
        if affinity_model is not None:
            assert parallel <= len(data_list)
            loader = DataLoader(data_list, batch_size=parallel)
            complex_graph_batch = next(iter(loader)).to(device)
            positions = complex_graph_batch['ligand'].pos

            assert affinity_data_list is not None
            complex_graph = affinity_data_list[0]
            N = complex_graph['ligand'].num_nodes
            complex_graph['ligand'].x = complex_graph['ligand'].x.repeat(parallel, 1)
            complex_graph['ligand'].edge_mask = complex_graph['ligand'].edge_mask.repeat(parallel)
            complex_graph['ligand', 'ligand'].edge_index = torch.cat(
                [N * i + complex_graph['ligand', 'ligand'].edge_index for i in range(parallel)], dim=1)
            complex_graph['ligand', 'ligand'].edge_attr = complex_graph['ligand', 'ligand'].edge_attr.repeat(parallel, 1)
            complex_graph['ligand'].pos = positions

            affinity_loader = DataLoader([complex_graph], batch_size=1)
            affinity_batch = next(iter(affinity_loader)).to(device)
            set_time(affinity_batch, 0, 0, 0, 0, 1, all_atoms, device, include_miscellaneous_atoms=include_miscellaneous_atoms)
            _, affinity = affinity_model(affinity_batch)
        else:
            affinity = None

    return affinity