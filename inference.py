import copy
import os
import torch

import time
from argparse import ArgumentParser, Namespace, FileType
from rdkit.Chem import RemoveHs
from functools import partial
import numpy as np
import pandas as pd
from rdkit import RDLogger
from rdkit.Chem import MolFromSmiles, AddHs
from torch_geometric.loader import DataLoader

from datasets.process_mols import read_molecule, generate_conformer, write_mol_with_coords
from datasets.pdbbind import PDBBind
from utils.diffusion_utils import t_to_sigma as t_to_sigma_compl, get_t_schedule
from utils.sampling import randomize_position, sampling
from utils.utils import get_model
from utils.visualise import PDBFile
from tqdm import tqdm

RDLogger.DisableLog('rdApp.*')
import yaml
parser = ArgumentParser()
parser.add_argument('--config', type=FileType(mode='r'), default=None)
parser.add_argument('--protein_ligand_csv', type=str, default=None, help='Path to a .csv file specifying the input as described in the README. If this is not None, it will be used instead of the --protein_path and --ligand parameters')
parser.add_argument('--protein_path', type=str, default='data/dummy_data/1a0q_protein.pdb', help='Path to the protein .pdb file')
parser.add_argument('--ligand', type=str, default='COc(cc1)ccc1C#N', help='Either a SMILES string or the path to a molecule file that rdkit can read')
parser.add_argument('--out_dir', type=str, default='results/user_inference', help='Directory where the outputs will be written to')
parser.add_argument('--esm_embeddings_path', type=str, default='data/esm2_output', help='If this is set then the LM embeddings at that path will be used for the receptor features')
parser.add_argument('--save_visualisation', action='store_true', default=False, help='Save a pdb file with all of the steps of the reverse diffusion')
parser.add_argument('--samples_per_complex', type=int, default=10, help='Number of samples to generate')

parser.add_argument('--model_dir', type=str, default='workdir/paper_score_model', help='Path to folder with trained score model and hyperparameters')
parser.add_argument('--ckpt', type=str, default='best_ema_inference_epoch_model.pt', help='Checkpoint to use for the score model')
parser.add_argument('--confidence_model_dir', type=str, default='workdir/paper_confidence_model', help='Path to folder with trained confidence model and hyperparameters')
parser.add_argument('--confidence_ckpt', type=str, default='best_model_epoch75.pt', help='Checkpoint to use for the confidence model')

parser.add_argument('--batch_size', type=int, default=32, help='')
parser.add_argument('--cache_path', type=str, default='data/cache', help='Folder from where to load/restore cached dataset')
parser.add_argument('--no_random', action='store_true', default=False, help='Use no randomness in reverse diffusion')
parser.add_argument('--no_final_step_noise', action='store_true', default=False, help='Use no noise in the final step of the reverse diffusion')
parser.add_argument('--ode', action='store_true', default=False, help='Use ODE formulation for inference')
parser.add_argument('--inference_steps', type=int, default=20, help='Number of denoising steps')
parser.add_argument('--num_workers', type=int, default=1, help='Number of workers for creating the dataset')
parser.add_argument('--sigma_schedule', type=str, default='expbeta', help='')
parser.add_argument('--actual_steps', type=int, default=None, help='Number of denoising steps that are actually performed')
parser.add_argument('--keep_local_structures', action='store_true', default=False, help='Keeps the local structure when specifying an input with 3D coordinates instead of generating them with RDKit')
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

os.makedirs(args.out_dir, exist_ok=True)

with open(f'{args.model_dir}/model_parameters.yml') as f:
    score_model_args = Namespace(**yaml.full_load(f))

if args.confidence_model_dir is not None:
    with open(f'{args.confidence_model_dir}/model_parameters.yml') as f:
        confidence_args = Namespace(**yaml.full_load(f))

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

if args.protein_ligand_csv is not None:
    df = pd.read_csv(args.protein_ligand_csv)
    protein_path_list = df['protein'].tolist()
    ligand_descriptions = df['ligand'].tolist()
else:
    protein_path_list = [args.protein_path]
    ligand_descriptions = [args.ligand]


test_dataset = PDBBind(transform=None, root='', protein_path_list=protein_path_list, ligand_descriptions=ligand_descriptions,
                       receptor_radius=score_model_args.receptor_radius, cache_path=args.cache_path,
                       remove_hs=score_model_args.remove_hs, max_lig_size=None,
                       c_alpha_max_neighbors=score_model_args.c_alpha_max_neighbors, matching=False, keep_original=False,
                       popsize=score_model_args.matching_popsize, maxiter=score_model_args.matching_maxiter,
                       all_atoms=score_model_args.all_atoms, atom_radius=score_model_args.atom_radius,
                       atom_max_neighbors=score_model_args.atom_max_neighbors,
                       esm_embeddings_path= args.esm_embeddings_path if score_model_args.esm_embeddings_path is not None else None,
                       require_ligand=True, num_workers=args.num_workers, keep_local_structures=args.keep_local_structures)
test_loader = DataLoader(dataset=test_dataset, batch_size=1, shuffle=False)

if args.confidence_model_dir is not None:
    if not (confidence_args.use_original_model_cache or confidence_args.transfer_weights): # if the confidence model uses the same type of data as the original model then we do not need this dataset and can just use the complexes
        print('HAPPENING | confidence model uses different type of graphs than the score model. Loading (or creating if not existing) the data for the confidence model now.')
        confidence_test_dataset = PDBBind(transform=None, root='', protein_path_list=protein_path_list,
                                         ligand_descriptions=ligand_descriptions, receptor_radius=confidence_args.receptor_radius,
                                         cache_path=args.cache_path, remove_hs=confidence_args.remove_hs, max_lig_size=None,
                                         c_alpha_max_neighbors=confidence_args.c_alpha_max_neighbors, matching=False, keep_original=False,
                                         popsize=confidence_args.matching_popsize, maxiter=confidence_args.matching_maxiter,
                                         all_atoms=confidence_args.all_atoms, atom_radius=confidence_args.atom_radius,
                                         atom_max_neighbors=confidence_args.atom_max_neighbors,
                                         esm_embeddings_path=args.esm_embeddings_path if confidence_args.esm_embeddings_path is not None else None,
                                         require_ligand=True, num_workers=args.num_workers)
        confidence_complex_dict = {d.name: d for d in confidence_test_dataset}

t_to_sigma = partial(t_to_sigma_compl, args=score_model_args)

model = get_model(score_model_args, device, t_to_sigma=t_to_sigma, no_parallel=True)
state_dict = torch.load(f'{args.model_dir}/{args.ckpt}', map_location=torch.device('cpu'))
model.load_state_dict(state_dict, strict=True)
model = model.to(device)
model.eval()

if args.confidence_model_dir is not None:
    if confidence_args.transfer_weights:
        with open(f'{confidence_args.original_model_dir}/model_parameters.yml') as f:
            confidence_model_args = Namespace(**yaml.full_load(f))
    else:
        confidence_model_args = confidence_args
    confidence_model = get_model(confidence_model_args, device, t_to_sigma=t_to_sigma, no_parallel=True, confidence_mode=True)
    state_dict = torch.load(f'{args.confidence_model_dir}/{args.confidence_ckpt}', map_location=torch.device('cpu'))
    confidence_model.load_state_dict(state_dict, strict=True)
    confidence_model = confidence_model.to(device)
    confidence_model.eval()
else:
    confidence_model = None
    confidence_args = None
    confidence_model_args = None

tr_schedule = get_t_schedule(inference_steps=args.inference_steps)
rot_schedule = tr_schedule
tor_schedule = tr_schedule
print('common t schedule', tr_schedule)

failures, skipped, confidences_list, names_list, run_times, min_self_distances_list = 0, 0, [], [], [], []
N = args.samples_per_complex
print('Size of test dataset: ', len(test_dataset))
for idx, orig_complex_graph in tqdm(enumerate(test_loader)):
    if confidence_model is not None and not (confidence_args.use_original_model_cache or confidence_args.transfer_weights) and orig_complex_graph.name[0] not in confidence_complex_dict.keys():
        skipped += 1
        print(f"HAPPENING | The confidence dataset did not contain {orig_complex_graph.name[0]}. We are skipping this complex.")
        continue
    try:
        data_list = [copy.deepcopy(orig_complex_graph) for _ in range(N)]
        randomize_position(data_list, score_model_args.no_torsion, args.no_random,score_model_args.tr_sigma_max)
        pdb = None
        lig = orig_complex_graph.mol[0]
        if args.save_visualisation:
            visualization_list = []
            for graph in data_list:
                pdb = PDBFile(lig)
                pdb.add(lig, 0, 0)
                pdb.add((orig_complex_graph['ligand'].pos + orig_complex_graph.original_center).detach().cpu(), 1, 0)
                pdb.add((graph['ligand'].pos + graph.original_center).detach().cpu(), part=1, order=1)
                visualization_list.append(pdb)
        else:
            visualization_list = None

        start_time = time.time()
        if confidence_model is not None and not (confidence_args.use_original_model_cache or confidence_args.transfer_weights):
            confidence_data_list = [copy.deepcopy(confidence_complex_dict[orig_complex_graph.name[0]]) for _ in range(N)]
        else:
            confidence_data_list = None

        data_list, confidence = sampling(data_list=data_list, model=model,
                                         inference_steps=args.actual_steps if args.actual_steps is not None else args.inference_steps,
                                         tr_schedule=tr_schedule, rot_schedule=rot_schedule, tor_schedule=tor_schedule,
                                         device=device, t_to_sigma=t_to_sigma, model_args=score_model_args, no_random=args.no_random,
                                         ode=args.ode, visualization_list=visualization_list, confidence_model=confidence_model,
                                         confidence_data_list=confidence_data_list, confidence_model_args=confidence_model_args,
                                         batch_size=args.batch_size, no_final_step_noise=args.no_final_step_noise)
        ligand_pos = np.asarray([complex_graph['ligand'].pos.cpu().numpy() + orig_complex_graph.original_center.cpu().numpy() for complex_graph in data_list])
        run_times.append(time.time() - start_time)

        if confidence is not None and isinstance(confidence_args.rmsd_classification_cutoff, list):
            confidence = confidence[:,0]
        if confidence is not None:
            confidence = confidence.cpu().numpy()
            re_order = np.argsort(confidence)[::-1]
            confidence = confidence[re_order]
            confidences_list.append(confidence)
            ligand_pos = ligand_pos[re_order]
        write_dir = f'{args.out_dir}/index{idx}_{data_list[0]["name"][0].replace("/","-")}'
        os.makedirs(write_dir, exist_ok=True)
        for rank, pos in enumerate(ligand_pos):
            mol_pred = copy.deepcopy(lig)
            if score_model_args.remove_hs: mol_pred = RemoveHs(mol_pred)
            if rank == 0: write_mol_with_coords(mol_pred, pos, os.path.join(write_dir, f'rank{rank+1}.sdf'))
            write_mol_with_coords(mol_pred, pos, os.path.join(write_dir, f'rank{rank+1}_confidence{confidence[rank]:.2f}.sdf'))
        self_distances = np.linalg.norm(ligand_pos[:, :, None, :] - ligand_pos[:, None, :, :], axis=-1)
        self_distances = np.where(np.eye(self_distances.shape[2]), np.inf, self_distances)
        min_self_distances_list.append(np.min(self_distances, axis=(1, 2)))

        if args.save_visualisation:
            if confidence is not None:
                for rank, batch_idx in enumerate(re_order):
                    visualization_list[batch_idx].write(os.path.join(write_dir, f'rank{rank+1}_reverseprocess.pdb'))
            else:
                for rank, batch_idx in enumerate(ligand_pos):
                    visualization_list[batch_idx].write(os.path.join(write_dir, f'rank{rank+1}_reverseprocess.pdb'))
        names_list.append(orig_complex_graph.name[0])
    except Exception as e:
        print("Failed on", orig_complex_graph["name"], e)
        failures += 1

print(f'Failed for {failures} complexes')
print(f'Skipped {skipped} complexes')

min_self_distances = np.array(min_self_distances_list)
confidences = np.array(confidences_list)
names = np.array(names_list)
run_times = np.array(run_times)
np.save(f'{args.out_dir}/min_self_distances.npy', min_self_distances)
np.save(f'{args.out_dir}/confidences.npy', confidences)
np.save(f'{args.out_dir}/run_times.npy', run_times)
np.save(f'{args.out_dir}/complex_names.npy', np.array(names))

print(f'Results are in {args.out_dir}')


