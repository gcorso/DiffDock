import copy
import os
import torch
from argparse import ArgumentParser, Namespace, BooleanOptionalAction
from rdkit.Chem import RemoveHs
from functools import partial
import numpy as np
import pandas as pd
from rdkit import RDLogger
from torch_geometric.loader import DataLoader

from datasets.process_mols import write_mol_with_coords
from utils.diffusion_utils import t_to_sigma as t_to_sigma_compl, get_t_schedule
from utils.inference_utils import InferenceDataset, set_nones
from utils.sampling import randomize_position, sampling
from utils.utils import get_model
from utils.visualise import PDBFile
from tqdm import tqdm

RDLogger.DisableLog('rdApp.*')
import yaml

class InferenceArgs():
    def __init__(self, obj):
        # protein_ligand_csv: Path to a .csv file specifying the input as described in the README. If this is not None, it will be used instead of the --protein_path, --protein_sequence and --ligand parameters
        self.protein_ligand_csv = obj.get('protein_ligand_csv', None)
        # complex_name: Name that the complex will be saved with
        self.complex_name = obj.get('complex_name', '1a0q')
        # protein_path: Path to the protein file
        self.protein_path = obj.get('protein_path', None)
        # protein_sequence: Sequence of the protein for ESMFold, this is ignored if --protein_path is not None
        self.protein_sequence = obj.get('protein_sequence', None)
        # ligand_description: Either a SMILES string or the path to a molecule file that rdkit can read
        self.ligand_description = obj.get('ligand_description', 'CCCCC(NC(=O)CCC(=O)O)P(=O)(O)OC1=CC=CC=C1')
        # out_dir: Directory where the outputs will be written to
        self.out_dir = obj.get('out_dir', 'results/user_inference')
        # save_visualisation: Save a pdb file with all of the steps of the reverse diffusion
        self.save_visualisation = obj.get('save_visualisation', False)
        # samples_per_complex: Number of samples to generate
        self.samples_per_complex = obj.get('samples_per_complex', 10)
        # model_dir: Path to folder with trained score model and hyperparameters
        self.model_dir = obj.get('model_dir', 'workdir/paper_score_model')
        # ckpt: Checkpoint to use for the score model
        self.ckpt = obj.get('ckpt', 'best_ema_inference_epoch_model.pt')
        # confidence_model_dir: Path to folder with trained confidence model and hyperparameters
        self.confidence_model_dir = obj.get('confidence_model_dir', 'workdir/paper_confidence_model')
        # confidence_ckpt: Checkpoint to use for the confidence model
        self.confidence_ckpt = obj.get('confidence_ckpt', 'best_model_epoch75.pt')
        # batch_size: 
        self.batch_size = obj.get('batch_size', 32)
        # no_final_step_noise: Use no noise in the final step of the reverse diffusion
        self.no_final_step_noise = obj.get('no_final_step_noise', False)
        # inference_steps: Number of denoising steps
        self.inference_steps = obj.get('inference_steps', 20)
        # actual_steps: Number of denoising steps that are actually performed
        self.actual_steps = obj.get('actual_steps', None)

        # DiffDock: import model args

        # CONIFER POINT PARAMS
        # add_hs: Add hydrogens to poses before writing out poses
        self.add_hs = obj.get('add_hs', False)
        # keep_src_3d: Whether or not to drop the 3D coordinates for 3D input
        self.keep_src_3d = obj.get('keep_src_3d', False)
        # keep_hs: Keep hydrogens during processing
        self.keep_hs = obj.get('keep_hs', False)

def RunInference(args):
    os.makedirs(args.out_dir, exist_ok=True)
    with open(f'{args.model_dir}/model_parameters.yml') as f:
        score_model_args = Namespace(**yaml.full_load(f))
    if args.confidence_model_dir is not None:
        with open(f'{args.confidence_model_dir}/model_parameters.yml') as f:
            confidence_args = Namespace(**yaml.full_load(f))

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    if args.protein_ligand_csv is not None:
        df = pd.read_csv(args.protein_ligand_csv)
        complex_name_list = set_nones(df['complex_name'].tolist())
        protein_path_list = set_nones(df['protein_path'].tolist())
        protein_sequence_list = set_nones(df['protein_sequence'].tolist())
        ligand_description_list = set_nones(df['ligand_description'].tolist())
    else:
        complex_name_list = [args.complex_name]
        protein_path_list = [args.protein_path]
        protein_sequence_list = [args.protein_sequence]
        ligand_description_list = [args.ligand_description]

    complex_name_list = [name if name is not None else f"complex_{i}" for i, name in enumerate(complex_name_list)]
    for name in complex_name_list:
        write_dir = f'{args.out_dir}/{name}'
        os.makedirs(write_dir, exist_ok=True)

    # CONIFER POINT
    # Defaults
    import json
    add_hs = args.add_hs
    keep_src_3d = args.keep_src_3d
    # Annoying params for removing Hs: defaults, general override via keep_hs, individual overrides
    # Remove H defaults
    remove_hs_score = score_model_args.remove_hs
    remove_hs_confidence = confidence_args.remove_hs if args.confidence_model_dir is not None else None
    remove_hs_output = score_model_args.remove_hs
        
    # General override
    if args.keep_hs:
        remove_hs_score = False
        remove_hs_confidence = False
        remove_hs_output = False

    print(f"Special vars: add_hs={json.dumps(add_hs)} keep_src_3d={json.dumps(keep_src_3d)} remove_hs_score={json.dumps(remove_hs_score)} remove_hs_confidence={json.dumps(remove_hs_confidence)} remove_hs_output={json.dumps(remove_hs_output)}")
    # DONE CONIFER POINT

    # preprocessing of complexes into geometric graphs
    test_dataset = InferenceDataset(out_dir=args.out_dir, complex_names=complex_name_list, protein_files=protein_path_list,
                                ligand_descriptions=ligand_description_list, protein_sequences=protein_sequence_list,
                                lm_embeddings=score_model_args.esm_embeddings_path is not None,
                                receptor_radius=score_model_args.receptor_radius, remove_hs=remove_hs_score,
                                c_alpha_max_neighbors=score_model_args.c_alpha_max_neighbors,
                                all_atoms=score_model_args.all_atoms, atom_radius=score_model_args.atom_radius,
                                atom_max_neighbors=score_model_args.atom_max_neighbors,
                                keep_src_3d=keep_src_3d,
                                )
    test_loader = DataLoader(dataset=test_dataset, batch_size=1, shuffle=False)

    if args.confidence_model_dir is not None and not confidence_args.use_original_model_cache:
        print('HAPPENING | confidence model uses different type of graphs than the score model. '
              'Loading (or creating if not existing) the data for the confidence model now.')
        confidence_test_dataset = \
            InferenceDataset(out_dir=args.out_dir, complex_names=complex_name_list, protein_files=protein_path_list,
                         ligand_descriptions=ligand_description_list, protein_sequences=protein_sequence_list,
                         lm_embeddings=confidence_args.esm_embeddings_path is not None,
                         receptor_radius=confidence_args.receptor_radius, remove_hs=remove_hs_confidence,
                         c_alpha_max_neighbors=confidence_args.c_alpha_max_neighbors,
                         all_atoms=confidence_args.all_atoms, atom_radius=confidence_args.atom_radius,
                         atom_max_neighbors=confidence_args.atom_max_neighbors,
                         precomputed_lm_embeddings=test_dataset.lm_embeddings,
                         keep_src_3d=keep_src_3d,
                         )
    else:
        confidence_test_dataset = None

    t_to_sigma = partial(t_to_sigma_compl, args=score_model_args)

    model = get_model(score_model_args, device, t_to_sigma=t_to_sigma, no_parallel=True)
    state_dict = torch.load(f'{args.model_dir}/{args.ckpt}', map_location=torch.device('cpu'))
    model.load_state_dict(state_dict, strict=True)
    model = model.to(device)
    model.eval()

    if args.confidence_model_dir is not None:
        confidence_model = get_model(confidence_args, device, t_to_sigma=t_to_sigma, no_parallel=True, confidence_mode=True)
        state_dict = torch.load(f'{args.confidence_model_dir}/{args.confidence_ckpt}', map_location=torch.device('cpu'))
        confidence_model.load_state_dict(state_dict, strict=True)
        confidence_model = confidence_model.to(device)
        confidence_model.eval()
    else:
        confidence_model = None
        confidence_args = None

    tr_schedule = get_t_schedule(inference_steps=args.inference_steps)

    failures, skipped = 0, 0
    N = args.samples_per_complex
    print('Size of test dataset: ', len(test_dataset))
    for idx, orig_complex_graph in tqdm(enumerate(test_loader)):
        if not orig_complex_graph.success[0]:
            skipped += 1
            print(f"HAPPENING | The test dataset did not contain {test_dataset.complex_names[idx]} for {test_dataset.ligand_descriptions[idx]} and {test_dataset.protein_files[idx]}. We are skipping this complex.")
            continue
        try:
            if confidence_test_dataset is not None:
                confidence_complex_graph = confidence_test_dataset[idx]
                if not confidence_complex_graph.success:
                    skipped += 1
                    print(f"HAPPENING | The confidence dataset did not contain {orig_complex_graph.name}. We are skipping this complex.")
                    continue
                confidence_data_list = [copy.deepcopy(confidence_complex_graph) for _ in range(N)]
            else:
                confidence_data_list = None
            data_list = [copy.deepcopy(orig_complex_graph) for _ in range(N)]
            randomize_position(data_list, score_model_args.no_torsion, False, score_model_args.tr_sigma_max)
            lig = orig_complex_graph.mol[0]

            # initialize visualisation
            pdb = None
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

            # run reverse diffusion
            data_list, confidence = sampling(data_list=data_list, model=model,
                                         inference_steps=args.actual_steps if args.actual_steps is not None else args.inference_steps,
                                         tr_schedule=tr_schedule, rot_schedule=tr_schedule, tor_schedule=tr_schedule,
                                         device=device, t_to_sigma=t_to_sigma, model_args=score_model_args,
                                         visualization_list=visualization_list, confidence_model=confidence_model,
                                         confidence_data_list=confidence_data_list, confidence_model_args=confidence_args,
                                         batch_size=args.batch_size, no_final_step_noise=args.no_final_step_noise)
            ligand_pos = np.asarray([complex_graph['ligand'].pos.cpu().numpy() + orig_complex_graph.original_center.cpu().numpy() for complex_graph in data_list])

            # reorder predictions based on confidence output
            if confidence is not None and isinstance(confidence_args.rmsd_classification_cutoff, list):
                confidence = confidence[:, 0]
            if confidence is not None:
                confidence = confidence.cpu().numpy()
                re_order = np.argsort(confidence)[::-1]
                confidence = confidence[re_order]
                ligand_pos = ligand_pos[re_order]

            # save predictions
            write_dir = f'{args.out_dir}/{complex_name_list[idx]}'
            for rank, pos in enumerate(ligand_pos):
                mol_pred = copy.deepcopy(lig)
                if remove_hs_output: mol_pred = RemoveHs(mol_pred)
                if rank == 0: write_mol_with_coords(mol_pred, pos, os.path.join(write_dir, f'rank{rank+1}.sdf'), add_hs=add_hs)
                write_mol_with_coords(mol_pred, pos, os.path.join(write_dir, f'rank{rank+1}_confidence{confidence[rank]:.2f}.sdf'), add_hs=add_hs)

            # save visualisation frames
            if args.save_visualisation:
                if confidence is not None:
                    for rank, batch_idx in enumerate(re_order):
                        visualization_list[batch_idx].write(os.path.join(write_dir, f'rank{rank+1}_reverseprocess.pdb'))
                else:
                    for rank, batch_idx in enumerate(ligand_pos):
                        visualization_list[batch_idx].write(os.path.join(write_dir, f'rank{rank+1}_reverseprocess.pdb'))

        except Exception as e:
            print("Failed on", orig_complex_graph["name"], e)
            failures += 1

        print(f'Failed for {failures} complexes')
        print(f'Skipped {skipped} complexes')
        print(f'Results are in {args.out_dir}')


