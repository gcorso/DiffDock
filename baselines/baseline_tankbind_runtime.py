# This file needs to be ran in the TANKBind repository together with baseline_run_tankbind_parallel.sh

import sys
import time
from multiprocessing import Pool


import copy
import warnings
from argparse import ArgumentParser

from rdkit.Chem import AllChem, RemoveHs

from feature_utils import save_cleaned_protein, read_mol
from generation_utils import get_LAS_distance_constraint_mask, get_info_pred_distance, write_with_new_coords
import logging
from torch_geometric.loader import DataLoader
from tqdm import tqdm  # pip install tqdm if fails.
from model import get_model
# from utils import *
import torch


from data import TankBind_prediction

import os
import numpy as np
import pandas as pd
import rdkit.Chem as Chem
from feature_utils import generate_sdf_from_smiles_using_rdkit
from feature_utils import get_protein_feature
from Bio.PDB import PDBParser
from feature_utils import extract_torchdrug_feature_from_mol


def read_strings_from_txt(path):
    # every line will be one element of the returned list
    with open(path) as file:
        lines = file.readlines()
        return [line.rstrip() for line in lines]


def read_molecule(molecule_file, sanitize=False, calc_charges=False, remove_hs=False):
    if molecule_file.endswith('.mol2'):
        mol = Chem.MolFromMol2File(molecule_file, sanitize=False, removeHs=False)
    elif molecule_file.endswith('.sdf'):
        supplier = Chem.SDMolSupplier(molecule_file, sanitize=False, removeHs=False)
        mol = supplier[0]
    elif molecule_file.endswith('.pdbqt'):
        with open(molecule_file) as file:
            pdbqt_data = file.readlines()
        pdb_block = ''
        for line in pdbqt_data:
            pdb_block += '{}\n'.format(line[:66])
        mol = Chem.MolFromPDBBlock(pdb_block, sanitize=False, removeHs=False)
    elif molecule_file.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(molecule_file, sanitize=False, removeHs=False)
    else:
        return ValueError('Expect the format of the molecule_file to be '
                          'one of .mol2, .sdf, .pdbqt and .pdb, got {}'.format(molecule_file))
    try:
        if sanitize or calc_charges:
            Chem.SanitizeMol(mol)

        if calc_charges:
            # Compute Gasteiger charges on the molecule.
            try:
                AllChem.ComputeGasteigerCharges(mol)
            except:
                warnings.warn('Unable to compute charges for the molecule.')

        if remove_hs:
            mol = Chem.RemoveHs(mol, sanitize=sanitize)
    except:
        return None

    return mol


def parallel_save_prediction(arguments):
    dataset, y_pred_list, chosen,rdkit_mol_path, result_folder, name = arguments
    for idx, line in chosen.iterrows():
        pocket_name = line['pocket_name']
        compound_name = line['compound_name']
        ligandName = compound_name.split("_")[1]
        dataset_index = line['dataset_index']
        coords = dataset[dataset_index].coords.to('cpu')
        protein_nodes_xyz = dataset[dataset_index].node_xyz.to('cpu')
        n_compound = coords.shape[0]
        n_protein = protein_nodes_xyz.shape[0]
        y_pred = y_pred_list[dataset_index].reshape(n_protein, n_compound).to('cpu')
        compound_pair_dis_constraint = torch.cdist(coords, coords)
        mol = Chem.MolFromMolFile(rdkit_mol_path)
        LAS_distance_constraint_mask = get_LAS_distance_constraint_mask(mol).bool()
        pred_dist_info = get_info_pred_distance(coords, y_pred, protein_nodes_xyz, compound_pair_dis_constraint,
                                                LAS_distance_constraint_mask=LAS_distance_constraint_mask,
                                                n_repeat=1, show_progress=False)

        toFile = f'{result_folder}/{name}_tankbind_chosen.sdf'
        new_coords = pred_dist_info.sort_values("loss")['coords'].iloc[0].astype(np.double)
        write_with_new_coords(mol, new_coords, toFile)

if __name__ == '__main__':
    tankbind_src_folder = "../tankbind"
    sys.path.insert(0, tankbind_src_folder)
    torch.set_num_threads(16)
    parser = ArgumentParser()
    parser.add_argument('--data_dir', type=str, default='/Users/hstark/projects/ligbind/data/PDBBind_processed', help='')
    parser.add_argument('--split_path', type=str, default='/Users/hstark/projects/ligbind/data/splits/timesplit_test', help='')
    parser.add_argument('--prank_path', type=str, default='/Users/hstark/projects/p2rank_2.3/prank', help='')
    parser.add_argument('--results_path', type=str, default='results/tankbind_results', help='')
    parser.add_argument('--skip_existing', action='store_true', default=False, help='')
    parser.add_argument('--skip_p2rank', action='store_true', default=False, help='')
    parser.add_argument('--skip_multiple_pocket_outputs', action='store_true', default=False, help='')
    parser.add_argument('--device', type=str, default='cpu', help='')
    parser.add_argument('--num_workers', type=int, default=1, help='')
    parser.add_argument('--parallel_id', type=int, default=0, help='')
    parser.add_argument('--parallel_tot', type=int, default=1, help='')
    args = parser.parse_args()

    device = args.device
    cache_path = "tankbind_cache"
    os.makedirs(cache_path, exist_ok=True)
    os.makedirs(args.results_path, exist_ok=True)



    logging.basicConfig(level=logging.INFO)
    model = get_model(0, logging, device)
    # re-dock model
    # modelFile = "../saved_models/re_dock.pt"
    # self-dock model
    modelFile = f"{tankbind_src_folder}/../saved_models/self_dock.pt"

    model.load_state_dict(torch.load(modelFile, map_location=device))
    _ = model.eval()
    batch_size = 5
    names = read_strings_from_txt(args.split_path)
    if args.parallel_tot > 1:
        size = len(names) // args.parallel_tot + 1
        names = names[args.parallel_id*size:(args.parallel_id+1)*size]
    rmsds = []

    forward_pass_time = []
    times_preprocess = []
    times_inference = []
    top_10_generation_time = []
    top_1_generation_time = []
    start_time = time.time()
    if not args.skip_p2rank:
        for name in names:
            if args.skip_existing and os.path.exists(f'{args.results_path}/{name}/{name}_tankbind_1.sdf'): continue
            print("Now processing: ", name)
            protein_path = f'{args.data_dir}/{name}/{name}_protein_processed.pdb'
            cleaned_protein_path = f"{cache_path}/{name}_protein_tankbind_cleaned.pdb"  # if you change this you also need to change below
            parser = PDBParser(QUIET=True)
            s = parser.get_structure(name, protein_path)
            c = s[0]
            clean_res_list, ligand_list = save_cleaned_protein(c, cleaned_protein_path)

        with open(f"{cache_path}/pdb_list_p2rank.txt", "w") as out:
            for name in names:
                out.write(f"{name}_protein_tankbind_cleaned.pdb\n")
        cmd = f"bash {args.prank_path} predict {cache_path}/pdb_list_p2rank.txt -o {cache_path}/p2rank -threads 4"
        os.system(cmd)
    times_preprocess.append(time.time() - start_time)
    p2_rank_time = time.time() - start_time




    list_to_parallelize = []
    for name in tqdm(names):
        single_preprocess_time = time.time()
        if args.skip_existing and os.path.exists(f'{args.results_path}/{name}/{name}_tankbind_1.sdf'): continue
        print("Now processing: ", name)
        protein_path = f'{args.data_dir}/{name}/{name}_protein_processed.pdb'
        ligand_path = f"{args.data_dir}/{name}/{name}_ligand.sdf"
        cleaned_protein_path = f"{cache_path}/{name}_protein_tankbind_cleaned.pdb"  # if you change this you also need to change below
        rdkit_mol_path = f"{cache_path}/{name}_rdkit_ligand.sdf"

        parser = PDBParser(QUIET=True)
        s = parser.get_structure(name, protein_path)
        c = s[0]
        clean_res_list, ligand_list = save_cleaned_protein(c, cleaned_protein_path)
        lig, _ = read_mol(f"{args.data_dir}/{name}/{name}_ligand.sdf", f"{args.data_dir}/{name}/{name}_ligand.mol2")

        lig = RemoveHs(lig)
        smiles = Chem.MolToSmiles(lig)
        generate_sdf_from_smiles_using_rdkit(smiles, rdkit_mol_path, shift_dis=0)

        parser = PDBParser(QUIET=True)
        s = parser.get_structure("x", cleaned_protein_path)
        res_list = list(s.get_residues())

        protein_dict = {}
        protein_dict[name] = get_protein_feature(res_list)
        compound_dict = {}

        mol = Chem.MolFromMolFile(rdkit_mol_path)
        compound_dict[name + f"_{name}" + "_rdkit"] = extract_torchdrug_feature_from_mol(mol, has_LAS_mask=True)

        info = []
        for compound_name in list(compound_dict.keys()):
            # use protein center as the block center.
            com = ",".join([str(a.round(3)) for a in protein_dict[name][0].mean(axis=0).numpy()])
            info.append([name, compound_name, "protein_center", com])

            p2rankFile = f"{cache_path}/p2rank/{name}_protein_tankbind_cleaned.pdb_predictions.csv"
            pocket = pd.read_csv(p2rankFile)
            pocket.columns = pocket.columns.str.strip()
            pocket_coms = pocket[['center_x', 'center_y', 'center_z']].values
            for ith_pocket, com in enumerate(pocket_coms):
                com = ",".join([str(a.round(3)) for a in com])
                info.append([name, compound_name, f"pocket_{ith_pocket + 1}", com])
        info = pd.DataFrame(info, columns=['protein_name', 'compound_name', 'pocket_name', 'pocket_com'])

        dataset_path = f"{cache_path}/{name}_dataset/"
        os.system(f"rm -r {dataset_path}")
        os.system(f"mkdir -p {dataset_path}")
        dataset = TankBind_prediction(dataset_path, data=info, protein_dict=protein_dict, compound_dict=compound_dict)

        # dataset = TankBind_prediction(dataset_path)
        times_preprocess.append(time.time() - single_preprocess_time)
        single_forward_pass_time = time.time()
        data_loader = DataLoader(dataset, batch_size=batch_size, follow_batch=['x', 'y', 'compound_pair'], shuffle=False,
                                 num_workers=0)
        affinity_pred_list = []
        y_pred_list = []
        for data in tqdm(data_loader):
            data = data.to(device)
            y_pred, affinity_pred = model(data)
            affinity_pred_list.append(affinity_pred.detach().cpu())
            for i in range(data.y_batch.max() + 1):
                y_pred_list.append((y_pred[data['y_batch'] == i]).detach().cpu())

        affinity_pred_list = torch.cat(affinity_pred_list)
        forward_pass_time.append(time.time() - single_forward_pass_time)
        output_info = copy.deepcopy(dataset.data)
        output_info['affinity'] = affinity_pred_list
        output_info['dataset_index'] = range(len(output_info))
        output_info_sorted = output_info.sort_values('affinity', ascending=False)


        result_folder = f'{args.results_path}/{name}'
        os.makedirs(result_folder, exist_ok=True)
        output_info_sorted.to_csv(f"{result_folder}/output_info_sorted_by_affinity.csv")

        if not args.skip_multiple_pocket_outputs:
            for idx, (dataframe_idx, line) in enumerate(copy.deepcopy(output_info_sorted).iterrows()):
                single_top10_generation_time = time.time()
                pocket_name = line['pocket_name']
                compound_name = line['compound_name']
                ligandName = compound_name.split("_")[1]
                coords = dataset[dataframe_idx].coords.to('cpu')
                protein_nodes_xyz = dataset[dataframe_idx].node_xyz.to('cpu')
                n_compound = coords.shape[0]
                n_protein = protein_nodes_xyz.shape[0]
                y_pred = y_pred_list[dataframe_idx].reshape(n_protein, n_compound).to('cpu')
                y = dataset[dataframe_idx].dis_map.reshape(n_protein, n_compound).to('cpu')
                compound_pair_dis_constraint = torch.cdist(coords, coords)
                mol = Chem.MolFromMolFile(rdkit_mol_path)
                LAS_distance_constraint_mask = get_LAS_distance_constraint_mask(mol).bool()
                pred_dist_info = get_info_pred_distance(coords, y_pred, protein_nodes_xyz, compound_pair_dis_constraint,
                                              LAS_distance_constraint_mask=LAS_distance_constraint_mask,
                                              n_repeat=1, show_progress=False)

                toFile = f'{result_folder}/{name}_tankbind_{idx}.sdf'
                new_coords = pred_dist_info.sort_values("loss")['coords'].iloc[0].astype(np.double)
                write_with_new_coords(mol, new_coords, toFile)
                if idx < 10:
                    top_10_generation_time.append(time.time() - single_top10_generation_time)
                if idx == 0:
                    top_1_generation_time.append(time.time() - single_top10_generation_time)

        output_info_chosen = copy.deepcopy(dataset.data)
        output_info_chosen['affinity'] = affinity_pred_list
        output_info_chosen['dataset_index'] = range(len(output_info_chosen))
        chosen = output_info_chosen.loc[
            output_info_chosen.groupby(['protein_name', 'compound_name'], sort=False)['affinity'].agg(
                'idxmax')].reset_index()

        list_to_parallelize.append((dataset, y_pred_list, chosen, rdkit_mol_path, result_folder, name))

    chosen_generation_start_time = time.time()
    if args.num_workers > 1:
        p = Pool(args.num_workers, maxtasksperchild=1)
        p.__enter__()
    with tqdm(total=len(list_to_parallelize), desc=f'running optimization {i}/{len(list_to_parallelize)}') as pbar:
        map_fn = p.imap_unordered if args.num_workers > 1 else map
        for t in map_fn(parallel_save_prediction, list_to_parallelize):
            pbar.update()
    if args.num_workers > 1: p.__exit__(None, None, None)
    chosen_generation_time = time.time() - chosen_generation_start_time
    """
        lig, _ = read_mol(f"{args.data_dir}/{name}/{name}_ligand.sdf", f"{args.data_dir}/{name}/{name}_ligand.mol2")
        sm = Chem.MolToSmiles(lig)
        m_order = list(lig.GetPropsAsDict(includePrivate=True, includeComputed=True)['_smilesAtomOutputOrder'])
        lig = Chem.RenumberAtoms(lig, m_order)
        lig = Chem.RemoveAllHs(lig)
        lig = RemoveHs(lig)
        true_ligand_pos = np.array(lig.GetConformer().GetPositions())
    
        toFile = f'{result_folder}/{name}_tankbind_chosen.sdf'
        mol_pred, _ = read_mol(toFile, None)
        sm = Chem.MolToSmiles(mol_pred)
        m_order = list(mol_pred.GetPropsAsDict(includePrivate=True, includeComputed=True)['_smilesAtomOutputOrder'])
        mol_pred = Chem.RenumberAtoms(mol_pred, m_order)
        mol_pred = RemoveHs(mol_pred)
        mol_pred_pos = np.array(mol_pred.GetConformer().GetPositions())
        rmsds.append(np.sqrt(((true_ligand_pos - mol_pred_pos) ** 2).sum(axis=1).mean(axis=0)))
        print(np.sqrt(((true_ligand_pos - mol_pred_pos) ** 2).sum(axis=1).mean(axis=0)))
    """
    forward_pass_time  = np.array(forward_pass_time).sum()
    times_preprocess  = np.array(times_preprocess).sum()
    times_inference  = np.array(times_inference).sum()
    top_10_generation_time  = np.array(top_10_generation_time).sum()
    top_1_generation_time  = np.array(top_1_generation_time).sum()

    rmsds = np.array(rmsds)

    print(f'forward_pass_time: {forward_pass_time}')
    print(f'times_preprocess: {times_preprocess}')
    print(f'times_inference: {times_inference}')
    print(f'top_10_generation_time: {top_10_generation_time}')
    print(f'top_1_generation_time: {top_1_generation_time}')
    print(f'chosen_generation_time: {chosen_generation_time}')
    print(f'rmsds_below_2: {(100 * (rmsds < 2).sum() / len(rmsds))}')
    print(f'p2rank Time: {p2_rank_time}')
    print(
        f'total_time: '
        f'{forward_pass_time + times_preprocess + times_inference + top_10_generation_time + top_1_generation_time + p2_rank_time}')

    with open(os.path.join(args.results_path, 'tankbind_log.log'), 'w') as file:
        file.write(f'forward_pass_time: {forward_pass_time}')
        file.write(f'times_preprocess: {times_preprocess}')
        file.write(f'times_inference: {times_inference}')
        file.write(f'top_10_generation_time: {top_10_generation_time}')
        file.write(f'top_1_generation_time: {top_1_generation_time}')
        file.write(f'rmsds_below_2: {(100 * (rmsds < 2).sum() / len(rmsds))}')
        file.write(f'p2rank Time: {p2_rank_time}')
        file.write(f'total_time: {forward_pass_time + times_preprocess + times_inference + top_10_generation_time + top_1_generation_time + p2_rank_time}')
