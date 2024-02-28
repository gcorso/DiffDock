import os
import subprocess

import numpy as np
from rdkit.Chem import AllChem, RemoveHs, RemoveAllHs

from datasets.process_mols import write_mol_with_coords, read_molecule
import re

from utils.utils import remove_all_hs


def read_gnina_metrics(gnina_sdf_path):
    with open(gnina_sdf_path, 'r') as f:
        pattern = re.compile(r'> <(.*?)>\n(.*?)\n')
        content = f.read()
        matches = pattern.findall(content)
        metrics = {k: v for k, v in matches}
    return metrics


def read_gnina_score(gnina_sdf_path):
    with open(gnina_sdf_path, 'r') as f:
        pattern = re.compile(r'> <CNNscore>\n(.*?)\n')
        content = f.read()
        matches = pattern.findall(content)
    return float(matches[0])


def invert_permutation(p):
    """Return an array s with which np.array_equal(arr[p][s], arr) is True.
    The array_like argument p must be some permutation of 0, 1, ..., len(p)-1.
    """
    p = np.asanyarray(p) # in case p is a tuple, etc.
    s = np.empty_like(p)
    s[p] = np.arange(p.size)
    return s


def get_gnina_poses(args, mol, pos, orig_center, name, folder, gnina_path, thread_id=0):
    #folder = "data/MOAD_new_test_processed" if args.split == 'test' else "data/MOAD_new_val_processed"
    out_dir = args.out_dir if hasattr(args, 'out_dir') else args.inference_out_dir
    rec_path = os.path.join(folder, name[:6] + '_protein_chain_removed.pdb')
    pred_lig_path = os.path.join(out_dir, f'pred_{name}_tid{thread_id}_lig.sdf')

    if not os.path.exists(os.path.dirname(pred_lig_path)):
        os.mkdir(os.path.dirname(pred_lig_path))
    
    print(f'Ligand path {pred_lig_path}')
    write_mol_with_coords(mol, pos + orig_center, pred_lig_path)
    gnina_pred_path = os.path.join(out_dir, f'gnina_{name}_tid{thread_id}_lig.sdf')
    
    gnina_logs_dir = os.path.join(out_dir, "gnina_logs")
    
    with open(os.path.join(gnina_logs_dir, f'{name}'), "w+") as f:
        if args.gnina_full_dock:
            return_code = subprocess.run(
                f'{gnina_path} -r {rec_path} -l "{pred_lig_path}" --autobox_ligand "{pred_lig_path}" -o "{gnina_pred_path}" --no_gpu --autobox_add {args.gnina_autobox_add}',
                shell=True, stdout=f, stderr=f)
        else:
            return_code = subprocess.run(
                f'{gnina_path} --receptor {rec_path} --ligand "{pred_lig_path}" --minimize -o "{gnina_pred_path}"',
                shell=True, stdout=f, stderr=f)

    # print(f'gnina return code: {return_code}')
    try:
        gnina_mol = RemoveAllHs(read_molecule(gnina_pred_path, remove_hs=True, sanitize=True))
        gnina_minimized_ligand_pos = np.array(gnina_mol.GetConformer(0).GetPositions())
        gnina_atoms = np.array([atom.GetSymbol() for atom in gnina_mol.GetAtoms()])
        gnina_filter_Hs = np.where(gnina_atoms != 'H')
        gnina_ligand_pos = gnina_minimized_ligand_pos[gnina_filter_Hs] - orig_center

        try:
            gnina_score = read_gnina_score(gnina_pred_path)
            if gnina_score is None:
                gnina_score = 0
        except Exception as e:
            print(f'Error reading gnina score: {e}')
            gnina_score = 0

    except Exception as e:
        print(f'Error when running gnina with {name} to minimize energy')
        print('Error:', e)
        print('Using score model output pos instead.')
        gnina_ligand_pos = pos
        gnina_mol = RemoveAllHs(mol)
        gnina_score = 0

    return gnina_ligand_pos, gnina_mol, gnina_score
