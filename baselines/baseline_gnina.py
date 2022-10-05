# small script to extract the ligand and save it in a separate file because GNINA will use the ligand position as
# initial pose
import os
import shutil
import subprocess
import sys

import time
from argparse import ArgumentParser, FileType
from datetime import datetime

import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb
from rdkit import Chem
from rdkit.Chem import AllChem, MolToPDBFile
from scipy.spatial.distance import cdist

from datasets.pdbbind import read_mol
from utils.utils import read_strings_from_txt

parser = ArgumentParser()
parser.add_argument('--data_dir', type=str, default='data/PDBBind_processed', help='')
parser.add_argument('--file_suffix', type=str, default='_baseline_ligand', help='Path to folder with trained model and hyperparameters')
parser.add_argument('--results_path', type=str, default='results/gnina_predictions', help='')
parser.add_argument('--complex_names_path', type=str, default='data/splits/timesplit_test', help='')
parser.add_argument('--seed_molecules_path', type=str, default=None, help='Use the molecules at seed molecule path as initialization and only search around them')
parser.add_argument('--seed_molecule_filename', type=str, default='equibind_corrected.sdf', help='Use the molecules at seed molecule path as initialization and only search around them')
parser.add_argument('--smina', action='store_true', default=False, help='')
parser.add_argument('--no_gpu', action='store_true', default=False, help='')
parser.add_argument('--exhaustiveness', type=int, default=8, help='')
parser.add_argument('--num_cpu', type=int, default=16, help='')
parser.add_argument('--pocket_mode', action='store_true', default=False, help='')
parser.add_argument('--pocket_cutoff', type=int, default=5, help='')
parser.add_argument('--num_modes', type=int, default=10, help='')
parser.add_argument('--autobox_add', type=int, default=4, help='')
parser.add_argument('--use_p2rank_pocket', action='store_true', default=False, help='')
parser.add_argument('--skip_p2rank', action='store_true', default=False, help='')
parser.add_argument('--prank_path', type=str, default='/Users/hstark/projects/p2rank_2.3/prank', help='')
parser.add_argument('--skip_existing', action='store_true', default=False, help='')





args = parser.parse_args()

class Logger(object):
    def __init__(self, logpath, syspart=sys.stdout):
        self.terminal = syspart
        self.log = open(logpath, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass


def log(*args):
    print(f'[{datetime.now()}]', *args)


# parameters
names = read_strings_from_txt(args.complex_names_path)

if os.path.exists(args.results_path) and not args.skip_existing:
    shutil.rmtree(args.results_path)
os.makedirs(args.results_path, exist_ok=True)
sys.stdout = Logger(logpath=f'{args.results_path}/gnina.log', syspart=sys.stdout)
sys.stderr = Logger(logpath=f'{args.results_path}/error.log', syspart=sys.stderr)

p2rank_cache_path = "results/.p2rank_cache"
if args.use_p2rank_pocket and not args.skip_p2rank:
    os.makedirs(p2rank_cache_path, exist_ok=True)
    pdb_files_cache = os.path.join(p2rank_cache_path,'pdb_files')
    os.makedirs(pdb_files_cache, exist_ok=True)
    with open(f"{p2rank_cache_path}/pdb_list_p2rank.txt", "w") as out:
        for name in names:
            shutil.copy(os.path.join(args.data_dir, name, f'{name}_protein_processed.pdb'), f'{pdb_files_cache}/{name}_protein_processed.pdb')
            out.write(os.path.join('pdb_files', f'{name}_protein_processed.pdb\n'))
    cmd = f"bash {args.prank_path} predict {p2rank_cache_path}/pdb_list_p2rank.txt -o {p2rank_cache_path}/p2rank_output -threads 4"
    os.system(cmd)


all_times = []
start_time = time.time()
for i, name in enumerate(names):
    os.makedirs(os.path.join(args.results_path, name), exist_ok=True)
    log('\n')
    log(f'complex {i} of {len(names)}')
    # call gnina to find binding pose
    rec_path = os.path.join(args.data_dir, name, f'{name}_protein_processed.pdb')
    prediction_output_name = os.path.join(args.results_path, name, f'{name}{args.file_suffix}.pdb')
    log_path = os.path.join(args.results_path, name, f'{name}{args.file_suffix}.log')
    if args.seed_molecules_path is not None: seed_mol_path = os.path.join(args.seed_molecules_path, name, f'{args.seed_molecule_filename}')
    if args.skip_existing and os.path.exists(prediction_output_name): continue

    if args.pocket_mode:
        mol = read_mol(args.data_dir, name, remove_hs=False)
        rec = PandasPdb().read_pdb(rec_path)
        rec_df = rec.get(s='c-alpha')
        rec_pos = rec_df[['x_coord', 'y_coord', 'z_coord']].to_numpy().squeeze().astype(np.float32)
        lig_pos = mol.GetConformer().GetPositions()
        d = cdist(rec_pos, lig_pos)
        label = np.any(d < args.pocket_cutoff, axis=1)

        if np.any(label):
            center_pocket = rec_pos[label].mean(axis=0)
        else:
            print("No pocket residue below minimum distance ", args.pocket_cutoff, "taking closest at", np.min(d))
            center_pocket = rec_pos[np.argmin(np.min(d, axis=1)[0])]
        radius_pocket = np.max(np.linalg.norm(lig_pos - center_pocket[None, :], axis=1))
        diameter_pocket = radius_pocket * 2
        center_x = center_pocket[0]
        size_x = diameter_pocket + 8
        center_y = center_pocket[1]
        size_y = diameter_pocket + 8
        center_z = center_pocket[2]
        size_z = diameter_pocket + 8


    mol_rdkit = read_mol(args.data_dir, name, remove_hs=False)
    single_time = time.time()

    mol_rdkit.RemoveAllConformers()
    ps = AllChem.ETKDGv2()
    id = AllChem.EmbedMolecule(mol_rdkit, ps)
    if id == -1:
        print('rdkit pos could not be generated without using random pos. using random pos now.')
        ps.useRandomCoords = True
        AllChem.EmbedMolecule(mol_rdkit, ps)
        AllChem.MMFFOptimizeMolecule(mol_rdkit, confId=0)
    rdkit_mol_path = os.path.join(args.data_dir, name, f'{name}_rdkit_ligand.pdb')
    MolToPDBFile(mol_rdkit, rdkit_mol_path)

    fallback_without_p2rank = False
    if args.use_p2rank_pocket:
        df = pd.read_csv(f'{p2rank_cache_path}/p2rank_output/{name}_protein_processed.pdb_predictions.csv')
        rdkit_lig_pos = mol_rdkit.GetConformer().GetPositions()
        diameter_pocket = np.max(cdist(rdkit_lig_pos, rdkit_lig_pos))
        size_x = diameter_pocket + args.autobox_add * 2
        size_y = diameter_pocket + args.autobox_add * 2
        size_z = diameter_pocket + args.autobox_add * 2
        if df.empty:
            fallback_without_p2rank = True
        else:
            center_x = df.iloc[0]['   center_x']
            center_y = df.iloc[0]['   center_y']
            center_z = df.iloc[0]['   center_z']



    log(f'processing {rec_path}')
    if not args.pocket_mode and not args.use_p2rank_pocket or fallback_without_p2rank:
        return_code = subprocess.run(
            f"gnina --receptor {rec_path} --ligand {rdkit_mol_path} --num_modes {args.num_modes} -o {prediction_output_name} {'--no_gpu' if args.no_gpu else ''} --autobox_ligand {rec_path if args.seed_molecules_path is None else seed_mol_path} --autobox_add {args.autobox_add} --log {log_path} --exhaustiveness {args.exhaustiveness} --cpu {args.num_cpu} {'--cnn_scoring none' if args.smina else ''}",
            shell=True)
    else:
        return_code = subprocess.run(
            f"gnina --receptor {rec_path} --ligand {rdkit_mol_path} --num_modes {args.num_modes} -o {prediction_output_name} {'--no_gpu' if args.no_gpu else ''} --log {log_path} --exhaustiveness {args.exhaustiveness} --cpu {args.num_cpu} {'--cnn_scoring none' if args.smina else ''} --center_x {center_x} --center_y {center_y} --center_z {center_z} --size_x {size_x} --size_y {size_y} --size_z {size_z}",
            shell=True)
    log(return_code)
    all_times.append(time.time() - single_time)

    log("single time: --- %s seconds ---" % (time.time() - single_time))
    log("time so far: --- %s seconds ---" % (time.time() - start_time))
    log('\n')
log(all_times)
log("--- %s seconds ---" % (time.time() - start_time))
