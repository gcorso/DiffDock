import binascii
import glob
import hashlib
import os
import pickle
from collections import defaultdict
from multiprocessing import Pool
import random
import copy

import numpy as np
import torch
from rdkit.Chem import MolToSmiles, MolFromSmiles, AddHs
from torch_geometric.data import Dataset, HeteroData
from torch_geometric.loader import DataLoader, DataListLoader
from torch_geometric.transforms import BaseTransform
from tqdm import tqdm

from datasets.process_mols import read_molecule, get_rec_graph, generate_conformer, \
    get_lig_graph_with_matching, extract_receptor_structure, parse_receptor, parse_pdb_from_path
from utils.diffusion_utils import modify_conformer, set_time
from utils.utils import read_strings_from_txt
from utils import so3, torus


class NoiseTransform(BaseTransform):
    def __init__(self, t_to_sigma, no_torsion, all_atom):
        self.t_to_sigma = t_to_sigma
        self.no_torsion = no_torsion
        self.all_atom = all_atom

    def __call__(self, data):
        t = np.random.uniform()
        t_tr, t_rot, t_tor = t, t, t
        return self.apply_noise(data, t_tr, t_rot, t_tor)

    def apply_noise(self, data, t_tr, t_rot, t_tor, tr_update = None, rot_update=None, torsion_updates=None):
        if not torch.is_tensor(data['ligand'].pos):
            data['ligand'].pos = random.choice(data['ligand'].pos)

        tr_sigma, rot_sigma, tor_sigma = self.t_to_sigma(t_tr, t_rot, t_tor)
        set_time(data, t_tr, t_rot, t_tor, 1, self.all_atom, device=None)

        tr_update = torch.normal(mean=0, std=tr_sigma, size=(1, 3)) if tr_update is None else tr_update
        rot_update = so3.sample_vec(eps=rot_sigma) if rot_update is None else rot_update
        torsion_updates = np.random.normal(loc=0.0, scale=tor_sigma, size=data['ligand'].edge_mask.sum()) if torsion_updates is None else torsion_updates
        torsion_updates = None if self.no_torsion else torsion_updates
        modify_conformer(data, tr_update, torch.from_numpy(rot_update).float(), torsion_updates)

        data.tr_score = -tr_update / tr_sigma ** 2
        data.rot_score = torch.from_numpy(so3.score_vec(vec=rot_update, eps=rot_sigma)).float().unsqueeze(0)
        data.tor_score = None if self.no_torsion else torch.from_numpy(torus.score(torsion_updates, tor_sigma)).float()
        data.tor_sigma_edge = None if self.no_torsion else np.ones(data['ligand'].edge_mask.sum()) * tor_sigma
        return data


class PDBBind(Dataset):
    def __init__(self, root, transform=None, cache_path='data/cache', split_path='data/', limit_complexes=0,
                 receptor_radius=30, num_workers=1, c_alpha_max_neighbors=None, popsize=15, maxiter=15,
                 matching=True, keep_original=False, max_lig_size=None, remove_hs=False, num_conformers=1, all_atoms=False,
                 atom_radius=5, atom_max_neighbors=None, esm_embeddings_path=None, require_ligand=False,
                 ligands_list=None, protein_path_list=None, ligand_descriptions=None, keep_local_structures=False):

        super(PDBBind, self).__init__(root, transform)
        self.pdbbind_dir = root
        self.max_lig_size = max_lig_size
        self.split_path = split_path
        self.limit_complexes = limit_complexes
        self.receptor_radius = receptor_radius
        self.num_workers = num_workers
        self.c_alpha_max_neighbors = c_alpha_max_neighbors
        self.remove_hs = remove_hs
        self.esm_embeddings_path = esm_embeddings_path
        self.require_ligand = require_ligand
        self.protein_path_list = protein_path_list
        self.ligand_descriptions = ligand_descriptions
        self.keep_local_structures = keep_local_structures
        if matching or protein_path_list is not None and ligand_descriptions is not None:
            cache_path += '_torsion'
        if all_atoms:
            cache_path += '_allatoms'
        self.full_cache_path = os.path.join(cache_path, f'limit{self.limit_complexes}'
                                                        f'_INDEX{os.path.splitext(os.path.basename(self.split_path))[0]}'
                                                        f'_maxLigSize{self.max_lig_size}_H{int(not self.remove_hs)}'
                                                        f'_recRad{self.receptor_radius}_recMax{self.c_alpha_max_neighbors}'
                                            + ('' if not all_atoms else f'_atomRad{atom_radius}_atomMax{atom_max_neighbors}')
                                            + ('' if not matching or num_conformers == 1 else f'_confs{num_conformers}')
                                            + ('' if self.esm_embeddings_path is None else f'_esmEmbeddings')
                                            + ('' if not keep_local_structures else f'_keptLocalStruct')
                                            + ('' if protein_path_list is None or ligand_descriptions is None else str(binascii.crc32(''.join(ligand_descriptions + protein_path_list).encode()))))
        self.popsize, self.maxiter = popsize, maxiter
        self.matching, self.keep_original = matching, keep_original
        self.num_conformers = num_conformers
        self.all_atoms = all_atoms
        self.atom_radius, self.atom_max_neighbors = atom_radius, atom_max_neighbors
        if not os.path.exists(os.path.join(self.full_cache_path, "heterographs.pkl"))\
                or (require_ligand and not os.path.exists(os.path.join(self.full_cache_path, "rdkit_ligands.pkl"))):
            os.makedirs(self.full_cache_path, exist_ok=True)
            if protein_path_list is None or ligand_descriptions is None:
                self.preprocessing()
            else:
                self.inference_preprocessing()

        print('loading data from memory: ', os.path.join(self.full_cache_path, "heterographs.pkl"))
        with open(os.path.join(self.full_cache_path, "heterographs.pkl"), 'rb') as f:
            self.complex_graphs = pickle.load(f)
        if require_ligand:
            with open(os.path.join(self.full_cache_path, "rdkit_ligands.pkl"), 'rb') as f:
                self.rdkit_ligands = pickle.load(f)

        print_statistics(self.complex_graphs)

    def len(self):
        return len(self.complex_graphs)

    def get(self, idx):
        if self.require_ligand:
            complex_graph = copy.deepcopy(self.complex_graphs[idx])
            complex_graph.mol = copy.deepcopy(self.rdkit_ligands[idx])
            return complex_graph
        else:
            return copy.deepcopy(self.complex_graphs[idx])

    def preprocessing(self):
        print(f'Processing complexes from [{self.split_path}] and saving it to [{self.full_cache_path}]')

        complex_names_all = read_strings_from_txt(self.split_path)
        if self.limit_complexes is not None and self.limit_complexes != 0:
            complex_names_all = complex_names_all[:self.limit_complexes]
        print(f'Loading {len(complex_names_all)} complexes.')

        if self.esm_embeddings_path is not None:
            id_to_embeddings = torch.load(self.esm_embeddings_path)
            chain_embeddings_dictlist = defaultdict(list)
            for key, embedding in id_to_embeddings.items():
                key_name = key.split('_')[0]
                if key_name in complex_names_all:
                    chain_embeddings_dictlist[key_name].append(embedding)
            lm_embeddings_chains_all = []
            for name in complex_names_all:
                lm_embeddings_chains_all.append(chain_embeddings_dictlist[name])
        else:
            lm_embeddings_chains_all = [None] * len(complex_names_all)

        if self.num_workers > 1:
            # running preprocessing in parallel on multiple workers and saving the progress every 1000 complexes
            for i in range(len(complex_names_all)//1000+1):
                if os.path.exists(os.path.join(self.full_cache_path, f"heterographs{i}.pkl")):
                    continue
                complex_names = complex_names_all[1000*i:1000*(i+1)]
                lm_embeddings_chains = lm_embeddings_chains_all[1000*i:1000*(i+1)]
                complex_graphs, rdkit_ligands = [], []
                if self.num_workers > 1:
                    p = Pool(self.num_workers, maxtasksperchild=1)
                    p.__enter__()
                with tqdm(total=len(complex_names), desc=f'loading complexes {i}/{len(complex_names_all)//1000+1}') as pbar:
                    map_fn = p.imap_unordered if self.num_workers > 1 else map
                    for t in map_fn(self.get_complex, zip(complex_names, lm_embeddings_chains, [None] * len(complex_names), [None] * len(complex_names))):
                        complex_graphs.extend(t[0])
                        rdkit_ligands.extend(t[1])
                        pbar.update()
                if self.num_workers > 1: p.__exit__(None, None, None)

                with open(os.path.join(self.full_cache_path, f"heterographs{i}.pkl"), 'wb') as f:
                    pickle.dump((complex_graphs), f)
                with open(os.path.join(self.full_cache_path, f"rdkit_ligands{i}.pkl"), 'wb') as f:
                    pickle.dump((rdkit_ligands), f)

            complex_graphs_all = []
            for i in range(len(complex_names_all)//1000+1):
                with open(os.path.join(self.full_cache_path, f"heterographs{i}.pkl"), 'rb') as f:
                    l = pickle.load(f)
                    complex_graphs_all.extend(l)
            with open(os.path.join(self.full_cache_path, f"heterographs.pkl"), 'wb') as f:
                pickle.dump((complex_graphs_all), f)

            rdkit_ligands_all = []
            for i in range(len(complex_names_all) // 1000 + 1):
                with open(os.path.join(self.full_cache_path, f"rdkit_ligands{i}.pkl"), 'rb') as f:
                    l = pickle.load(f)
                    rdkit_ligands_all.extend(l)
            with open(os.path.join(self.full_cache_path, f"rdkit_ligands.pkl"), 'wb') as f:
                pickle.dump((rdkit_ligands_all), f)
        else:
            complex_graphs, rdkit_ligands = [], []
            with tqdm(total=len(complex_names_all), desc='loading complexes') as pbar:
                for t in map(self.get_complex, zip(complex_names_all, lm_embeddings_chains_all, [None] * len(complex_names_all), [None] * len(complex_names_all))):
                    complex_graphs.extend(t[0])
                    rdkit_ligands.extend(t[1])
                    pbar.update()
            with open(os.path.join(self.full_cache_path, "heterographs.pkl"), 'wb') as f:
                pickle.dump((complex_graphs), f)
            with open(os.path.join(self.full_cache_path, "rdkit_ligands.pkl"), 'wb') as f:
                pickle.dump((rdkit_ligands), f)

    def inference_preprocessing(self):
        ligands_list = []
        print('Reading molecules and generating local structures with RDKit')
        for ligand_description in tqdm(self.ligand_descriptions):
            mol = MolFromSmiles(ligand_description)  # check if it is a smiles or a path
            if mol is not None:
                mol = AddHs(mol)
                generate_conformer(mol)
                ligands_list.append(mol)
            else:
                mol = read_molecule(ligand_description, remove_hs=False, sanitize=True)
                if not self.keep_local_structures:
                    mol.RemoveAllConformers()
                    mol = AddHs(mol)
                    generate_conformer(mol)
                ligands_list.append(mol)

        if self.esm_embeddings_path is not None:
            print('Reading language model embeddings.')
            lm_embeddings_chains_all = []
            if not os.path.exists(self.esm_embeddings_path): raise Exception('ESM embeddings path does not exist: ',self.esm_embeddings_path)
            for protein_path in self.protein_path_list:
                embeddings_paths = sorted(glob.glob(os.path.join(self.esm_embeddings_path, os.path.basename(protein_path)) + '*'))
                lm_embeddings_chains = []
                for embeddings_path in embeddings_paths:
                    lm_embeddings_chains.append(torch.load(embeddings_path)['representations'][33])
                lm_embeddings_chains_all.append(lm_embeddings_chains)
        else:
            lm_embeddings_chains_all = [None] * len(self.protein_path_list)

        print('Generating graphs for ligands and proteins')
        if self.num_workers > 1:
            # running preprocessing in parallel on multiple workers and saving the progress every 1000 complexes
            for i in range(len(self.protein_path_list)//1000+1):
                if os.path.exists(os.path.join(self.full_cache_path, f"heterographs{i}.pkl")):
                    continue
                protein_paths_chunk = self.protein_path_list[1000*i:1000*(i+1)]
                ligand_description_chunk = self.ligand_descriptions[1000*i:1000*(i+1)]
                ligands_chunk = ligands_list[1000 * i:1000 * (i + 1)]
                lm_embeddings_chains = lm_embeddings_chains_all[1000*i:1000*(i+1)]
                complex_graphs, rdkit_ligands = [], []
                if self.num_workers > 1:
                    p = Pool(self.num_workers, maxtasksperchild=1)
                    p.__enter__()
                with tqdm(total=len(protein_paths_chunk), desc=f'loading complexes {i}/{len(protein_paths_chunk)//1000+1}') as pbar:
                    map_fn = p.imap_unordered if self.num_workers > 1 else map
                    for t in map_fn(self.get_complex, zip(protein_paths_chunk, lm_embeddings_chains, ligands_chunk,ligand_description_chunk)):
                        complex_graphs.extend(t[0])
                        rdkit_ligands.extend(t[1])
                        pbar.update()
                if self.num_workers > 1: p.__exit__(None, None, None)

                with open(os.path.join(self.full_cache_path, f"heterographs{i}.pkl"), 'wb') as f:
                    pickle.dump((complex_graphs), f)
                with open(os.path.join(self.full_cache_path, f"rdkit_ligands{i}.pkl"), 'wb') as f:
                    pickle.dump((rdkit_ligands), f)

            complex_graphs_all = []
            for i in range(len(self.protein_path_list)//1000+1):
                with open(os.path.join(self.full_cache_path, f"heterographs{i}.pkl"), 'rb') as f:
                    l = pickle.load(f)
                    complex_graphs_all.extend(l)
            with open(os.path.join(self.full_cache_path, f"heterographs.pkl"), 'wb') as f:
                pickle.dump((complex_graphs_all), f)

            rdkit_ligands_all = []
            for i in range(len(self.protein_path_list) // 1000 + 1):
                with open(os.path.join(self.full_cache_path, f"rdkit_ligands{i}.pkl"), 'rb') as f:
                    l = pickle.load(f)
                    rdkit_ligands_all.extend(l)
            with open(os.path.join(self.full_cache_path, f"rdkit_ligands.pkl"), 'wb') as f:
                pickle.dump((rdkit_ligands_all), f)
        else:
            complex_graphs, rdkit_ligands = [], []
            with tqdm(total=len(self.protein_path_list), desc='loading complexes') as pbar:
                for t in map(self.get_complex, zip(self.protein_path_list, lm_embeddings_chains_all, ligands_list, self.ligand_descriptions)):
                    complex_graphs.extend(t[0])
                    rdkit_ligands.extend(t[1])
                    pbar.update()
            with open(os.path.join(self.full_cache_path, "heterographs.pkl"), 'wb') as f:
                pickle.dump((complex_graphs), f)
            with open(os.path.join(self.full_cache_path, "rdkit_ligands.pkl"), 'wb') as f:
                pickle.dump((rdkit_ligands), f)

    def get_complex(self, par):
        name, lm_embedding_chains, ligand, ligand_description = par
        if not os.path.exists(os.path.join(self.pdbbind_dir, name)) and ligand is None:
            print("Folder not found", name)
            return [], []

        if ligand is not None:
            rec_model = parse_pdb_from_path(name)
            name = f'{name}____{ligand_description}'
            ligs = [ligand]
        else:
            try:
                rec_model = parse_receptor(name, self.pdbbind_dir)
            except Exception as e:
                print(f'Skipping {name} because of the error:')
                print(e)
                return [], []

            ligs = read_mols(self.pdbbind_dir, name, remove_hs=False)
        complex_graphs = []
        for i, lig in enumerate(ligs):
            if self.max_lig_size is not None and lig.GetNumHeavyAtoms() > self.max_lig_size:
                print(f'Ligand with {lig.GetNumHeavyAtoms()} heavy atoms is larger than max_lig_size {self.max_lig_size}. Not including {name} in preprocessed data.')
                continue
            complex_graph = HeteroData()
            complex_graph['name'] = name
            try:
                get_lig_graph_with_matching(lig, complex_graph, self.popsize, self.maxiter, self.matching, self.keep_original,
                                            self.num_conformers, remove_hs=self.remove_hs)
                rec, rec_coords, c_alpha_coords, n_coords, c_coords, lm_embeddings = extract_receptor_structure(copy.deepcopy(rec_model), lig, lm_embedding_chains=lm_embedding_chains)
                if lm_embeddings is not None and len(c_alpha_coords) != len(lm_embeddings):
                    print(f'LM embeddings for complex {name} did not have the right length for the protein. Skipping {name}.')
                    continue

                get_rec_graph(rec, rec_coords, c_alpha_coords, n_coords, c_coords, complex_graph, rec_radius=self.receptor_radius,
                              c_alpha_max_neighbors=self.c_alpha_max_neighbors, all_atoms=self.all_atoms,
                              atom_radius=self.atom_radius, atom_max_neighbors=self.atom_max_neighbors, remove_hs=self.remove_hs, lm_embeddings=lm_embeddings)

            except Exception as e:
                print(f'Skipping {name} because of the error:')
                print(e)
                continue

            protein_center = torch.mean(complex_graph['receptor'].pos, dim=0, keepdim=True)
            complex_graph['receptor'].pos -= protein_center
            if self.all_atoms:
                complex_graph['atom'].pos -= protein_center

            if (not self.matching) or self.num_conformers == 1:
                complex_graph['ligand'].pos -= protein_center
            else:
                for p in complex_graph['ligand'].pos:
                    p -= protein_center

            complex_graph.original_center = protein_center
            complex_graphs.append(complex_graph)
        return complex_graphs, ligs


def print_statistics(complex_graphs):
    statistics = ([], [], [], [])

    for complex_graph in complex_graphs:
        lig_pos = complex_graph['ligand'].pos if torch.is_tensor(complex_graph['ligand'].pos) else complex_graph['ligand'].pos[0]
        radius_protein = torch.max(torch.linalg.vector_norm(complex_graph['receptor'].pos, dim=1))
        molecule_center = torch.mean(lig_pos, dim=0)
        radius_molecule = torch.max(
            torch.linalg.vector_norm(lig_pos - molecule_center.unsqueeze(0), dim=1))
        distance_center = torch.linalg.vector_norm(molecule_center)
        statistics[0].append(radius_protein)
        statistics[1].append(radius_molecule)
        statistics[2].append(distance_center)
        if "rmsd_matching" in complex_graph:
            statistics[3].append(complex_graph.rmsd_matching)
        else:
            statistics[3].append(0)

    name = ['radius protein', 'radius molecule', 'distance protein-mol', 'rmsd matching']
    print('Number of complexes: ', len(complex_graphs))
    for i in range(4):
        array = np.asarray(statistics[i])
        print(f"{name[i]}: mean {np.mean(array)}, std {np.std(array)}, max {np.max(array)}")


def construct_loader(args, t_to_sigma):
    transform = NoiseTransform(t_to_sigma=t_to_sigma, no_torsion=args.no_torsion,
                               all_atom=args.all_atoms)

    common_args = {'transform': transform, 'root': args.data_dir, 'limit_complexes': args.limit_complexes,
                   'receptor_radius': args.receptor_radius,
                   'c_alpha_max_neighbors': args.c_alpha_max_neighbors,
                   'remove_hs': args.remove_hs, 'max_lig_size': args.max_lig_size,
                   'matching': not args.no_torsion, 'popsize': args.matching_popsize, 'maxiter': args.matching_maxiter,
                   'num_workers': args.num_workers, 'all_atoms': args.all_atoms,
                   'atom_radius': args.atom_radius, 'atom_max_neighbors': args.atom_max_neighbors,
                   'esm_embeddings_path': args.esm_embeddings_path}

    train_dataset = PDBBind(cache_path=args.cache_path, split_path=args.split_train, keep_original=True,
                            num_conformers=args.num_conformers, **common_args)
    val_dataset = PDBBind(cache_path=args.cache_path, split_path=args.split_val, keep_original=True, **common_args)

    loader_class = DataListLoader if torch.cuda.is_available() else DataLoader
    train_loader = loader_class(dataset=train_dataset, batch_size=args.batch_size, num_workers=args.num_dataloader_workers, shuffle=True, pin_memory=args.pin_memory)
    val_loader = loader_class(dataset=val_dataset, batch_size=args.batch_size, num_workers=args.num_dataloader_workers, shuffle=True, pin_memory=args.pin_memory)

    return train_loader, val_loader


def read_mol(pdbbind_dir, name, remove_hs=False):
    lig = read_molecule(os.path.join(pdbbind_dir, name, f'{name}_ligand.sdf'), remove_hs=remove_hs, sanitize=True)
    if lig is None:  # read mol2 file if sdf file cannot be sanitized
        lig = read_molecule(os.path.join(pdbbind_dir, name, f'{name}_ligand.mol2'), remove_hs=remove_hs, sanitize=True)
    return lig


def read_mols(pdbbind_dir, name, remove_hs=False):
    ligs = []
    for file in os.listdir(os.path.join(pdbbind_dir, name)):
        if file.endswith(".sdf") and 'rdkit' not in file:
            lig = read_molecule(os.path.join(pdbbind_dir, name, file), remove_hs=remove_hs, sanitize=True)
            if lig is None and os.path.exists(os.path.join(pdbbind_dir, name, file[:-4] + ".mol2")):  # read mol2 file if sdf file cannot be sanitized
                print('Using the .sdf file failed. We found a .mol2 file instead and are trying to use that.')
                lig = read_molecule(os.path.join(pdbbind_dir, name, file[:-4] + ".mol2"), remove_hs=remove_hs, sanitize=True)
            if lig is not None:
                ligs.append(lig)
    return ligs