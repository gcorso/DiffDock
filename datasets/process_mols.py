import copy
import os
import warnings

import numpy as np
import scipy.spatial as spa
import torch
from Bio.PDB import PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from rdkit import Chem
from rdkit.Chem.rdchem import BondType as BT
from rdkit.Chem import AllChem, GetPeriodicTable, RemoveHs
from rdkit.Geometry import Point3D
from scipy import spatial
from scipy.special import softmax
from torch_cluster import radius_graph


import torch.nn.functional as F

from datasets.conformer_matching import get_torsion_angles, optimize_rotatable_bonds
from utils.torsion import get_transformation_mask


biopython_parser = PDBParser()
periodic_table = GetPeriodicTable()
allowable_features = {
    'possible_atomic_num_list': list(range(1, 119)) + ['misc'],
    'possible_chirality_list': [
        'CHI_UNSPECIFIED',
        'CHI_TETRAHEDRAL_CW',
        'CHI_TETRAHEDRAL_CCW',
        'CHI_OTHER'
    ],
    'possible_degree_list': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 'misc'],
    'possible_numring_list': [0, 1, 2, 3, 4, 5, 6, 'misc'],
    'possible_implicit_valence_list': [0, 1, 2, 3, 4, 5, 6, 'misc'],
    'possible_formal_charge_list': [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 'misc'],
    'possible_numH_list': [0, 1, 2, 3, 4, 5, 6, 7, 8, 'misc'],
    'possible_number_radical_e_list': [0, 1, 2, 3, 4, 'misc'],
    'possible_hybridization_list': [
        'SP', 'SP2', 'SP3', 'SP3D', 'SP3D2', 'misc'
    ],
    'possible_is_aromatic_list': [False, True],
    'possible_is_in_ring3_list': [False, True],
    'possible_is_in_ring4_list': [False, True],
    'possible_is_in_ring5_list': [False, True],
    'possible_is_in_ring6_list': [False, True],
    'possible_is_in_ring7_list': [False, True],
    'possible_is_in_ring8_list': [False, True],
    'possible_amino_acids': ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET',
                             'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HIP', 'HIE', 'TPO', 'HID', 'LEV', 'MEU',
                             'PTR', 'GLV', 'CYT', 'SEP', 'HIZ', 'CYM', 'GLM', 'ASQ', 'TYS', 'CYX', 'GLZ', 'misc'],
    'possible_atom_type_2': ['C*', 'CA', 'CB', 'CD', 'CE', 'CG', 'CH', 'CZ', 'N*', 'ND', 'NE', 'NH', 'NZ', 'O*', 'OD',
                             'OE', 'OG', 'OH', 'OX', 'S*', 'SD', 'SG', 'misc'],
    'possible_atom_type_3': ['C', 'CA', 'CB', 'CD', 'CD1', 'CD2', 'CE', 'CE1', 'CE2', 'CE3', 'CG', 'CG1', 'CG2', 'CH2',
                             'CZ', 'CZ2', 'CZ3', 'N', 'ND1', 'ND2', 'NE', 'NE1', 'NE2', 'NH1', 'NH2', 'NZ', 'O', 'OD1',
                             'OD2', 'OE1', 'OE2', 'OG', 'OG1', 'OH', 'OXT', 'SD', 'SG', 'misc'],
}
bonds = {BT.SINGLE: 0, BT.DOUBLE: 1, BT.TRIPLE: 2, BT.AROMATIC: 3}

lig_feature_dims = (list(map(len, [
    allowable_features['possible_atomic_num_list'],
    allowable_features['possible_chirality_list'],
    allowable_features['possible_degree_list'],
    allowable_features['possible_formal_charge_list'],
    allowable_features['possible_implicit_valence_list'],
    allowable_features['possible_numH_list'],
    allowable_features['possible_number_radical_e_list'],
    allowable_features['possible_hybridization_list'],
    allowable_features['possible_is_aromatic_list'],
    allowable_features['possible_numring_list'],
    allowable_features['possible_is_in_ring3_list'],
    allowable_features['possible_is_in_ring4_list'],
    allowable_features['possible_is_in_ring5_list'],
    allowable_features['possible_is_in_ring6_list'],
    allowable_features['possible_is_in_ring7_list'],
    allowable_features['possible_is_in_ring8_list'],
])), 0)  # number of scalar features

rec_atom_feature_dims = (list(map(len, [
    allowable_features['possible_amino_acids'],
    allowable_features['possible_atomic_num_list'],
    allowable_features['possible_atom_type_2'],
    allowable_features['possible_atom_type_3'],
])), 0)

rec_residue_feature_dims = (list(map(len, [
    allowable_features['possible_amino_acids']
])), 0)


def lig_atom_featurizer(mol):
    ringinfo = mol.GetRingInfo()
    atom_features_list = []
    for idx, atom in enumerate(mol.GetAtoms()):
        atom_features_list.append([
            safe_index(allowable_features['possible_atomic_num_list'], atom.GetAtomicNum()),
            allowable_features['possible_chirality_list'].index(str(atom.GetChiralTag())),
            safe_index(allowable_features['possible_degree_list'], atom.GetTotalDegree()),
            safe_index(allowable_features['possible_formal_charge_list'], atom.GetFormalCharge()),
            safe_index(allowable_features['possible_implicit_valence_list'], atom.GetImplicitValence()),
            safe_index(allowable_features['possible_numH_list'], atom.GetTotalNumHs()),
            safe_index(allowable_features['possible_number_radical_e_list'], atom.GetNumRadicalElectrons()),
            safe_index(allowable_features['possible_hybridization_list'], str(atom.GetHybridization())),
            allowable_features['possible_is_aromatic_list'].index(atom.GetIsAromatic()),
            safe_index(allowable_features['possible_numring_list'], ringinfo.NumAtomRings(idx)),
            allowable_features['possible_is_in_ring3_list'].index(ringinfo.IsAtomInRingOfSize(idx, 3)),
            allowable_features['possible_is_in_ring4_list'].index(ringinfo.IsAtomInRingOfSize(idx, 4)),
            allowable_features['possible_is_in_ring5_list'].index(ringinfo.IsAtomInRingOfSize(idx, 5)),
            allowable_features['possible_is_in_ring6_list'].index(ringinfo.IsAtomInRingOfSize(idx, 6)),
            allowable_features['possible_is_in_ring7_list'].index(ringinfo.IsAtomInRingOfSize(idx, 7)),
            allowable_features['possible_is_in_ring8_list'].index(ringinfo.IsAtomInRingOfSize(idx, 8)),
        ])

    return torch.tensor(atom_features_list)


def rec_residue_featurizer(rec):
    feature_list = []
    for residue in rec.get_residues():
        feature_list.append([safe_index(allowable_features['possible_amino_acids'], residue.get_resname())])
    return torch.tensor(feature_list, dtype=torch.float32)  # (N_res, 1)


def safe_index(l, e):
    """ Return index of element e in list l. If e is not present, return the last index """
    try:
        return l.index(e)
    except:
        return len(l) - 1



def parse_receptor(pdbid, pdbbind_dir):
    rec = parsePDB(pdbid, pdbbind_dir)
    return rec


def parsePDB(pdbid, pdbbind_dir):
    rec_path = os.path.join(pdbbind_dir, pdbid, f'{pdbid}_protein_processed.pdb')
    return parse_pdb_from_path(rec_path)

def parse_pdb_from_path(path):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=PDBConstructionWarning)
        structure = biopython_parser.get_structure('random_id', path)
        rec = structure[0]
    return rec


def extract_receptor_structure(rec, lig, lm_embedding_chains=None):
    conf = lig.GetConformer()
    lig_coords = conf.GetPositions()
    min_distances = []
    coords = []
    c_alpha_coords = []
    n_coords = []
    c_coords = []
    valid_chain_ids = []
    lengths = []
    for i, chain in enumerate(rec):
        chain_coords = []  # num_residues, num_atoms, 3
        chain_c_alpha_coords = []
        chain_n_coords = []
        chain_c_coords = []
        count = 0
        invalid_res_ids = []
        for res_idx, residue in enumerate(chain):
            if residue.get_resname() == 'HOH':
                invalid_res_ids.append(residue.get_id())
                continue
            residue_coords = []
            c_alpha, n, c = None, None, None
            for atom in residue:
                if atom.name == 'CA':
                    c_alpha = list(atom.get_vector())
                if atom.name == 'N':
                    n = list(atom.get_vector())
                if atom.name == 'C':
                    c = list(atom.get_vector())
                residue_coords.append(list(atom.get_vector()))

            if c_alpha != None and n != None and c != None:
                # only append residue if it is an amino acid and not some weird molecule that is part of the complex
                chain_c_alpha_coords.append(c_alpha)
                chain_n_coords.append(n)
                chain_c_coords.append(c)
                chain_coords.append(np.array(residue_coords))
                count += 1
            else:
                invalid_res_ids.append(residue.get_id())
        for res_id in invalid_res_ids:
            chain.detach_child(res_id)
        if len(chain_coords) > 0:
            all_chain_coords = np.concatenate(chain_coords, axis=0)
            distances = spatial.distance.cdist(lig_coords, all_chain_coords)
            min_distance = distances.min()
        else:
            min_distance = np.inf

        min_distances.append(min_distance)
        lengths.append(count)
        coords.append(chain_coords)
        c_alpha_coords.append(np.array(chain_c_alpha_coords))
        n_coords.append(np.array(chain_n_coords))
        c_coords.append(np.array(chain_c_coords))
        if not count == 0: valid_chain_ids.append(chain.get_id())

    min_distances = np.array(min_distances)
    if len(valid_chain_ids) == 0:
        valid_chain_ids.append(np.argmin(min_distances))
    valid_coords = []
    valid_c_alpha_coords = []
    valid_n_coords = []
    valid_c_coords = []
    valid_lengths = []
    invalid_chain_ids = []
    valid_lm_embeddings = []
    for i, chain in enumerate(rec):
        if chain.get_id() in valid_chain_ids:
            valid_coords.append(coords[i])
            valid_c_alpha_coords.append(c_alpha_coords[i])
            if lm_embedding_chains is not None:
                if i >= len(lm_embedding_chains):
                    raise ValueError('Encountered valid chain id that was not present in the LM embeddings')
                valid_lm_embeddings.append(lm_embedding_chains[i])
            valid_n_coords.append(n_coords[i])
            valid_c_coords.append(c_coords[i])
            valid_lengths.append(lengths[i])
        else:
            invalid_chain_ids.append(chain.get_id())
    coords = [item for sublist in valid_coords for item in sublist]  # list with n_residues arrays: [n_atoms, 3]

    c_alpha_coords = np.concatenate(valid_c_alpha_coords, axis=0)  # [n_residues, 3]
    n_coords = np.concatenate(valid_n_coords, axis=0)  # [n_residues, 3]
    c_coords = np.concatenate(valid_c_coords, axis=0)  # [n_residues, 3]
    lm_embeddings = np.concatenate(valid_lm_embeddings, axis=0) if lm_embedding_chains is not None else None
    for invalid_id in invalid_chain_ids:
        rec.detach_child(invalid_id)

    assert len(c_alpha_coords) == len(n_coords)
    assert len(c_alpha_coords) == len(c_coords)
    assert sum(valid_lengths) == len(c_alpha_coords)
    return rec, coords, c_alpha_coords, n_coords, c_coords, lm_embeddings


def get_lig_graph(mol, complex_graph):
    lig_coords = torch.from_numpy(mol.GetConformer().GetPositions()).float()
    atom_feats = lig_atom_featurizer(mol)

    row, col, edge_type = [], [], []
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        row += [start, end]
        col += [end, start]
        edge_type += 2 * [bonds[bond.GetBondType()]] if bond.GetBondType() != BT.UNSPECIFIED else [0, 0]

    edge_index = torch.tensor([row, col], dtype=torch.long)
    edge_type = torch.tensor(edge_type, dtype=torch.long)
    edge_attr = F.one_hot(edge_type, num_classes=len(bonds)).to(torch.float)

    complex_graph['ligand'].x = atom_feats
    complex_graph['ligand'].pos = lig_coords
    complex_graph['ligand', 'lig_bond', 'ligand'].edge_index = edge_index
    complex_graph['ligand', 'lig_bond', 'ligand'].edge_attr = edge_attr
    return

def generate_conformer(mol):
    ps = AllChem.ETKDGv2()
    id = AllChem.EmbedMolecule(mol, ps)
    if id == -1:
        print('rdkit coords could not be generated without using random coords. using random coords now.')
        ps.useRandomCoords = True
        AllChem.EmbedMolecule(mol, ps)
        AllChem.MMFFOptimizeMolecule(mol, confId=0)
    # else:
    #    AllChem.MMFFOptimizeMolecule(mol_rdkit, confId=0)

def get_lig_graph_with_matching(mol_, complex_graph, popsize, maxiter, matching, keep_original, num_conformers, remove_hs):
    if matching:
        mol_maybe_noh = copy.deepcopy(mol_)
        if remove_hs:
            mol_maybe_noh = RemoveHs(mol_maybe_noh, sanitize=True)
        if keep_original:
            complex_graph['ligand'].orig_pos = mol_maybe_noh.GetConformer().GetPositions()

        rotable_bonds = get_torsion_angles(mol_maybe_noh)
        if not rotable_bonds: print("no_rotable_bonds but still using it")

        for i in range(num_conformers):
            mol_rdkit = copy.deepcopy(mol_)

            mol_rdkit.RemoveAllConformers()
            mol_rdkit = AllChem.AddHs(mol_rdkit)
            generate_conformer(mol_rdkit)
            if remove_hs:
                mol_rdkit = RemoveHs(mol_rdkit, sanitize=True)
            mol = copy.deepcopy(mol_maybe_noh)
            if rotable_bonds:
                optimize_rotatable_bonds(mol_rdkit, mol, rotable_bonds, popsize=popsize, maxiter=maxiter)
            mol.AddConformer(mol_rdkit.GetConformer())
            rms_list = []
            AllChem.AlignMolConformers(mol, RMSlist=rms_list)
            mol_rdkit.RemoveAllConformers()
            mol_rdkit.AddConformer(mol.GetConformers()[1])

            if i == 0:
                complex_graph.rmsd_matching = rms_list[0]
                get_lig_graph(mol_rdkit, complex_graph)
            else:
                if torch.is_tensor(complex_graph['ligand'].pos):
                    complex_graph['ligand'].pos = [complex_graph['ligand'].pos]
                complex_graph['ligand'].pos.append(torch.from_numpy(mol_rdkit.GetConformer().GetPositions()).float())

    else:  # no matching
        complex_graph.rmsd_matching = 0
        if remove_hs: mol_ = RemoveHs(mol_)
        get_lig_graph(mol_, complex_graph)

    edge_mask, mask_rotate = get_transformation_mask(complex_graph)
    complex_graph['ligand'].edge_mask = torch.tensor(edge_mask)
    complex_graph['ligand'].mask_rotate = mask_rotate

    return


def get_calpha_graph(rec, c_alpha_coords, n_coords, c_coords, complex_graph, cutoff=20, max_neighbor=None, lm_embeddings=None):
    n_rel_pos = n_coords - c_alpha_coords
    c_rel_pos = c_coords - c_alpha_coords
    num_residues = len(c_alpha_coords)
    if num_residues <= 1:
        raise ValueError(f"rec contains only 1 residue!")

    # Build the k-NN graph
    distances = spa.distance.cdist(c_alpha_coords, c_alpha_coords)
    src_list = []
    dst_list = []
    mean_norm_list = []
    for i in range(num_residues):
        dst = list(np.where(distances[i, :] < cutoff)[0])
        dst.remove(i)
        if max_neighbor != None and len(dst) > max_neighbor:
            dst = list(np.argsort(distances[i, :]))[1: max_neighbor + 1]
        if len(dst) == 0:
            dst = list(np.argsort(distances[i, :]))[1:2]  # choose second because first is i itself
            print(f'The c_alpha_cutoff {cutoff} was too small for one c_alpha such that it had no neighbors. '
                  f'So we connected it to the closest other c_alpha')
        assert i not in dst
        src = [i] * len(dst)
        src_list.extend(src)
        dst_list.extend(dst)
        valid_dist = list(distances[i, dst])
        valid_dist_np = distances[i, dst]
        sigma = np.array([1., 2., 5., 10., 30.]).reshape((-1, 1))
        weights = softmax(- valid_dist_np.reshape((1, -1)) ** 2 / sigma, axis=1)  # (sigma_num, neigh_num)
        assert weights[0].sum() > 1 - 1e-2 and weights[0].sum() < 1.01
        diff_vecs = c_alpha_coords[src, :] - c_alpha_coords[dst, :]  # (neigh_num, 3)
        mean_vec = weights.dot(diff_vecs)  # (sigma_num, 3)
        denominator = weights.dot(np.linalg.norm(diff_vecs, axis=1))  # (sigma_num,)
        mean_vec_ratio_norm = np.linalg.norm(mean_vec, axis=1) / denominator  # (sigma_num,)
        mean_norm_list.append(mean_vec_ratio_norm)
    assert len(src_list) == len(dst_list)

    node_feat = rec_residue_featurizer(rec)
    mu_r_norm = torch.from_numpy(np.array(mean_norm_list).astype(np.float32))
    side_chain_vecs = torch.from_numpy(
        np.concatenate([np.expand_dims(n_rel_pos, axis=1), np.expand_dims(c_rel_pos, axis=1)], axis=1))

    complex_graph['receptor'].x = torch.cat([node_feat, torch.tensor(lm_embeddings)], axis=1) if lm_embeddings is not None else node_feat
    complex_graph['receptor'].pos = torch.from_numpy(c_alpha_coords).float()
    complex_graph['receptor'].mu_r_norm = mu_r_norm
    complex_graph['receptor'].side_chain_vecs = side_chain_vecs.float()
    complex_graph['receptor', 'rec_contact', 'receptor'].edge_index = torch.from_numpy(np.asarray([src_list, dst_list]))

    return


def rec_atom_featurizer(rec):
    atom_feats = []
    for i, atom in enumerate(rec.get_atoms()):
        atom_name, element = atom.name, atom.element
        if element == 'CD':
            element = 'C'
        assert not element == ''
        try:
            atomic_num = periodic_table.GetAtomicNumber(element)
        except:
            atomic_num = -1
        atom_feat = [safe_index(allowable_features['possible_amino_acids'], atom.get_parent().get_resname()),
                     safe_index(allowable_features['possible_atomic_num_list'], atomic_num),
                     safe_index(allowable_features['possible_atom_type_2'], (atom_name + '*')[:2]),
                     safe_index(allowable_features['possible_atom_type_3'], atom_name)]
        atom_feats.append(atom_feat)

    return atom_feats


def get_rec_graph(rec, rec_coords, c_alpha_coords, n_coords, c_coords, complex_graph, rec_radius, c_alpha_max_neighbors=None, all_atoms=False,
                  atom_radius=5, atom_max_neighbors=None, remove_hs=False, lm_embeddings=None):
    if all_atoms:
        return get_fullrec_graph(rec, rec_coords, c_alpha_coords, n_coords, c_coords, complex_graph,
                                 c_alpha_cutoff=rec_radius, c_alpha_max_neighbors=c_alpha_max_neighbors,
                                 atom_cutoff=atom_radius, atom_max_neighbors=atom_max_neighbors, remove_hs=remove_hs,lm_embeddings=lm_embeddings)
    else:
        return get_calpha_graph(rec, c_alpha_coords, n_coords, c_coords, complex_graph, rec_radius, c_alpha_max_neighbors,lm_embeddings=lm_embeddings)


def get_fullrec_graph(rec, rec_coords, c_alpha_coords, n_coords, c_coords, complex_graph, c_alpha_cutoff=20,
                      c_alpha_max_neighbors=None, atom_cutoff=5, atom_max_neighbors=None, remove_hs=False, lm_embeddings=None):
    # builds the receptor graph with both residues and atoms

    n_rel_pos = n_coords - c_alpha_coords
    c_rel_pos = c_coords - c_alpha_coords
    num_residues = len(c_alpha_coords)
    if num_residues <= 1:
        raise ValueError(f"rec contains only 1 residue!")

    # Build the k-NN graph of residues
    distances = spa.distance.cdist(c_alpha_coords, c_alpha_coords)
    src_list = []
    dst_list = []
    mean_norm_list = []
    for i in range(num_residues):
        dst = list(np.where(distances[i, :] < c_alpha_cutoff)[0])
        dst.remove(i)
        if c_alpha_max_neighbors != None and len(dst) > c_alpha_max_neighbors:
            dst = list(np.argsort(distances[i, :]))[1: c_alpha_max_neighbors + 1]
        if len(dst) == 0:
            dst = list(np.argsort(distances[i, :]))[1:2]  # choose second because first is i itself
            print(f'The c_alpha_cutoff {c_alpha_cutoff} was too small for one c_alpha such that it had no neighbors. '
                  f'So we connected it to the closest other c_alpha')
        assert i not in dst
        src = [i] * len(dst)
        src_list.extend(src)
        dst_list.extend(dst)
        valid_dist = list(distances[i, dst])
        valid_dist_np = distances[i, dst]
        sigma = np.array([1., 2., 5., 10., 30.]).reshape((-1, 1))
        weights = softmax(- valid_dist_np.reshape((1, -1)) ** 2 / sigma, axis=1)  # (sigma_num, neigh_num)
        assert 1 - 1e-2 < weights[0].sum() < 1.01
        diff_vecs = c_alpha_coords[src, :] - c_alpha_coords[dst, :]  # (neigh_num, 3)
        mean_vec = weights.dot(diff_vecs)  # (sigma_num, 3)
        denominator = weights.dot(np.linalg.norm(diff_vecs, axis=1))  # (sigma_num,)
        mean_vec_ratio_norm = np.linalg.norm(mean_vec, axis=1) / denominator  # (sigma_num,)
        mean_norm_list.append(mean_vec_ratio_norm)
    assert len(src_list) == len(dst_list)

    node_feat = rec_residue_featurizer(rec)
    mu_r_norm = torch.from_numpy(np.array(mean_norm_list).astype(np.float32))
    side_chain_vecs = torch.from_numpy(
        np.concatenate([np.expand_dims(n_rel_pos, axis=1), np.expand_dims(c_rel_pos, axis=1)], axis=1))

    complex_graph['receptor'].x = torch.cat([node_feat, torch.tensor(lm_embeddings)], axis=1) if lm_embeddings is not None else node_feat
    complex_graph['receptor'].pos = torch.from_numpy(c_alpha_coords).float()
    complex_graph['receptor'].mu_r_norm = mu_r_norm
    complex_graph['receptor'].side_chain_vecs = side_chain_vecs.float()
    complex_graph['receptor', 'rec_contact', 'receptor'].edge_index = torch.from_numpy(np.asarray([src_list, dst_list]))

    src_c_alpha_idx = np.concatenate([np.asarray([i]*len(l)) for i, l in enumerate(rec_coords)])
    atom_feat = torch.from_numpy(np.asarray(rec_atom_featurizer(rec)))
    atom_coords = torch.from_numpy(np.concatenate(rec_coords, axis=0)).float()

    if remove_hs:
        not_hs = (atom_feat[:, 1] != 0)
        src_c_alpha_idx = src_c_alpha_idx[not_hs]
        atom_feat = atom_feat[not_hs]
        atom_coords = atom_coords[not_hs]

    atoms_edge_index = radius_graph(atom_coords, atom_cutoff, max_num_neighbors=atom_max_neighbors if atom_max_neighbors else 1000)
    atom_res_edge_index = torch.from_numpy(np.asarray([np.arange(len(atom_feat)), src_c_alpha_idx])).long()

    complex_graph['atom'].x = atom_feat
    complex_graph['atom'].pos = atom_coords
    complex_graph['atom', 'atom_contact', 'atom'].edge_index = atoms_edge_index
    complex_graph['atom', 'atom_rec_contact', 'receptor'].edge_index = atom_res_edge_index

    return

def write_mol_with_coords(mol, new_coords, path, add_hs=False):
    w = Chem.SDWriter(path)
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        x,y,z = new_coords.astype(np.double)[i]
        conf.SetAtomPosition(i,Point3D(x,y,z))
    if add_hs:
        mol = Chem.AddHs(mol, addCoords=True)
    w.write(mol)
    w.close()

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
        raise ValueError('Expect the format of the molecule_file to be '
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
    except Exception as e:
        print(e)
        print("RDKit was unable to read the molecule.")
        return None

    return mol


def read_sdf_or_mol2(sdf_fileName, mol2_fileName):

    mol = Chem.MolFromMolFile(sdf_fileName, sanitize=False)
    problem = False
    try:
        Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)
    except Exception as e:
        problem = True
    if problem:
        mol = Chem.MolFromMol2File(mol2_fileName, sanitize=False)
        try:
            Chem.SanitizeMol(mol)
            mol = Chem.RemoveHs(mol)
            problem = False
        except Exception as e:
            problem = True

    return mol, problem
