import copy
import warnings
import numpy as np
import torch
from Bio.PDB import PDBParser
from rdkit import Chem
from rdkit.Chem.rdchem import BondType as BT
from rdkit.Chem import AllChem, GetPeriodicTable, RemoveHs
from rdkit.Geometry import Point3D
from torch import cdist
from torch_cluster import knn_graph
import prody as pr

import torch.nn.functional as F

from datasets.conformer_matching import get_torsion_angles, optimize_rotatable_bonds
from datasets.constants import aa_short2long, atom_order, three_to_one
from datasets.parse_chi import get_chi_angles, get_coords, aa_idx2aa_short, get_onehot_sequence
from utils.torsion import get_transformation_mask
from utils.logging_utils import get_logger


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
        chiral_tag = str(atom.GetChiralTag())
        if chiral_tag  in ['CHI_SQUAREPLANAR', 'CHI_TRIGONALBIPYRAMIDAL', 'CHI_OCTAHEDRAL']:
            chiral_tag = 'CHI_OTHER'

        atom_features_list.append([
            safe_index(allowable_features['possible_atomic_num_list'], atom.GetAtomicNum()),
            allowable_features['possible_chirality_list'].index(str(chiral_tag)),
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
            #g_charge if not np.isnan(g_charge) and not np.isinf(g_charge) else 0.
        ])
    return torch.tensor(atom_features_list)


def safe_index(l, e):
    """ Return index of element e in list l. If e is not present, return the last index """
    try:
        return l.index(e)
    except:
        return len(l) - 1


def moad_extract_receptor_structure(path, complex_graph, neighbor_cutoff=20, max_neighbors=None, sequences_to_embeddings=None,
                                    knn_only_graph=False, lm_embeddings=None, all_atoms=False, atom_cutoff=None, atom_max_neighbors=None):
    # load the entire pdb file
    pdb = pr.parsePDB(path)
    seq = pdb.ca.getSequence()
    coords = get_coords(pdb)
    one_hot = get_onehot_sequence(seq)

    chain_ids = np.zeros(len(one_hot))
    res_chain_ids = pdb.ca.getChids()
    res_seg_ids = pdb.ca.getSegnames()
    res_chain_ids = np.asarray([s + c for s, c in zip(res_seg_ids, res_chain_ids)])
    ids = np.unique(res_chain_ids)
    sequences = []
    lm_embeddings = lm_embeddings if sequences_to_embeddings is None else []

    for i, id in enumerate(ids):
        chain_ids[res_chain_ids == id] = i

        s = np.argmax(one_hot[res_chain_ids == id], axis=1)
        s = ''.join([aa_idx2aa_short[aa_idx] for aa_idx in s])
        sequences.append(s)
        if sequences_to_embeddings is not None:
            lm_embeddings.append(sequences_to_embeddings[s])

    complex_graph['receptor'].sequence = sequences
    complex_graph['receptor'].chain_ids = torch.from_numpy(np.asarray(chain_ids)).long()

    new_extract_receptor_structure(seq, coords, complex_graph, neighbor_cutoff=neighbor_cutoff, max_neighbors=max_neighbors,
                                   lm_embeddings=lm_embeddings, knn_only_graph=knn_only_graph, all_atoms=all_atoms,
                                   atom_cutoff=atom_cutoff, atom_max_neighbors=atom_max_neighbors)


def new_extract_receptor_structure(seq, all_coords, complex_graph, neighbor_cutoff=20, max_neighbors=None, lm_embeddings=None,
                                   knn_only_graph=False, all_atoms=False, atom_cutoff=None, atom_max_neighbors=None):
    chi_angles, one_hot = get_chi_angles(all_coords, seq, return_onehot=True)
    n_rel_pos, c_rel_pos = all_coords[:, 0, :] - all_coords[:, 1, :], all_coords[:, 2, :] - all_coords[:, 1, :]
    side_chain_vecs = torch.from_numpy(np.concatenate([chi_angles / 360, n_rel_pos, c_rel_pos], axis=1))

    # Build the k-NN graph
    coords = torch.tensor(all_coords[:, 1, :], dtype=torch.float)
    if len(coords) > 3000:
        raise ValueError(f'The receptor is too large {len(coords)}')
    if knn_only_graph:
        edge_index = knn_graph(coords, k=max_neighbors if max_neighbors else 32)
    else:
        distances = cdist(coords, coords)
        src_list = []
        dst_list = []
        for i in range(len(coords)):
            dst = list(np.where(distances[i, :] < neighbor_cutoff)[0])
            dst.remove(i)
            max_neighbors = max_neighbors if max_neighbors else 1000
            if max_neighbors != None and len(dst) > max_neighbors:
                dst = list(np.argsort(distances[i, :]))[1: max_neighbors + 1]
            if len(dst) == 0:
                dst = list(np.argsort(distances[i, :]))[1:2]  # choose second because first is i itself
                print(
                    f'The cutoff {neighbor_cutoff} was too small for one atom such that it had no neighbors. '
                    f'So we connected it to the closest other atom')
            assert i not in dst
            src = [i] * len(dst)
            src_list.extend(src)
            dst_list.extend(dst)
        edge_index = torch.from_numpy(np.asarray([dst_list, src_list]))

    res_names_list = [aa_short2long[seq[i]] if seq[i] in aa_short2long else 'misc' for i in range(len(seq))]
    feature_list = [[safe_index(allowable_features['possible_amino_acids'], res)] for res in res_names_list]
    node_feat = torch.tensor(feature_list, dtype=torch.float32)

    lm_embeddings = torch.tensor(np.concatenate(lm_embeddings, axis=0)) if lm_embeddings is not None else None
    complex_graph['receptor'].x = torch.cat([node_feat, lm_embeddings], axis=1) if lm_embeddings is not None else node_feat
    complex_graph['receptor'].pos = coords
    complex_graph['receptor'].side_chain_vecs = side_chain_vecs.float()
    complex_graph['receptor', 'rec_contact', 'receptor'].edge_index = edge_index
    if all_atoms:
        atom_coords = all_coords.reshape(-1, 3)
        atom_coords = torch.from_numpy(atom_coords[~np.any(np.isnan(atom_coords), axis=1)]).float()

        if knn_only_graph:
            atoms_edge_index = knn_graph(atom_coords, k=atom_max_neighbors if atom_max_neighbors else 1000)
        else:
            atoms_distances = cdist(atom_coords, atom_coords)
            atom_src_list = []
            atom_dst_list = []
            for i in range(len(atom_coords)):
                dst = list(np.where(atoms_distances[i, :] < atom_cutoff)[0])
                dst.remove(i)
                max_neighbors = atom_max_neighbors if atom_max_neighbors else 1000
                if max_neighbors != None and len(dst) > max_neighbors:
                    dst = list(np.argsort(atoms_distances[i, :]))[1: max_neighbors + 1]
                if len(dst) == 0:
                    dst = list(np.argsort(atoms_distances[i, :]))[1:2]  # choose second because first is i itself
                    print(
                        f'The atom_cutoff {atom_cutoff} was too small for one atom such that it had no neighbors. '
                        f'So we connected it to the closest other atom')
                assert i not in dst
                src = [i] * len(dst)
                atom_src_list.extend(src)
                atom_dst_list.extend(dst)
            atoms_edge_index = torch.from_numpy(np.asarray([atom_dst_list, atom_src_list]))
        
        feats = [get_moad_atom_feats(res, all_coords[i]) for i, res in enumerate(seq)]
        atom_feat = torch.from_numpy(np.concatenate(feats, axis=0)).float()
        c_alpha_idx = np.concatenate([np.zeros(len(f)) + i for i, f in enumerate(feats)])
        np_array = np.stack([np.arange(len(atom_feat)), c_alpha_idx])
        atom_res_edge_index = torch.from_numpy(np_array).long()
        complex_graph['atom'].x = atom_feat
        complex_graph['atom'].pos = atom_coords
        assert len(complex_graph['atom'].x) == len(complex_graph['atom'].pos)
        complex_graph['atom', 'atom_contact', 'atom'].edge_index = atoms_edge_index
        complex_graph['atom', 'atom_rec_contact', 'receptor'].edge_index = atom_res_edge_index

    return


def get_moad_atom_feats(res, coords):
    feats = []
    res_long = aa_short2long[res]
    res_order = atom_order[res]
    for i, c in enumerate(coords):
        if np.any(np.isnan(c)):
            continue
        atom_feats = []
        if res == '-':
            atom_feats = [safe_index(allowable_features['possible_amino_acids'], 'misc'),
                     safe_index(allowable_features['possible_atomic_num_list'], 'misc'),
                     safe_index(allowable_features['possible_atom_type_2'], 'misc'),
                     safe_index(allowable_features['possible_atom_type_3'], 'misc')]
        else:
            atom_feats.append(safe_index(allowable_features['possible_amino_acids'], res_long))
            if i >= len(res_order):
                atom_feats.extend([safe_index(allowable_features['possible_atomic_num_list'], 'misc'),
                                   safe_index(allowable_features['possible_atom_type_2'], 'misc'),
                                   safe_index(allowable_features['possible_atom_type_3'], 'misc')])
            else:
                atom_name = res_order[i]
                try:
                    atomic_num = periodic_table.GetAtomicNumber(atom_name[:1])
                except:
                    print("element", res_order[i][:1], 'not found')
                    atomic_num = -1

                atom_feats.extend([safe_index(allowable_features['possible_atomic_num_list'], atomic_num),
                                   safe_index(allowable_features['possible_atom_type_2'], (atom_name + '*')[:2]),
                                   safe_index(allowable_features['possible_atom_type_3'], atom_name)])
        feats.append(atom_feats)
    feats = np.asarray(feats)
    return feats


def get_lig_graph(mol, complex_graph):
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
    complex_graph['ligand', 'lig_bond', 'ligand'].edge_index = edge_index
    complex_graph['ligand', 'lig_bond', 'ligand'].edge_attr = edge_attr

    if mol.GetNumConformers() > 0:
        lig_coords = torch.from_numpy(mol.GetConformer().GetPositions()).float()
        complex_graph['ligand'].pos = lig_coords

    return


def generate_conformer(mol):
    ps = AllChem.ETKDGv2()
    failures, id = 0, -1
    while failures < 3 and id == -1:
        if failures > 0:
            get_logger().debug(f'rdkit coords could not be generated. trying again {failures}.')
        id = AllChem.EmbedMolecule(mol, ps)
        failures += 1
    if id == -1:
        get_logger().info('rdkit coords could not be generated without using random coords. using random coords now.')
        ps.useRandomCoords = True
        AllChem.EmbedMolecule(mol, ps)
        AllChem.MMFFOptimizeMolecule(mol, confId=0)
        return True
    #else:
    #    AllChem.MMFFOptimizeMolecule(mol, confId=0)
    return False


def get_lig_graph_with_matching(mol_, complex_graph, popsize, maxiter, matching, keep_original, num_conformers, remove_hs, tries=10, skip_matching=False):
    if matching:
        mol_maybe_noh = copy.deepcopy(mol_)
        if remove_hs:
            mol_maybe_noh = RemoveHs(mol_maybe_noh, sanitize=True)
            mol_maybe_noh = AllChem.RemoveAllHs(mol_maybe_noh)
        if keep_original:
            positions = []
            for conf in mol_maybe_noh.GetConformers():
                positions.append(conf.GetPositions())
            complex_graph['ligand'].orig_pos = np.asarray(positions) if len(positions) > 1 else positions[0]

        # rotatable_bonds = get_torsion_angles(mol_maybe_noh)
        _tmp = copy.deepcopy(mol_)
        if remove_hs:
            _tmp = RemoveHs(_tmp, sanitize=True)
        _tmp = AllChem.RemoveAllHs(_tmp)
        rotatable_bonds = get_torsion_angles(_tmp)

        for i in range(num_conformers):
            mols, rmsds = [], []
            for _ in range(tries):
                mol_rdkit = copy.deepcopy(mol_)

                mol_rdkit.RemoveAllConformers()
                mol_rdkit = AllChem.AddHs(mol_rdkit)
                generate_conformer(mol_rdkit)
                if remove_hs:
                    mol_rdkit = RemoveHs(mol_rdkit, sanitize=True)
                mol_rdkit = AllChem.RemoveAllHs(mol_rdkit)
                mol = AllChem.RemoveAllHs(copy.deepcopy(mol_maybe_noh))
                if rotatable_bonds and not skip_matching:
                    optimize_rotatable_bonds(mol_rdkit, mol, rotatable_bonds, popsize=popsize, maxiter=maxiter)
                mol.AddConformer(mol_rdkit.GetConformer())
                rms_list = []
                AllChem.AlignMolConformers(mol, RMSlist=rms_list)
                mol_rdkit.RemoveAllConformers()
                mol_rdkit.AddConformer(mol.GetConformers()[1])
                mols.append(mol_rdkit)
                rmsds.append(rms_list[0])

            # select molecule with lowest rmsd
            #print("mean std min max", np.mean(rmsds), np.std(rmsds), np.min(rmsds), np.max(rmsds))
            mol_rdkit = mols[np.argmin(rmsds)]
            if i == 0:
                complex_graph.rmsd_matching = min(rmsds)
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


def get_rec_misc_atom_feat(bio_atom=None, atom_name=None, element=None, get_misc_features=False):
    if get_misc_features:
        return [safe_index(allowable_features['possible_amino_acids'], 'misc'),
                 safe_index(allowable_features['possible_atomic_num_list'], 'misc'),
                 safe_index(allowable_features['possible_atom_type_2'], 'misc'),
                 safe_index(allowable_features['possible_atom_type_3'], 'misc')]
    if atom_name is not None:
        atom_name = atom_name
    else:
        atom_name = bio_atom.name
    if element is not None:
        element = element
    else:
        element = bio_atom.element
    if element == 'CD':
        element = 'C'
    assert not element == ''
    try:
        atomic_num = periodic_table.GetAtomicNumber(element.lower().capitalize())
    except:
        atomic_num = -1

    atom_feat = [safe_index(allowable_features['possible_amino_acids'], bio_atom.get_parent().get_resname()),
                 safe_index(allowable_features['possible_atomic_num_list'], atomic_num),
                 safe_index(allowable_features['possible_atom_type_2'], (atom_name + '*')[:2]),
                 safe_index(allowable_features['possible_atom_type_3'], atom_name)]
    return atom_feat


def write_mol_with_coords(mol, new_coords, path):
    w = Chem.SDWriter(path)
    conf = mol.GetConformer()
    for i in range(mol.GetNumAtoms()):
        x,y,z = new_coords.astype(np.double)[i]
        conf.SetAtomPosition(i,Point3D(x,y,z))
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
        # Print stacktrace
        import traceback
        msg = traceback.format_exc()
        get_logger().warning(f"Failed to process molecule: {molecule_file}\n{msg}")
        return None

    return mol
