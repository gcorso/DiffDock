import copy, time
import numpy as np
from collections import defaultdict
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolTransforms
from rdkit import Geometry
import networkx as nx
from scipy.optimize import differential_evolution

RDLogger.DisableLog('rdApp.*')

"""
    Conformer matching routines from Torsional Diffusion
"""

def GetDihedral(conf, atom_idx):
    return rdMolTransforms.GetDihedralRad(conf, atom_idx[0], atom_idx[1], atom_idx[2], atom_idx[3])


def SetDihedral(conf, atom_idx, new_vale):
    rdMolTransforms.SetDihedralRad(conf, atom_idx[0], atom_idx[1], atom_idx[2], atom_idx[3], new_vale)


def apply_changes(mol, values, rotatable_bonds, conf_id):
    opt_mol = copy.copy(mol)
    [SetDihedral(opt_mol.GetConformer(conf_id), rotatable_bonds[r], values[r]) for r in range(len(rotatable_bonds))]
    return opt_mol


def optimize_rotatable_bonds(mol, true_mol, rotatable_bonds, probe_id=-1, ref_id=-1, seed=0, popsize=15, maxiter=500,
                             mutation=(0.5, 1), recombination=0.8):
    opt = OptimizeConformer(mol, true_mol, rotatable_bonds, seed=seed, probe_id=probe_id, ref_id=ref_id)
    max_bound = [np.pi] * len(opt.rotatable_bonds)
    min_bound = [-np.pi] * len(opt.rotatable_bonds)
    bounds = (min_bound, max_bound)
    bounds = list(zip(bounds[0], bounds[1]))

    # Optimize conformations
    result = differential_evolution(opt.score_conformation, bounds,
                                    maxiter=maxiter, popsize=popsize,
                                    mutation=mutation, recombination=recombination, disp=False, seed=seed)
    opt_mol = apply_changes(opt.mol, result['x'], opt.rotatable_bonds, conf_id=probe_id)

    return opt_mol


class OptimizeConformer:
    def __init__(self, mol, true_mol, rotatable_bonds, probe_id=-1, ref_id=-1, seed=None):
        super(OptimizeConformer, self).__init__()
        if seed:
            np.random.seed(seed)
        self.rotatable_bonds = rotatable_bonds
        self.mol = mol
        self.true_mol = true_mol
        self.probe_id = probe_id
        self.ref_id = ref_id

    def score_conformation(self, values):
        for i, r in enumerate(self.rotatable_bonds):
            SetDihedral(self.mol.GetConformer(self.probe_id), r, values[i])
        return AllChem.AlignMol(self.mol, self.true_mol, self.probe_id, self.ref_id)


def get_torsion_angles(mol):
    torsions_list = []
    G = nx.Graph()
    for i, atom in enumerate(mol.GetAtoms()):
        G.add_node(i)
    nodes = set(G.nodes())
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        G.add_edge(start, end)
    for e in G.edges():
        G2 = copy.deepcopy(G)
        G2.remove_edge(*e)
        if nx.is_connected(G2): continue
        l = list(sorted(nx.connected_components(G2), key=len)[0])
        if len(l) < 2: continue
        n0 = list(G2.neighbors(e[0]))
        n1 = list(G2.neighbors(e[1]))
        torsions_list.append(
            (n0[0], e[0], e[1], n1[0])
        )
    return torsions_list

