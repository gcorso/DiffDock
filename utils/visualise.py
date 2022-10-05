from rdkit.Chem.rdmolfiles import MolToPDBBlock, MolToPDBFile
import rdkit.Chem 
from rdkit import Geometry
from collections import defaultdict
import copy
import numpy as np
import torch

    
class PDBFile:
    def __init__(self, mol):
        self.parts = defaultdict(dict)
        self.mol = copy.deepcopy(mol)
        [self.mol.RemoveConformer(j) for j in range(mol.GetNumConformers()) if j]        
    def add(self, coords, order, part=0, repeat=1):
        if type(coords) in [rdkit.Chem.Mol, rdkit.Chem.RWMol]:
            block = MolToPDBBlock(coords).split('\n')[:-2]
            self.parts[part][order] = {'block': block, 'repeat': repeat}
            return
        elif type(coords) is np.ndarray:
            coords = coords.astype(np.float64)
        elif type(coords) is torch.Tensor:
            coords = coords.double().numpy()
        for i in range(coords.shape[0]):
            self.mol.GetConformer(0).SetAtomPosition(i, Geometry.Point3D(coords[i, 0], coords[i, 1], coords[i, 2]))
        block = MolToPDBBlock(self.mol).split('\n')[:-2]
        self.parts[part][order] = {'block': block, 'repeat': repeat}
        
    def write(self, path=None, limit_parts=None):
        is_first = True
        str_ = ''
        for part in sorted(self.parts.keys()):
            if limit_parts and part >= limit_parts:
                break
            part = self.parts[part]
            keys_positive = sorted(filter(lambda x: x >=0, part.keys()))
            keys_negative = sorted(filter(lambda x: x < 0, part.keys()))
            keys = list(keys_positive) + list(keys_negative)
            for key in keys:
                block = part[key]['block']
                times = part[key]['repeat']
                for _ in range(times):
                    if not is_first:
                        block = [line for line in block if 'CONECT' not in line]
                    is_first = False
                    str_ += 'MODEL\n'
                    str_ += '\n'.join(block)
                    str_ += '\nENDMDL\n'
        if not path:
            return str_
        with open(path, 'w') as f:
            f.write(str_)