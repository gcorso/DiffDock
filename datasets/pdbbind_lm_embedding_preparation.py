import os
from argparse import FileType, ArgumentParser

import numpy as np
from Bio.PDB import PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

parser = ArgumentParser()
parser.add_argument('--data_dir', type=str, default='data/PDBBind_processed', help='')
parser.add_argument('--chain_cutoff', type=int, default=10, help='')
parser.add_argument('--out_file', type=str, default="data/pdbbind_sequences.fasta")
args = parser.parse_args()

cutoff = args.chain_cutoff
data_dir = args.data_dir
names = os.listdir(data_dir)
#%%
from Bio import SeqIO
biopython_parser = PDBParser()

three_to_one = {'ALA':	'A',
'ARG':	'R',
'ASN':	'N',
'ASP':	'D',
'CYS':	'C',
'GLN':	'Q',
'GLU':	'E',
'GLY':	'G',
'HIS':	'H',
'ILE':	'I',
'LEU':	'L',
'LYS':	'K',
'MET':	'M',
'MSE':  'M', # this is almost the same AA as MET. The sulfur is just replaced by Selen
'PHE':	'F',
'PRO':	'P',
'PYL':	'O',
'SER':	'S',
'SEC':	'U',
'THR':	'T',
'TRP':	'W',
'TYR':	'Y',
'VAL':	'V',
'ASX':	'B',
'GLX':	'Z',
'XAA':	'X',
'XLE':	'J'}

sequences = []
ids = []
for name in tqdm(names):
    if name == '.DS_Store': continue
    if os.path.exists(os.path.join(data_dir, name, f'{name}_protein_processed.pdb')):
        rec_path = os.path.join(data_dir, name, f'{name}_protein_processed.pdb')
    elif os.path.exists(os.path.join(data_dir, name, f'{name}_protein.pdb')):
        rec_path = os.path.join(data_dir, name, f'{name}_protein.pdb')
    else:
        continue
    if cutoff > 10:
        rec_path = os.path.join(data_dir, name, f'{name}_protein_obabel_reduce.pdb')
        if not os.path.exists(rec_path):
            rec_path = os.path.join(data_dir, name, f'{name}_protein.pdb')
    structure = biopython_parser.get_structure('random_id', rec_path)
    structure = structure[0]
    for i, chain in enumerate(structure):
        seq = ''
        for res_idx, residue in enumerate(chain):
            if residue.get_resname() == 'HOH':
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
            if c_alpha != None and n != None and c != None:  # only append residue if it is an amino acid and not
                try:
                    seq += three_to_one[residue.get_resname()]
                except Exception as e:
                    seq += '-'
                    print("encountered unknown AA: ", residue.get_resname(), ' in the complex ', name, '. Replacing it with a dash - .')
        sequences.append(seq)
        ids.append(f'{name}_chain_{i}')
records = []
for (index, seq) in zip(ids,sequences):
    record = SeqRecord(Seq(seq), str(index))
    record.description = ''
    records.append(record)
SeqIO.write(records, args.out_file, "fasta")


