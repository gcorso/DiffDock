
import os
import pickle
from argparse import ArgumentParser

import torch
from tqdm import tqdm


parser = ArgumentParser()
parser.add_argument('--esm_embeddings_path', type=str, default='data/BindingMOAD_2020_ab_processed_biounit/moad_sequences_new', help='')
parser.add_argument('--output_path', type=str, default='data/BindingMOAD_2020_ab_processed_biounit/moad_sequences_new.pt', help='')
args = parser.parse_args()

dic = {}

# read text file with all sequences
with open('data/pdb_2021aug02/sequences_to_id.fasta') as f:
    lines = f.readlines()

# read sequences
with open('data/pdb_2021aug02/useful_sequences.pkl', 'rb') as f:
    sequences = pickle.load(f)

ids = set()
dict_seq_id = {seq[:-1]: str(id) for id, seq in enumerate(lines)}

for i, seq in tqdm(enumerate(sequences)):
    ids.add(dict_seq_id[seq])
    if i == 20000: break

print("total", len(ids), "out of", len(os.listdir(args.esm_embeddings_path)))

available = set([filename.split('.')[0] for filename in os.listdir(args.esm_embeddings_path)])
final = available.intersection(ids)

for idp in tqdm(final):
    dic[idp] = torch.load(os.path.join(args.esm_embeddings_path, idp+'.pt'))['representations'][33]
torch.save(dic,args.output_path)