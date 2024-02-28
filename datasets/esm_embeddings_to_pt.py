import os
from argparse import ArgumentParser

import torch
from tqdm import tqdm


parser = ArgumentParser()
parser.add_argument('--esm_embeddings_path', type=str, default='data/BindingMOAD_2020_ab_processed_biounit/moad_sequences_new', help='')
parser.add_argument('--output_path', type=str, default='data/BindingMOAD_2020_ab_processed_biounit/moad_sequences_new.pt', help='')
args = parser.parse_args()

dict = {}
for filename in tqdm(os.listdir(args.esm_embeddings_path)):
    dict[filename.split('.')[0]] = torch.load(os.path.join(args.esm_embeddings_path,filename))['representations'][33]
torch.save(dict,args.output_path)