
import os
from argparse import ArgumentParser

import torch
from tqdm import tqdm


parser = ArgumentParser()
parser.add_argument('--esm_embeddings_path', type=str, default='data/embeddings_output', help='')
parser.add_argument('--output_path', type=str, default='data/esm2_3billion_embeddings.pt', help='')
args = parser.parse_args()

dict = {}
for filename in tqdm(os.listdir(args.esm_embeddings_path)):
    dict[filename.split('.')[0]] = torch.load(os.path.join(args.esm_embeddings_path,filename))['representations'][33]
torch.save(dict,args.output_path)