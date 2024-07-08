import os
from esm import FastaBatchedDataset, pretrained
from rdkit.Chem import AddHs, MolFromSmiles
from torch_geometric.data import Dataset, HeteroData
import numpy as np
import torch
import prody as pr
import esm

from datasets.process_mols import generate_conformer, read_molecule, get_lig_graph_with_matching, moad_extract_receptor_structure
from datasets.parse_chi import aa_idx2aa_short, get_onehot_sequence


def get_sequences_from_pdbfile(file_path):
    sequence = None

    pdb = pr.parsePDB(file_path)
    seq = pdb.ca.getSequence()
    one_hot = get_onehot_sequence(seq)

    chain_ids = np.zeros(len(one_hot))
    res_chain_ids = pdb.ca.getChids()
    res_seg_ids = pdb.ca.getSegnames()
    res_chain_ids = np.asarray([s + c for s, c in zip(res_seg_ids, res_chain_ids)])
    ids = np.unique(res_chain_ids)

    for i, id in enumerate(ids):
        chain_ids[res_chain_ids == id] = i

        s_temp = np.argmax(one_hot[res_chain_ids == id], axis=1)
        s = ''.join([aa_idx2aa_short[aa_idx] for aa_idx in s_temp])

        if sequence is None:
            sequence = s
        else:
            sequence += (":" + s)

    return sequence


def set_nones(l):
    return [s if str(s) != 'nan' else None for s in l]


def get_sequences(protein_files, protein_sequences):
    new_sequences = []
    for i in range(len(protein_files)):
        if protein_files[i] is not None:
            new_sequences.append(get_sequences_from_pdbfile(protein_files[i]))
        else:
            new_sequences.append(protein_sequences[i])
    return new_sequences


def compute_ESM_embeddings(model, alphabet, labels, sequences):
    # settings used
    toks_per_batch = 4096
    repr_layers = [33]
    include = "per_tok"
    truncation_seq_length = 1022

    dataset = FastaBatchedDataset(labels, sequences)
    batches = dataset.get_batch_indices(toks_per_batch, extra_toks_per_seq=1)
    data_loader = torch.utils.data.DataLoader(
        dataset, collate_fn=alphabet.get_batch_converter(truncation_seq_length), batch_sampler=batches
    )

    assert all(-(model.num_layers + 1) <= i <= model.num_layers for i in repr_layers)
    repr_layers = [(i + model.num_layers + 1) % (model.num_layers + 1) for i in repr_layers]
    embeddings = {}

    with torch.no_grad():
        for batch_idx, (labels, strs, toks) in enumerate(data_loader):
            print(f"Processing {batch_idx + 1} of {len(batches)} batches ({toks.size(0)} sequences)")
            if torch.cuda.is_available():
                toks = toks.to(device="cuda", non_blocking=True)

            out = model(toks, repr_layers=repr_layers, return_contacts=False)
            representations = {layer: t.to(device="cpu") for layer, t in out["representations"].items()}

            for i, label in enumerate(labels):
                truncate_len = min(truncation_seq_length, len(strs[i]))
                embeddings[label] = representations[33][i, 1: truncate_len + 1].clone()
    return embeddings


def generate_ESM_structure(model, filename, sequence):
    model.set_chunk_size(256)
    chunk_size = 256
    output = None

    while output is None:
        try:
            with torch.no_grad():
                output = model.infer_pdb(sequence)

            with open(filename, "w") as f:
                f.write(output)
                print("saved", filename)
        except RuntimeError as e:
            if 'out of memory' in str(e):
                print('| WARNING: ran out of memory on chunk_size', chunk_size)
                for p in model.parameters():
                    if p.grad is not None:
                        del p.grad  # free some memory
                torch.cuda.empty_cache()
                chunk_size = chunk_size // 2
                if chunk_size > 2:
                    model.set_chunk_size(chunk_size)
                else:
                    print("Not enough memory for ESMFold")
                    break
            else:
                raise e
    return output is not None


class InferenceDataset(Dataset):
    def __init__(self, out_dir, complex_names, protein_files, ligand_descriptions, protein_sequences, lm_embeddings,
                 receptor_radius=30, c_alpha_max_neighbors=None, precomputed_lm_embeddings=None,
                 remove_hs=False, all_atoms=False, atom_radius=5, atom_max_neighbors=None, knn_only_graph=False):

        super(InferenceDataset, self).__init__()
        self.receptor_radius = receptor_radius
        self.c_alpha_max_neighbors = c_alpha_max_neighbors
        self.remove_hs = remove_hs
        self.all_atoms = all_atoms
        self.atom_radius, self.atom_max_neighbors = atom_radius, atom_max_neighbors
        self.knn_only_graph = knn_only_graph

        self.complex_names = complex_names
        self.protein_files = protein_files
        self.ligand_descriptions = ligand_descriptions
        self.protein_sequences = protein_sequences

        # generate LM embeddings
        if lm_embeddings and (precomputed_lm_embeddings is None or precomputed_lm_embeddings[0] is None):
            print("Generating ESM language model embeddings")
            model_location = "esm2_t33_650M_UR50D"
            model, alphabet = pretrained.load_model_and_alphabet(model_location)
            model.eval()
            if torch.cuda.is_available():
                model = model.cuda()

            protein_sequences = get_sequences(protein_files, protein_sequences)
            labels, sequences = [], []
            for i in range(len(protein_sequences)):
                s = protein_sequences[i].split(':')
                sequences.extend(s)
                labels.extend([complex_names[i] + '_chain_' + str(j) for j in range(len(s))])

            lm_embeddings = compute_ESM_embeddings(model, alphabet, labels, sequences)

            self.lm_embeddings = []
            for i in range(len(protein_sequences)):
                s = protein_sequences[i].split(':')
                self.lm_embeddings.append([lm_embeddings[f'{complex_names[i]}_chain_{j}'] for j in range(len(s))])

        elif not lm_embeddings:
            self.lm_embeddings = [None] * len(self.complex_names)

        else:
            self.lm_embeddings = precomputed_lm_embeddings

        # generate structures with ESMFold
        if None in protein_files:
            print("generating missing structures with ESMFold")
            model = esm.pretrained.esmfold_v1()
            model = model.eval().cuda()

            for i in range(len(protein_files)):
                if protein_files[i] is None:
                    self.protein_files[i] = f"{out_dir}/{complex_names[i]}/{complex_names[i]}_esmfold.pdb"
                    if not os.path.exists(self.protein_files[i]):
                        print("generating", self.protein_files[i])
                        generate_ESM_structure(model, self.protein_files[i], protein_sequences[i])

    def len(self):
        return len(self.complex_names)

    def get(self, idx):

        name, protein_file, ligand_description, lm_embedding = \
            self.complex_names[idx], self.protein_files[idx], self.ligand_descriptions[idx], self.lm_embeddings[idx]

        # build the pytorch geometric heterogeneous graph
        complex_graph = HeteroData()
        complex_graph['name'] = name

        # parse the ligand, either from file or smile
        try:
            mol = MolFromSmiles(ligand_description)  # check if it is a smiles or a path

            if mol is not None:
                mol = AddHs(mol)
                generate_conformer(mol)
            else:
                mol = read_molecule(ligand_description, remove_hs=False, sanitize=True)
                if mol is None:
                    raise Exception('RDKit could not read the molecule ', ligand_description)
                mol.RemoveAllConformers()
                mol = AddHs(mol)
                generate_conformer(mol)
        except Exception as e:
            print('Failed to read molecule ', ligand_description, ' We are skipping it. The reason is the exception: ', e)
            complex_graph['success'] = False
            return complex_graph

        try:
            # parse the receptor from the pdb file
            get_lig_graph_with_matching(mol, complex_graph, popsize=None, maxiter=None, matching=False, keep_original=False,
                                        num_conformers=1, remove_hs=self.remove_hs)

            moad_extract_receptor_structure(
                path=os.path.join(protein_file),
                complex_graph=complex_graph,
                neighbor_cutoff=self.receptor_radius,
                max_neighbors=self.c_alpha_max_neighbors,
                lm_embeddings=lm_embedding,
                knn_only_graph=self.knn_only_graph,
                all_atoms=self.all_atoms,
                atom_cutoff=self.atom_radius,
                atom_max_neighbors=self.atom_max_neighbors)

        except Exception as e:
            print(f'Skipping {name} because of the error:')
            print(e)
            complex_graph['success'] = False
            return complex_graph

        protein_center = torch.mean(complex_graph['receptor'].pos, dim=0, keepdim=True)
        complex_graph['receptor'].pos -= protein_center
        if self.all_atoms:
            complex_graph['atom'].pos -= protein_center

        ligand_center = torch.mean(complex_graph['ligand'].pos, dim=0, keepdim=True)
        complex_graph['ligand'].pos -= ligand_center

        complex_graph.original_center = protein_center
        complex_graph.mol = mol
        complex_graph['success'] = True
        return complex_graph
