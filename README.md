# DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking
[![Open in HuggingFace](https://huggingface.co/datasets/huggingface/badges/raw/main/open-in-hf-spaces-sm.svg)](https://huggingface.co/spaces/reginabarzilaygroup/DiffDock-Web)


![Alt Text](overview.png)

### [Original paper on arXiv](https://arxiv.org/abs/2210.01776)

Implementation of DiffDock, state-of-the-art method for molecular docking, by Gabriele Corso*, Hannes Stark*, Bowen Jing*, Regina Barzilay and Tommi Jaakkola.
This repository contains code and instructions to run the method. Since 2024, Jacob Silterra has been leading the effort to maintain and improve the code.
If you have any question, feel free to open an issue or reach out to us: [gcorso@mit.edu](gcorso@mit.edu) and [silterra@mit.edu](silterra@mit.edu).

**Update February 2024:** We have released DiffDock-L, a new version of DiffDock that provides a significant improvement in performance and generalization capacity (see the description of the new version in [our new paper](https://arxiv.org/abs/2402.18396)). By default the repository now runs the new model, please use GitHub commit history to run the original DiffDock model. Further we now provide instructions for Docker and to set up your own local UI interface.

You can also try out the model on [Hugging Face Spaces](https://huggingface.co/spaces/reginabarzilaygroup/DiffDock-Web).



<details><summary><b>Citation</b></summary>

If you use this code or the models in your research, please cite the following paper:

```bibtex
    @inproceedings{corso2023diffdock,
        title={DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking}, 
        author = {Corso, Gabriele and Stärk, Hannes and Jing, Bowen and Barzilay, Regina and Jaakkola, Tommi},
        booktitle={International Conference on Learning Representations (ICLR)},
        year={2023}
    }
```

If you use the latest version, DiffDock-L, please also cite the following paper:

```bibtex
    @inproceedings{corso2024discovery,
        title={Deep Confident Steps to New Pockets: Strategies for Docking Generalization},
        author={Corso, Gabriele and Deng, Arthur and Polizzi, Nicholas and Barzilay, Regina and Jaakkola, Tommi},
        booktitle={International Conference on Learning Representations (ICLR)},
        year={2024}
    }

```
</details>

<details open><summary><b>Table of contents</b></summary>

- [Usage](#usage)
  - [Quick Start](#quickstart)
  - [Setup Environment](#environment)
  - [Docking Prediction](#inference)
- [FAQ](#faq)
- [Datasets](#datasets)
- [Replicate results](#replicate)
- [Citations](#citations)
- [License](#license)
- [Acknowledgements](#acknowledgements)
</details>



## Usage  <a name="usage"></a>

### Quick Start  <a name="quickstart"></a>

You can directly try out the model without the need of installing anything through [Hugging Face Spaces](https://huggingface.co/spaces/reginabarzilaygroup/DiffDock-Web). Credit for the current HF interface goes to Jacob Silterra and for the previous version to Simon Duerr.

### Setup Environment  <a name="environment"></a>

We will set up the environment using [Anaconda](https://docs.anaconda.com/anaconda/install/index.html). Clone the
current repo

    git clone https://github.com/gcorso/DiffDock.git

To set up an appropriate environment, navigate to the root of the repository and run the following commands:

    conda env create --file environment.yml
    conda activate diffdock

See [conda documentation](https://conda.io/projects/conda/en/latest/commands/env/create.html) for more information.

### Using a Docker container

A Dockerfile is provided for building a container:

    docker build . -f Dockerfile -t diffdock

Alternatively, you can use a pre-built container to run the code.
First, download the container from Docker Hub:

    docker pull rbgcsail/diffdock

Check if you have a GPU available

    docker run --rm --gpus all nvidia/cuda:11.7.1-devel-ubuntu22.04 nvidia-smi

Then, run the container:

    docker run -it --gpus all --entrypoint /bin/bash rbgcsail/diffdock 

If you don't have a GPU, run (it will be significantly slower):
    
    docker run -it --entrypoint /bin/bash rbgcsail/diffdock 

Inside the container

    micromamba activate diffdock

You can now run the code as described below.

### Docking Prediction  <a name="inference"></a>

We support multiple input formats depending on whether you only want to make predictions for a single complex or for many at once.\
The protein inputs need to be `.pdb` files or sequences that will be folded with ESMFold. The ligand input can either be a SMILES string or a filetype that RDKit can read like `.sdf` or `.mol2`.

For a single complex: specify the protein with `--protein_path protein.pdb` or `--protein_sequence GIQSYCTPPYSVLQDPPQPVV` and the ligand with `--ligand ligand.sdf` or `--ligand "COc(cc1)ccc1C#N"`

For many complexes: create a csv file with paths to proteins and ligand files or SMILES. It contains as columns `complex_name` (name used to save predictions, can be left empty), `protein_path` (path to `.pdb` file, if empty uses sequence), `ligand_description` (SMILE or file path)  and `protein_sequence` (to fold with ESMFold in case the protein_path is empty).
An example .csv is at `data/protein_ligand_example.csv` and you would use it with `--protein_ligand_csv protein_ligand_example.csv`.

And you are ready to run inference:

    python -m inference --config default_inference_args.yaml  --protein_ligand_csv data/protein_ligand_example.csv --out_dir results/user_predictions_small 

When providing the `.pdb` files you can run DiffDock also on CPU, however, if possible, we recommend using a GPU as the model runs significantly faster. Note that the first time you run DiffDock on a device the program will precompute and store in cache look-up tables for SO(2) and SO(3) distributions (typically takes a couple of minutes), this won't be repeated in following runs.  

## Graphical UI  <a name="gui"></a>

We provide a simple graphical user interface to run DiffDock on a single complex. To use it, run the following command:

    python app/main.py

and navigate to http://localhost:7860 in your browser. 


## FAQ  <a name="faq"></a>

<details>
<summary><b>How to interpret the DiffDock output confidence score?</b> </summary>
It can be hard to interpret and compare confidence score of different complexes or different protein conformations, however, here a rough guideline that we typically use (c is the confidence score of the top pose):

 - c  > 0 high confidence
 - -1.5 < c < 0 moderate confidence
 - c < -1.5 low confidence

This is assuming the complex is similar to what DiffDock saw in the training set i.e. a not too large drug-like molecule bound to medium size protein (1 or 2 chains) in a conformation that is similar to the bound one (e.g. if it comes from an homologue crystal structure). If you are dealing with a large ligand, a large protein complex and/or an app/unbound protein conformation you should shift these intervals down.
</details>


<details>
<summary><b>Does DiffDock predict the binding affinity of the ligand to the protein?</b> </summary>
No, DiffDock does not predict the binding affinity of the ligand to the protein. It predicts the 3D structure of the complex and it outputs a confidence score. This latter is a measure of the quality of the prediction, i.e. the model's confidence in its prediction of the binding structure. Several of our collaborators have seen this to have some correlation with binding affinity (intuitively if a ligand does not bind there will be no good pose), but it is not a direct measure of it. 

We are working on better affinity prediction models, but in the meantime we recommend combining DiffDock's prediction with other tools such as docking function (e.g. GNINA), MM/GBSA or absolute binding free energy calculations. For this we recommend to first relax the DiffDock's structure predictions with the tool/force field used for the affinity prediction.
</details>

<details>
<summary><b>Can I use DiffDock for protein-protein or protein-nucleic acid interactions?</b> </summary>
While the program might not throw and error when fed with a large biomolecules as input, the model has only been designed, trained and tested for small molecule docking to proteins. Therefore, DiffDock is only likely to be able to deal with small peptides and nucleic acids as ligands, we do not recommend using DiffDock for the interactions of larger biomolecules. For other interactions we recommend looking at [DiffDock-PP](https://github.com/ketatam/DiffDock-PP) (rigid protein-protein interactions), [AlphaFold-Multimer](https://github.com/google-deepmind/alphafold) (flexible protein-protein interactions) or [RoseTTAFold2NA](https://github.com/uw-ipd/RoseTTAFold2NA) (protein-nucleic acid interactions).
</details>

## Datasets  <a name="datasets"></a>

The files in `data` contain the splits used for the various datasets. Below instructions for how to download each of the different datasets used for training and evaluation:

 - **PDBBind:** download the processed complexes from [zenodo](https://zenodo.org/record/6408497), unzip the directory and place it into `data` such that you have the path `data/PDBBind_processed`.
 - **BindingMOAD:** download the processed complexes from [zenodo](https://zenodo.org/records/10656052) under `BindingMOAD_2020_processed.tar`, unzip the directory and place it into `data` such that you have the path `data/BindingMOAD_2020_processed`.
 - **DockGen:** to evaluate the performance of `DiffDock-L` with this repository you should use directly the data from BindingMOAD above. For other purposes you can download exclusively the complexes of the DockGen benchmark already processed (e.g. chain cutoff) from [zenodo](https://zenodo.org/records/10656052) downloading the `DockGen.tar` file.
 - **PoseBusters:** download the processed complexes from [zenodo](https://zenodo.org/records/8278563).
 - **van der Mers:** the protein structures used for the van der Mers data augmentation strategy were downloaded [here](https://files.ipd.uw.edu/pub/training_sets/pdb_2021aug02.tar.gz).


## Replicate results  <a name="replicate"></a>

If you are interested in replicating the results of the original DiffDock paper please checkout to the following commit:

    git checkout v1.0

Otherwise download the data and place it as described in the "Dataset" section above.

### Generate the ESM2 embeddings for the proteins
To avoid having to compute ESM embeddings every time we evaluate on a dataset we first cache them and then run the evaluation script. Here the instructions for generating these for PDBBind but it also applies similarly to the other benchmarks. First run the following command to save the list of ESM embeddings:

    python datasets/esm_embedding_preparation.py

Use the generated file `data/pdbbind_sequences.fasta` to generate the ESM2 language model embeddings using the library https://github.com/facebookresearch/esm by installing their repository and executing the following in their repository:

    python scripts/extract.py esm2_t33_650M_UR50D pdbbind_sequences.fasta embeddings_output --repr_layers 33 --include per_tok --truncation_seq_length 4096

This generates the `embeddings_output` directory which you have to copy into the `data` folder of our repository to have `data/embeddings_output`.
Then run the command:

    python datasets/esm_embeddings_to_pt.py

### Run DiffDock-L

For PDBBind: 

    python -m evaluate  --config default_inference_args.yaml --split_path data/splits/timesplit_test --split_path data/splits/timesplit_test --batch_size 10 --esm_embeddings_path data/esm2_embeddings.pt --data_dir data/PDBBind_processed/ --tqdm --split test --chain_cutoff 10 --dataset pdbbind

For DockGen:

    python -m evaluate  --config default_inference_args.yaml --dataset moad --data_dir data/BindingMOAD_2020_processed --unroll_clusters --tqdm --split test --esm_embeddings_path data/moad_esm2_embeddings.pt --min_ligand_size 2 --moad_esm_embeddings_sequences_path data/moad_sequences_to_id.fasta --chain_cutoff 10 --batch_size 10 

For PoseBusters:

    python -m evaluate  --config default_inference_args.yaml --data_dir data/posebusters_benchmark_set --tqdm --dataset posebusters --split_path data/splits/posebusters_benchmark_set_ids.txt --esm_embeddings_path data/posebusters_ESM.pt --chain_cutoff 10 --batch_size 10 --protein_file protein --ligand_file ligands 

To additionally save the .sdf files of the generated molecules, add the flag `--save_visualisation`.

Note: the notebook `data/apo_alignment.ipynb` contains the code used to align the ESMFold-generated apo-structures to the holo-structures.

## Citations <a name="citations"></a>
If you use this code or the models in your research, please cite the following paper:

```bibtex
@inproceedings{corso2023diffdock,
    title={DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking}, 
    author = {Corso, Gabriele and Stärk, Hannes and Jing, Bowen and Barzilay, Regina and Jaakkola, Tommi},
    booktitle={International Conference on Learning Representations (ICLR)},
    year={2023}
}
```

If you use the latest version of our model, DiffDock-L, please also cite the following paper:

```bibtex
@inproceedings{corso2024discovery,
    title={Deep Confident Steps to New Pockets: Strategies for Docking Generalization},
    author={Corso, Gabriele and Deng, Arthur and Polizzi, Nicholas and Barzilay, Regina and Jaakkola, Tommi},
    booktitle={International Conference on Learning Representations (ICLR)},
    year={2024}
}
```

## License <a name="license"></a>
The code and model weights are released under MIT license. See the [LICENSE](LICENSE) file for details.

Components of the code of the [spyrmsd](spyrmsd) package by Rocco Meli (also MIT license) were integrated in the repo.

## Acknowledgements <a name="acknowledgements"></a>
We sincerely thank:
* Jacob Silterra for his help with the publishing and deployment of the code.
* Arthur Deng, Nicholas Polizzi and Ben Fry for their critical contributions to part of the code in this repository. 
* Wei Lu and Rachel Wu for pointing out some issues with the code.
