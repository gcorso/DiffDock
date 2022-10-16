#!/bin/bash

# edit here if you put diffdock into another folder
DIFF_DOCK_HOME=/repo/DiffDock/
ESM_HOME=/repo/esm
ESM_PRETRAINED=/mnt/db/esm/

source activate diffdock

set -e

usage() {
        echo ""
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters:"
        # edited by Yinying
        echo "-p <protein_pdb>          Protein in PDB format."
        echo "-l <ligand>               Ligand in SMILE string or in sdf/mol2 format"
        echo "-L <ligand_name>          Ligand name"
        echo "-o <save_dir>             Where to save the results"
        echo "-j <nproc>                number of parallel worker"
        echo "-n <nstruct>              number of predictions"
        echo "-s <inference_steps>      number of predictions"
        echo ""
        exit 1
}

while getopts ":p:l:L:o:j:n:s:" i; do
        case "${i}" in

        p)
                protein_pdb=$OPTARG
        ;;
        l)
                ligand=$OPTARG
        ;;
        L)
                ligand_name=$OPTARG
        ;;
        o)
                save_dir=$OPTARG
        ;;
        j)
                nproc=$OPTARG
        ;;
        n)
                nstruct=$OPTARG
        ;;
        s)
                inference_steps=$OPTARG
        ;;
        *)
                echo Unknown argument!
                usage
        ;;

        esac
done

if [[ ! -f "$protein_pdb" && "$ligand" == "" ]];then
  usage
fi

pdb_basename=$(basename $protein_pdb)
instance=${pdb_basename%.pdb}

#protein_pdb=$(readlink -f $protein_pdb)


if [[ -f "$ligand" ]];then
  ligand=$(readlink -f $ligand)
else
  ligand="\"${ligand}\""
fi

if [[ "$save_dir" == "" ]];then
  save_dir=$PWD/save
fi

if [[ ! -p "$save_dir" ]];then
  mkdir -p $save_dir
fi

if [[ "$nstruct" == "" ]];then
  nstruct=20
fi

if [[ "$nproc" == "" ]];then
  nproc=1
fi

if [[ "$inference_steps" == "" ]];then
  inference_steps=100
fi

if [[ "$ligand_name" == "" ]];then
  ligand_name='UNK'
fi


# prepare esm2 embedding
cmd="python ${DIFF_DOCK_HOME}/datasets/esm_embedding_preparation.py \
  --protein_path ${protein_pdb}\
  --out_file ${save_dir}/${instance}_prepared_for_esm.fasta"

echo "$cmd"
eval "$cmd"

# run esm2 embedding
cmd="python ${ESM_HOME}/scripts/extract.py \
  ${ESM_PRETRAINED}/esm2_t33_650M_UR50D.pt \
  ${save_dir}/${instance}_prepared_for_esm.fasta \
  ${save_dir}/${instance}_esm2_output \
  --repr_layers 33 \
  --include per_tok"

echo "$cmd"
eval "$cmd"

# run diffdock
cmd="python ${DIFF_DOCK_HOME}/inference.py \
  --protein_path ${protein_pdb} \
  --ligand ${ligand} \
  --out_dir ${save_dir}/${instance} \
  --inference_steps ${inference_steps} \
  --samples_per_complex ${nstruct} \
  --batch_size 10 \
  --num_workers ${nproc} \
  --esm_embeddings_path ${save_dir}/${instance}_esm2_output \
  --save_visualisation"

echo "$cmd"
eval "$cmd"

run_name=$(date +'%Y-%m-%d_%H%M%S')

# relocate run name as ${save_dir}/${instance}/${run_name}
mv ${save_dir}/${instance}/index* ${save_dir}/${instance}/${run_name}


# run visualization w/ PyMOL
cmd="python ${DIFF_DOCK_HOME}/analysis/analysis_pymol.py \
  --protein_pdb ${protein_pdb} \
  --ligand_name ${ligand_name} \
  --ligand_dir ${save_dir}/${instance}/${run_name} \
  --rank 1 ${nstruct} \
  --suffix diffdock_${run_name}_${nstruct}_${inference_steps}"

echo "$cmd"
eval "$cmd"