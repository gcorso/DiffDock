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
        echo "-p <protein_pdb>  Protein in PDB format."
        echo "-l <ligand>       Ligand in SMILE string or in sdf/mol2 format"
        echo "-o <save_dir>     Where to save the results"
        echo "-j <nproc>     Where to save the results"
        echo ""
        exit 1
}

while getopts ":p:l:o:j:" i; do
        case "${i}" in

        p)
                protein_pdb=$OPTARG
        ;;
        l)
                ligand=$OPTARG
        ;;
        o)
                save_dir=$OPTARG
        ;;
        j)
                nproc=$OPTARG
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

protein_pdb=$(readlink -f $protein_pdb)


if [[ -f "$ligand" ]];then
  ligand=$(readlink -f $ligand)
fi

if [[ "$save_dir" == "" ]];then
  save_dir=$PWD/save
fi

if [[ ! -p "$save_dir" ]];then
  mkdir -p $save_dir
fi

if [[ "$nproc" == "" ]];then
  nproc=1
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

cmd="python ${DIFF_DOCK_HOME}/inference.py \
  --protein_path ${protein_pdb} \
  --out_dir ${save_dir}/user_predictions_small \
  --inference_steps 20 \
  --samples_per_complex 40 \
  --batch_size 10 \
  --num_workers ${nproc} \
  --esm_embeddings_path ${save_dir}/${instance}_esm2_output "

echo "$cmd"
eval "$cmd"

# evaluate
cmd="python ${DIFF_DOCK_HOME}evaluate_files.py \
  --results_path ${save_dir}/user_predictions_small \
  --file_to_exclude rank1.sdf
	"


cmd="python  ${DIFF_DOCK_HOME}/evaluate.py \
  --model_dir workdir/paper_score_model \
  --ckpt best_ema_inference_epoch_model.pt \
  --confidence_ckpt best_model_epoch75.pt \
  --confidence_model_dir workdir/paper_confidence_model \
  --run_name DiffDockInference \
  --inference_steps 20 \
  --split_path data/splits/timesplit_test \
  --samples_per_complex 40 --batch_size 10 \
  --save_visualisation
"