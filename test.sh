# run within container
# install esm submodule 
#TODO #66 should not be necessary with latest update to dockerfile
cd /home/ubuntu/diffdock
rm -r esm
git clone https://github.com/facebookresearch/esm.git
cd esm
pip install -e .
cd ..

# prep ligand and protein for inference by generating ESM embedding - has to run every cycle
python datasets/esm_embedding_preparation.py --protein_path test/test.pdb --out_file data/prepared_for_esm.fasta
HOME=esm/model_weights python esm/scripts/extract.py esm2_t33_650M_UR50D data/prepared_for_esm.fasta data/esm2_output --repr_layers 33 --include per_tok

# run simplified inference
python -m inference --protein_path test/test.pdb --ligand test/test.sdf --out_dir /outputs --inference_steps 20 --samples_per_complex 40 --batch_size 10 --actual_steps 18 --no_final_step_noise