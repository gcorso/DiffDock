# prep ligand and protein for inference by generating ESM embedding - has to run every cycle
python datasets/esm_embedding_preparation.py --protein_ligand_csv data/protein_ligand_example_csv.csv --out_file data/prepared_for_esm.fasta
HOME=esm/model_weights python esm/scripts/extract.py esm2_t33_650M_UR50D data/prepared_for_esm.fasta data/esm2_output --repr_layers 33 --include per_tok

# run simplified inference
python -m inference --protein_path /home/ubuntu/test/test.pdb --ligand /home/ubuntu/test/test/test.sdf --out_dir home/ubuntu --inference_steps 20 --samples_per_complex 40 --batch_size 10 --actual_steps 18 --no_final_step_noise