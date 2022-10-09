#python prep.py 
python datasets/esm_embedding_preparation.py \
            --protein_ligand_csv /home/tmp/input_protein_ligand.csv         \
            --out_file /home/data/prepared_for_esm.fasta 
python /app/DiffDock/esm/scripts/extract.py esm2_t33_650M_UR50D \
                                            ./data/prepared_for_esm.fasta \
                                            ./data/esm2_output \
                                            --repr_layers 33 \
                                            --include per_tok
python -m inference \
            --protein_ligand_csv ./tmp/input_protein_ligand.csv \
            --out_dir /home/results/user_predictions_7nei_5PET \
            --inference_steps 100 \
            --samples_per_complex 40 \
            --batch_size 10