for i in $(seq 0 15); do
 	python baseline_tankbind_runtime.py --parallel_id $i --parallel_tot 16 --prank_path /data/rsg/nlp/hstark/TankBind/packages/p2rank_2.3/prank --data_dir /data/rsg/nlp/hstark/ligbind/data/PDBBind_processed --split_path /data/rsg/nlp/hstark/ligbind/data/splits/timesplit_test --results_path /data/rsg/nlp/hstark/ligbind/results/tankbind_16_worker_runtime --device cpu --skip_p2rank --num_workers 1 --skip_multiple_pocket_outputs &
done
wait

