import os
from argparse import ArgumentParser

import pandas as pd
import plotly.express as px
import numpy as np
import scipy

from utils.utils import read_strings_from_txt

parser = ArgumentParser()


parser.add_argument('--data_dir', type=str, default='data/PDBBind_processed', help='')
parser.add_argument('--results_path', type=str, default='inference_out_dir_not_specified/TEST_top40_epoch75_FILTER_restart_cacheNewRestart_big_ema_ESM2emb_tr34_WITH_fixedSamples28_id1_FILTERFROM_temp_restart_ema_ESM2emb_tr34', help='')
parser.add_argument('--gnina_results_path', type=str, default='results/gnina_rosetta13', help='')
parser.add_argument('--smina_results_path', type=str, default='results/smina_rosetta13', help='')
parser.add_argument('--glide_results_path', type=str, default='results/glide', help='')
parser.add_argument('--qvinaw_results_path', type=str, default='results/qvinaw', help='')
parser.add_argument('--tankbind_results_path', type=str, default='results/tankbind_top5', help='')
parser.add_argument('--equibind_results_path', type=str, default='results/equibind_paper', help='')
parser.add_argument('--no_rec_overlap', action='store_true', default=False, help='')
args = parser.parse_args()



min_cross_distances = np.load(f'{args.results_path}/min_cross_distances.npy')
#min_self_distances = np.load(f'{args.results_path}/min_self_distances.npy')
base_min_cross_distances = np.load(f'{args.results_path}/base_min_cross_distances.npy')
rmsds = np.load(f'{args.results_path}/rmsds.npy')
centroid_distances = np.load(f'{args.results_path}/centroid_distances.npy')
confidences = np.load(f'{args.results_path}/confidences.npy')
#complex_names = np.load(f'{args.results_path}/complex_names.npy')
complex_names = read_strings_from_txt('data/splits/timesplit_test')
if args.no_rec_overlap:
    names_no_rec_overlap = read_strings_from_txt(f'data/splits/timesplit_test_no_rec_overlap')
    without_rec_overlap_list = []
    for name in complex_names:
        if name in names_no_rec_overlap:
            without_rec_overlap_list.append(1)
        else:
            without_rec_overlap_list.append(0)
    without_rec_overlap = np.array(without_rec_overlap_list, dtype=bool)
    rmsds = np.array(rmsds)[without_rec_overlap]
    #min_self_distances = np.array(min_self_distances)[without_rec_overlap]
    centroid_distances = np.array(centroid_distances)[without_rec_overlap]
    confidences = np.array(confidences)[without_rec_overlap]
    min_cross_distances = np.array(min_cross_distances)[without_rec_overlap]
    base_min_cross_distances = np.array(base_min_cross_distances)[without_rec_overlap]
    complex_names = names_no_rec_overlap




N = rmsds.shape[1]
performance_metrics = {
    'steric_clash_fraction': (100 * (min_cross_distances < 0.4).sum() / len(min_cross_distances) / N).__round__(2),
    'mean_rmsd': rmsds.mean(),
    'rmsds_below_2': (100 * (rmsds < 2).sum() / len(rmsds) / N),
    'rmsds_below_5': (100 * (rmsds < 5).sum() / len(rmsds) / N),
    'rmsds_percentile_25': np.percentile(rmsds, 25).round(2),
    'rmsds_percentile_50': np.percentile(rmsds, 50).round(2),
    'rmsds_percentile_75': np.percentile(rmsds, 75).round(2),

    'mean_centroid': centroid_distances.mean().__round__(2),
    'centroid_below_2': (100 * (centroid_distances < 2).sum() / len(centroid_distances) / N).__round__(2),
    'centroid_below_5': (100 * (centroid_distances < 5).sum() / len(centroid_distances) / N).__round__(2),
    'centroid_percentile_25': np.percentile(centroid_distances, 25).round(2),
    'centroid_percentile_50': np.percentile(centroid_distances, 50).round(2),
    'centroid_percentile_75': np.percentile(centroid_distances, 75).round(2),
}

if N >= 5:
    top5_rmsds = np.min(rmsds[:, :5], axis=1)
    top5_centroid_distances = centroid_distances[np.arange(rmsds.shape[0])[:, None], np.argsort(rmsds[:, :5], axis=1)][ :, 0]
    top5_min_cross_distances = min_cross_distances[ np.arange(rmsds.shape[0])[:, None], np.argsort(rmsds[:, :5], axis=1)][:, 0]
    performance_metrics.update({
        'top5_steric_clash_fraction': (100 * (top5_min_cross_distances < 0.4).sum() / len(top5_min_cross_distances)).__round__(2),
        'top5_rmsds_below_2': (100 * (top5_rmsds < 2).sum() / len(top5_rmsds)).__round__(2),
        'top5_rmsds_below_5': (100 * (top5_rmsds < 5).sum() / len(top5_rmsds)).__round__(2),
        'top5_rmsds_percentile_25': np.percentile(top5_rmsds, 25).round(2),
        'top5_rmsds_percentile_50': np.percentile(top5_rmsds, 50).round(2),
        'top5_rmsds_percentile_75': np.percentile(top5_rmsds, 75).round(2),

        'top5_centroid_below_2': (100 * (top5_centroid_distances < 2).sum() / len(top5_centroid_distances)).__round__(2),
        'top5_centroid_below_5': (100 * (top5_centroid_distances < 5).sum() / len(top5_centroid_distances)).__round__(2),
        'top5_centroid_percentile_25': np.percentile(top5_centroid_distances, 25).round(2),
        'top5_centroid_percentile_50': np.percentile(top5_centroid_distances, 50).round(2),
        'top5_centroid_percentile_75': np.percentile(top5_centroid_distances, 75).round(2),
    })

if N >= 10:
    top10_rmsds = np.min(rmsds[:, :10], axis=1)
    top10_centroid_distances = centroid_distances[np.arange(rmsds.shape[0])[:, None], np.argsort(rmsds[:, :10], axis=1)][:, 0]
    top10_min_cross_distances = min_cross_distances[np.arange(rmsds.shape[0])[:, None], np.argsort(rmsds[:, :10], axis=1)][:, 0]
    performance_metrics.update({
        'top10_steric_clash_fraction': (100 * (top10_min_cross_distances < 0.4).sum() / len(top10_min_cross_distances)).__round__(2),
        'top10_rmsds_below_2': (100 * (top10_rmsds < 2).sum() / len(top10_rmsds)).__round__(2),
        'top10_rmsds_below_5': (100 * (top10_rmsds < 5).sum() / len(top10_rmsds)).__round__(2),
        'top10_rmsds_percentile_25': np.percentile(top10_rmsds, 25).round(2),
        'top10_rmsds_percentile_50': np.percentile(top10_rmsds, 50).round(2),
        'top10_rmsds_percentile_75': np.percentile(top10_rmsds, 75).round(2),

        'top10_centroid_below_2': (100 * (top10_centroid_distances < 2).sum() / len(top10_centroid_distances)).__round__(2),
        'top10_centroid_below_5': (100 * (top10_centroid_distances < 5).sum() / len(top10_centroid_distances)).__round__(2),
        'top10_centroid_percentile_25': np.percentile(top10_centroid_distances, 25).round(2),
        'top10_centroid_percentile_50': np.percentile(top10_centroid_distances, 50).round(2),
        'top10_centroid_percentile_75': np.percentile(top10_centroid_distances, 75).round(2),
    })


confidence_ordering = np.argsort(confidences,axis=1)[:,::-1]
filtered_rmsds = rmsds[np.arange(rmsds.shape[0])[:,None],confidence_ordering][:,0]
filtered_centroid_distances = centroid_distances[np.arange(rmsds.shape[0])[:,None],confidence_ordering][:,0]
filtered_min_cross_distances = min_cross_distances[np.arange(rmsds.shape[0])[:, None], confidence_ordering][:, 0]
performance_metrics.update({
    'filtered_steric_clash_fraction': (100 * (filtered_min_cross_distances < 0.4).sum() / len(filtered_min_cross_distances)).__round__(2),
    'filtered_rmsds_below_2': (100 * (filtered_rmsds < 2).sum() / len(filtered_rmsds)).__round__(2),
    'filtered_rmsds_below_5': (100 * (filtered_rmsds < 5).sum() / len(filtered_rmsds)).__round__(2),
    'filtered_rmsds_percentile_25': np.percentile(filtered_rmsds, 25).round(2),
    'filtered_rmsds_percentile_50': np.percentile(filtered_rmsds, 50).round(2),
    'filtered_rmsds_percentile_75': np.percentile(filtered_rmsds, 75).round(2),

    'filtered_centroid_below_2': (100 * (filtered_centroid_distances < 2).sum() / len(filtered_centroid_distances)).__round__(2),
    'filtered_centroid_below_5': (100 * (filtered_centroid_distances < 5).sum() / len(filtered_centroid_distances)).__round__(2),
    'filtered_centroid_percentile_25': np.percentile(filtered_centroid_distances, 25).round(2),
    'filtered_centroid_percentile_50': np.percentile(filtered_centroid_distances, 50).round(2),
    'filtered_centroid_percentile_75': np.percentile(filtered_centroid_distances, 75).round(2),
})

if N >= 5:
    top5_filtered_rmsds = np.min(rmsds[np.arange(rmsds.shape[0])[:,None],confidence_ordering][:,:5], axis=1)
    top5_filtered_centroid_distances = centroid_distances[np.arange(rmsds.shape[0])[:,None],confidence_ordering][:,:5][ np.arange(rmsds.shape[0])[:, None], np.argsort(rmsds[np.arange(rmsds.shape[0])[:,None],confidence_ordering][:, :5], axis=1)][:, 0]
    top5_filtered_min_cross_distances = min_cross_distances[np.arange(rmsds.shape[0])[:, None], confidence_ordering][:, :5][ np.arange(rmsds.shape[0])[:, None], np.argsort(rmsds[np.arange(rmsds.shape[0])[:,None],confidence_ordering][:, :5], axis=1)][:, 0]
    performance_metrics.update({
        'top5_filtered_steric_clash_fraction': (100 * (top5_filtered_min_cross_distances < 0.4).sum() / len(top5_filtered_min_cross_distances)).__round__(2),
        'top5_filtered_rmsds_below_2': (100 * (top5_filtered_rmsds < 2).sum() / len(top5_filtered_rmsds)).__round__(2),
        'top5_filtered_rmsds_below_5': (100 * (top5_filtered_rmsds < 5).sum() / len(top5_filtered_rmsds)).__round__(2),
        'top5_filtered_rmsds_percentile_25': np.percentile(top5_filtered_rmsds, 25).round(2),
        'top5_filtered_rmsds_percentile_50': np.percentile(top5_filtered_rmsds, 50).round(2),
        'top5_filtered_rmsds_percentile_75': np.percentile(top5_filtered_rmsds, 75).round(2),

        'top5_filtered_centroid_below_2': (100 * (top5_filtered_centroid_distances < 2).sum() / len(top5_filtered_centroid_distances)).__round__(2),
        'top5_filtered_centroid_below_5': (100 * (top5_filtered_centroid_distances < 5).sum() / len(top5_filtered_centroid_distances)).__round__(2),
        'top5_filtered_centroid_percentile_25': np.percentile(top5_filtered_centroid_distances, 25).round(2),
        'top5_filtered_centroid_percentile_50': np.percentile(top5_filtered_centroid_distances, 50).round(2),
        'top5_filtered_centroid_percentile_75': np.percentile(top5_filtered_centroid_distances, 75).round(2),
    })
if N >= 10:
    top10_filtered_rmsds = np.min(rmsds[np.arange(rmsds.shape[0])[:,None],confidence_ordering][:,:10], axis=1)
    top10_filtered_centroid_distances = centroid_distances[np.arange(rmsds.shape[0])[:,None],confidence_ordering][:,:10][ np.arange(rmsds.shape[0])[:, None], np.argsort(rmsds[np.arange(rmsds.shape[0])[:,None],confidence_ordering][:, :10], axis=1)][:, 0]
    top10_filtered_min_cross_distances = min_cross_distances[np.arange(rmsds.shape[0])[:, None], confidence_ordering][:, :10][ np.arange(rmsds.shape[0])[:, None], np.argsort(rmsds[np.arange(rmsds.shape[0])[:,None],confidence_ordering][:, :10], axis=1)][:, 0]
    performance_metrics.update({
        'top10_filtered_steric_clash_fraction': (100 * (top10_filtered_min_cross_distances < 0.4).sum() / len(top10_filtered_min_cross_distances)).__round__(2),
        'top10_filtered_rmsds_below_2': (100 * (top10_filtered_rmsds < 2).sum() / len(top10_filtered_rmsds)).__round__(2),
        'top10_filtered_rmsds_below_5': (100 * (top10_filtered_rmsds < 5).sum() / len(top10_filtered_rmsds)).__round__(2),
        'top10_filtered_rmsds_percentile_25': np.percentile(top10_filtered_rmsds, 25).round(2),
        'top10_filtered_rmsds_percentile_50': np.percentile(top10_filtered_rmsds, 50).round(2),
        'top10_filtered_rmsds_percentile_75': np.percentile(top10_filtered_rmsds, 75).round(2),

        'top10_filtered_centroid_below_2': (100 * (top10_filtered_centroid_distances < 2).sum() / len(top10_filtered_centroid_distances)).__round__(2),
        'top10_filtered_centroid_below_5': (100 * (top10_filtered_centroid_distances < 5).sum() / len(top10_filtered_centroid_distances)).__round__(2),
        'top10_filtered_centroid_percentile_25': np.percentile(top10_filtered_centroid_distances, 25).round(2),
        'top10_filtered_centroid_percentile_50': np.percentile(top10_filtered_centroid_distances, 50).round(2),
        'top10_filtered_centroid_percentile_75': np.percentile(top10_filtered_centroid_distances, 75).round(2),
    })

reverse_confidence_ordering = np.argsort(confidences,axis=1)
reverse_filtered_rmsds = rmsds[np.arange(rmsds.shape[0])[:, None], reverse_confidence_ordering][:, 0]
reverse_filtered_centroid_distances = centroid_distances[np.arange(rmsds.shape[0])[:, None], reverse_confidence_ordering][:, 0]
reverse_filtered_min_cross_distances = min_cross_distances[np.arange(rmsds.shape[0])[:, None], reverse_confidence_ordering][:, 0]
performance_metrics.update({
    'reversefiltered_steric_clash_fraction': (100 * (reverse_filtered_min_cross_distances < 0.4).sum() / len(reverse_filtered_min_cross_distances)).__round__(2),
    'reversefiltered_rmsds_below_2': (100 * (reverse_filtered_rmsds < 2).sum() / len(reverse_filtered_rmsds)).__round__(2),
    'reversefiltered_rmsds_below_5': (100 * (reverse_filtered_rmsds < 5).sum() / len(reverse_filtered_rmsds)).__round__(2),
    'reversefiltered_rmsds_percentile_25': np.percentile(reverse_filtered_rmsds, 25).round(2),
    'reversefiltered_rmsds_percentile_50': np.percentile(reverse_filtered_rmsds, 50).round(2),
    'reversefiltered_rmsds_percentile_75': np.percentile(reverse_filtered_rmsds, 75).round(2),

    'reversefiltered_centroid_below_2': (100 * (reverse_filtered_centroid_distances < 2).sum() / len(reverse_filtered_centroid_distances)).__round__(2),
    'reversefiltered_centroid_below_5': (100 * (reverse_filtered_centroid_distances < 5).sum() / len(reverse_filtered_centroid_distances)).__round__(2),
    'reversefiltered_centroid_percentile_25': np.percentile(reverse_filtered_centroid_distances, 25).round(2),
    'reversefiltered_centroid_percentile_50': np.percentile(reverse_filtered_centroid_distances, 50).round(2),
    'reversefiltered_centroid_percentile_75': np.percentile(reverse_filtered_centroid_distances, 75).round(2),
})

if N >= 5:
    top5_reverse_filtered_rmsds = np.min(rmsds[np.arange(rmsds.shape[0])[:, None], reverse_confidence_ordering][:, :5], axis=1)
    top5_reverse_filtered_centroid_distances = np.min(centroid_distances[np.arange(rmsds.shape[0])[:, None], reverse_confidence_ordering][:, :5], axis=1)
    top5_reverse_filtered_min_cross_distances = np.max(min_cross_distances[np.arange(rmsds.shape[0])[:, None], reverse_confidence_ordering][:, :5], axis=1)
    performance_metrics.update({
        'top5_reverse_filtered_steric_clash_fraction': (100 * (top5_reverse_filtered_min_cross_distances < 0.4).sum() / len(top5_reverse_filtered_min_cross_distances)).__round__(2),
        'top5_reversefiltered_rmsds_below_2': (100 * (top5_reverse_filtered_rmsds < 2).sum() / len(top5_reverse_filtered_rmsds)).__round__(2),
        'top5_reversefiltered_rmsds_below_5': (100 * (top5_reverse_filtered_rmsds < 5).sum() / len(top5_reverse_filtered_rmsds)).__round__(2),
        'top5_reversefiltered_rmsds_percentile_25': np.percentile(top5_reverse_filtered_rmsds, 25).round(2),
        'top5_reversefiltered_rmsds_percentile_50': np.percentile(top5_reverse_filtered_rmsds, 50).round(2),
        'top5_reversefiltered_rmsds_percentile_75': np.percentile(top5_reverse_filtered_rmsds, 75).round(2),

        'top5_reversefiltered_centroid_below_2': (100 * (top5_reverse_filtered_centroid_distances < 2).sum() / len(top5_reverse_filtered_centroid_distances)).__round__(2),
        'top5_reversefiltered_centroid_below_5': (100 * (top5_reverse_filtered_centroid_distances < 5).sum() / len(top5_reverse_filtered_centroid_distances)).__round__(2),
        'top5_reversefiltered_centroid_percentile_25': np.percentile(top5_reverse_filtered_centroid_distances, 25).round(2),
        'top5_reversefiltered_centroid_percentile_50': np.percentile(top5_reverse_filtered_centroid_distances, 50).round(2),
        'top5_reversefiltered_centroid_percentile_75': np.percentile(top5_reverse_filtered_centroid_distances, 75).round(2),
    })

if N >= 10:
    top10_reverse_filtered_rmsds = np.min(rmsds[np.arange(rmsds.shape[0])[:, None], reverse_confidence_ordering][:, :10], axis=1)
    top10_reverse_filtered_centroid_distances = np.min(centroid_distances[np.arange(rmsds.shape[0])[:, None], reverse_confidence_ordering][:, :10], axis=1)
    top10_reverse_filtered_min_cross_distances = np.max(min_cross_distances[np.arange(rmsds.shape[0])[:, None], reverse_confidence_ordering][:, :10], axis=1)
    performance_metrics.update({
        'top10_reverse_filtered_steric_clash_fraction': (100 * (top10_reverse_filtered_min_cross_distances < 0.4).sum() / len(top10_reverse_filtered_min_cross_distances)).__round__(2),
        'top10_reversefiltered_rmsds_below_2': (100 * (top10_reverse_filtered_rmsds < 2).sum() / len(top10_reverse_filtered_rmsds)).__round__(2),
        'top10_reversefiltered_rmsds_below_5': (100 * (top10_reverse_filtered_rmsds < 5).sum() / len(top10_reverse_filtered_rmsds)).__round__(2),
        'top10_reversefiltered_rmsds_percentile_25': np.percentile(top10_reverse_filtered_rmsds, 25).round(2),
        'top10_reversefiltered_rmsds_percentile_50': np.percentile(top10_reverse_filtered_rmsds, 50).round(2),
        'top10_reversefiltered_rmsds_percentile_75': np.percentile(top10_reverse_filtered_rmsds, 75).round(2),

        'top10_reversefiltered_centroid_below_2': (100 * (top10_reverse_filtered_centroid_distances < 2).sum() / len(top10_reverse_filtered_centroid_distances)).__round__(2),
        'top10_reversefiltered_centroid_below_5': (100 * (top10_reverse_filtered_centroid_distances < 5).sum() / len(top10_reverse_filtered_centroid_distances)).__round__(2),
        'top10_reversefiltered_centroid_percentile_25': np.percentile(top10_reverse_filtered_centroid_distances, 25).round(2),
        'top10_reversefiltered_centroid_percentile_50': np.percentile(top10_reverse_filtered_centroid_distances, 50).round(2),
        'top10_reversefiltered_centroid_percentile_75': np.percentile(top10_reverse_filtered_centroid_distances, 75).round(2),
    })

filtered_confidences = confidences[np.arange(confidences.shape[0])[:,None],confidence_ordering][:,0]

confident_mask = filtered_confidences > 0
confident_rmsds = filtered_rmsds[confident_mask]
confident_centroid_distances = filtered_centroid_distances[confident_mask]
confident_min_cross_distances = filtered_min_cross_distances[confident_mask]

performance_metrics.update({
    'fraction_confident_predictions': (100 * len(confident_rmsds) / len(rmsds)).__round__(2),
    'confident_steric_clash_fraction': (100 * (confident_min_cross_distances < 0.4).sum() / len(confident_min_cross_distances)).__round__(2),
    'confident_rmsds_below_2': (100 * (confident_rmsds < 2).sum() / len(confident_rmsds)).__round__(2),
    'confident_rmsds_below_5': (100 * (confident_rmsds < 5).sum() / len(confident_rmsds)).__round__(2),
    'confident_rmsds_percentile_25': np.percentile(confident_rmsds, 25).round(2),
    'confident_rmsds_percentile_50': np.percentile(confident_rmsds, 50).round(2),
    'confident_rmsds_percentile_75': np.percentile(confident_rmsds, 75).round(2),

    'confident_centroid_below_2': (100 * (confident_centroid_distances < 2).sum() / len(confident_centroid_distances)).__round__(2),
    'confident_centroid_below_5': (100 * (confident_centroid_distances < 5).sum() / len(confident_centroid_distances)).__round__(2),
    'confident_centroid_percentile_25': np.percentile(confident_centroid_distances, 25).round(2),
    'confident_centroid_percentile_50': np.percentile(confident_centroid_distances, 50).round(2),
    'confident_centroid_percentile_75': np.percentile(confident_centroid_distances, 75).round(2),
})

for k in performance_metrics:
    print(k, performance_metrics[k])

fraction_dataset_rmsds_below_2 = []
perfect_calibration = []
no_calibration = []
for dataset_percentage in range(100):
    dataset_percentage += 1
    dataset_fraction = (dataset_percentage)/100
    num_samples = round(len(rmsds)*dataset_fraction)
    per_complex_confidence_ordering = np.argsort(filtered_confidences)[::-1]
    confident_complexes_rmsds = filtered_rmsds[per_complex_confidence_ordering][:num_samples]
    confident_complexes_centroid_distances = filtered_centroid_distances[per_complex_confidence_ordering][:num_samples]
    confident_complexes_min_cross_distances = filtered_min_cross_distances[per_complex_confidence_ordering][:num_samples]
    confident_complexes_metrics = {
        'fraction_confident_complexes_predictions': (100 * len(confident_complexes_rmsds) / len(rmsds)).__round__(2),
        'confident_complexes_steric_clash_fraction': (100 * (confident_complexes_min_cross_distances < 0.4).sum() / len(confident_complexes_min_cross_distances)).__round__(2),
        'confident_complexes_rmsds_below_2': (100 * (confident_complexes_rmsds < 2).sum() / len(confident_complexes_rmsds)).__round__(2),
        'confident_complexes_rmsds_below_5': (100 * (confident_complexes_rmsds < 5).sum() / len(confident_complexes_rmsds)).__round__(2),
        'confident_complexes_rmsds_percentile_25': np.percentile(confident_complexes_rmsds, 25).round(2),
        'confident_complexes_rmsds_percentile_50': np.percentile(confident_complexes_rmsds, 50).round(2),
        'confident_complexes_rmsds_percentile_75': np.percentile(confident_complexes_rmsds, 75).round(2),

        'confident_complexes_centroid_below_2': (100 * (confident_complexes_centroid_distances < 2).sum() / len(confident_complexes_centroid_distances)).__round__(2),
        'confident_complexes_centroid_below_5': (100 * (confident_complexes_centroid_distances < 5).sum() / len(confident_complexes_centroid_distances)).__round__(2),
        'confident_complexes_centroid_percentile_25': np.percentile(confident_complexes_centroid_distances, 25).round(2),
        'confident_complexes_centroid_percentile_50': np.percentile(confident_complexes_centroid_distances, 50).round(2),
        'confident_complexes_centroid_percentile_75': np.percentile(confident_complexes_centroid_distances, 75).round(2),
    }
    fraction_dataset_rmsds_below_2.append(confident_complexes_metrics['confident_complexes_rmsds_below_2'])
    perfect_calibration.append((100 * (np.sort(filtered_rmsds)[:num_samples] < 2).sum() / len(confident_complexes_rmsds)).__round__(2))
    no_calibration.append(performance_metrics['filtered_rmsds_below_2'])
    #print('percentage: ',dataset_percentage)
    #print(confident_complexes_metrics['confident_complexes_rmsds_below_2'])

print(scipy.stats.spearmanr(filtered_rmsds, filtered_confidences))
df = {'conf': filtered_confidences, 'rmsd': filtered_rmsds}
fig = px.scatter(df, x='rmsd',y='conf').update_layout(
    xaxis_title="Percentage of datapoints that may be abstained", yaxis_title="Percentage of predictions with RMSD < 2A"
)
fig.update_layout(margin={'l': 0, 'r': 0, 't': 20, 'b': 100}, plot_bgcolor='white',
                paper_bgcolor='white', legend_title_text='', legend_title_font_size=1,
                legend=dict(yanchor="bottom", y=0.1, xanchor="right", x=0.99, font=dict(size=17), ),
                )
fig.update_xaxes(showgrid=True, gridcolor='lightgrey',title_font=dict(size=19),mirror=True,ticks='outside',showline=True,)
fig.update_yaxes(showgrid=True, gridcolor='lightgrey',title_font=dict(size=19),mirror=True,ticks='outside',showline=True,)
fig.show()

df = {'Confidence Model': reversed(fraction_dataset_rmsds_below_2),'No Calibration': reversed(no_calibration),'Perfect Calibration': reversed(perfect_calibration),}
fig = px.line(df, y=list(df.keys())).update_layout(
    xaxis_title="Percentage of datapoints that may be abstained", yaxis_title="Percentage of predictions with RMSD < 2A"
)
fig.update_yaxes(range = [0,103])
fig.update_layout(margin={'l': 0, 'r': 0, 't': 20, 'b': 100}, plot_bgcolor='white',
                paper_bgcolor='white', legend_title_text='', legend_title_font_size=1,
                legend=dict(yanchor="bottom", y=0.1, xanchor="right", x=0.99, font=dict(size=17), ),
                )
fig.update_xaxes(showgrid=True, gridcolor='lightgrey',title_font=dict(size=19),mirror=True,ticks='outside',showline=True,)
fig.update_yaxes(showgrid=True, gridcolor='lightgrey',title_font=dict(size=19),mirror=True,ticks='outside',showline=True,)
fig.write_image('results/confidence_calibration.pdf')
fig.show()

def filter_by_names(method_names, method_array, names_to_keep):
    output_array = []
    output_names = []
    for method_name, array_element in zip(method_names,method_array):
        if method_name in names_to_keep:
            output_array.append(array_element)
            output_names.append(method_name)
    return np.array(output_array), np.array(output_names)

qvinaw_rmsds = np.load(os.path.join(args.qvinaw_results_path, 'rmsds.npy'))
qvinaw_names = np.load(os.path.join(args.qvinaw_results_path, 'names.npy'))
qvinaw_rmsds, qvinaw_names = filter_by_names(qvinaw_names, qvinaw_rmsds, complex_names)
qvinaw_rmsds = np.concatenate([qvinaw_rmsds, np.random.choice(qvinaw_rmsds, size=len(complex_names) - len(qvinaw_rmsds))])

glide_rmsds = np.load(os.path.join(args.glide_results_path, 'rmsds.npy'))
glide_names = np.load(os.path.join(args.glide_results_path, 'names.npy')).tolist()
glide_rmsds, glide_names = filter_by_names(glide_names, glide_rmsds, complex_names)
glide_rmsds = np.concatenate([glide_rmsds, np.random.choice(glide_rmsds, size=len(complex_names) - len(glide_rmsds))])

smina_rmsds = np.load(os.path.join(args.smina_results_path, 'rmsds.npy'))[:,0]
smina_names = np.load(os.path.join(args.smina_results_path, 'names.npy'))
smina_rmsds, smina_names = filter_by_names(smina_names, smina_rmsds, complex_names)
smina_rmsds = np.concatenate([smina_rmsds, np.random.choice(smina_rmsds, size=len(complex_names) - len(smina_rmsds))])

gnina_rmsds = np.load(os.path.join(args.gnina_results_path, 'rmsds.npy'))[:,0]
gnina_names = np.load(os.path.join(args.gnina_results_path, 'names.npy'))
gnina_rmsds, gnina_names = filter_by_names(gnina_names, gnina_rmsds, complex_names)
gnina_rmsds = np.concatenate([gnina_rmsds, np.random.choice(gnina_rmsds, size=len(complex_names) - len(gnina_rmsds))])

tankbind_rmsds = np.load(os.path.join(args.tankbind_results_path, 'rmsds.npy'))[:,0]
tankbind_names = np.load(os.path.join(args.tankbind_results_path, 'names.npy'))
tankbind_rmsds, tankbind_names = filter_by_names(tankbind_names, tankbind_rmsds, complex_names)

equibind_rmsds = np.load(os.path.join(args.equibind_results_path, 'rmsds.npy'))
equibind_names = np.load(os.path.join(args.equibind_results_path, 'names.npy'))
equibind_rmsds, equibind_names = filter_by_names(equibind_names, equibind_rmsds, complex_names)


df = {'DiffDock': filtered_rmsds, 'GLIDE': glide_rmsds, 'GNINA': gnina_rmsds, 'SMINA': smina_rmsds, 'QVinaW':qvinaw_rmsds, 'TANKBind': tankbind_rmsds, 'EquiBind': equibind_rmsds}
fig = px.ecdf(df, range_x=[0, 5], range_y=[0.001, 0.75],  width=600, height=400)
fig.add_vline(x=2, annotation_text='', annotation_font_size=20, annotation_position="top right",
              line_dash='dash', line_color='firebrick', annotation_font_color='firebrick')
fig.update_xaxes(title=f'RMSD (Å)')
fig.update_yaxes(title=f'Fraction with lower RMSD')
fig.update_layout(autosize=False, margin={'l': 65, 'r': 5, 't': 5, 'b': 60}, plot_bgcolor='white',
                  paper_bgcolor='white', legend_title_text='', legend_title_font_size=18,
                  legend=dict(yanchor="top", y=0.995, xanchor="left", x=0.02, font=dict(size=18, color='black'), ), )
fig.update_xaxes(showgrid=True, gridcolor='lightgrey',title_font=dict(size=23, color='black'),mirror=True,ticks='outside',showline=True, linewidth=1, linecolor='black', tickfont = dict(size = 18, color='black'))
fig.update_yaxes(showgrid=True, gridcolor='lightgrey',title_font=dict(size=23, color='black'),mirror=True,ticks='outside',showline=True, linewidth=1, linecolor='black', tickfont = dict(size = 18, color='black'))
fig.update_traces(line=dict(width=3))
fig.write_image('results/rmsds_nooverlap.pdf')
fig.show()