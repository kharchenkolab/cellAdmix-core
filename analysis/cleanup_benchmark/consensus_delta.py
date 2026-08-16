"""Majority-vote molecule-action consensus: for each dataset/variant/method,
vote per molecule across the three seeds' corrections and write the
consensus-removed molecules aggregated to (gene, cell_id) counts as CSV.gz,
for evaluation in R (which lacks a parquet reader).
"""
import json
import gzip
import os
import sys
import numpy as np
import pyarrow.parquet as pq

CONFIGS = {
    'pancreas': ('examples/xenium_pancreas_membrane_377_full/out', ['membrane', 'bridge']),
    'breast_crop': ('examples/xenium_breast_membrane_5k_full/out_medium_crop', ['membrane', 'bridge']),
    'nsclc': ('examples/cosmx_nsclc_giotto/out', ['bridge']),
}
OUT_DIR = 'analysis/cleanup_benchmark/results'

for dataset, (out, methods) in CONFIGS.items():
    run1 = f'{out}/runs/bench_seed1_ls_nmf'
    mol = pq.read_table(f'{run1}/molecules.parquet',
        columns=['obs_id', 'gene_idx', 'cell_idx']).to_pandas()
    genes = json.load(open(f'{run1}/run.json'))['genes']
    cells = pq.read_table(f'{run1}/cells.parquet',
        columns=['cell_idx', 'cell_id']).to_pandas()
    cell_id_by_idx = dict(zip(cells['cell_idx'], cells['cell_id'].astype(str)))
    obs = mol['obs_id'].to_numpy()
    max_obs = int(obs.max()) + 1
    for variant in ['ls_nmf', 'invsqrt_kl']:
        for method in methods:
            votes = np.zeros(max_obs, dtype=np.int8)
            n_seeds = 0
            for s in [1, 2, 3]:
                cdir = (f'{out}/runs/bench_seed{s}_{variant}/corrected/'
                        f'cmp_{method}_{variant}_s{s}')
                if not os.path.isdir(cdir):
                    print(f'missing: {cdir}', flush=True)
                    continue
                kept = pq.read_table(f'{cdir}/molecules.parquet',
                    columns=['obs_id'])['obs_id'].to_numpy()
                removed = np.ones(max_obs, dtype=bool)
                removed[kept] = False
                votes += removed
                n_seeds += 1
            if n_seeds < 2:
                continue
            for thr, tag in [(2, 'molcons'), (1, 'moluni')]:
                rm = votes[obs] >= thr
                sub = mol.loc[rm]
                agg = sub.groupby(['gene_idx', 'cell_idx']).size().reset_index(name='n')
                path = f'{OUT_DIR}/{dataset}_{tag}_{variant}_{method}.csv.gz'
                with gzip.open(path, 'wt') as f:
                    f.write('gene,cell_id,n\n')
                    for gi, ci, n in agg.itertuples(index=False):
                        f.write(f'{genes[gi]},{cell_id_by_idx[ci]},{n}\n')
                print(f'{dataset} {variant}/{method} vote>={thr}: removes {int(rm.sum())} '
                      f'molecules -> {path}', flush=True)
print('DELTAS DONE')
