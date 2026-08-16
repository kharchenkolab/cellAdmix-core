"""Vote-threshold sweep over N seed corrections (pancreas): write
consensus-removal deltas for a range of vote thresholds, for the R evaluator.
"""
import json
import gzip
import os
import numpy as np
import pyarrow.parquet as pq

OUT = 'examples/xenium_pancreas_membrane_377_full/out'
RESULTS = 'analysis/cleanup_benchmark/results'
SEEDS = list(range(1, 11))
THRESHOLDS = [1, 2, 3, 5, 7, 10]

run1 = f'{OUT}/runs/bench_seed1_ls_nmf'
mol = pq.read_table(f'{run1}/molecules.parquet',
    columns=['obs_id', 'gene_idx', 'cell_idx']).to_pandas()
genes = json.load(open(f'{run1}/run.json'))['genes']
cells = pq.read_table(f'{run1}/cells.parquet',
    columns=['cell_idx', 'cell_id']).to_pandas()
cell_id_by_idx = dict(zip(cells['cell_idx'], cells['cell_id'].astype(str)))
obs = mol['obs_id'].to_numpy()
max_obs = int(obs.max()) + 1

for variant in ['ls_nmf', 'invsqrt_kl']:
    for method in ['membrane', 'bridge']:
        votes = np.zeros(max_obs, dtype=np.int16)
        n_seeds = 0
        for s in SEEDS:
            cdir = (f'{OUT}/runs/bench_seed{s}_{variant}/corrected/'
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
        print(f'{variant}/{method}: {n_seeds} seeds pooled', flush=True)
        for thr in THRESHOLDS:
            if thr > n_seeds:
                continue
            rm = votes[obs] >= thr
            sub = mol.loc[rm]
            agg = sub.groupby(['gene_idx', 'cell_idx']).size().reset_index(name='n')
            path = f'{RESULTS}/pancreas_vote{thr}of{n_seeds}_{variant}_{method}.csv.gz'
            with gzip.open(path, 'wt') as f:
                f.write('gene,cell_id,n\n')
                for gi, ci, n in agg.itertuples(index=False):
                    f.write(f'{genes[gi]},{cell_id_by_idx[ci]},{n}\n')
            print(f'  vote>={thr}: removes {int(rm.sum())} -> {path}', flush=True)
print('SWEEP DELTAS DONE')
