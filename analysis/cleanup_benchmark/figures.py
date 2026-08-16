"""Render the figures for docs/benchmarks.md from benchmark result files."""
import json
import os
import re
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = 'analysis/cleanup_benchmark/results'
OUT = 'docs/figures'
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'figure.dpi': 150})

BIN_ORDER = ['0', '1', '2', '3+']
short = lambda s: (s.replace(' epithelial', ' ep.').replace('Ductal/tumor', 'Ductal/tum.')
                    .replace('Fibroblast / CAF', 'Fibro/CAF').replace('Mural / pericyte', 'Mural'))

# ---------- Figure 1: the neighbor benchmark ----------
bins = pd.read_csv(f'{R}/pancreas_fig_bins.csv')
bins['pair'] = bins['source'] + ' -> ' + bins['target']
pairs_df = pd.read_csv(f'{R}/pancreas_pairs.csv')
bare = pairs_df[pairs_df['arm'] == 'ls_nmf/membrane/s1'].copy()

show_pairs = ['Exocrine epithelial -> Endothelial', 'Ductal/tumor epithelial -> Immune',
              'Endocrine -> Endothelial', 'Fibroblast / CAF -> Immune']
fig = plt.figure(figsize=(11, 6.2))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.5], hspace=0.45)
for i, pair in enumerate(show_pairs):
    ax = fig.add_subplot(gs[0, i])
    d = bins[bins['pair'] == pair].set_index('bin').reindex(BIN_ORDER)
    x = np.arange(4)
    ax.plot(x, d['rate_after'] * 1e3, 'o-', color='#2980b9', label='after cleanup')
    ax.plot(x, d['rate_before'] * 1e3, 's--', color='#c0392b', mfc='none',
        label='before cleanup')
    ax.axhline(d['rate_before'].iloc[0] * 1e3, color='grey', lw=0.7, ls=':')
    ax.set_xticks(x, BIN_ORDER)
    ax.set_xlabel('source-type neighbors')
    ax.set_title(short(pair), fontsize=8.5)
    if i == 0:
        ax.set_ylabel('source-marker rate\n(per 1000 molecules)')
        ax.legend(frameon=False, fontsize=8)

ax = fig.add_subplot(gs[1, :])
b = bare.dropna(subset=['power_strict']).sort_values('power_strict', ascending=False)
x = np.arange(len(b))
bars = ax.bar(x, b['power_strict'], color='#34495e', alpha=0.85)
ax2 = ax.twinx()
ax2.plot(x, b['excess_strict'], 'o', color='#e67e22', ms=4, alpha=0.8)
ax2.set_yscale('log')
ax2.set_ylabel('estimated leaked molecules', color='#e67e22')
ax2.tick_params(axis='y', colors='#e67e22')
ax.set_xticks(x, [short(p) for p in b['pair']], rotation=60, ha='right', fontsize=6.5)
ax.set_ylabel('estimated sensitivity (strict tier)')
ax.set_ylim(0, 1.05)
ax.set_title('all detected pairs under a standard single-fit cleanup (bare ls-NMF, membrane scoring)',
    fontsize=9)
fig.savefig(f'{OUT}/benchmark_fig1.png', bbox_inches='tight')
print('fig1 done')

# ---------- Figure 2: stochasticity ----------
ov = json.load(open(f'{R}/pancreas_removal_overlap.json'))
fig, axes = plt.subplots(1, 3, figsize=(11, 3.1),
    gridspec_kw={'width_ratios': [1, 1, 1.4]})
ax = axes[0]
x = np.arange(1, 11)
ax.bar(x - 0.2, np.array(ov['invsqrt_kl']['n_removed']) / 1e6, width=0.4,
    color='#8e44ad', label='invsqrt KL-NMF')
ax.bar(x + 0.2, np.array(ov['ls_nmf']['n_removed']) / 1e6, width=0.4,
    color='#16a085', label='ls-NMF')
ax.set_xlabel('random seed'); ax.set_ylabel('molecules removed (millions)')
ax.set_xticks(range(1, 11))
ax.legend(frameon=False, fontsize=8)
ax.set_title('total removal by seed')

ax = axes[1]
J = np.array(ov['invsqrt_kl']['jaccard'])
im = ax.imshow(J, vmin=0, vmax=1, cmap='viridis')
ax.set_xticks(range(10), range(1, 11), fontsize=7)
ax.set_yticks(range(10), range(1, 11), fontsize=7)
ax.set_xlabel('seed'); ax.set_ylabel('seed')
ax.set_title('removal-set overlap (Jaccard),\ninvsqrt KL-NMF')
plt.colorbar(im, ax=ax, fraction=0.045)

ax = axes[2]
pp = pairs_df[(pairs_df['variant'] == 'invsqrt_kl') & (pairs_df['method'] == 'membrane')
              & pairs_df['arm'].str.match(r'invsqrt_kl/membrane/s\d')]
piv = pp.pivot_table(index='pair', columns='seed', values='power_strict')
piv = piv.dropna().sort_values(1)
im = ax.imshow(piv.values, aspect='auto', cmap='RdYlBu', vmin=0, vmax=1)
ax.set_xticks(range(piv.shape[1]), [f's{int(c)}' for c in piv.columns])
ax.set_yticks(range(len(piv)), [short(p) for p in piv.index], fontsize=5)
ax.set_title('per-pair sensitivity by seed,\ninvsqrt KL-NMF / membrane')
plt.colorbar(im, ax=ax, fraction=0.03, label='sensitivity')
fig.tight_layout()
fig.savefig(f'{OUT}/benchmark_fig2.png', bbox_inches='tight')
print('fig2 done')

# ---------- Figure 3: vote-threshold ROC ----------
sweep_log = sys.argv[1] if len(sys.argv) > 1 else None
if sweep_log and os.path.exists(sweep_log):
    rows = []
    pat = re.compile(r'VOTE(\d+)OF(\d+) (\S+) (\S+)/(\S+): sens_strict=([\d.]+) '
                     r'sens_broad=([\d.]+) fpr_native=([\d.]+) spec_worst_pair=([\d.]+)')
    for line in open(sweep_log):
        m = pat.search(line)
        if m:
            rows.append(dict(thr=int(m[1]), n=int(m[2]), dataset=m[3], variant=m[4],
                method=m[5], sens=float(m[6]), fpr=float(m[7]), worst=float(m[8])))
    sw = pd.DataFrame(rows)
    sw.to_csv(f'{R}/vote_sweep_all.csv', index=False)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    colors = {'pancreas': '#c0392b', 'breast_crop': '#2980b9', 'nsclc': '#27ae60'}
    styles = {'membrane': '-', 'bridge': '--'}
    ax = axes[0]
    for (dsname, variant, method), g in sw[sw['variant'] == 'invsqrt_kl'].groupby(
            ['dataset', 'variant', 'method']):
        g = g.sort_values('thr')
        ax.plot(g['fpr'] * 100, g['sens'], styles[method], marker='o', ms=3,
            color=colors[dsname], label=f'{dsname} {method}')
        for _, r in g.iterrows():
            if r['thr'] in (1, 3, 5, 10):
                ax.annotate(f"≥{int(r['thr'])}", (r['fpr'] * 100, r['sens']),
                    fontsize=6, xytext=(3, 2), textcoords='offset points')
    ax.set_xlabel('estimated FPR, native stratum (%)')
    ax.set_ylabel('estimated sensitivity (strict tier)')
    ax.set_title('invsqrt KL-NMF: vote-threshold ROC')
    ax.legend(frameon=False, fontsize=7)
    ax = axes[1]
    for (dsname, variant, method), g in sw[sw['variant'] == 'ls_nmf'].groupby(
            ['dataset', 'variant', 'method']):
        g = g.sort_values('thr')
        ax.plot(g['fpr'] * 100, g['sens'], styles[method], marker='o', ms=3,
            color=colors[dsname], label=f'{dsname} {method}')
    ax.set_xlabel('estimated FPR, native stratum (%)')
    ax.set_title('ls-NMF: vote-threshold ROC')
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(f'{OUT}/benchmark_fig3.png', bbox_inches='tight')
    print('fig3 done')

# ---------- Figure 4: cross-dataset summary ----------
if os.path.exists(f'{R}/vote_sweep_all.csv'):
    sw = pd.read_csv(f'{R}/vote_sweep_all.csv')
    cards = []
    for dsname in ['pancreas', 'breast_crop', 'nsclc']:
        pares = pd.read_csv(f'{R}/{dsname}_pairs.csv')
        pares['cell'] = pares['variant'] + '/' + pares['method']
        for (variant, method), g in pares.groupby(['variant', 'method']):
            if variant not in ('ls_nmf', 'invsqrt_kl'):
                continue
            for kind, pat_ in [('bare', r'/s\d$'), ('rule-consensus', r'/cons_s\d$')]:
                gg = g[g['arm'].str.contains(pat_, regex=True)]
                for seed, gs in gg.groupby('seed'):
                    ok = gs.dropna(subset=['power_strict'])
                    cards.append(dict(dataset=dsname, variant=variant, method=method,
                        arm=kind, seed=seed,
                        sens=(ok['power_strict'] * ok['excess_strict']).sum() /
                             ok['excess_strict'].sum()))
    card = pd.DataFrame(cards)
    vote = sw[sw['thr'] == 3].rename(columns={'sens': 'sens3'})
    fig, ax = plt.subplots(figsize=(9.5, 3.6))
    cells = [(d, m) for d in ['pancreas', 'breast_crop', 'nsclc']
             for m in (['membrane', 'bridge'] if d != 'nsclc' else ['bridge'])]
    xbase = np.arange(len(cells)) * 2.6
    offs = {'ls_nmf': -0.5, 'invsqrt_kl': 0.5}
    for vi, variant in enumerate(['ls_nmf', 'invsqrt_kl']):
        for ci, (d, m) in enumerate(cells):
            x0 = xbase[ci] + offs[variant]
            b = card[(card['dataset'] == d) & (card['method'] == m) &
                     (card['variant'] == variant) & (card['arm'] == 'bare')]
            ax.scatter(np.full(len(b), x0 - 0.18), b['sens'], s=14, c='#95a5a6',
                label='single fits (seeds)' if vi == 0 and ci == 0 else None)
            rc = card[(card['dataset'] == d) & (card['method'] == m) &
                      (card['variant'] == variant) & (card['arm'] == 'rule-consensus')]
            ax.scatter(np.full(len(rc), x0), rc['sens'], s=14, c='#e67e22', marker='s',
                label='rule consensus (seeds)' if vi == 0 and ci == 0 else None)
            v = vote[(vote['dataset'] == d) & (vote['method'] == m) &
                     (vote['variant'] == variant)]
            if len(v):
                ax.scatter([x0 + 0.18], v['sens3'], s=70, c='#c0392b', marker='*',
                    label='molecule vote (3 of 10)' if vi == 0 and ci == 0 else None)
            ax.annotate('ls' if variant == 'ls_nmf' else 'inv',
                (x0, -0.06), ha='center', fontsize=7, annotation_clip=False)
    ax.set_xticks(xbase, [f"{d.replace('_crop','')}\n{m}" for d, m in cells], fontsize=8)
    ax.set_ylabel('estimated sensitivity (strict tier)')
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False, fontsize=8, loc='lower left')
    ax.set_title('cleanup across datasets, scoring methods, variants, and correction strategies')
    fig.tight_layout()
    fig.savefig(f'{OUT}/benchmark_fig4.png', bbox_inches='tight')
    print('fig4 done')

print('FIGURES DONE')
