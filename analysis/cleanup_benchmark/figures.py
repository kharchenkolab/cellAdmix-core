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
fig = plt.figure(figsize=(11, 11.2))
gs = fig.add_gridspec(3, 4, height_ratios=[1, 1.6, 1.45], hspace=0.75)
for i, pair in enumerate(show_pairs):
    ax = fig.add_subplot(gs[0, i])
    if i == 0:
        ax.annotate('a', (-0.32, 1.12), xycoords='axes fraction',
            fontsize=13, fontweight='bold')
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
        ax.set_ylabel('admixture-marker rate $\\hat{\\rho}_B$\n(per 1000 molecules)')
        ax.legend(frameon=False, fontsize=8)

# panel b: extrapolated per-pair admixture rates r as a source x target map
rates = pd.read_csv(f'{R}/pancreas_pair_rates.csv')
rate_types = sorted(set(rates['source']) | set(rates['target']))
mat = np.full((len(rate_types), len(rate_types)), np.nan)
for _, row in rates.iterrows():
    mat[rate_types.index(row['source']), rate_types.index(row['target'])] = row['rate'] * 100
ax = fig.add_subplot(gs[1, 1:3])
ax.annotate('b', (-0.75, 1.05), xycoords='axes fraction',
    fontsize=13, fontweight='bold')
masked = np.ma.masked_invalid(mat)
im = ax.imshow(masked, cmap='Reds', vmin=0)
for i in range(len(rate_types)):
    for j in range(len(rate_types)):
        if not np.isnan(mat[i, j]):
            im_color = 'white' if mat[i, j] > 0.6 * np.nanmax(mat) else 'black'
            ax.text(j, i, f'{mat[i, j]:.1f}', ha='center', va='center',
                fontsize=7.5, color=im_color)
ax.set_xticks(range(len(rate_types)), [short(t) for t in rate_types],
    rotation=45, ha='right', fontsize=7.5)
ax.set_yticks(range(len(rate_types)), [short(t) for t in rate_types], fontsize=7.5)
ax.set_xlabel('target cell type')
ax.set_ylabel('source cell type')
ax.set_title('estimated admixture rate $\\hat{r}_{S \\to T}$\n(% of target-type molecules)', fontsize=9)
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

ax = fig.add_subplot(gs[2, :])
ax.annotate('c', (-0.065, 1.05), xycoords='axes fraction',
    fontsize=13, fontweight='bold')
b = bare.dropna(subset=['power_strict']).sort_values('power_strict', ascending=False)
admixed = rates.assign(pair=rates['source'] + ' -> ' + rates['target']) \
    .set_index('pair')['admixed_molecules']
x = np.arange(len(b))
bars = ax.bar(x, b['power_strict'], color='#34495e', alpha=0.85)
ax2 = ax.twinx()
ax2.plot(x, b['pair'].map(admixed), 'o', color='#e67e22', ms=4, alpha=0.8)
ax2.set_yscale('log')
ax2.set_ylabel('estimated admixed molecules $\\hat{A}_{S \\to T}$', color='#e67e22')
ax2.tick_params(axis='y', colors='#e67e22')
ax.set_xticks(x, [short(p) for p in b['pair']], rotation=60, ha='right', fontsize=6.5)
ax.set_ylabel('estimated sensitivity (strict tier)')
ax.set_ylim(0, 1.05)
ax.set_title('all detected pairs under a standard single-fit cleanup (bare ls-NMF, membrane scoring)',
    fontsize=9)
fig.savefig(f'{OUT}/benchmark_fig1.png', bbox_inches='tight')
print('fig1 done')

# ---------- Figure 2: stochasticity ----------
fo = pd.read_csv(f'{R}/pancreas_factor_overlap.csv')
ov = json.load(open(f'{R}/pancreas_removal_overlap.json'))
dec_path = f'{R}/pancreas_decisions.csv'
have_dec = os.path.exists(dec_path)
vcolors = {'ls_nmf': '#16a085', 'invsqrt_kl': '#8e44ad'}
vnames = {'ls_nmf': 'ls-NMF', 'invsqrt_kl': 'invsqrt\nKL-NMF'}

def dot_panel(ax, data, letter, ylab, title):
    for k, variant in enumerate(['ls_nmf', 'invsqrt_kl']):
        v = np.asarray(data[variant])
        x = k + (np.random.default_rng(0).uniform(-0.14, 0.14, len(v)))
        ax.scatter(x, v, s=9, color=vcolors[variant], alpha=0.6)
        ax.hlines(np.median(v), k - 0.25, k + 0.25, color='black', lw=1.8)
    ax.set_xticks([0, 1], [vnames['ls_nmf'], vnames['invsqrt_kl']])
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel(ylab)
    ax.set_title(title, fontsize=9)
    ax.annotate(letter, (-0.28, 1.05), xycoords='axes fraction',
        fontsize=13, fontweight='bold')

fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.0))
dot_panel(axes[0],
    {v: fo[fo['variant'] == v]['corr'].values for v in ['ls_nmf', 'invsqrt_kl']},
    'a', 'gene-ownership correlation\nof matched factors', 'factor variability')
if have_dec:
    dec = pd.read_csv(dec_path)
    dec = dec[dec['method'] == 'membrane']
    dd = {}
    for variant in ['ls_nmf', 'invsqrt_kl']:
        sets = {s_: set(map(tuple, g[['source', 'target']].values))
                for s_, g in dec[dec['variant'] == variant].groupby('seed')}
        vals = []
        seeds = sorted(sets)
        for i in range(len(seeds)):
            for j in range(i + 1, len(seeds)):
                a, b = sets[seeds[i]], sets[seeds[j]]
                vals.append(len(a & b) / max(len(a | b), 1))
        dd[variant] = vals
    dot_panel(axes[1], dd, 'b', 'overlap of removal decisions\n(Jaccard)',
        'scoring variability')
J = {v: np.array(ov[v]['jaccard'])[np.triu_indices(10, 1)] for v in ['ls_nmf', 'invsqrt_kl']}
dot_panel(axes[2], J, 'c', 'overlap of removed molecules\n(Jaccard)', 'net effect')
fig.tight_layout()
fig.savefig(f'{OUT}/benchmark_fig2.png', bbox_inches='tight')
print('fig2 done')

# ---------- Figure 3: vote-threshold trade-off ----------
sweep_log = sys.argv[1] if len(sys.argv) > 1 else None
if sweep_log and os.path.exists(sweep_log):
    rows = []
    pat = re.compile(r'VOTE(\d+)OF(\d+) (\S+) (\S+)/(\S+): sens_strict=([\d.]+) '
                     r'sens_broad=([\d.]+) fpr_native=([\d.]+) spec_worst_pair=([\d.]+)')
    for line in open(sweep_log):
        m = pat.search(line)
        if m:
            rows.append(dict(thr=int(m[1]), n=int(m[2]), dataset=m[3], variant=m[4],
                method=m[5], sens=float(m[6]), fpr=float(m[8]), worst=float(m[9])))
    sw = pd.DataFrame(rows)
    sw.to_csv(f'{R}/vote_sweep_all.csv', index=False)
    combos = [('pancreas', 'invsqrt_kl', 'membrane', '#c0392b', '-', 'pancreas membrane (invsqrt)'),
              ('breast_crop', 'invsqrt_kl', 'membrane', '#2980b9', '-', 'breast membrane (invsqrt)'),
              ('pancreas', 'ls_nmf', 'bridge', '#c0392b', '--', 'pancreas bridge (ls)'),
              ('breast_crop', 'ls_nmf', 'bridge', '#2980b9', '--', 'breast bridge (ls)'),
              ('nsclc', 'ls_nmf', 'bridge', '#27ae60', '--', 'NSCLC bridge (ls)')]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.2))
    for ax in axes:
        ax.axvline(3, color='grey', lw=0.8, alpha=0.6)
        ax.set_xticks([1, 2, 3, 5, 7, 10])
        ax.set_xlabel('votes required (of 10 runs)')
    for dsname, variant, method, color, style, lab in combos:
        g = sw[(sw['dataset'] == dsname) & (sw['variant'] == variant) &
               (sw['method'] == method)].sort_values('thr')
        axes[0].plot(g['thr'], g['sens'], style, marker='o', ms=3.5, color=color, label=lab)
        axes[1].plot(g['thr'], np.maximum(g['fpr'] * 100, 0.005), style, marker='o',
            ms=3.5, color=color, label=lab)
    axes[0].set_ylabel('estimated sensitivity (strict tier)')
    axes[0].set_ylim(0, 1)
    axes[0].annotate('a', (-0.2, 1.04), xycoords='axes fraction', fontsize=13, fontweight='bold')
    axes[0].legend(frameon=False, fontsize=7)
    axes[1].set_yscale('log')
    axes[1].set_yticks([0.01, 0.1, 1, 10], ['0.01', '0.1', '1', '10'])
    axes[1].set_ylabel('own-marker false-removal rate (%)')
    axes[1].annotate('b', (-0.22, 1.04), xycoords='axes fraction', fontsize=13, fontweight='bold')
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
