# cellAdmix Scoring Methods

This note describes the current cellAdmix admixture scoring methods: bridge
scoring, membrane scoring, and coherence scoring. All three methods operate
after a fitted cellAdmix run has assigned molecules to NMF factors. The scores
are then used to decide which factor molecules should be treated as potential
admixture in a target cell type.

## Motivation

The NMF stage identifies molecular factors. Some factors correspond to native
cell states. Others may correspond to molecules that entered a segmented target
cell from a neighboring, out-of-plane, or otherwise incorrectly segmented source
cell. The scoring stage asks a separate spatial question:

For factor $f$, does target cell type $T$ show spatial evidence that factor $f$
behaves like admixture from source cell type $S$?

The source identity of a factor and the target evidence for cleanup are related
but distinct:

- Source identity asks which annotated cell type best explains the factor.
- Target evidence asks where that factor appears in an unexpected spatial
  pattern.

Bridge and membrane scoring produce direct target/source/factor summaries.
Coherence scoring separates source/factor enrichment from target-cell-local
spatial coherence.

## Common Notation

Let:

- $c$ denote a cell.
- $a$ denote an ordered target cell.
- $b$ denote an ordered source cell.
- $T(c)$ denote the annotated cell type of cell $c$.
- $M_c$ denote the set of molecules assigned to cell $c$.
- $x_i$ denote the spatial coordinate of molecule $i$.
- $g_i$ denote the gene of molecule $i$.
- $z_i$ denote the fitted factor label of molecule $i$.
- $m_i$ denote the factor assignment margin for molecule $i$, when available.
- $f$ denote a factor.
- $K$ denote the number of fitted factors.

The number of molecules in cell $c$ assigned to factor $f$ is:

$$
N_{cf} = \left|\{i \in M_c : z_i = f\}\right|.
$$

The within-cell factor fraction is:

$$
p_{cf} = \frac{N_{cf}}{\max(|M_c|, 1)}.
$$

Each method produces scores for ordered triples: `target_cell_type`,
`source_cell_type`, and `factor`.

Most summaries report:

- `mean_score`: average evidence score across selected cells or cell pairs.
- `p_value`: one-sided statistical evidence that observed scores exceed a null
  or zero baseline.
- `neg_log10_p`: $-\log_{10}(p)$, where $p$ is the reported `p_value`.

The final cleanup rule for a target cell type is usually based on a threshold
such as `p_value < p_thresh`, plus method-specific gates.

## Bridge Score

### Intuition

Bridge scoring looks for factor molecules in a target cell that physically
touch or cross toward a neighboring source cell. It is meant to capture a
classic segmentation-boundary problem: a target cell contains a factor, and
factor-positive molecules sit near the molecular interface with a candidate
source cell.

The bridge score combines two ideas:

- The target and source cells differ in how much they contain the factor.
- Factor-positive target molecules are located close enough to source-cell
  molecules that local molecular neighborhoods cross the target/source
  boundary.

It is therefore conservative. It needs an explicit target/source cell pair and
enough local molecular contact between them.

### Candidate Pairs

Bridge scoring first builds ordered candidate cell pairs $(a, b)$, where $a$ is
the target and $b$ is the source.

Two candidate discovery modes are implemented:

- `molecule_global`: use molecule-level spatial kNN to find cells that contact
  through nearby transcripts.
- `cell_center`: use cell-centroid kNN to nominate nearby cells, then perform
  exact pair scoring locally.

For each ordered candidate pair, bridge scoring constructs a local spatial kNN
index from the molecules in both cells.

### Per-Pair Formula

For ordered pair $(a,b)$ and factor $f$, the target factor count is:

$$
N_{af} = \left|\{i \in M_a : z_i = f\}\right|.
$$

The pair is skipped for factor $f$ if $N_{af}$ is below
`min_factor_molecules`.

For each target molecule $i \in M_a$, query its `crossing_k` nearest neighbors
in the two-cell molecule set $M_a \cup M_b$. Define:

$$
r_i^{ab} =
\begin{cases}
1, & \text{if at least one neighbor of } i \text{ belongs to source cell } b,\\
0, & \text{otherwise.}
\end{cases}
$$

Then:

$$
C_{abf} = \sum_{\{i \in M_a : z_i = f\}} r_i^{ab},
$$

$$
C_{ab} = \sum_{i \in M_a} r_i^{ab},
$$

$$
q_{abf} = \frac{C_{abf}}{\max(N_{af}, 1)}.
$$

Here, $C_{abf}$ is the number of factor-$f$ target molecules whose local
neighborhood crosses into source cell $b$, $C_{ab}$ is stored as
`total_crossing_count`, and $q_{abf}$ is stored as `crossing_fraction`.

The bridge score is:

$$
s_{abf}^{\mathrm{bridge}} =
\left|p_{bf} - p_{af}\right| q_{abf}.
$$

The current implementation uses an absolute abundance difference. Directionality
enters through the ordered crossing term $q_{abf}$, which is computed from
target-cell molecules crossing toward source cell $b$.

### Statistical Summary

For each target/source cell-type pair and factor:

$$
G_{TSf} =
\{(a,b) : T(a)=T,\ T(b)=S,\ (a,b) \text{ was selected for factor } f\}.
$$

Groups with fewer than `min_pairs` are skipped. If the group has more than
`max_cells_per_type_pair` rows, a bounded subset is selected.

The summary mean score is:

$$
\bar{s}_{TSf}^{\mathrm{bridge}} =
\frac{1}{|G_{TSf}|}
\sum_{(a,b) \in G_{TSf}} s_{abf}^{\mathrm{bridge}}.
$$

If null computation is enabled, bridge scoring computes matched null scores for
the same selected cell pairs. The null procedure attempts to move factor
assignments or compare against replacement assignments while preserving the
same local pair structure. The summary p-value is a one-sided paired Wilcoxon
test:

$$
H_1: s_{abf}^{\mathrm{bridge}} > s_{abf}^{\mathrm{bridge,null}}.
$$

The output includes `mean_null_score`, `p_value`, and `neg_log10_p`.

## Membrane Score

### Intuition

Membrane scoring uses a stain image, typically a membrane-like Xenium morphology
channel, to ask whether factor-assigned target molecules lie along a
source-facing membrane boundary. It is designed for cases where the source cell
is visible or partially visible and the target cell contains source-like factor
molecules near the physical boundary between cells.

For each target/source cell pair and factor, the method compares stain signal
along rays from factor molecules to the target-cell centroid against control
rays from other target-cell molecules at similar radial distance. This controls
approximately for the fact that peripheral molecules see more membrane signal
than central molecules.

### Candidate Pairs

Membrane scoring uses cell centroid candidates. For each active target cell, it
queries nearby source cells by centroid distance. Candidate pairs may be
filtered by a halo around estimated cell radii:

$$
d(c_a, c_b) \le r_a + r_b + h,
$$

where $c_a$ and $c_b$ are cell centroids, $r_a$ and $r_b$ are estimated cell
radii, and $h$ is the halo distance.

Candidates are grouped by ordered target/source cell type and capped by
`candidate_pairs_per_type_pair`.

### Image Preprocessing

The stain image is read in physical coordinates. Image intensities are sampled
by bilinear interpolation. For score computation, the image is usually robustly
normalized to $[0,1]$ using low and high quantiles.

Let $I(x)$ be the normalized stain intensity at coordinate $x$. For a segment
from $u$ to $v$, define:

$$
L(u,v) = \max_{\ell \in [u,v]} I(\ell).
$$

In code, this is the maximum sampled line signal along the molecule-to-centroid
segment.

### Per-Molecule and Per-Pair Formula

For ordered pair $(a,b)$, let $c_a$ be the target-cell centroid and $c_b$ be
the source-cell centroid. For target-cell molecule $i$, define a source-facing
directional weight $w_i$. The implementation computes $w_i$ from the alignment
between the molecule direction from $c_a$ and the source-cell direction from
$c_a$ to $c_b$. Molecules on the source-facing side of the target cell get
higher weight.

For factor molecule $i \in M_a$ with $z_i=f$, choose a control molecule
$j \in M_a$ at a similar radial distance from the target centroid:

$$
\left| d(x_j,c_a) - d(x_i,c_a) \right| \le \tau.
$$

The tolerance is based on:

$$
\tau = \max(r_{\max,a} \alpha_{\mathrm{control}}, \epsilon).
$$

The per-molecule membrane score is:

$$
u_i =
\log
\left(
\frac{L(x_i,c_a)\max(w_i,10^{-6}) + \epsilon}
{L(x_j,c_a)\max(w_j,10^{-6}) + \epsilon}
\right).
$$

The pair/factor score is the mean over scored factor molecules:

$$
s_{abf}^{\mathrm{membrane}} =
\frac{1}{|\{i \in M_a : z_i=f\}|}
\sum_{\{i \in M_a : z_i=f\}} u_i.
$$

The pair output also reports:

$$
F_{abf}^{+} =
\frac{|\{i : u_i > 0\}|}{|\{i : u_i \text{ scored}\}|},
\qquad
\bar{w}_{abf} = \frac{1}{n}\sum_i w_i.
$$

These are stored as `fraction_positive` and `mean_directional_weight`.
Pairs/factors with fewer than `min_factor_molecules` target molecules are
skipped.

### Statistical Summary

For each target/source cell-type pair and factor:

$$
G_{TSf} =
\{(a,b) : T(a)=T,\ T(b)=S,\ s_{abf}^{\mathrm{membrane}} \text{ exists}\}.
$$

The highest-scoring rows are kept up to `max_cells_per_type_pair`; groups with
fewer than `min_pairs` are skipped.

The summary reports:

$$
\bar{s}_{TSf}^{\mathrm{membrane}} =
\frac{1}{|G_{TSf}|}
\sum_{(a,b)\in G_{TSf}} s_{abf}^{\mathrm{membrane}},
$$

$$
q_{TSf}^{0.75,\mathrm{membrane}} =
Q_{0.75}
\left(
\{s_{abf}^{\mathrm{membrane}} : (a,b)\in G_{TSf}\}
\right),
$$

$$
F_{TSf}^{+} =
\frac{
|\{(a,b)\in G_{TSf}: s_{abf}^{\mathrm{membrane}}>0\}|
}{
|G_{TSf}|
}.
$$

These are stored as `mean_score`, `q75_score`, and
`fraction_positive_pairs`. The p-value is a one-sided Wilcoxon signed-rank test
against zero:

$$
H_1: s_{abf}^{\mathrm{membrane}} > 0.
$$

Thus membrane scoring asks whether source-facing factor molecules see more
membrane signal than radial control molecules in the same target cells.

## Factor Source Score

The factor source score is not itself a cleanup score. It asks whether the gene
loading profile of a factor resembles marker expression for a particular
annotated cell type. It is useful as a source prior for bridge and membrane
scoring.

For each annotated cell type $S$, a marker weight vector $w_S$ is built from
the annotated cells' expression profiles: per gene, the one-vs-rest margin
(the type's mean expression minus the best other type's, floored at zero and
at `min_marker_logfc`, optionally restricted to the top
`top_markers_per_type` genes), raised to `specificity_power`. The factor's
non-negative loading vector $a_f$ (raised to `factor_power`) is then compared
to each type's weight vector by cosine similarity:

$$
\mathrm{score}_{Sf} =
\frac{a_f \cdot w_S}{\lVert a_f \rVert \, \lVert w_S \rVert}.
$$

The implementation reports, per factor × cell type, the `score`, its `rank`
among types, the factor's `best_score` and `second_score`, their difference
as `margin`, and an `is_best` flag. The source prior uses the best-scoring
cell type per factor, with `margin` measuring how decisive the call is.

## Coherence Score

### Intuition

Coherence scoring does not require an explicit source cell. Instead, it asks
whether source-like factor molecules form locally coherent patches inside a
target cell type. This helps when the source may be out of plane, partially
visible, missing from segmentation, or otherwise not well represented by a
source cell pair.

Coherence has two distinct components:

- A source/factor prior based on how enriched a factor is in each annotated
  source cell type.
- A target-cell-local spatial coherence score based on molecular kNN structure
  inside each target cell.

If a membrane image is supplied, same-factor molecular graph edges can be
downweighted when they cross strong membrane signal. This does not make
coherence a membrane score; it prevents the local molecular graph from
connecting patches across membrane barriers too easily.

### Source/Factor Enrichment

For source cell type $S$ and factor $f$, define:

$$
A_{Sf} = |\{i : T(c_i)=S,\ z_i=f\}|,
$$

$$
A_{S\cdot} = |\{i : T(c_i)=S\}|,
\qquad
A_{\cdot f} = |\{i : z_i=f\}|,
\qquad
A_{\cdot\cdot} = |\{i\}|.
$$

With pseudocount $\epsilon$, the implementation computes:

$$
u_{Sf} =
\frac{A_{Sf}+\epsilon}{A_{S\cdot}+\epsilon K},
$$

$$
v_{Sf} =
\frac{A_{\cdot f}-A_{Sf}+\epsilon}
{A_{\cdot\cdot}-A_{S\cdot}+\epsilon K},
$$

$$
e_{Sf} =
\log\left(\frac{u_{Sf}}{v_{Sf}}\right),
$$

$$
\pi_{Sf} =
\frac{A_{Sf}+\epsilon}
{A_{\cdot f}+\epsilon C},
$$

where $C$ is the number of annotated cell types. These are stored as
`type_fraction`, `other_fraction`, `source_log_enrichment`, and
`source_probability`.

By default, the R annotation layer uses `source_probability` to call the source
cell type for each factor. Target cleanup evidence is then evaluated
source-conditionally: the relevant row is `(target type, inferred source type,
factor)`.

### Within-Cell Molecular Graph

For each target cell $a$, coherence builds a molecule kNN graph using the
molecules in $M_a$. Let:

$$
N_i = \text{neighbors of molecule } i,
\qquad
w_{ij} = \text{edge weight from molecule } i \text{ to neighbor } j.
$$

If no membrane image is used, edge weights are based on distance or are equal
weights. If membrane barrier mode is used, edge weights are multiplied by:

$$
\exp(-\alpha_{\mathrm{membrane}} B_{ij}),
$$

where $B_{ij}$ is the maximum normalized membrane signal sampled along the
segment from $x_i$ to $x_j$.

### Per-Molecule Coherence Score

For target cell $a$, factor $f$, and source type $S$, consider the selected
molecules:

$$
Q_{af} = \{i \in M_a : z_i=f\}.
$$

The row is skipped if $|Q_{af}|$ is below `min_factor_molecules`.

For each selected molecule $i \in Q_{af}$, the raw source/factor score is:

$$
r_i = e_{Sf} + \beta_{\mathrm{margin}} m_i.
$$

The local same-factor support is:

$$
h_i =
\frac{
\sum_{\{j \in N_i : j \in Q_{af}\}} w_{ij}\max(r_j,0)
}{
\max\left(\sum_{j \in N_i} w_{ij}, \epsilon\right)
}.
$$

The final molecule score is:

$$
u_i = r_i + \lambda_{\mathrm{coherence}} h_i.
$$

The per-cell mean coherence score is:

$$
s_{aSf}^{\mathrm{mean}} =
\frac{1}{|Q_{af}|}
\sum_{i \in Q_{af}} u_i.
$$

The output also records:

$$
\bar{r}_{aSf} =
\frac{1}{|Q_{af}|}\sum_{i \in Q_{af}} r_i,
$$

$$
\bar{h}_{aSf} =
\frac{1}{|Q_{af}|}\sum_{i \in Q_{af}} h_i,
$$

$$
A_{aSf}^{\mathrm{active}} =
\frac{|\{i \in Q_{af}: u_i>\theta\}|}{|Q_{af}|}.
$$

These are stored as `mean_raw_score`, `mean_coherence_score`, and
`active_fraction`.

### Patch Score

The patch score is intended to identify a compact local group of source-like
factor molecules rather than just a diffuse cell-wide mean.

For each selected seed molecule $i \in Q_{af}$, define a one-hop selected
neighborhood:

$$
P_i =
\{i\}
\cup
\{j \in N_i : j \in Q_{af},\ w_{ij} \ge w_{\min}^{\mathrm{patch}}\}.
$$

The local patch mean is:

$$
\mu_i^{\mathrm{patch}} =
\frac{1}{|P_i|}
\sum_{j \in P_i} u_j.
$$

The patch size weight is:

$$
\omega_i^{\mathrm{patch}} =
\sqrt{\frac{|P_i|}{\max(|Q_{af}|,1)}}.
$$

The per-cell patch score is:

$$
s_{aSf}^{\mathrm{patch}} =
\max_{i \in Q_{af}}
\mu_i^{\mathrm{patch}}\omega_i^{\mathrm{patch}}.
$$

The implementation also computes `largest_patch_fraction` as a diagnostic: it
is the size of the largest connected component of selected factor molecules,
using graph edges with $w_{ij} \ge w_{\min}^{\mathrm{patch}}$, divided by
$|Q_{af}|$. This is not currently the primary patch score because connected
components can percolate through dense molecular graphs.

### Null Model and Summary

For each target cell, source type, and factor, coherence samples matched null
subsets with the same number of molecules as $Q_{af}$. By default, null subsets
avoid the observed factor molecules when enough other molecules are available;
otherwise they sample from all molecules in the target cell.

For each null subset, coherence recomputes the mean and patch scores, producing
cell-level null expectations $\bar{s}_{aSf}^{\mathrm{mean,null}}$ and
$\bar{s}_{aSf}^{\mathrm{patch,null}}$.

The per-cell deltas are:

$$
\Delta_{aSf}^{\mathrm{mean}} =
s_{aSf}^{\mathrm{mean}} -
\bar{s}_{aSf}^{\mathrm{mean,null}},
$$

$$
\Delta_{aSf}^{\mathrm{patch}} =
s_{aSf}^{\mathrm{patch}} -
\bar{s}_{aSf}^{\mathrm{patch,null}}.
$$

For each target/source cell-type pair and factor:

$$
G_{TSf} =
\{a : T(a)=T \text{ and score exists for source } S \text{ and factor } f\}.
$$

The summary reports `mean_score`, `mean_null_score`, `mean_delta_score`,
`mean_patch_score`, `mean_null_patch_score`, `mean_delta_patch_score`, and
`mean_largest_patch_fraction`.

The mean p-value tests:

$$
H_1:
s_{aSf}^{\mathrm{mean}} >
\bar{s}_{aSf}^{\mathrm{mean,null}}.
$$

The patch p-value tests:

$$
H_1:
s_{aSf}^{\mathrm{patch}} >
\bar{s}_{aSf}^{\mathrm{patch,null}}.
$$

The R annotation layer currently uses patch p-values by default when available,
equivalent to `score_mode = "patch"`. It can also use `score_mode = "mean"` or
`score_mode = "combined"`, where `combined` requires both mean and patch
evidence.

## Source Calls and Cleanup Calls

The plotting/correction layer collapses target/source/factor summaries into
factor-level source calls and target cleanup calls.

For coherence:

- Source calls default to the source cell type with highest
  `source_probability`.
- Target cleanup calls are source-conditioned: for factor $f$, if source $S^*$
  is inferred, target evidence is evaluated from row $(T,S^*,f)$.
- Cleanup calls can be gated by `p_thresh`, `min_delta_score`,
  `min_delta_patch_score`, `min_largest_patch_fraction`, `min_cells`, and
  `min_molecules`.

For bridge and membrane:

- Source calls are inferred from target/source/factor score summaries during
  annotation collapse.
- Cleanup calls are based on target evidence passing the p-value threshold,
  excluding the inferred source cell type itself.

This distinction matters: source identity should be interpreted as a factor
annotation, while cleanup calls are target-specific evidence.

## Native-Factor Check

Correction rules generated from any score pass a false-positive control by
default (`native_check = TRUE` in `rules()`): a counterfactual test of whether
the factor is native to the target cell type rather than admixed into it. A
factor can be genuinely expressed by a target type that is biologically
related to the inferred source — for example a malignant-cell factor in
epithelial cells — in which case removal would delete native signal.

For a rule removing factor $f$ from target type $T$ with inferred source type
$S$, target cells are stratified by source exposure: cell $a$ is
**source-distant** when none of its $k$ nearest neighbor cells (default
$k = 15$, cell centroids) are of type $S$, and **source-exposed** otherwise.
Let $x_{af}$ denote the cell-level fraction of factor $f$ in cell $a$, and let
$D$ and $E$ denote the source-distant and source-exposed cells of type $T$.
The rule is flagged as a likely false positive (`keep = FALSE`) when any of
the following holds:

- **Native persistence:** $\mathrm{median}_{a \in D}(x_{af}) > 0.1$ — the
  factor persists in target cells with no source contact.
- **Cross-type outlier:** the fraction of source-distant $T$ cells with
  $x_{af} > 0.05$ exceeds $Q_3 + 1.5\,\mathrm{IQR}$ of the same rate computed
  over all cell types except the source (the target's own rate is included
  in the reference), and exceeds $0.1$ — the target expresses the
  factor distinctly even in isolation.
- **No exposure gradient:**
  $\mathrm{mean}_{a \in E}(x_{af}) - \mathrm{mean}_{a \in D}(x_{af}) < 0$ —
  source proximity does not explain the signal.
- **Unverifiable:** the source or target type is absent from the annotation,
  fewer than two source-distant or source-exposed cells exist, or the factor
  is not expressed in either group. Rules that cannot be verified are excluded
  to be safe.

Flagged rules stay in the returned rule table with the reason in the
`native_check` column and the supporting evidence
(`native_distant_median`, `native_exposure_gradient`,
`native_distant_expr_frac`, group sizes); correction skips rules with
`keep = FALSE`. The source-exposure stratification uses the same
cell-neighborhood construction as the neighbor-enrichment score: enrichment
reads the positive direction (exposed target cells carry more of the factor),
while this check reads the negative control (unexposed cells still carry it).

## Strengths and Weaknesses

### Bridge Score

Strengths:

- Directly uses molecule-level adjacency between target and source cells.
- Provides high-specificity evidence for segmentation-boundary transfer.
- Does not require stain images.
- Uses a paired null model for the same local cell-pair structure.

Weaknesses:

- Requires an explicit visible or segmented source cell.
- Can miss out-of-plane sources, partial cells, or missing source cells.
- Can be conservative when molecule crossing is sparse.
- Runtime can become high if too many target/source pairs are scored.

### Membrane Score

Strengths:

- Uses independent imaging evidence.
- Can detect source-facing boundary effects even when molecule contacts are
  sparse.
- Radial controls reduce simple center-versus-periphery artifacts.
- Often works well on Xenium datasets with membrane-like morphology channels.

Weaknesses:

- Requires a useful stain image and correct image registration.
- Can be less informative when source cells are out of plane or not visible.
- Can confuse true membrane enrichment with other peripheral biological
  compartments if controls are insufficient.

### Coherence Score

Strengths:

- Does not require an explicit source cell.
- Can detect source-like molecular patches inside target cells.
- Can optionally incorporate membrane barriers without requiring a visible
  source pair.
- Separates source prior and target evidence.

Weaknesses:

- Depends more strongly on factor source priors.
- Patch statistics and null models are less direct than bridge or membrane
  scoring.
- Can be sensitive to molecule density, cell size, and local graph parameters.
- Current null models are approximate and may need dataset-specific validation.
