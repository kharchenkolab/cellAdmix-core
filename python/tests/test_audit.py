"""Planted-contamination tests for the admixture-audit metric internals."""

import numpy as np
import pandas as pd
from scipy import sparse

from celladmix.audit import (
    BIN_EDGES,
    BIN_LABELS,
    _bin_rates,
    _excess,
    _marker_pool,
    _power,
    _pseudobulk,
)


def make_scenario():
    genes = [f"smk{i}" for i in range(4)] + [f"tmk{i}" for i in range(4)] + ["shared1"]
    n_t, n_s = 200, 100
    cells = [f"T{i}" for i in range(n_t)] + [f"S{i}" for i in range(n_s)]
    types = np.array(["T"] * n_t + ["S"] * n_s)
    exposure = np.array([i % 5 for i in range(n_t)] + [0] * n_s)
    counts = np.zeros((len(genes), len(cells)))
    counts[4:8, types == "T"] = np.array([50, 40, 30, 20])[:, None]
    counts[0:4, types == "T"] = 2
    counts[8, :] = 30
    counts[0:4, types == "S"] = 60
    planted = np.zeros_like(counts)
    planted[0:4, :n_t] = 8 * exposure[:n_t][None, :]
    return genes, cells, types, exposure, counts, planted


def test_marker_pool_and_strict_tier():
    genes, cells, types, exposure, counts, planted = make_scenario()
    before = sparse.csc_matrix(counts + planted)
    cell_index = pd.Index(cells)
    t_cells = pd.Index([c for c, t in zip(cells, types) if t == "T"])
    ref = t_cells[exposure[: len(t_cells)] == 0]
    profiles = np.column_stack([
        _pseudobulk(before, cell_index, pd.Index([c for c, t in zip(cells, types) if t == ty]))
        for ty in ["S", "T"]
    ])
    baseline = _pseudobulk(before, cell_index, ref)
    pool = _marker_pool(profiles, ["S", "T"], "S", np.nan_to_num(baseline), n_pool=4)
    assert sorted(np.array(genes)[pool]) == [f"smk{i}" for i in range(4)]


def test_excess_detection_and_power():
    genes, cells, types, exposure, counts, planted = make_scenario()
    cell_index = pd.Index(cells)
    t_mask = types == "T"
    t_cols = np.flatnonzero(t_mask)
    bins = pd.cut(exposure[t_mask], BIN_EDGES, labels=BIN_LABELS).astype(str)
    pool = np.arange(4)

    def rates_of(mat):
        m = sparse.csc_matrix(mat)
        totals = np.asarray(sparse.csc_matrix(counts + planted).sum(axis=0)).ravel()
        mk = np.asarray(m[pool][:, t_cols].sum(axis=0)).ravel()
        return _bin_rates(mk, totals[t_cols], np.asarray(bins))

    rb = rates_of(counts + planted)
    assert (rb["rate_lo"] <= rb["rate"]).all() and (rb["rate"] <= rb["rate_hi"]).all()
    excess, p = _excess(rb)
    assert p < 1e-10
    assert abs(excess - planted.sum()) / planted.sum() < 0.15

    for corrected, expected in [
        (counts, 1.0),
        (counts + planted, 0.0),
        (counts + planted / 2, 0.5),
    ]:
        assert abs(_power(rb, rates_of(corrected)) - expected) < 0.06
