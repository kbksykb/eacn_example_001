"""
RareScore — label-free detectability for rare-subpopulation loss in single-cell batch integration.

Three complementary channels, late-fused per candidate cluster:

    1. mutual_knn_purity        — CompBio baseline (this file).
    2. procrustes_displacement  — TumorBio suggestion (this file; bootstrap non-reproducibility).
    3. ot_density_flow          — ML channel (stub; ML agent owns formal version).

Ensemble: Fisher's method on per-channel p-values from permutation null.
Primary output:  per-candidate LossProbability and LossRate@k.

Design goals:
    - No external labels required.
    - CPU-usable on pancreas (~40k cells).
    - GPU-accelerable on HLCA / Kang 2024 via faiss-gpu for kNN.
    - Deterministic under a fixed seed.

Not yet implemented: OT channel (awaiting ML formal spec), protection-module gradient hooks.
"""

from __future__ import annotations

import dataclasses
import logging
from typing import Literal

import numpy as np
import scipy.sparse as sp

logger = logging.getLogger(__name__)


@dataclasses.dataclass
class RareScoreConfig:
    k: int = 30                             # kNN size
    n_permutations: int = 200               # permutation null
    n_bootstrap: int = 100                  # batch-resample bootstrap
    rare_abundance_threshold: float = 0.01  # candidates with ≤ 1% abundance
    purity_tau: float = 0.5                 # post-integration purity floor
    seed: int = 0
    device: Literal["cpu", "gpu"] = "cpu"


@dataclasses.dataclass
class CandidateReport:
    candidate_id: int
    abundance: float
    purity_pre: float
    purity_post: float
    nn_identity_pre: int
    nn_identity_post: int
    procrustes_displacement: float
    bootstrap_fraction_stable: float
    loss_probability: float


def _mutual_knn(
    emb: np.ndarray, k: int, seed: int
) -> sp.csr_matrix:
    """
    Mutual-kNN graph. Returns symmetric boolean CSR.
    Uses sklearn NearestNeighbors on CPU; swap for faiss-gpu on large data.
    """
    from sklearn.neighbors import NearestNeighbors

    rng = np.random.default_rng(seed)
    n = emb.shape[0]
    nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean", n_jobs=-1)
    nn.fit(emb)
    _, idx = nn.kneighbors(emb)
    idx = idx[:, 1:]  # drop self

    rows = np.repeat(np.arange(n), k)
    cols = idx.reshape(-1)
    data = np.ones_like(cols, dtype=np.uint8)
    directed = sp.csr_matrix((data, (rows, cols)), shape=(n, n))
    mutual = directed.multiply(directed.T).tocsr()
    return mutual.astype(bool)


def _cluster_purity(mutual: sp.csr_matrix, labels: np.ndarray) -> np.ndarray:
    """
    Per-cell mutual-kNN purity: fraction of mutual neighbours with the same cluster label.
    For a cell with 0 mutual neighbours, purity = 0 (it is evidence of isolation/absorption).
    """
    purity = np.zeros(mutual.shape[0], dtype=np.float32)
    labels = labels.astype(np.int32)
    for i in range(mutual.shape[0]):
        nbrs = mutual[i].indices
        if len(nbrs) == 0:
            continue
        purity[i] = (labels[nbrs] == labels[i]).mean()
    return purity


def _nn_majority_label(
    mutual: sp.csr_matrix, labels: np.ndarray
) -> np.ndarray:
    """
    For each cell, plurality-vote of mutual-kNN labels.
    Returns the cell's own label if no mutual neighbours.
    """
    out = labels.copy().astype(np.int32)
    for i in range(mutual.shape[0]):
        nbrs = mutual[i].indices
        if len(nbrs) == 0:
            continue
        vals, counts = np.unique(labels[nbrs], return_counts=True)
        out[i] = vals[np.argmax(counts)]
    return out


def _procrustes(pre: np.ndarray, post: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Orthogonal Procrustes pre→post in a shared normalized frame.
    Returns (pre_in_shared_frame, post_aligned_in_shared_frame, mean_residual).
    Per-cell displacement is computed in this shared frame so channels are comparable.
    """
    pre_c = pre - pre.mean(axis=0, keepdims=True)
    post_c = post - post.mean(axis=0, keepdims=True)
    s_pre = np.linalg.norm(pre_c) + 1e-12
    s_post = np.linalg.norm(post_c) + 1e-12
    pre_n = pre_c / s_pre
    post_n = post_c / s_post
    u, _, vt = np.linalg.svd(post_n.T @ pre_n, full_matrices=False)
    r = (u @ vt).T  # maps post→pre
    aligned = post_n @ r
    residual = float(np.linalg.norm(pre_n - aligned) / np.sqrt(pre.shape[0]))
    return pre_n, aligned, residual


def _candidate_displacement(
    pre_frame: np.ndarray, post_aligned: np.ndarray, mask: np.ndarray
) -> float:
    """Mean per-cell residual in the shared normalized frame."""
    if mask.sum() == 0:
        return 0.0
    d = np.linalg.norm(pre_frame[mask] - post_aligned[mask], axis=1)
    return float(d.mean())


def score(
    emb_pre: np.ndarray,
    emb_post: np.ndarray,
    candidate_labels: np.ndarray,
    batch_labels: np.ndarray,
    cfg: RareScoreConfig | None = None,
) -> list[CandidateReport]:
    """
    Compute RareScore per rare candidate cluster.

    emb_pre, emb_post: shape (n_cells, d). Must be row-aligned (same cells).
    candidate_labels: pre-integration over-clustering (int).
    batch_labels: batch id per cell (int).

    Returns one CandidateReport per candidate with abundance ≤ cfg.rare_abundance_threshold.
    """
    cfg = cfg or RareScoreConfig()
    n = emb_pre.shape[0]
    assert emb_post.shape[0] == n
    assert candidate_labels.shape == (n,)
    assert batch_labels.shape == (n,)

    logger.info("RareScore: %d cells, %d candidate clusters", n, len(np.unique(candidate_labels)))

    mutual_pre = _mutual_knn(emb_pre, cfg.k, cfg.seed)
    mutual_post = _mutual_knn(emb_post, cfg.k, cfg.seed)

    purity_pre = _cluster_purity(mutual_pre, candidate_labels)
    purity_post = _cluster_purity(mutual_post, candidate_labels)
    majority_post = _nn_majority_label(mutual_post, candidate_labels)

    post_aligned_pre_frame, post_aligned, _ = _procrustes(emb_pre, emb_post)

    rng = np.random.default_rng(cfg.seed)
    unique = np.unique(candidate_labels)
    reports: list[CandidateReport] = []

    for cid in unique:
        mask = candidate_labels == cid
        abundance = float(mask.mean())
        if abundance > cfg.rare_abundance_threshold:
            continue

        c_purity_pre = float(purity_pre[mask].mean()) if mask.any() else 0.0
        c_purity_post = float(purity_post[mask].mean()) if mask.any() else 0.0
        c_nn_pre = int(cid)
        # Post-integration identity: plurality of the candidate's cells'
        # post-integration kNN-majority label.
        if mask.any():
            vals, counts = np.unique(majority_post[mask], return_counts=True)
            c_nn_post = int(vals[np.argmax(counts)])
        else:
            c_nn_post = int(cid)

        displacement = _candidate_displacement(post_aligned_pre_frame, post_aligned, mask)

        stable = 0
        for _ in range(cfg.n_bootstrap):
            boot_batches = rng.choice(
                np.unique(batch_labels),
                size=np.unique(batch_labels).size,
                replace=True,
            )
            boot_mask = np.isin(batch_labels, boot_batches) & mask
            if boot_mask.sum() < 3:
                continue
            vals, counts = np.unique(majority_post[boot_mask], return_counts=True)
            if vals[np.argmax(counts)] == cid:
                stable += 1
        bootstrap_stable = stable / max(cfg.n_bootstrap, 1)

        loss_prob = _combine_channels(
            purity_pre=c_purity_pre,
            purity_post=c_purity_post,
            nn_shift=(c_nn_pre != c_nn_post),
            displacement=displacement,
            bootstrap_stable=bootstrap_stable,
        )

        reports.append(
            CandidateReport(
                candidate_id=int(cid),
                abundance=abundance,
                purity_pre=c_purity_pre,
                purity_post=c_purity_post,
                nn_identity_pre=c_nn_pre,
                nn_identity_post=c_nn_post,
                procrustes_displacement=displacement,
                bootstrap_fraction_stable=bootstrap_stable,
                loss_probability=loss_prob,
            )
        )

    reports.sort(key=lambda r: r.loss_probability, reverse=True)
    return reports


def _combine_channels(
    purity_pre: float,
    purity_post: float,
    nn_shift: bool,
    displacement: float,
    bootstrap_stable: float,
) -> float:
    """
    Heuristic ensemble of the implemented channels, returning a 0..1 loss probability.
    Formal Fisher/Stouffer version pending OT channel integration.
    """
    purity_drop = max(0.0, purity_pre - purity_post)
    displacement_norm = 1.0 - np.exp(-displacement * 500.0)  # displacement is in normalized frame (~1e-3 scale)
    instability = 1.0 - bootstrap_stable
    shift = 1.0 if nn_shift else 0.0

    # Weighted heuristic; weights to be re-fit on pancreas epsilon hold-out.
    loss = (
        0.35 * purity_drop
        + 0.25 * shift
        + 0.20 * displacement_norm
        + 0.20 * instability
    )
    return float(np.clip(loss, 0.0, 1.0))


def loss_rate_at_k(reports: list[CandidateReport], k: int, tau: float = 0.5) -> float:
    """LossRate@k — fraction of top-k rarest candidates called as lost."""
    top_k = sorted(reports, key=lambda r: r.abundance)[:k]
    if not top_k:
        return 0.0
    return sum(r.loss_probability >= tau for r in top_k) / len(top_k)
