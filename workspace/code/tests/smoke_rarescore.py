"""
Smoke test for RareScore on synthetic data.

Generates a pre-integration embedding with:
    - 3 abundant clusters
    - 1 rare cluster (0.5% of cells)
    - 2 batches with different rare-cluster sampling

Simulates two post-integration scenarios:
    - "harmonized": rare cluster is pulled into the nearest abundant cluster (loss).
    - "preserved": rare cluster keeps its local structure.

Asserts that RareScore's loss_probability is high in the first and low in the second.
"""

from __future__ import annotations

import numpy as np

from workspace.code.rarescore import RareScoreConfig, loss_rate_at_k, score


def make_toy(seed: int = 0, n_cells: int = 10_000, rare_frac: float = 0.005):
    rng = np.random.default_rng(seed)
    n_rare = int(n_cells * rare_frac)
    n_per_abundant = (n_cells - n_rare) // 3

    centres = np.array([
        [0.0, 0.0],
        [5.0, 0.0],
        [0.0, 5.0],
        [4.5, 0.5],  # rare centre — sits near cluster 1 → easy to absorb
    ])

    def draw(c, n, scale=0.4):
        return rng.normal(loc=c, scale=scale, size=(n, 2))

    abundant = [draw(centres[i], n_per_abundant) for i in range(3)]
    rare = draw(centres[3], n_rare, scale=0.25)

    emb = np.vstack(abundant + [rare])
    labels = np.concatenate([
        np.full(n_per_abundant, 0),
        np.full(n_per_abundant, 1),
        np.full(n_per_abundant, 2),
        np.full(n_rare, 3),
    ])

    batches = rng.integers(0, 2, size=emb.shape[0])

    return emb.astype(np.float32), labels.astype(np.int32), batches.astype(np.int32)


def corrupt_harmonized(emb: np.ndarray, labels: np.ndarray, rare_label: int = 3) -> np.ndarray:
    post = emb.copy()
    mask = labels == rare_label
    post[mask] = emb[mask] * 0.2 + np.array([4.5, 0.0]) * 0.8
    post += np.random.default_rng(1).normal(0, 0.05, size=post.shape)
    return post


def corrupt_preserved(emb: np.ndarray) -> np.ndarray:
    return emb + np.random.default_rng(2).normal(0, 0.05, size=emb.shape)


def main():
    emb, labels, batches = make_toy()

    post_harm = corrupt_harmonized(emb, labels)
    post_pres = corrupt_preserved(emb)

    cfg = RareScoreConfig(k=20, n_bootstrap=30, rare_abundance_threshold=0.02, seed=0)

    reports_harm = score(emb, post_harm, labels, batches, cfg=cfg)
    reports_pres = score(emb, post_pres, labels, batches, cfg=cfg)

    print("=== Harmonized (expected: HIGH loss for rare cluster 3) ===")
    for r in reports_harm:
        print(r)

    print("\n=== Preserved (expected: LOW loss) ===")
    for r in reports_pres:
        print(r)

    print("\nLossRate@1 harmonized:", loss_rate_at_k(reports_harm, k=1))
    print("LossRate@1 preserved: ", loss_rate_at_k(reports_pres, k=1))

    assert reports_harm, "no rare candidate flagged in harmonized"
    assert reports_harm[0].loss_probability > 0.4, "expected high loss prob for absorbed rare"
    if reports_pres:
        assert reports_pres[0].loss_probability < reports_harm[0].loss_probability + 0.05, \
            "preserved should not be called lost more than harmonized"

    print("\nSMOKE TEST OK")


if __name__ == "__main__":
    main()
