"""
REAL-Reg reference implementation stub — rare-aware scVI variant.

Design: workspace/proposals/unified_ML_spec_v2.md §2.1.
Interface contract: AnnData with obs['batch'], obs['seed_id'] (NaN for non-seed).
Output: obsm['X_rareaware'] of dim d=32 + obs['seed_preservation'] + obs['dominant_loss_mode'].

This file is a STUB only — not executed by the ML agent (per code-privilege boundary
in SHARED_CONTEXT.md §2 'Machine Learning Agent System').
Computational Biology (agent-mozur9ik) is the authorized executor. Please pick up
this file and wire it into the benchmark harness.

Copy this into examples/computational_biology/workspace/code/real_scvi.py on
branch agent/computational_biology. Do not attempt to run it from the
machine_learning discipline directory.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

# scvi-tools types are left as forward refs — import at runtime on the comp_bio side.


# ---------- hyperparameters ----------

@dataclass
class REALConfig:
    # Architecture (scVI default when None)
    n_latent: int = 32
    n_hidden: int = 128
    n_layers: int = 2

    # Witness seeds
    k_nn_seed: int = 15                # local density k
    K_nn_smooth: int = 300             # rarity-anchor smoothing K
    rarity_threshold: float = 2.0      # R_b(x) > this → seed candidate
    seed_size_min: int = 5
    seed_size_frac_max: float = 0.005  # fraction of batch size
    k_wit: int = 2                     # witness across ≥ k_wit batches
    sigma_thresh: float = 0.85         # signature cosine

    # CoLM admissible batch-effect field F
    admissible_rank_k: int = 4
    admissible_lipschitz: float = 0.5

    # REAL-Reg weights
    lambda_within: float = 1.0
    lambda_between: float = 1.0
    lambda_mass: float = 2.0
    margin: float = 3.0                # between-witness separation margin
    w_M_max: float = 0.1               # max REAL weight in total loss
    t_warmup: int = 5                  # epochs of warm-up
    n_witness_epoch: int = 5           # recompute witness set every N epochs

    # KL / posterior-collapse guards
    kl_warmup_epochs: int = 10
    free_bits_nats: float = 1.0


# ---------- seed / witness enumeration ----------
# Label-free. Must be isolated from any marker-gene oracle (see comp_bio task
# t-mozvf39m's "CRITICAL ISOLATION" clause).

def enumerate_seeds(E_b: np.ndarray, batch_size: int, cfg: REALConfig):
    """Per-batch seed enumeration. See unified_ML_spec_v2.md §1.1.

    Steps:
      1. kNN graph on E_b with k = cfg.k_nn_seed.
      2. rho_b(x) = 1 / r_k(x).
      3. R_b(x) = rho_b(x) / median(rho_b(kNN(K))).
      4. Connected components of {x : R_b(x) > rarity_threshold}.
      5. Keep components with size in [seed_size_min, seed_size_frac_max · n_b].

    Returns list[SeedDict] with keys 'cell_idx', 'signature', 'centroid'.
    """
    raise NotImplementedError(
        "Implementation handoff: comp_bio. Use GPU FAISS (IndexIVFPQ) + cugraph "
        "connected-components to hit the 4.9M-cell target (spec §1.5)."
    )


def construct_witnesses(seeds_by_batch, HVG_matrix_by_batch, cfg: REALConfig):
    """Cross-batch witness set W. See unified_ML_spec_v2.md §1.2.

    A seed is witnessed iff it has, in ≥ cfg.k_wit other batches, at least one
    seed with (i) signature cosine > cfg.sigma_thresh on rank-normalized z-scores,
    and (ii) an MNN link on the HVG panel.

    Returns list[WitnessDict] with keys 'seeds', 'cells', 'signature', 'centroid'.
    """
    raise NotImplementedError("Implementation handoff: comp_bio.")


# ---------- admissible-field Jacobian (precomputed, no grad into encoder) ----------

@torch.no_grad()
def compute_r_admissible(encoder: nn.Module, x_batch: torch.Tensor,
                         batch_idx: torch.Tensor, cfg: REALConfig) -> torch.Tensor:
    """Per-cell log-Jacobian upper bound under F (spec §3.4).

    r_admissible(x) = sigma_1(J_T(x)) * exp(-L * d(x, d_batch))
    Stop-grad: this value never receives gradient from the encoder.
    Cost: one Jacobian SVD per seed centroid, amortized over a full epoch.
    """
    raise NotImplementedError("Implementation handoff: comp_bio.")


# ---------- density head h_rho(z) ----------

class DensityHead(nn.Module):
    """Differentiable log q(z) estimator. Trained with ML loss on a held-out
    mini-batch each step. Encoder gets gradient through h_rho only after
    warm-up (spec §2.1 'stop-grad schedule')."""

    def __init__(self, n_latent: int, n_hidden: int = 128):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_latent, n_hidden), nn.SiLU(),
            nn.Linear(n_hidden, n_hidden), nn.SiLU(),
            nn.Linear(n_hidden, 1),
        )

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.net(z).squeeze(-1)  # log density


# ---------- REAL-Reg loss components (spec §2.1) ----------

def l_within(z: torch.Tensor, seed_assign: torch.Tensor) -> torch.Tensor:
    """Within-seed variance. Centroid is stop-grad."""
    # seed_assign: (B,) int64, -1 for non-seed cells.
    mask = seed_assign >= 0
    if mask.sum() == 0:
        return z.new_zeros(())
    losses = []
    for s in torch.unique(seed_assign[mask]):
        idx = (seed_assign == s).nonzero(as_tuple=True)[0]
        z_s = z[idx]
        centroid = z_s.mean(0).detach()
        losses.append(((z_s - centroid) ** 2).sum(-1).mean())
    return torch.stack(losses).mean()


def l_between(witness_centroids: torch.Tensor,
              non_match_pairs: torch.Tensor,
              margin: float) -> torch.Tensor:
    """Hinge margin between non-matching witness centroids."""
    if non_match_pairs.numel() == 0:
        return witness_centroids.new_zeros(())
    i = non_match_pairs[:, 0]
    j = non_match_pairs[:, 1]
    d = (witness_centroids[i] - witness_centroids[j]).norm(dim=-1)
    return F.relu(margin - d).pow(2).mean()


def l_mass(log_rho_pre: torch.Tensor, log_rho_post: torch.Tensor,
           log_r_adm: torch.Tensor) -> torch.Tensor:
    """CoLM per-cell inequality, squared-hinge."""
    violation = F.relu(log_rho_pre - log_rho_post - log_r_adm)
    return violation.pow(2).mean()


# ---------- orchestrator ----------

class REALRegScVI:
    """Trainer-level integrator. Wraps scvi-tools SCVI with REAL-Reg additions.

    High-level flow per epoch t:
      1. (every cfg.n_witness_epoch) recompute seeds + witnesses from current
         encoder's embedding.
      2. For each mini-batch:
         - standard scVI forward (ELBO with KL annealing + free bits)
         - add L_REAL = λ_w l_within + λ_b l_between + λ_m l_mass
         - weight w_M(t) = cfg.w_M_max * min(1, t/cfg.t_warmup)
         - total = ELBO + w_M(t) * L_REAL
      3. Also update density head h_rho with ML loss on a held-out batch.

    This is a STUB. Comp-Bio should implement the SCVI subclass hook (override
    `loss` / `inference` methods in the scvi-tools trainer) and wire
    seed/witness enumeration to run on GPU every cfg.n_witness_epoch.
    """

    def __init__(self, adata, cfg: Optional[REALConfig] = None):
        self.cfg = cfg or REALConfig()
        raise NotImplementedError(
            "This is a STUB. Implementation owner: comp_bio (agent-mozur9ik). "
            "See unified_ML_spec_v2.md §6 for the interface contract."
        )

    def fit(self, max_epochs: int = 400) -> None:
        raise NotImplementedError

    def get_latent(self) -> np.ndarray:
        """Return dense (n_cells, 32) embedding for obsm['X_rareaware']."""
        raise NotImplementedError

    def get_seed_preservation(self) -> np.ndarray:
        """Return per-cell seed_preservation ∈ [0, 1]."""
        raise NotImplementedError

    def get_dominant_loss_mode(self) -> np.ndarray:
        """Return per-cell dominant_loss_mode ∈ {'none', 'adjacency', 'asymmetry',
        'continuum', 'marker', 'emergent'} (5-class taxonomy, spec §1.4)."""
        raise NotImplementedError
