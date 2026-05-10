"""
REAL — Rare-subpopulation Erasure under Alignment, Label-free.

This package implements the channels that make up the REAL detectability
framework (team-unified name for DS's LRS + CompBio's RareScore + TumorBio's
Procrustes signature, see data_science/workspace/05_naming_unification.md).

Channels currently implemented here:
    channels.score_purity_and_procrustes  — mutual-kNN purity (P_n) + Procrustes
                                             displacement + bootstrap non-
                                             reproducibility.
    channels.loss_rate_at_k               — aggregate LossRate@k metric.

ML owns the CoLM channel (density-flow matching + admissible field F) and will
add `channels_colm.py`; DS owns motif (W) construction in lrs_framework and will
add `channels_mass_topology.py` for P_d / P_t / B.

The aggregate L(f; w) score is computed by combine.late_fuse() once all channels
return per-motif p-values.
"""

from .channels import (
    CandidateReport,
    RareScoreConfig as ChannelConfig,
    loss_rate_at_k,
    score as score_purity_and_procrustes,
)

__all__ = [
    "CandidateReport",
    "ChannelConfig",
    "loss_rate_at_k",
    "score_purity_and_procrustes",
]
