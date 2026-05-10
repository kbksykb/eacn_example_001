"""
Dataset fetchers for the REAL / RareShield pilots.

All datasets are cached to /mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/data/
if that path is writable, else to workspace/data_cache/.

Staged for the CNS submission:
    - pancreas_baron_muraro_segerstolpe_wang (Tier A unit test, epsilon ground-truth)
    - hlca_core (Tier A, ionocyte ground-truth)
    - kang2024_pancancer (Tier B scale headline)
    - pbmc_10k_3batch (synthetic injection substrate)
"""

from __future__ import annotations

import logging
import os
import pathlib
from typing import Any

logger = logging.getLogger(__name__)

DEFAULT_SHARED = "/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/data"
FALLBACK_LOCAL = pathlib.Path(__file__).resolve().parents[2] / "data_cache"


def shared_root() -> pathlib.Path:
    try:
        p = pathlib.Path(DEFAULT_SHARED)
        p.mkdir(parents=True, exist_ok=True)
        probe = p / ".writable_probe"
        probe.write_text("ok")
        probe.unlink()
        return p
    except Exception as exc:
        logger.warning("shared path not writable (%s), falling back to %s", exc, FALLBACK_LOCAL)
        FALLBACK_LOCAL.mkdir(parents=True, exist_ok=True)
        return FALLBACK_LOCAL


def fetch_pancreas_scIB():
    """
    scIB pancreas benchmark (Baron + Muraro + Segerstolpe + Wang + Xin),
    hosted on figshare by the scIB team. Rare labels include epsilon, activated-stellate,
    and Schwann.

    Source: Luecken et al. 2022 Nat. Methods scIB benchmark datasets.
    URL: https://figshare.com/ndownloader/files/24539828  (human_pancreas_norm_complexBatch.h5ad)
    """
    import anndata as ad
    import urllib.request

    root = shared_root()
    dest = root / "pancreas_scib.h5ad"
    if dest.exists():
        logger.info("pancreas already cached at %s", dest)
        return ad.read_h5ad(dest)

    url = "https://figshare.com/ndownloader/files/24539828"
    logger.info("downloading pancreas scIB benchmark → %s", dest)
    tmp = dest.with_suffix(".tmp")
    urllib.request.urlretrieve(url, tmp)
    tmp.rename(dest)
    return ad.read_h5ad(dest)


def fetch_hlca_core():
    """
    HLCA core (Sikkema 2023). Primary ionocyte validation dataset.
    Download via cellxgene-census to get the core-only subset (~600k cells).
    """
    try:
        import cellxgene_census
    except ImportError as exc:
        raise RuntimeError(
            "cellxgene-census required: pip install cellxgene-census"
        ) from exc
    import anndata as ad

    root = shared_root()
    dest = root / "hlca_core.h5ad"
    if dest.exists():
        logger.info("HLCA core already cached at %s", dest)
        return ad.read_h5ad(dest)

    logger.info("downloading HLCA core via cellxgene-census → %s", dest)
    with cellxgene_census.open_soma(census_version="stable") as census:
        adata = cellxgene_census.get_anndata(
            census=census,
            organism="Homo sapiens",
            measurement_name="RNA",
            obs_value_filter=(
                "collection_name == 'Human Lung Cell Atlas (HLCA)' "
                "and is_primary_data == True "
                "and tissue_general == 'lung'"
            ),
        )
    adata.write_h5ad(dest)
    return adata


def fetch_kang2024_subset(cancer_types: list[str] | None = None, max_cells: int = 500_000):
    """
    Kang 2024 pan-cancer immune atlas (4.9M cells). For v1 pilot, return a
    stratified subsample by cancer type, capped at max_cells.
    """
    raise NotImplementedError(
        "Kang 2024 access path pending. Most likely CELLxGENE hosted; "
        "once confirmed, wire identically to fetch_hlca_core()."
    )


def mask_rare_label(adata, cell_type_key: str, rare_label: str) -> Any:
    """Set obs['heldout_rare']=True for rare, and NaN the cell_type label. Returns modified copy."""
    import numpy as np
    a = adata.copy()
    mask = a.obs[cell_type_key].astype(str).str.lower().str.strip() == rare_label.lower().strip()
    a.obs["heldout_rare"] = mask.values
    # Store the true label under a dunder key we'll only reveal ex post.
    a.obs["__hidden_label__"] = a.obs[cell_type_key].astype(str).copy()
    new = a.obs[cell_type_key].astype(object).copy()
    new[mask] = np.nan
    a.obs[cell_type_key] = new.astype("category")
    return a


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logger.info("fetching pancreas...")
    a = fetch_pancreas_scIB()
    logger.info("pancreas: %s cells x %s genes", a.n_obs, a.n_vars)
    logger.info("obs keys: %s", list(a.obs.columns))
    if "celltype" in a.obs.columns:
        key = "celltype"
    elif "cell_type" in a.obs.columns:
        key = "cell_type"
    else:
        key = None
    if key:
        logger.info("cell types (top): %s", a.obs[key].value_counts().head(30).to_dict())
    if "tech" in a.obs.columns:
        logger.info("tech/batch (top): %s", a.obs["tech"].value_counts().head(10).to_dict())
    elif "batch" in a.obs.columns:
        logger.info("batch (top): %s", a.obs["batch"].value_counts().head(10).to_dict())
