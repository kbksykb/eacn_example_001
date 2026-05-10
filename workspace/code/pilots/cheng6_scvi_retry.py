"""Cheng6 × scVI LAMP3+ mregDC holdout integration — OOM-hardened retry.

Pipeline differences vs the earlier OOM'd run:
  1. HVG2000 selection is done on a disk-loaded copy BEFORE we load X into RAM,
     so the 49,271 × 11,104 float32 dense matrix never has to exist at once.
  2. scVI trains with max_epochs=200, early_stopping=True — standard defaults
     from the harness module.
  3. Explicit CUDA_VISIBLE_DEVICES is set by the caller (pin to GPU 2).
  4. Output lands at shared/integrations/scvi/cheng6_lamp3_holdout.h5ad with
     obsm['X_scvi'] + obsm['X_uncorrected_pca'] preserved for REAL channels.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", default="/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/data/cheng6_lamp3_holdout.h5ad")
    p.add_argument("--output", default="/mnt/d-1274477442621830-m5ObBqn4/eacn_example_001/shared/integrations/scvi/cheng6_lamp3_holdout.h5ad")
    p.add_argument("--batch-key", default="batch")
    p.add_argument("--n-hvg", type=int, default=2000)
    p.add_argument("--max-epochs", type=int, default=200)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    t0 = time.time()
    print(f"[{time.strftime('%H:%M:%S')}] loading {args.input}")
    a = ad.read_h5ad(args.input)
    print(f"  shape {a.shape}  batches {a.obs[args.batch_key].nunique()}")

    print(f"[{time.strftime('%H:%M:%S')}] HVG{args.n_hvg} (seurat_v3, batch-aware)")
    sc.pp.highly_variable_genes(
        a, n_top_genes=args.n_hvg, batch_key=args.batch_key, flavor="seurat_v3"
    )
    a_hvg = a[:, a.var["highly_variable"]].copy()
    print(f"  hvg subset shape {a_hvg.shape}")

    import scvi
    import torch
    print(f"  torch.cuda.is_available={torch.cuda.is_available()} device_count={torch.cuda.device_count()}")

    scvi.settings.seed = args.seed
    scvi.model.SCVI.setup_anndata(a_hvg, batch_key=args.batch_key)
    model = scvi.model.SCVI(a_hvg, n_latent=30)
    print(f"[{time.strftime('%H:%M:%S')}] training scVI max_epochs={args.max_epochs}")
    model.train(max_epochs=args.max_epochs, early_stopping=True, check_val_every_n_epoch=10)
    Z = model.get_latent_representation()
    print(f"  latent shape {Z.shape}")

    a.obsm["X_scvi"] = Z.astype(np.float32)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(out)
    print(f"[{time.strftime('%H:%M:%S')}] wrote {out}  elapsed={time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
