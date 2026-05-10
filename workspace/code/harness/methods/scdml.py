"""
scDML / scDREAMER / scCRAFT / Seurat-RPCA placeholders.

Each is a separate installable package; to keep the harness loadable without
optional deps, these wrappers lazy-import and raise a clear error if the package
is missing. CompBio will install each in sequence; until then, calls to these
methods raise ImportError and the harness records the error in the runtime JSON
instead of crashing the whole batch.
"""

from __future__ import annotations
import numpy as np


def run(adata, batch_key: str, seed: int) -> np.ndarray:
    raise NotImplementedError(
        "scDML wrapper not yet installed. `pip install scDML` then replace this stub."
    )
