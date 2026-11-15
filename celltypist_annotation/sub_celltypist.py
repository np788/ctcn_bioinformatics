#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc
import celltypist


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate cells with CellTypist."
    )
    parser.add_argument("--sample", required=True)
    parser.add_argument("--h5", required=True, type=Path)
    parser.add_argument("--barcodes", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--model", default="Immune_All_Low.pkl")
    return parser.parse_args()


def load_qc_barcodes(path: Path) -> list[str]:
    barcodes = pd.read_csv(path, header=None, sep="\t")[0].astype(str).tolist()
    if not barcodes:
        raise ValueError(f"No barcodes found in {path}")
    return barcodes


def resolve_barcodes(sample_id: str, adata_barcodes: pd.Index, qc_names: list[str]) -> list[str]:
    adata_set = set(map(str, adata_barcodes))

    exact = [b for b in qc_names if b in adata_set]
    if len(exact) == len(qc_names):
        return exact

    prefix = f"{sample_id}_"
    stripped = [b[len(prefix):] if b.startswith(prefix) else b for b in qc_names]
    missing = sorted(set(stripped) - adata_set)

    return stripped


def main() -> None:
    args = parse_args()

    if not args.h5.exists():
        raise FileNotFoundError(args.h5)
    if not args.barcodes.exists():
        raise FileNotFoundError(args.barcodes)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    # Keep choices from last CellTypist fl
    adata = sc.read_10x_h5(args.h5, gex_only=True)
    adata.var_names_make_unique()

    qc_names = load_qc_barcodes(args.barcodes)
    raw_barcodes = resolve_barcodes(args.sample, adata.obs_names, qc_names)
    adata = adata[raw_barcodes].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    result = celltypist.annotate(
        adata,
        model=args.model,
        majority_voting=True,
    )

    predictions = result.predicted_labels.copy()
    predictions.index = qc_names
    predictions.index.name = "barcode"
    predictions.to_csv(args.output, sep="\t")

    print(f"Annotated {len(predictions):,} cells")
    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
