"""Quick stats for data/Showcase/participant.parquet.

Usage: python parquet_stats.py [path/to/participant.parquet]
"""
import sys
from pathlib import Path

import pandas as pd

# Default to the repo's parquet regardless of CWD (playground/ is one level
# under the repo root).
_DEFAULT = Path(__file__).resolve().parent.parent / "data" / "Showcase" / "participant.parquet"

path = sys.argv[1] if len(sys.argv) > 1 else str(_DEFAULT)
df = pd.read_parquet(path)

print(f"file: {path}")
print(f"samples (rows): {len(df):,}")
print(f"columns: {df.shape[1]}")

# Cohort split is encoded in the EID range (1kg: 1,000,000+, sgdp: 2,000,000+).
if "eid" in df.columns:
    src = pd.cut(df["eid"], [0, 2_000_000, 9_999_999], labels=["1kg", "sgdp"], right=False)
    print("\nby cohort (from EID range):")
    print(src.value_counts().to_string())

def dist(col, label):
    if col in df.columns:
        print(f"\n{label} ({col}):")
        print(df[col].value_counts(dropna=False).sort_index().to_string())

dist("p31", "sex (1=male, 2=female)")
dist("p21000_i0", "ethnic background")
dist("p54_i0", "assessment centre")

# Availability flags
flags = {"p22418": "array", "p22828": "imputed", "p24051": "WGS gVCF"}
present = {c: l for c, l in flags.items() if c in df.columns}
if present:
    print("\ndata-availability flags (count == 1):")
    for c, l in present.items():
        print(f"  {l:10s} {c}: {(df[c] == 1).sum():,}")

# GEUVADIS expression-PC coverage (non-null)
pc_cols = [c for c in df.columns if c.startswith("geuvadis_pc")]
if pc_cols:
    n = df[pc_cols[0]].notna().sum()
    print(f"\nGEUVADIS expression PCs ({len(pc_cols)} cols): {n:,} samples have scores, "
          f"{len(df) - n:,} null")

# Columns that are entirely null (the placeholder clinical fields)
all_null = [c for c in df.columns if df[c].isna().all()]
print(f"\nall-null columns ({len(all_null)}): {', '.join(all_null) if all_null else 'none'}")
