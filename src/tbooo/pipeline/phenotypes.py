"""Build synthetic UKB-mirrored phenotype Parquet table.

Populated fields (all others are null):
    eid
    p31                  sex (1=male, 2=female)
    p21000_i0            ethnic background (UKB Data-Coding 1001)
    p22006               in White British ancestry subset (0/1)
    p22020               used in PCA calculation (0/1)
    p22000               genotyping array batch code
    p22418               genotype calls available (1)
    p22828               imputed genotypes available (1)
    p23149               WGS CRAM available (1)
    p54_i0               assessment centre (synthetic code)

Output:
    data/Showcase/participant.parquet
"""

from __future__ import annotations

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from tbooo.config import Config
from tbooo.utils import ensure_dirs, log

# UKB Field 21000 (ethnic background) — Data-Coding 1001
_ETHNIC_BG: dict[str, int] = {
    # 1KGP superpopulations
    "EUR": 1,   # White
    "AFR": 4,   # Black or Black British
    "EAS": 3,   # Asian or Asian British
    "SAS": 3,   # Asian or Asian British
    "AMR": 2,   # Mixed
    # SGDP regions
    "West Eurasia": 1,
    "Africa": 4,
    "East Asia": 3,
    "South Asia": 3,
    "Central Asia / Siberia": 6,   # Other ethnic group
    "Oceania": 6,
    "Native Americas": 2,
    "Unknown": -1,
}

# Synthetic assessment centre codes (UKB Field 54)
_CENTRE: dict[str, int] = {
    "GBR": 11010,   # Leeds
    "EUR": 11020,   # Nottingham
    "AFR": 11021,   # Bristol
    "EAS": 11022,   # Hounslow
    "SAS": 11023,   # Croydon
    "AMR": 11024,   # Birmingham
}

# Null phenotype columns included so downstream scripts don't break
_NULL_FIELDS: list[str] = [
    "p21022",        # Age at recruitment
    "p53_i0",        # Date of attending assessment centre
    "p21001_i0",     # BMI
    "p4079_i0",      # Diastolic blood pressure
    "p4080_i0",      # Systolic blood pressure
    "p41270_a0",     # Diagnoses – ICD10 (first)
    "p20002_i0_a0",  # Non-cancer illness codes
    "p20003_i0_a0",  # Treatment/medication codes
    "p40001_i0",     # Cause of death – ICD10
]


def build_phenotype_table(cfg: Config) -> None:
    ensure_dirs(cfg.showcase_dir())

    frames: list[pd.DataFrame] = []

    kg_path = cfg.metadata_dir() / "eid_map_1kg.tsv"
    if kg_path.exists():
        frames.append(_build_1kg_rows(pd.read_csv(kg_path, sep="\t")))

    sgdp_path = cfg.metadata_dir() / "eid_map_sgdp.tsv"
    if sgdp_path.exists():
        frames.append(_build_sgdp_rows(pd.read_csv(sgdp_path, sep="\t")))

    if not frames:
        raise RuntimeError("No EID maps found. Run `tbooo map eids` first.")

    df = pd.concat(frames, ignore_index=True)

    # Merge expression PCs if already computed
    pcs_path = cfg.metadata_dir() / "geuvadis_expression_pcs.tsv"
    if pcs_path.exists():
        pcs = pd.read_csv(pcs_path, sep="\t")
        df = df.merge(pcs, on="eid", how="left")
        n_matched = pcs["eid"].isin(df["eid"]).sum()
        log(f"  merged GEUVADIS expression PCs: {n_matched}/{len(df)} rows")
    else:
        log("  GEUVADIS PCs not found; run `tbooo map geuvadis` to add expression scores")

    # Add null columns for commonly referenced clinical fields
    for col in _NULL_FIELDS:
        df[col] = pd.NA

    out = cfg.showcase_dir() / "participant.parquet"
    pq.write_table(pa.Table.from_pandas(df), str(out), compression="snappy")
    log(f"Phenotype table written → {out} ({len(df)} rows, {len(df.columns)} columns)")


def _build_1kg_rows(df: pd.DataFrame) -> pd.DataFrame:
    rows = pd.DataFrame()
    rows["eid"] = df["eid"].astype(int)
    rows["p31"] = df["sex"].astype(int)
    rows["p21000_i0"] = df["super_pop"].map(_ETHNIC_BG).fillna(-1).astype(int)
    rows["p22006"] = (
        (df["super_pop"] == "EUR") & (df["pop"] == "GBR")
    ).astype(int)
    rows["p22020"] = 1          # all Phase 3 unrelated samples used in PCA
    rows["p22000"] = df["batch"].fillna(0).astype(int)
    rows["p22418"] = 1          # array data available
    rows["p22828"] = 1          # imputed data available
    rows["p23149"] = 1          # WGS CRAM available
    rows["p54_i0"] = df.apply(
        lambda r: _CENTRE.get(r["pop"], _CENTRE.get(r["super_pop"], 11020)), axis=1
    )
    return rows


def _build_sgdp_rows(df: pd.DataFrame) -> pd.DataFrame:
    rows = pd.DataFrame()
    rows["eid"] = df["eid"].astype(int)
    rows["p31"] = df["sex"].fillna(0).astype(int)
    rows["p21000_i0"] = df["region"].map(_ETHNIC_BG).fillna(-1).astype(int)
    rows["p22006"] = 0          # SGDP samples are not White British
    rows["p22020"] = 0          # SGDP samples excluded from PCA
    rows["p22000"] = 0          # no batch code
    rows["p22418"] = 0          # no array data
    rows["p22828"] = 0          # no imputed data
    rows["p23149"] = 1          # WGS CRAM available
    rows["p54_i0"] = df["region"].map({
        "West Eurasia": 11020,
        "Africa": 11021,
        "East Asia": 11022,
        "South Asia": 11023,
        "Central Asia / Siberia": 11024,
        "Oceania": 11024,
        "Native Americas": 11024,
    }).fillna(11024).astype(int)
    return rows
