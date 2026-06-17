"""Assign synthetic 7-digit EIDs to all 1KGP and SGDP samples.

Outputs:
    data/metadata/eid_map_1kg.tsv   columns: eid, sample_id, pop, super_pop, sex, source
    data/metadata/eid_map_sgdp.tsv  columns: eid, ena_accession, population, region, sex, source
    data/metadata/vcf_sample_rename_1kg.txt   old_id → eid  (one per line, space-separated)
    data/metadata/vcf_sample_rename_sgdp.txt  old_id → eid
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from tbooo.config import Config
from tbooo.integrity import remove, table_ok
from tbooo.utils import ensure_dirs, log

# 1KGP superpop → batch code (used in FAM column 6)
_BATCH_MAP = {"EUR": 1, "AFR": 2, "EAS": 3, "SAS": 4, "AMR": 5}


def assign_eids(cfg: Config) -> None:
    ensure_dirs(cfg.metadata_dir())
    _assign_1kg(cfg)
    _assign_sgdp(cfg)
    log("EID assignment complete.")


# ── 1KGP ─────────────────────────────────────────────────────────────────────

def _assign_1kg(cfg: Config) -> None:
    out_map = cfg.metadata_dir() / "eid_map_1kg.tsv"
    out_rename = cfg.metadata_dir() / "vcf_sample_rename_1kg.txt"

    if table_ok(out_map, min_lines=2) and table_ok(out_rename, min_lines=1):
        log(f"  skip (valid): {out_map.name} + {out_rename.name}")
        return
    remove(out_map, out_rename)  # clear any partial write before rebuilding

    # Try NYGC panel first (3,202 samples); fall back to Phase 3 panel (2,504)
    nygc_panel = cfg.kg_raw_dir() / "20130606_g1k_3202_samples_ped_population.txt"
    phase3_panel = cfg.kg_raw_dir() / f"integrated_call_samples_v3.{cfg.kg_phase3_release_date}.ALL.panel"

    if nygc_panel.exists():
        log("  reading NYGC 30x sample panel (3,202 samples)…")
        panel = _read_nygc_panel(nygc_panel)
    elif phase3_panel.exists():
        log("  reading Phase 3 sample panel (2,504 samples)…")
        panel = _read_phase3_panel(phase3_panel)
    else:
        raise FileNotFoundError(
            f"No 1KGP sample panel found. Run `tbooo download 1kg` first.\n"
            f"Expected: {nygc_panel} or {phase3_panel}"
        )

    panel = panel.reset_index(drop=True)
    panel.insert(0, "eid", range(cfg.kg_eid_start, cfg.kg_eid_start + len(panel)))
    panel["source"] = "1kg"
    panel["batch"] = panel["super_pop"].map(_BATCH_MAP).fillna(0).astype(int)

    panel.to_csv(out_map, sep="\t", index=False)
    log(f"  wrote {out_map} ({len(panel)} samples)")

    # rename file: "NA12878 1000001" (original_id eid)
    with open(out_rename, "w") as f:
        for _, row in panel.iterrows():
            f.write(f"{row['sample_id']}\t{row['eid']}\n")
    log(f"  wrote {out_rename}")


def _read_phase3_panel(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    # columns: sample, pop, super_pop, gender
    df = df.rename(columns={"sample": "sample_id", "gender": "sex_label"})
    df["sex"] = df["sex_label"].map({"male": 1, "female": 2, "unknown": 0}).fillna(0).astype(int)
    return df[["sample_id", "pop", "super_pop", "sex"]]


def _read_nygc_panel(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", engine="python")

    # 20130606_g1k_3202_samples_ped_population.txt is PED+pop format:
    # FamilyID  SampleID  FatherID  MotherID  Sex  Phenotype  Population  Superpopulation
    # Detect by checking for PED-style columns.
    cols_lower = [c.lower() for c in df.columns]
    is_ped = "fatherid" in cols_lower or "father_id" in cols_lower or (
        len(df.columns) >= 6 and cols_lower[2] in ("fatherid", "father_id", "father")
    )

    if is_ped:
        # Normalise column names regardless of capitalisation
        rename = {}
        for c in df.columns:
            lc = c.lower()
            if lc in ("sampleid", "sample_id", "individualid", "individual_id") or (
                "sample" in lc and "id" in lc
            ):
                rename[c] = "sample_id"
            elif lc in ("sex", "gender"):
                rename[c] = "sex_raw"
            elif lc in ("population", "pop"):
                rename[c] = "pop"
            elif "super" in lc:
                rename[c] = "super_pop"
        df = df.rename(columns=rename)
        # PED sex encoding: 1=male, 2=female, 0=unknown
        if "sex_raw" in df.columns:
            df["sex"] = pd.to_numeric(df["sex_raw"], errors="coerce").fillna(0).astype(int)
        else:
            df["sex"] = 0
    else:
        # Generic fallback for TSV panels with labelled columns
        col_map = {}
        for c in df.columns:
            lc = c.lower()
            if "sample" in lc:
                col_map[c] = "sample_id"
            elif lc in ("population", "pop"):
                col_map[c] = "pop"
            elif "super" in lc:
                col_map[c] = "super_pop"
            elif lc in ("sex", "gender"):
                col_map[c] = "sex_label"
        df = df.rename(columns=col_map)
        if "sex_label" in df.columns:
            df["sex"] = df["sex_label"].map(
                {"male": 1, "female": 2, "Male": 1, "Female": 2, "1": 1, "2": 2}
            ).fillna(0).astype(int)
        else:
            df["sex"] = 0

    for col in ("pop", "super_pop"):
        if col not in df.columns:
            df[col] = "UNK"

    return df[["sample_id", "pop", "super_pop", "sex"]]


# ── SGDP ─────────────────────────────────────────────────────────────────────

def _assign_sgdp(cfg: Config) -> None:
    out_map = cfg.metadata_dir() / "eid_map_sgdp.tsv"
    out_rename = cfg.metadata_dir() / "vcf_sample_rename_sgdp.txt"

    if table_ok(out_map, min_lines=2) and table_ok(out_rename, min_lines=1):
        log(f"  skip (valid): {out_map.name} + {out_rename.name}")
        return
    remove(out_map, out_rename)  # clear any partial write before rebuilding

    samples_tsv = cfg.sgdp_raw_dir() / "sgdp_samples.tsv"
    if not samples_tsv.exists():
        log("  WARNING: SGDP metadata not found; skipping SGDP EID assignment.")
        log(f"  Run `tbooo download sgdp` to fetch {samples_tsv}")
        return

    df = pd.read_csv(samples_tsv, sep="\t")
    if df.empty:
        log("  WARNING: SGDP metadata file is empty; skipping.")
        return

    df = df.reset_index(drop=True)
    df.insert(0, "eid", range(cfg.sgdp_eid_start, cfg.sgdp_eid_start + len(df)))
    df["source"] = "sgdp"
    df.to_csv(out_map, sep="\t", index=False)
    log(f"  wrote {out_map} ({len(df)} samples)")

    with open(out_rename, "w") as f:
        for _, row in df.iterrows():
            f.write(f"{row['ena_accession']}\t{row['eid']}\n")
    log(f"  wrote {out_rename}")
