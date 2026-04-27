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
from tbooo.utils import ensure_dirs, log

# 1KGP superpop → batch code (used in FAM column 6)
_BATCH_MAP = {"EUR": 1, "AFR": 2, "EAS": 3, "SAS": 4, "AMR": 5}

# SGDP region → continent label (for phenotype mapping)
_SGDP_REGION_LABELS = {
    "WestEurasia": "West Eurasia",
    "Africa": "Africa",
    "EastAsia": "East Asia",
    "SouthAsia": "South Asia",
    "CentralAsiaSiberia": "Central Asia / Siberia",
    "Oceania": "Oceania",
    "America": "Native Americas",
}


def assign_eids(cfg: Config) -> None:
    ensure_dirs(cfg.metadata_dir())
    _assign_1kg(cfg)
    _assign_sgdp(cfg)
    log("EID assignment complete.")


# ── 1KGP ─────────────────────────────────────────────────────────────────────

def _assign_1kg(cfg: Config) -> None:
    out_map = cfg.metadata_dir() / "eid_map_1kg.tsv"
    out_rename = cfg.metadata_dir() / "vcf_sample_rename_1kg.txt"

    # Try NYGC panel first (3,202 samples); fall back to Phase 3 panel (2,504)
    nygc_panel = cfg.kg_raw_dir() / "20201028_3202_samples_5_subpopulations.tsv"
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
    df = pd.read_csv(path, sep="\t")
    # NYGC panel columns vary — handle common layouts
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

    pointers = cfg.sgdp_raw_dir() / "ena.ftp.pointers.txt"
    if not pointers.exists():
        log("  WARNING: SGDP pointers file not found; skipping SGDP EID assignment.")
        log(f"  Run `tbooo download sgdp` to fetch {pointers}")
        return

    rows: list[dict] = []
    seen: set[str] = set()
    with open(pointers) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            sample, population, url = parts[0], parts[1], parts[2]
            if not url.endswith(".cram"):
                continue
            if sample in seen:
                continue
            seen.add(sample)
            region = _infer_region(population)
            sex = _infer_sex_from_url(url)
            rows.append({"ena_accession": sample, "population": population,
                         "region": region, "sex": sex})

    df = pd.DataFrame(rows).reset_index(drop=True)
    df.insert(0, "eid", range(cfg.sgdp_eid_start, cfg.sgdp_eid_start + len(df)))
    df["source"] = "sgdp"
    df.to_csv(out_map, sep="\t", index=False)
    log(f"  wrote {out_map} ({len(df)} samples)")

    with open(out_rename, "w") as f:
        for _, row in df.iterrows():
            f.write(f"{row['ena_accession']}\t{row['eid']}\n")
    log(f"  wrote {out_rename}")


def _infer_region(population: str) -> str:
    # Best-effort region label; can be improved with a full SGDP metadata table
    mapping = {
        "Greek": "WestEurasia", "French": "WestEurasia", "English": "WestEurasia",
        "Yoruba": "Africa", "Mandinka": "Africa", "Zulu": "Africa",
        "Han": "EastAsia", "Japanese": "EastAsia", "Vietnamese": "EastAsia",
        "Bengali": "SouthAsia", "Punjabi": "SouthAsia", "Tamil": "SouthAsia",
        "Papuan": "Oceania", "Australian": "Oceania",
        "Maya": "America", "Quechua": "America",
        "Buryat": "CentralAsiaSiberia", "Kazakh": "CentralAsiaSiberia",
    }
    for key, region in mapping.items():
        if key.lower() in population.lower():
            return _SGDP_REGION_LABELS.get(region, region)
    return "Unknown"


def _infer_sex_from_url(url: str) -> int:
    # SGDP CRAMs sometimes encode sex in the filename; fall back to 0 (unknown)
    name = url.lower()
    if "_male_" in name and "_female_" not in name:
        return 1
    if "_female_" in name:
        return 2
    return 0
