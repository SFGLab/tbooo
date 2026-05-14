"""EDA: sex/region demographics and GEUVADIS RNA PCA visualisations.

Outputs (data/metadata/eda/):
    sex_distribution.png
    region_distribution.png
    rna_pca_variance.png       scree + cumulative (requires map geuvadis)
    rna_pca_pc1_pc2.png        PC1 vs PC2 coloured by super-pop and sex (requires map geuvadis)
    eda_summary.tsv            counts table
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from tbooo.config import Config
from tbooo.utils import ensure_dirs, log

_SUPERPOP_COLORS = {
    "EUR": "#4C72B0",
    "AFR": "#DD8452",
    "EAS": "#55A868",
    "SAS": "#C44E52",
    "AMR": "#8172B3",
}
_SEX_COLORS = {"Male": "#4C72B0", "Female": "#DD8452", "Unknown": "#929292"}
_SEX_MAP = {1: "Male", 2: "Female", 0: "Unknown"}


def run_eda(cfg: Config) -> None:
    """Generate EDA plots + summary for demographics and GEUVADIS RNA PCA."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = cfg.metadata_dir() / "eda"
    ensure_dirs(out_dir)

    parquet_path = cfg.showcase_dir() / "participant.parquet"
    if not parquet_path.exists():
        raise FileNotFoundError(
            f"participant.parquet not found: {parquet_path}\n"
            "Run `tbooo map phenotypes` first."
        )
    part = pd.read_parquet(parquet_path)

    map1kg_path = cfg.metadata_dir() / "eid_map_1kg.tsv"
    map_sgdp_path = cfg.metadata_dir() / "eid_map_sgdp.tsv"
    map1kg = pd.read_csv(map1kg_path, sep="\t") if map1kg_path.exists() else None
    map_sgdp = pd.read_csv(map_sgdp_path, sep="\t") if map_sgdp_path.exists() else None

    _plot_sex(part, out_dir, plt)
    _plot_regions(map1kg, map_sgdp, out_dir, plt)

    variance_path = cfg.metadata_dir() / "geuvadis_pca_variance.tsv"
    variance = pd.read_csv(variance_path, sep="\t") if variance_path.exists() else None
    if variance is not None:
        _plot_scree(variance, out_dir, plt)
    else:
        log("  skip scree: geuvadis_pca_variance.tsv not found — run `tbooo map geuvadis`")

    pc_cols = [c for c in part.columns if c.startswith("geuvadis_pc")]
    if pc_cols and map1kg is not None:
        _plot_pca_scatter(part, pc_cols, map1kg, variance, out_dir, plt)
    else:
        log("  skip PCA scatter: geuvadis_pc* columns or eid_map_1kg.tsv not found")

    _write_summary(part, map1kg, map_sgdp, out_dir)
    log(f"EDA complete → {out_dir}/")


# ── helpers ───────────────────────────────────────────────────────────────────

def _style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=9)


def _plot_sex(part: pd.DataFrame, out_dir: Path, plt) -> None:
    counts = (
        part["p31"].map(_SEX_MAP)
        .value_counts()
        .reindex(["Male", "Female", "Unknown"])
        .fillna(0).astype(int)
    )
    fig, ax = plt.subplots(figsize=(4, 3.5))
    colors = [_SEX_COLORS[k] for k in counts.index]
    bars = ax.bar(counts.index, counts.values, color=colors, width=0.5, edgecolor="white")
    ax.bar_label(bars, padding=3, fontsize=9)
    ax.set_ylabel("Samples")
    ax.set_title("Sex distribution", fontsize=11, pad=8)
    _style(ax)
    fig.tight_layout()
    fig.savefig(out_dir / "sex_distribution.png", dpi=150)
    plt.close(fig)
    log("  sex_distribution.png")


def _plot_regions(map1kg, map_sgdp, out_dir: Path, plt) -> None:
    panels = []
    if map1kg is not None:
        panels.append(("1KGP: super-population", map1kg["super_pop"].value_counts(),
                        lambda k: _SUPERPOP_COLORS.get(k, "#888888"), False))
    if map_sgdp is not None:
        panels.append(("SGDP: region", map_sgdp["region"].value_counts(),
                        lambda _: "#7FB3D3", True))

    if not panels:
        log("  skip regions: no EID maps found — run `tbooo map eids`")
        return

    fig, axes = plt.subplots(1, len(panels), figsize=(4.5 * len(panels), 4))
    if len(panels) == 1:
        axes = [axes]

    for ax, (title, counts, color_fn, rotate) in zip(axes, panels):
        colors = [color_fn(k) for k in counts.index]
        bars = ax.bar(counts.index, counts.values, color=colors, width=0.5, edgecolor="white")
        ax.bar_label(bars, padding=3, fontsize=8)
        ax.set_ylabel("Samples")
        ax.set_title(title, fontsize=11, pad=8)
        if rotate:
            ax.tick_params(axis="x", rotation=35)
        _style(ax)

    fig.tight_layout()
    fig.savefig(out_dir / "region_distribution.png", dpi=150)
    plt.close(fig)
    log("  region_distribution.png")


def _plot_scree(variance: pd.DataFrame, out_dir: Path, plt) -> None:
    pcs = variance["pc"].str.replace("geuvadis_", "").str.upper()
    x = np.arange(len(pcs))

    fig, ax1 = plt.subplots(figsize=(5.5, 4))
    ax1.bar(x, variance["variance_explained"] * 100, color="#4C72B0", alpha=0.85, label="Per-PC")
    ax1.set_ylabel("Variance explained (%)")
    ax1.set_xlabel("Principal component")
    ax1.set_xticks(x)
    ax1.set_xticklabels(pcs, fontsize=8)
    ax1.spines[["top"]].set_visible(False)

    ax2 = ax1.twinx()
    ax2.plot(x, variance["cumulative_variance"] * 100, color="#DD8452",
             marker="o", markersize=5, linewidth=1.5, label="Cumulative")
    ax2.set_ylabel("Cumulative variance (%)")
    ax2.set_ylim(0, 105)
    ax2.spines[["top"]].set_visible(False)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, fontsize=9, loc="lower right")
    ax1.set_title("GEUVADIS RNA PCA — variance explained", fontsize=11, pad=8)
    fig.tight_layout()
    fig.savefig(out_dir / "rna_pca_variance.png", dpi=150)
    plt.close(fig)
    log("  rna_pca_variance.png")


def _plot_pca_scatter(part: pd.DataFrame, pc_cols: list[str],
                      map1kg: pd.DataFrame, variance, out_dir: Path, plt) -> None:
    # Build variance-label lookup: "PC1 (5.2%)"
    var_pct: dict[str, str] = {}
    if variance is not None:
        for _, row in variance.iterrows():
            label = row["pc"].replace("geuvadis_", "").upper()
            var_pct[row["pc"]] = f"{label} ({row['variance_explained']*100:.1f}%)"

    def _ax_label(col: str) -> str:
        return var_pct.get(col, col.replace("geuvadis_", "").upper())

    # Pairs to show: PC1v2, PC1v3, PC2v3
    n = len(pc_cols)
    pairs = [(pc_cols[0], pc_cols[1]),
             (pc_cols[0], pc_cols[min(2, n-1)]),
             (pc_cols[1], pc_cols[min(2, n-1)])]

    needed = [c for pair in pairs for c in pair]
    df = (
        part[["eid", "p31"] + list(dict.fromkeys(needed))]
        .dropna(subset=[pc_cols[0]])
        .copy()
    )
    df["sex_label"] = df["p31"].map(_SEX_MAP).fillna("Unknown")
    df = df.merge(map1kg[["eid", "super_pop"]], on="eid", how="left")

    colorings = [
        ("super-population", "super_pop",
         lambda sp: _SUPERPOP_COLORS.get(sp, "#888888"),
         sorted(df["super_pop"].dropna().unique())),
        ("sex", "sex_label",
         lambda s: _SEX_COLORS.get(s, "#929292"),
         [k for k in ("Male", "Female", "Unknown") if k in df["sex_label"].values]),
    ]

    fig, axes = plt.subplots(len(colorings), len(pairs),
                             figsize=(5 * len(pairs), 4.5 * len(colorings)))

    for row_idx, (title, color_col, color_fn, groups) in enumerate(colorings):
        for col_idx, (xc, yc) in enumerate(pairs):
            ax = axes[row_idx][col_idx]
            for grp_val in groups:
                mask = df[color_col] == grp_val
                ax.scatter(df.loc[mask, xc], df.loc[mask, yc],
                           c=color_fn(grp_val), s=22, alpha=0.75,
                           linewidths=0, label=grp_val)
            ax.set_xlabel(_ax_label(xc), fontsize=9)
            ax.set_ylabel(_ax_label(yc), fontsize=9)
            if col_idx == 0:
                ax.set_title(f"RNA PCA — {title}", fontsize=10, pad=6)
            legend = ax.legend(fontsize=7, markerscale=1.8, frameon=False,
                               loc="best", ncol=1)
            _style(ax)

    fig.tight_layout(h_pad=3, w_pad=2)
    fig.savefig(out_dir / "rna_pca_grid.png", dpi=150)
    plt.close(fig)
    log("  rna_pca_grid.png")


def _write_summary(part: pd.DataFrame, map1kg, map_sgdp, out_dir: Path) -> None:
    rows: list[dict] = []

    for code, label in _SEX_MAP.items():
        rows.append({"category": "sex", "group": label, "count": int((part["p31"] == code).sum())})

    if map1kg is not None:
        for sp, n in map1kg["super_pop"].value_counts().items():
            rows.append({"category": "1kg_super_pop", "group": sp, "count": int(n)})
        for pop, n in map1kg["pop"].value_counts().items():
            rows.append({"category": "1kg_pop", "group": pop, "count": int(n)})

    if map_sgdp is not None:
        for region, n in map_sgdp["region"].value_counts().items():
            rows.append({"category": "sgdp_region", "group": region, "count": int(n)})

    pd.DataFrame(rows).to_csv(out_dir / "eda_summary.tsv", sep="\t", index=False)
    log(f"  eda_summary.tsv")
