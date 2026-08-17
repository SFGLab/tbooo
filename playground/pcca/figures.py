"""Figures for the PCCA replication."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

POP_COLORS = {
    "GBR": "#4C72B0",
    "FIN": "#DD8452",
    "TSI": "#55A868",
    "YRI": "#C44E52",
    "CEU": "#8172B3",
}
POP_ORDER = ["GBR", "FIN", "TSI", "YRI", "CEU"]
METHOD_COLORS = {"pca": "#4C72B0", "pcca": "#DD8452"}

# Continental groups, for cohorts too diverse to colour by population: MAGE spans
# 26 populations, which no categorical palette renders legibly.
CONTINENTAL_COLORS = {
    "AFR": "#C44E52",
    "EUR": "#4C72B0",
    "EAS": "#55A868",
    "SAS": "#8172B3",
    "AMR": "#DD8452",
    "unknown": "#BBBBBB",
}
CONTINENTAL_ORDER = ["AFR", "EUR", "EAS", "SAS", "AMR", "unknown"]


def scatter_by_group(
    coords: np.ndarray,
    labels,
    title: str,
    xlabel: str,
    ylabel: str,
    path: Path,
    palette: dict[str, str] | None = None,
    order: list[str] | None = None,
) -> Path:
    """Scatter coloured by an arbitrary categorical label.

    The population-specific helper below is a thin wrapper on this; cohorts with
    many populations colour by continental group instead.
    """
    labels = pd.Series(labels).astype(str).reset_index(drop=True)
    if coords.shape[0] != len(labels):
        raise ValueError(
            f"coords has {coords.shape[0]} rows but labels has {len(labels)}"
        )
    palette = palette or CONTINENTAL_COLORS
    groups = order or [g for g in CONTINENTAL_ORDER if g in set(labels)]
    if not groups:
        groups = sorted(labels.unique())

    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for group in groups:
        mask = (labels == group).to_numpy()
        if not mask.any():
            continue
        ax.scatter(coords[mask, 0], coords[mask, 1],
                   c=palette.get(group, "#888888"), s=20, alpha=0.8,
                   linewidths=0, label=group)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=11, pad=8)
    ax.legend(fontsize=8, frameon=False, markerscale=1.5)
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _style(ax) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=9)


def scatter_by_population(
    coords: np.ndarray,
    cohort: pd.DataFrame,
    title: str,
    xlabel: str,
    ylabel: str,
    path: Path,
) -> Path:
    """Two-column scatter coloured by 1KGP population."""
    if coords.shape[0] != len(cohort):
        raise ValueError(
            f"coords has {coords.shape[0]} rows but cohort has {len(cohort)}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for pop in POP_ORDER:
        mask = (cohort["pop"] == pop).to_numpy()
        if not mask.any():
            continue
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=POP_COLORS[pop],
            s=24,
            alpha=0.8,
            linewidths=0,
            label=pop,
        )
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=11, pad=8)
    ax.legend(fontsize=8, frameon=False, markerscale=1.5)
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def qq_plot(pvalue_sets: dict[str, np.ndarray], path: Path) -> Path:
    """Observed vs expected -log10(p) under the uniform null.

    The paper's claim is that the PCCA covariates reduce inflation at the high
    end, so the interesting region is the top right against the diagonal.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5, 4.5))
    x_limit = 0.0
    for name, pvals in pvalue_sets.items():
        pvals = np.asarray(pvals)
        pvals = pvals[np.isfinite(pvals) & (pvals > 0)]
        observed = -np.log10(np.sort(pvals))[::-1]
        expected = -np.log10((np.arange(len(observed), 0, -1) - 0.5) / len(observed))
        ax.scatter(
            expected,
            observed,
            s=6,
            alpha=0.6,
            linewidths=0,
            label=name,
            c=METHOD_COLORS.get(name),
        )
        x_limit = max(x_limit, expected.max())
    # The null caps expected at -log10(1/n), so squaring the axes would waste
    # most of the panel. Bound x by the null and let y follow the observed tail.
    ax.plot([0, x_limit], [0, x_limit], color="#666666", linewidth=1, linestyle="--")
    ax.set_xlim(0, x_limit * 1.05)
    ax.set_xlabel("expected $-\\log_{10}(p)$", fontsize=9)
    ax.set_ylabel("observed $-\\log_{10}(p)$", fontsize=9)
    ax.set_title("Panel B — eQTL p-value inflation", fontsize=11, pad=8)
    ax.legend(fontsize=8, frameon=False)
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def egene_curve(scans: dict[str, pd.DataFrame], alphas, path: Path) -> Path:
    """Number of eGenes as a function of the significance cutoff."""
    path.parent.mkdir(parents=True, exist_ok=True)
    alphas = np.asarray(alphas)
    fig, ax = plt.subplots(figsize=(5, 4.5))
    for name, scan in scans.items():
        counts = [(scan["best_p"] < a).sum() for a in alphas]
        ax.plot(
            alphas,
            counts,
            marker="o",
            markersize=4,
            linewidth=1.5,
            label=name,
            color=METHOD_COLORS.get(name),
        )
    ax.set_xscale("log")
    ax.set_xlabel("significance cutoff", fontsize=9)
    ax.set_ylabel("genes with an eQTL", fontsize=9)
    ax.set_title("Panel C — eGenes vs cutoff", fontsize=11, pad=8)
    ax.legend(fontsize=8, frameon=False)
    _style(ax)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def egene_venn(sets_by_method: dict[str, set], alpha: float, path: Path):
    """Two-set Venn of eGene identities, drawn with plain circles.

    matplotlib-venn is not installed and this is a two-set diagram, so the
    circles are drawn directly rather than adding a dependency.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    (name_a, set_a), (name_b, set_b) = list(sets_by_method.items())
    counts = {
        f"{name_a}_only": len(set_a - set_b),
        "shared": len(set_a & set_b),
        f"{name_b}_only": len(set_b - set_a),
    }

    fig, ax = plt.subplots(figsize=(5.5, 4))
    for x, colour, label in (
        (-0.35, METHOD_COLORS.get(name_a, "#4C72B0"), name_a),
        (0.35, METHOD_COLORS.get(name_b, "#DD8452"), name_b),
    ):
        ax.add_patch(plt.Circle((x, 0), 0.75, alpha=0.45, color=colour, linewidth=0))
        ax.text(x * 2.0, -1.05, label, ha="center", fontsize=9)
    ax.text(-0.75, 0, str(counts[f"{name_a}_only"]), ha="center", va="center", fontsize=10)
    ax.text(0.0, 0, str(counts["shared"]), ha="center", va="center", fontsize=10)
    ax.text(0.75, 0, str(counts[f"{name_b}_only"]), ha="center", va="center", fontsize=10)
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.3, 1.1)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(f"Panel D — eGenes at $\\alpha$ = {alpha:g}", fontsize=11, pad=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path, counts
