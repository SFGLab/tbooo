import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from pcca.figures import (  # noqa: E402
    POP_COLORS,
    egene_curve,
    egene_venn,
    qq_plot,
    scatter_by_population,
)


@pytest.fixture
def cohort():
    return pd.DataFrame(
        {
            "sample_id": [f"S{i}" for i in range(8)],
            "pop": ["GBR", "FIN", "TSI", "YRI"] * 2,
        }
    )


def test_every_geuvadis_population_has_a_colour():
    assert set(POP_COLORS) >= {"GBR", "FIN", "TSI", "YRI", "CEU"}


def test_scatter_writes_a_png(tmp_path, cohort):
    coords = np.arange(16, dtype=float).reshape(8, 2)
    out = scatter_by_population(coords, cohort, "t", "x", "y", tmp_path / "f.png")
    assert out.exists() and out.stat().st_size > 1000


def test_scatter_rejects_mismatched_lengths(tmp_path, cohort):
    """Row misalignment between coords and cohort is the silent-failure mode."""
    coords = np.zeros((3, 2))
    with pytest.raises(ValueError, match="rows"):
        scatter_by_population(coords, cohort, "t", "x", "y", tmp_path / "f.png")


def test_qq_plot_writes_a_png(tmp_path):
    rng = np.random.default_rng(0)
    sets = {"pca": rng.random(500), "pcca": rng.random(500)}
    out = qq_plot(sets, tmp_path / "qq.png")
    assert out.exists() and out.stat().st_size > 1000


def test_egene_curve_writes_a_png(tmp_path):
    scans = {
        "pca": pd.DataFrame({"best_p": np.logspace(-12, -1, 200)}),
        "pcca": pd.DataFrame({"best_p": np.logspace(-11, -1, 200)}),
    }
    out = egene_curve(scans, np.logspace(-9, -3, 7), tmp_path / "curve.png")
    assert out.exists() and out.stat().st_size > 1000


def test_egene_venn_writes_a_png_and_reports_counts(tmp_path):
    sets = {"pca": {"A", "B", "C"}, "pcca": {"B", "C", "D"}}
    out, counts = egene_venn(sets, alpha=1e-6, path=tmp_path / "venn.png")
    assert out.exists()
    assert counts == {"pca_only": 1, "shared": 2, "pcca_only": 1}
