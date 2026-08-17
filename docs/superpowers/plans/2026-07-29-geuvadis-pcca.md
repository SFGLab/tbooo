# GEUVADIS PCCA Replication Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replicate Brown, Bray & Pachter (PLOS Genetics 2018) on GEUVADIS — show that plain PCA on expression fails to cluster by population while PCCA (expression PCs coupled to genotype PCs via CCA) succeeds — and produce the four-panel figure plus a three-paper methodology document.

**Architecture:** Testable analysis functions live in a `playground/pcca/` package; `playground/geuvadis_pcca.ipynb` orchestrates them and emits figures. The shipped `tbooo` pipeline is touched in exactly one place: the `geuvadis` stage is corrected to read raw RPKM instead of PEER-residualized expression.

**Tech Stack:** Python 3 (`/opt/anaconda3/bin/python`), numpy 1.26, pandas 1.5.3, scipy 1.17, scikit-learn 1.8, matplotlib, pytest 9.0.3, plink2 (to be installed).

## Global Constraints

- **NEVER run `git commit`.** The user commits. No task in this plan ends with a commit step. This is deliberate and overrides the writing-plans skill default.
- **Python is `/opt/anaconda3/bin/python`.** The system `python3` has no pandas and will fail.
- **Disk headroom is ~22 GiB on a 98%-full volume.** Every genotype step processes one chromosome at a time and deletes decompressed intermediates before starting the next. Never decompress a whole-cohort callset at once.
- **Spec:** `docs/superpowers/specs/2026-07-29-geuvadis-pcca-design.md`. Read it before starting.
- **All coordinates are GRCh37.** GD660 uses GENCODE v12 on GRCh37; plink2 `all_phase3` is build 37. Do not mix in GRCh38 annotations.
- **Reference values to hit:** 343 individuals; ≈14k genes after filtering; first two canonical correlations 0.963 and 0.766; eGenes at α=1e-6 of 2,818 (PCA) vs 2,732 (PCCA).
- **Run tests with:** `/opt/anaconda3/bin/python -m pytest tests/ -v`

## File Structure

| File | Responsibility |
|---|---|
| `playground/pcca/__init__.py` | Package marker |
| `playground/pcca/paths.py` | Canonical data paths and download URLs |
| `playground/pcca/metadata.py` | SDRF parsing, cohort selection, covariate table |
| `playground/pcca/expression.py` | GD660 loading, gene filtering, normalization, confounder regression |
| `playground/pcca/core.py` | PCA, orthonormalization, CCA, leave-one-out projection |
| `playground/pcca/genotype.py` | plink2 driver: per-chromosome filter, merge, PCA, dosage export |
| `playground/pcca/eqtl.py` | Vectorized cis-eQTL scan |
| `playground/pcca/figures.py` | The five figures |
| `playground/geuvadis_pcca.ipynb` | Orchestration and narrative |
| `tests/conftest.py` | Puts `playground/` on `sys.path` |
| `tests/test_*.py` | One test module per source module |
| `src/tbooo/download/geuvadis.py` | Modified: also fetch GD660 and the SDRF |
| `src/tbooo/pipeline/geuvadis.py` | Modified: read GD660 raw RPKM |
| `docs/4_pcca_methods.md` | Three-paper methodology comparison |

---

### Task 1: Package skeleton, paths, and SDRF cohort selection

This is the task that proves the cohort is reproducible. Everything downstream keys off it.

**Files:**
- Create: `playground/pcca/__init__.py`, `playground/pcca/paths.py`, `playground/pcca/metadata.py`
- Create: `tests/conftest.py`, `tests/test_metadata.py`

**Interfaces:**
- Produces: `load_sdrf(path) -> pd.DataFrame` with columns `sample_id, assay_name, lab, sex, ancestry, pop, has_phase1`, one row per individual (462 rows). `replication_cohort(sdrf) -> pd.DataFrame` (343 rows). `production_cohort(sdrf) -> pd.DataFrame` (462 rows). `build_confounders(meta) -> np.ndarray` of shape `(n, 8)`.

- [ ] **Step 1: Download the SDRF**

```bash
cd /Users/nk2/Desktop/biobank/tbooo
curl -o data/raw/geuvadis/E-GEUV-1.sdrf.txt \
  "https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1/E-GEUV-1.sdrf.txt"
wc -l data/raw/geuvadis/E-GEUV-1.sdrf.txt
```

Expected: `925 data/raw/geuvadis/E-GEUV-1.sdrf.txt` (1 header + 924 rows; 2 rows per individual for paired-end FASTQs).

- [ ] **Step 2: Create the test bootstrap**

`tests/conftest.py`:

```python
"""Put playground/ on sys.path so tests can `import pcca`."""
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "playground"))
```

- [ ] **Step 3: Write the failing test**

`tests/test_metadata.py`:

```python
import numpy as np
import pytest

from pcca import paths
from pcca.metadata import (
    build_confounders,
    load_sdrf,
    production_cohort,
    replication_cohort,
)

pytestmark = pytest.mark.skipif(
    not paths.SDRF.exists(), reason="SDRF not downloaded; see Task 1 Step 1"
)


@pytest.fixture(scope="module")
def sdrf():
    return load_sdrf(paths.SDRF)


def test_one_row_per_individual(sdrf):
    assert len(sdrf) == 462
    assert sdrf["sample_id"].is_unique


def test_seven_labs(sdrf):
    assert sorted(sdrf["lab"].unique()) == ["1", "2", "3", "4", "5", "6", "7"]


def test_population_codes(sdrf):
    counts = sdrf["pop"].value_counts().to_dict()
    assert counts == {"FIN": 95, "GBR": 94, "TSI": 93, "CEU": 91, "YRI": 89}


def test_replication_cohort_is_the_papers_343(sdrf):
    """Non-CEU individuals with Phase 1 genotypes reproduce the paper's cohort exactly."""
    cohort = replication_cohort(sdrf)
    assert len(cohort) == 343
    assert cohort["pop"].value_counts().to_dict() == {
        "FIN": 89, "TSI": 92, "GBR": 85, "YRI": 77,
    }
    assert "CEU" not in set(cohort["pop"])


def test_production_cohort_is_all_462(sdrf):
    assert len(production_cohort(sdrf)) == 462


def test_confounders_shape_and_rank(sdrf):
    """Intercept + 6 lab dummies (7 labs, first dropped) + sex = 8 columns."""
    cohort = replication_cohort(sdrf)
    conf = build_confounders(cohort)
    assert conf.shape == (343, 8)
    assert np.linalg.matrix_rank(conf) == 8
    assert np.allclose(conf[:, 0], 1.0)
```

- [ ] **Step 4: Run the test to verify it fails**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_metadata.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pcca'`

- [ ] **Step 5: Create the package and paths module**

`playground/pcca/__init__.py`:

```python
"""GEUVADIS PCCA replication helpers.

Analysis code lives here rather than inline in the notebook so it can be tested.
See docs/superpowers/specs/2026-07-29-geuvadis-pcca-design.md
"""
```

`playground/pcca/paths.py`:

```python
"""Canonical data locations and source URLs for the PCCA replication."""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

DATA = REPO_ROOT / "data"
GEUVADIS_RAW = DATA / "raw" / "geuvadis"
METADATA = DATA / "metadata"
FIGURES = METADATA / "eda"
GENOTYPE = DATA / "raw" / "geuvadis" / "genotype"

SDRF = GEUVADIS_RAW / "E-GEUV-1.sdrf.txt"
GD660_GENE_RPKM = GEUVADIS_RAW / "GD660.GeneQuantRPKM.txt.gz"
GD660_TR_RPKM = GEUVADIS_RAW / "GD660.TrQuantRPKM.txt.gz"

_EBI = (
    "https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1"
)
SDRF_URL = f"{_EBI}/E-GEUV-1.sdrf.txt"
GD660_GENE_RPKM_URL = f"{_EBI}/analysis_results/GD660.GeneQuantRPKM.txt.gz"
GD660_TR_RPKM_URL = f"{_EBI}/analysis_results/GD660.TrQuantRPKM.txt.gz"

# GEUVADIS also hosts the paper's own Phase 1 genotypes under E-GEUV-1/genotypes/,
# but they total ~50 GB and are deliberately not used here. See the Deviations
# section of the spec: Phase 3 is the genotype source for this replication.

# plink2 prebuilt 1000 Genomes Phase 3 (build 37, 2504 samples). The per-chromosome
# files are hosted on Dropbox behind opaque hashes that change between releases, so
# the URLs are scraped from the resources page at runtime rather than hardcoded.
# Files are named chr{N}_phase3.{pgen,pvar}.zst plus phase3_orig.psam.
PLINK2_RESOURCES_URL = "https://www.cog-genomics.org/plink/2.0/resources"
PLINK2_MAC_ARM64_URL = (
    "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_mac_arm64_20260504.zip"
)
```

- [ ] **Step 6: Implement metadata.py**

`playground/pcca/metadata.py`:

```python
"""Cohort definition from the GEUVADIS SDRF.

The SDRF is the only source of the 7-lab assignment: the GD660 column headers
look like `HG00096.1.M_111124_6`, but all 660 runs share the prefix `M`, which
denotes the sequencing platform rather than the lab.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

ANCESTRY_TO_POP = {
    "British": "GBR",
    "Finnish": "FIN",
    "Tuscan": "TSI",
    "Utah": "CEU",
    "Yoruba": "YRI",
}


def load_sdrf(path) -> pd.DataFrame:
    """One row per individual, with the covariates the replication needs.

    The raw SDRF has two rows per individual (one per paired-end FASTQ);
    they are identical in every column used here.
    """
    raw = pd.read_csv(path, sep="\t", dtype=str)
    df = pd.DataFrame(
        {
            "sample_id": raw["Source Name"],
            "assay_name": raw["Assay Name"],
            "lab": raw["Characteristics[laboratory]"],
            "sex": raw["Characteristics[sex]"],
            "ancestry": raw["Characteristics[ancestry category]"],
            "has_phase1": raw["Comment[1000g Phase1 Genotypes]"] == "1",
        }
    )
    df["pop"] = df["ancestry"].map(ANCESTRY_TO_POP)
    if df["pop"].isna().any():
        unknown = sorted(df.loc[df["pop"].isna(), "ancestry"].unique())
        raise ValueError(f"Unmapped ancestry categories: {unknown}")
    return (
        df.drop_duplicates(subset="sample_id")
        .sort_values("sample_id")
        .reset_index(drop=True)
    )


def replication_cohort(sdrf: pd.DataFrame) -> pd.DataFrame:
    """The paper's exact 343 individuals: non-CEU, with Phase 1 genotypes.

    CEU is excluded in the paper for biased expression patterns.
    """
    mask = (sdrf["pop"] != "CEU") & sdrf["has_phase1"]
    return sdrf[mask].reset_index(drop=True)


def production_cohort(sdrf: pd.DataFrame) -> pd.DataFrame:
    """All 462 individuals including CEU — the input to the later fusion work."""
    return sdrf.reset_index(drop=True)


def build_confounders(meta: pd.DataFrame) -> np.ndarray:
    """Design matrix of intercept + 6 lab dummies + sex, shape (n, 8).

    One lab dummy is dropped to keep the matrix full rank alongside the intercept.
    """
    labs = sorted(meta["lab"].unique())
    dummies = np.column_stack(
        [(meta["lab"] == lab).to_numpy(dtype=float) for lab in labs[1:]]
    )
    sex = (meta["sex"] == "male").to_numpy(dtype=float)
    intercept = np.ones(len(meta))
    return np.column_stack([intercept, dummies, sex])
```

- [ ] **Step 7: Run the tests to verify they pass**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_metadata.py -v`
Expected: PASS, 6 passed.

If `test_replication_cohort_is_the_papers_343` fails, stop and investigate before continuing — every downstream number depends on the cohort being right.

---

### Task 2: GD660 loading and gene filtering

**Files:**
- Create: `playground/pcca/expression.py`, `tests/test_expression.py`
- Modify: `src/tbooo/download/geuvadis.py`

**Interfaces:**
- Consumes: `pcca.paths`, `pcca.metadata.load_sdrf`
- Produces: `load_gd660(path) -> (ann: pd.DataFrame, expr: pd.DataFrame)` where `ann` has columns `TargetID, Gene_Symbol, Chr, Coord` and `expr` is genes × 660 runs. `keep_gene(chrom: str, coord: int) -> bool`. `filter_genes(ann, expr, min_mean_rpkm=0.1) -> (ann, expr)`.

- [ ] **Step 1: Add GD660 and the SDRF to the downloader**

Modify `src/tbooo/download/geuvadis.py`. Replace the module docstring and the two URL constants with a table of files, and rewrite `download_expression` to loop:

```python
"""Download GEUVADIS gene expression data and sample metadata.

Downloads:
  - GD660 raw gene-level RPKM (660 runs, incl. technical replicates). This is the
    input to all population-structure work: it is unnormalized, so it retains the
    structure that the PEER-residualized GD462 file removes.
  - GD462 PEER-residualized RPKM. Retained because it is the correct input for
    eQTL work that wants PEER correction; no pipeline stage reads it.
  - E-GEUV-1 SDRF: the only source of the 7-lab assignment, plus sex, ancestry
    and the Phase 1 genotype availability flag.

Outputs:
    data/raw/geuvadis/GD660.GeneQuantRPKM.txt.gz
    data/raw/geuvadis/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz
    data/raw/geuvadis/E-GEUV-1.sdrf.txt
"""

from __future__ import annotations

from tbooo.config import Config
from tbooo.utils import ensure_dirs, log, wget_download

_EBI = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/GEUV/E-GEUV-1"

_FILES = {
    "GD660.GeneQuantRPKM.txt.gz": f"{_EBI}/analysis_results/GD660.GeneQuantRPKM.txt.gz",
    "GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz": (
        f"{_EBI}/analysis_results/GD462.GeneQuantRPKM.50FN.samplename.resk10.txt.gz"
    ),
    "E-GEUV-1.sdrf.txt": f"{_EBI}/E-GEUV-1.sdrf.txt",
}


def download_expression(cfg: Config) -> None:
    """Download the GEUVADIS quantification matrices and sample metadata."""
    out_dir = cfg.geuvadis_raw_dir()
    ensure_dirs(out_dir)
    for name, url in _FILES.items():
        dest = out_dir / name
        if dest.exists():
            log(f"  skip (exists): {name}")
            continue
        log(f"Downloading {name}…")
        wget_download(url, dest, tool_wget=cfg.tools.wget)
        log(f"  wrote {dest}")
```

- [ ] **Step 2: Run the downloader**

```bash
cd /Users/nk2/Desktop/biobank/tbooo
/opt/anaconda3/bin/tbooo download geuvadis
ls -la data/raw/geuvadis/
```

Expected: `GD660.GeneQuantRPKM.txt.gz` present (~250 MB). The SDRF and GD462 file are skipped as already present.

- [ ] **Step 3: Write the failing test**

`tests/test_expression.py`:

```python
import numpy as np
import pandas as pd
import pytest

from pcca import paths
from pcca.expression import filter_genes, keep_gene, load_gd660


def test_keep_gene_accepts_autosomes():
    assert keep_gene("1", 1_000_000)
    assert keep_gene("22", 500_000)


def test_keep_gene_rejects_sex_and_mito():
    assert not keep_gene("X", 1_000_000)
    assert not keep_gene("Y", 1_000_000)
    assert not keep_gene("MT", 1_000)


def test_keep_gene_rejects_mhc_window_on_chr6_only():
    assert not keep_gene("6", 30_000_000)     # inside MHC
    assert keep_gene("6", 28_999_999)         # just before
    assert keep_gene("6", 33_000_001)         # just after
    assert keep_gene("5", 30_000_000)         # same coord, different chromosome


def test_filter_genes_drops_low_expression_and_bad_locations():
    ann = pd.DataFrame(
        {
            "TargetID": ["keep", "low", "mhc", "chrX"],
            "Gene_Symbol": ["keep", "low", "mhc", "chrX"],
            "Chr": ["1", "1", "6", "X"],
            "Coord": [1_000_000, 2_000_000, 30_000_000, 1_000_000],
        }
    )
    expr = pd.DataFrame(
        np.array([[5.0, 5.0], [0.01, 0.01], [5.0, 5.0], [5.0, 5.0]]),
        columns=["s1", "s2"],
    )
    ann_f, expr_f = filter_genes(ann, expr, min_mean_rpkm=0.1)
    assert list(ann_f["TargetID"]) == ["keep"]
    assert expr_f.shape == (1, 2)


@pytest.mark.skipif(
    not paths.GD660_GENE_RPKM.exists(), reason="GD660 not downloaded; see Task 2 Step 2"
)
def test_load_gd660_real_file():
    ann, expr = load_gd660(paths.GD660_GENE_RPKM)
    assert expr.shape[1] == 660
    assert len(ann) == len(expr) == 23_722
    assert list(ann.columns) == ["TargetID", "Gene_Symbol", "Chr", "Coord"]
    assert expr.columns[0] == "HG00096.1.M_111124_6"


@pytest.mark.skipif(
    not paths.GD660_GENE_RPKM.exists(), reason="GD660 not downloaded; see Task 2 Step 2"
)
def test_filter_genes_real_file_yields_about_14k():
    ann, expr = load_gd660(paths.GD660_GENE_RPKM)
    ann_f, _ = filter_genes(ann, expr)
    assert 12_000 <= len(ann_f) <= 16_000, f"got {len(ann_f)} genes, paper reports 14,070"
```

- [ ] **Step 4: Run the test to verify it fails**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_expression.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pcca.expression'`

- [ ] **Step 5: Implement the loading and filtering half of expression.py**

`playground/pcca/expression.py`:

```python
"""GD660 expression matrix construction.

Pipeline: load raw RPKM -> filter genes -> quantile normalize -> standardize
-> regress out lab and sex. The result is the expression side of PCCA.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

ANNOTATION_COLUMNS = ["TargetID", "Gene_Symbol", "Chr", "Coord"]
AUTOSOMES = {str(i) for i in range(1, 23)}
MHC_START, MHC_END = 29_000_000, 33_000_000
MIN_MEAN_RPKM = 0.1


def load_gd660(path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (annotation, expression) with expression as genes x 660 runs."""
    df = pd.read_csv(path, sep="\t", compression="gzip")
    ann = df[ANNOTATION_COLUMNS].copy()
    ann["Chr"] = ann["Chr"].astype(str)
    ann["Coord"] = ann["Coord"].astype(int)
    expr = df.drop(columns=ANNOTATION_COLUMNS)
    return ann.reset_index(drop=True), expr.reset_index(drop=True)


def keep_gene(chrom: str, coord: int) -> bool:
    """Autosomal and outside the MHC, matching _keep_gene in pachterlab/PCCA.

    GD660 gives only a start coordinate, so MHC membership is tested on the start
    site rather than on interval overlap.
    """
    if chrom not in AUTOSOMES:
        return False
    if chrom == "6" and MHC_START <= coord <= MHC_END:
        return False
    return True


def filter_genes(
    ann: pd.DataFrame, expr: pd.DataFrame, min_mean_rpkm: float = MIN_MEAN_RPKM
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drop weakly expressed, non-autosomal, and MHC genes."""
    expressed = expr.mean(axis=1).to_numpy() >= min_mean_rpkm
    located = np.array(
        [keep_gene(c, p) for c, p in zip(ann["Chr"], ann["Coord"])], dtype=bool
    )
    keep = expressed & located
    return (
        ann[keep].reset_index(drop=True),
        expr[keep].reset_index(drop=True),
    )
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_expression.py -v`
Expected: PASS, 6 passed. Note the actual gene count printed by the last test; the paper's figure is 14,070.

---

### Task 3: Normalization and confounder regression

**Files:**
- Modify: `playground/pcca/expression.py`
- Modify: `tests/test_expression.py`

**Interfaces:**
- Consumes: `filter_genes` from Task 2, `build_confounders` from Task 1
- Produces: `quantile_normalize(mat: np.ndarray) -> np.ndarray` (genes × samples in, same shape out). `standardize_genes(mat) -> np.ndarray`. `regress_out(X, C) -> np.ndarray` (X is samples × genes). `build_expression_matrix(ann, expr, cohort) -> (ann, np.ndarray)` returning **samples × genes**.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_expression.py`:

```python
from pcca.expression import (
    build_expression_matrix,
    quantile_normalize,
    regress_out,
    standardize_genes,
)


def test_quantile_normalize_gives_columns_identical_distributions():
    mat = np.array([[5.0, 4.0], [2.0, 1.0], [3.0, 9.0]])
    out = quantile_normalize(mat)
    assert np.allclose(np.sort(out[:, 0]), np.sort(out[:, 1]))


def test_quantile_normalize_preserves_within_column_ranks():
    mat = np.array([[5.0, 4.0], [2.0, 1.0], [3.0, 9.0]])
    out = quantile_normalize(mat)
    for j in range(mat.shape[1]):
        assert list(np.argsort(mat[:, j])) == list(np.argsort(out[:, j]))


def test_quantile_normalize_preserves_shape():
    rng = np.random.default_rng(0)
    mat = rng.random((50, 7))
    assert quantile_normalize(mat).shape == (50, 7)


def test_standardize_genes_zero_mean_unit_variance_per_row():
    rng = np.random.default_rng(1)
    mat = rng.normal(loc=10.0, scale=3.0, size=(20, 100))
    out = standardize_genes(mat)
    assert np.allclose(out.mean(axis=1), 0.0, atol=1e-12)
    assert np.allclose(out.std(axis=1), 1.0, atol=1e-12)


def test_standardize_genes_handles_constant_rows():
    mat = np.array([[3.0, 3.0, 3.0], [1.0, 2.0, 3.0]])
    out = standardize_genes(mat)
    assert np.all(np.isfinite(out))
    assert np.allclose(out[0], 0.0)


def test_regress_out_removes_the_confounder_signal():
    """A gene that is a pure function of the confounder becomes ~zero."""
    rng = np.random.default_rng(2)
    n = 100
    group = rng.integers(0, 2, size=n).astype(float)
    C = np.column_stack([np.ones(n), group])
    gene_driven_by_group = 5.0 * group
    gene_independent = rng.normal(size=n)
    X = np.column_stack([gene_driven_by_group, gene_independent])
    resid = regress_out(X, C)
    assert np.allclose(resid[:, 0], 0.0, atol=1e-10)
    assert np.corrcoef(resid[:, 1], gene_independent)[0, 1] > 0.9


def test_regress_out_leaves_residuals_orthogonal_to_confounders():
    rng = np.random.default_rng(3)
    n = 60
    C = np.column_stack([np.ones(n), rng.normal(size=n)])
    X = rng.normal(size=(n, 10))
    resid = regress_out(X, C)
    assert np.allclose(C.T @ resid, 0.0, atol=1e-9)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_expression.py -v -k "quantile or standardize or regress"`
Expected: FAIL — `ImportError: cannot import name 'quantile_normalize'`

- [ ] **Step 3: Implement the normalization functions**

Append to `playground/pcca/expression.py`:

```python
def quantile_normalize(mat: np.ndarray) -> np.ndarray:
    """Force every column (sample) to share the same distribution.

    Equivalent to the reference implementation's pandas version, but done in
    numpy to avoid stack/unstack row-ordering surprises. Ties break by first
    occurrence, matching rank(method='first').
    """
    order = np.argsort(mat, axis=0, kind="stable")
    ranks = np.argsort(order, axis=0, kind="stable")
    reference = np.sort(mat, axis=0).mean(axis=1)
    return reference[ranks]


def standardize_genes(mat: np.ndarray) -> np.ndarray:
    """Center and scale each gene (row) to mean 0, variance 1.

    Constant genes would divide by zero; they are left at zero instead, which is
    the correct centered value and carries no signal either way.
    """
    mean = mat.mean(axis=1, keepdims=True)
    sd = mat.std(axis=1, keepdims=True)
    sd = np.where(sd == 0, 1.0, sd)
    return (mat - mean) / sd


def regress_out(X: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Residuals of samples x genes matrix X after regressing on design C."""
    beta, *_ = np.linalg.lstsq(C, X, rcond=None)
    return X - C @ beta
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_expression.py -v -k "quantile or standardize or regress"`
Expected: PASS, 7 passed.

- [ ] **Step 5: Write the failing test for the end-to-end builder**

Append to `tests/test_expression.py`:

```python
@pytest.mark.skipif(
    not (paths.GD660_GENE_RPKM.exists() and paths.SDRF.exists()),
    reason="GEUVADIS inputs not downloaded",
)
def test_build_expression_matrix_replication_cohort():
    from pcca.metadata import load_sdrf, replication_cohort

    cohort = replication_cohort(load_sdrf(paths.SDRF))
    ann, expr = load_gd660(paths.GD660_GENE_RPKM)
    ann_f, X = build_expression_matrix(ann, expr, cohort)

    assert X.shape[0] == 343, "one row per individual in the paper's cohort"
    assert X.shape[1] == len(ann_f)
    assert np.all(np.isfinite(X))
    # Residuals are orthogonal to the confounder design by construction.
    from pcca.metadata import build_confounders

    assert np.allclose(build_confounders(cohort).T @ X, 0.0, atol=1e-8)
```

- [ ] **Step 6: Run to verify it fails**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_expression.py -v -k build_expression_matrix`
Expected: FAIL — `ImportError: cannot import name 'build_expression_matrix'`

- [ ] **Step 7: Implement build_expression_matrix**

Append to `playground/pcca/expression.py`:

```python
def build_expression_matrix(
    ann: pd.DataFrame, expr: pd.DataFrame, cohort: pd.DataFrame
) -> tuple[pd.DataFrame, np.ndarray]:
    """Full expression side of PCCA for one cohort.

    Selects the cohort's nonredundant runs (the SDRF lists exactly one run per
    individual, so the 198 technical replicates in GD660 are simply not selected),
    filters genes, quantile normalizes, standardizes, and regresses out lab and sex.

    Returns the surviving annotation and a samples x genes residual matrix whose
    row order matches `cohort`.
    """
    from pcca.metadata import build_confounders

    missing = set(cohort["assay_name"]) - set(expr.columns)
    if missing:
        raise KeyError(f"{len(missing)} cohort assay names absent from GD660, e.g. {sorted(missing)[:3]}")

    selected = expr[list(cohort["assay_name"])]
    ann_f, expr_f = filter_genes(ann, selected)

    mat = quantile_normalize(expr_f.to_numpy(dtype=float))
    mat = standardize_genes(mat)

    X = mat.T  # samples x genes
    X = regress_out(X, build_confounders(cohort))
    return ann_f, X
```

- [ ] **Step 8: Run the whole expression suite**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_expression.py -v`
Expected: PASS, 14 passed.

---

### Task 4: PCCA core — PCA, CCA, and the leave-one-out Gram trick

The correctness of the whole replication rests on this module. The orthonormality requirement in Step 5 is the single easiest thing to get silently wrong.

**Files:**
- Create: `playground/pcca/core.py`, `tests/test_core.py`

**Interfaces:**
- Produces: `pca(arr, dim) -> (u, singular_values)` where `arr` is samples × features and `u` is samples × dim with orthonormal columns. `pca_from_gram(gram, dim) -> (u, singular_values)`. `orthonormalize(mat) -> np.ndarray`. `cca(u_x, u_y, dim) -> (l_x, l_y, rho)`. `loo_projection(X, gram_y, dim_x, dim_y, dim_z) -> np.ndarray` of shape (n, dim_z).

- [ ] **Step 1: Write the failing tests**

`tests/test_core.py`:

```python
import numpy as np
import pytest

from pcca.core import cca, loo_projection, orthonormalize, pca, pca_from_gram


def test_pca_returns_orthonormal_columns():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(50, 200))
    u, _ = pca(X, dim=5)
    assert u.shape == (50, 5)
    assert np.allclose(u.T @ u, np.eye(5), atol=1e-10)


def test_pca_recovers_a_planted_low_rank_structure():
    """Two well-separated clusters should be visible in the first component."""
    rng = np.random.default_rng(1)
    n = 60
    signal = np.concatenate([np.ones(n // 2), -np.ones(n // 2)])
    loadings = rng.normal(size=100)
    X = np.outer(signal, loadings) + 0.01 * rng.normal(size=(n, 100))
    u, _ = pca(X, dim=1)
    assert abs(np.corrcoef(u[:, 0], signal)[0, 1]) > 0.99


def test_pca_from_gram_matches_pca():
    """The Gram-matrix route is exact, which is what makes leave-one-out cheap."""
    rng = np.random.default_rng(2)
    X = rng.normal(size=(40, 300))
    u_direct, s_direct = pca(X, dim=6)
    u_gram, s_gram = pca_from_gram(X @ X.T, dim=6)
    assert np.allclose(np.abs(u_direct), np.abs(u_gram), atol=1e-8)
    assert np.allclose(s_direct, s_gram, atol=1e-8)


def test_pca_from_gram_on_a_submatrix_equals_leaving_that_sample_out():
    """This equivalence is the entire justification for the LOO shortcut."""
    rng = np.random.default_rng(3)
    X = rng.normal(size=(30, 200))
    gram = X @ X.T
    hold_out = 7
    idx = np.array([i for i in range(30) if i != hold_out])
    u_direct, _ = pca(X[idx], dim=4)
    u_sub, _ = pca_from_gram(gram[np.ix_(idx, idx)], dim=4)
    assert np.allclose(np.abs(u_direct), np.abs(u_sub), atol=1e-8)


def test_orthonormalize_fixes_non_unit_columns():
    """plink2 .eigenvec columns are orthogonal but not unit-norm."""
    rng = np.random.default_rng(4)
    q, _ = np.linalg.qr(rng.normal(size=(20, 3)))
    scaled = q * np.array([5.0, 0.2, 13.0])
    assert not np.allclose(scaled.T @ scaled, np.eye(3))
    out = orthonormalize(scaled)
    assert np.allclose(out.T @ out, np.eye(3), atol=1e-10)


def test_cca_of_identical_subspaces_gives_unit_correlations():
    rng = np.random.default_rng(5)
    u, _ = np.linalg.qr(rng.normal(size=(40, 5)))
    _, _, rho = cca(u, u, dim=3)
    assert np.allclose(rho, 1.0, atol=1e-10)


def test_cca_of_orthogonal_subspaces_gives_zero_correlations():
    rng = np.random.default_rng(6)
    q, _ = np.linalg.qr(rng.normal(size=(40, 10)))
    _, _, rho = cca(q[:, :5], q[:, 5:], dim=3)
    assert np.allclose(rho, 0.0, atol=1e-10)


def test_cca_recovers_a_shared_latent_variable():
    """The point of PCCA: a latent shared by both views appears in both projections."""
    rng = np.random.default_rng(7)
    n = 80
    latent = rng.normal(size=n)
    x_view = np.column_stack([latent + 0.05 * rng.normal(size=n), rng.normal(size=n)])
    y_view = np.column_stack([latent + 0.05 * rng.normal(size=n), rng.normal(size=n)])
    u_x, _ = pca(x_view, dim=2)
    u_y, _ = pca(y_view, dim=2)
    l_x, l_y, rho = cca(u_x, u_y, dim=1)
    assert rho[0] > 0.95
    assert abs(np.corrcoef(l_x[:, 0], latent)[0, 1]) > 0.9


def test_cca_correlations_are_descending_and_bounded():
    rng = np.random.default_rng(8)
    u_x, _ = np.linalg.qr(rng.normal(size=(50, 6)))
    u_y, _ = np.linalg.qr(rng.normal(size=(50, 4)))
    _, _, rho = cca(u_x, u_y, dim=4)
    assert np.all(rho <= 1.0 + 1e-12) and np.all(rho >= -1e-12)
    assert np.all(np.diff(rho) <= 1e-12)


def test_loo_projection_shape_and_finiteness():
    rng = np.random.default_rng(9)
    n = 25
    latent = rng.normal(size=n)
    X = np.outer(latent, rng.normal(size=40)) + 0.1 * rng.normal(size=(n, 40))
    Y = np.outer(latent, rng.normal(size=60)) + 0.1 * rng.normal(size=(n, 60))
    proj = loo_projection(X, Y @ Y.T, dim_x=5, dim_y=3, dim_z=2)
    assert proj.shape == (n, 2)
    assert np.all(np.isfinite(proj))
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_core.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pcca.core'`

- [ ] **Step 3: Implement core.py**

`playground/pcca/core.py`:

```python
"""PCA and CCA primitives, ported from pachterlab/PCCA `pcca/util.py`.

PCCA is PCA on each view followed by CCA between the two sets of principal
components. The paper uses 30 expression PCs, 5 genotype PCs, and 2 canonical
components.
"""

from __future__ import annotations

import numpy as np


def pca(arr: np.ndarray, dim: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Left singular vectors of a samples x features matrix.

    Computed via the sample-by-sample Gram matrix, which is cheap when features
    vastly outnumber samples (14k genes or 6M SNPs against a few hundred people).
    Equivalent to PCA on a centered and scaled matrix, which is what the callers
    supply.
    """
    if dim is None:
        dim = arr.shape[0]
    return pca_from_gram(arr @ arr.T, dim)


def pca_from_gram(gram: np.ndarray, dim: int) -> tuple[np.ndarray, np.ndarray]:
    """PCA directly from a precomputed Gram matrix.

    Leave-one-out PCA is then just the same call on a submatrix, which turns the
    cross-validation from hours of 6M-SNP decompositions into seconds.
    """
    u, eigenvalues, _ = np.linalg.svd(gram)
    return u[:, :dim], np.sqrt(np.maximum(eigenvalues[:dim], 0.0))


def orthonormalize(mat: np.ndarray) -> np.ndarray:
    """Return a matrix with the same column space but orthonormal columns.

    Required for anything read from plink2 `.eigenvec`: those columns are
    orthogonal but scaled by the singular values. Feeding them to `cca` unfixed
    silently yields wrong canonical correlations rather than an error.
    """
    q, _ = np.linalg.qr(mat)
    return q


def cca(u_x: np.ndarray, u_y: np.ndarray, dim: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Canonical correlation analysis between two orthonormal bases.

    With orthonormal inputs, the canonical correlations are exactly the singular
    values of u_x.T @ u_y.

    Returns (l_x, l_y, rho): the two projections and the canonical correlations.
    """
    for name, mat in (("u_x", u_x), ("u_y", u_y)):
        identity = np.eye(mat.shape[1])
        if not np.allclose(mat.T @ mat, identity, atol=1e-6):
            raise ValueError(
                f"{name} does not have orthonormal columns; pass it through "
                "orthonormalize() first (plink2 .eigenvec columns are not unit-norm)."
            )
    u_c, rho, v_c = np.linalg.svd(u_x.T @ u_y, full_matrices=False)
    l_x = u_x @ u_c[:, :dim]
    l_y = u_y @ v_c.T[:, :dim]
    return l_x, l_y, rho[:dim]


def loo_projection(
    X: np.ndarray, gram_y: np.ndarray, dim_x: int, dim_y: int, dim_z: int
) -> np.ndarray:
    """Leave-one-out PCCA coordinates for every sample.

    For each individual: drop them, recompute both views' principal components on
    the remainder, re-run CCA, and project the held-out individual onto the
    resulting expression-side canonical basis. If individuals land near their own
    population despite never contributing to the fit, the structure is real rather
    than an artefact of overfitting.

    Args:
        X: samples x genes expression residuals.
        gram_y: samples x samples genotype Gram matrix.
        dim_x, dim_y, dim_z: expression PCs, genotype PCs, canonical components.

    Returns:
        (n, dim_z) array of held-out projections, row order matching X.
    """
    n = X.shape[0]
    gram_x = X @ X.T
    out = np.empty((n, dim_z))
    for i in range(n):
        idx = np.delete(np.arange(n), i)
        u_x, _ = pca_from_gram(gram_x[np.ix_(idx, idx)], dim_x)
        u_y, _ = pca_from_gram(gram_y[np.ix_(idx, idx)], dim_y)
        l_x, _, _ = cca(u_x, u_y, dim_z)

        # Map the held-out sample into the training expression space, then onto
        # the canonical basis: its gene vector is scored against the training
        # samples' contributions to each canonical direction.
        weights = X[idx] @ X[i]
        denominator = np.sum(weights)
        out[i] = (weights @ l_x) / (denominator if denominator != 0 else 1.0)
    return out
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_core.py -v`
Expected: PASS, 10 passed.

---

### Task 5: Genotype pipeline

Produces the 5 genotype PCs that anchor the CCA, and the dosage matrices the eQTL scan needs.

**Disk budget:** ~3.5 GB downloaded across all 22 chromosomes (chr1 is 272 MB, chr22 is 52 MB), each deleted after filtering. Peak transient usage ~1.5 GB; the retained filtered callset is ~500 MB.

**Files:**
- Create: `playground/pcca/genotype.py`, `tests/test_genotype.py`
- Modify: `config.yaml` (record the plink2 path)

**Interfaces:**
- Consumes: `pcca.metadata.replication_cohort`
- Produces:
  - `fetch_phase3_urls() -> dict[str, str]`
  - `write_keep_file(cohort: pd.DataFrame, path: Path) -> Path`
  - `filter_chromosome(chrom: int, keep_file: Path, work_dir: Path, out_dir: Path, urls: dict) -> Path` (returns the output prefix)
  - `merge_and_pca(out_dir: Path, prefix_glob: str = "geuv_chr*", n_pcs: int = 5) -> Path` (returns the `.eigenvec` path)
  - `load_eigenvec(path: Path, cohort: pd.DataFrame) -> np.ndarray` of shape (n, 5), rows ordered to match `cohort`
  - `export_dosages(prefix: Path, out_prefix: Path) -> Path`
  - `load_dosages(raw_path: Path) -> (variant_ids: np.ndarray, dosages: np.ndarray)` — a **2-tuple**; positions are derived by the caller from the variant IDs

- [ ] **Step 1: Install plink2**

```bash
cd /private/tmp/claude-502/-Users-nk2-Desktop-biobank-tbooo/d633d09a-778d-47a8-8769-75cb77f0c79a/scratchpad
curl -LO "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_mac_arm64_20260504.zip"
unzip -o plink2_mac_arm64_20260504.zip
mkdir -p ~/bin && mv plink2 ~/bin/plink2 && chmod +x ~/bin/plink2
~/bin/plink2 --version
```

Expected: a version string like `PLINK v2.0.0-a.7`. If macOS Gatekeeper blocks it, run `xattr -d com.apple.quarantine ~/bin/plink2` and retry.

Then set the path in `config.yaml` under `tools:` — change the line `plink2: plink2` to:

```yaml
  plink2: ~/bin/plink2
```

- [ ] **Step 2: Write the failing tests**

`tests/test_genotype.py`:

```python
import numpy as np
import pandas as pd
import pytest

from pcca.genotype import fetch_phase3_urls, load_dosages, load_eigenvec, write_keep_file


def test_write_keep_file_uses_double_id_columns(tmp_path):
    """plink2 --keep with a .psam using IID-only still needs two columns here."""
    cohort = pd.DataFrame({"sample_id": ["HG00096", "NA12878"]})
    out = tmp_path / "keep.txt"
    write_keep_file(cohort, out)
    lines = out.read_text().strip().split("\n")
    assert lines == ["HG00096\tHG00096", "NA12878\tNA12878"]


def test_load_eigenvec_reorders_rows_to_match_cohort(tmp_path):
    """plink2 emits samples in .psam order; the CCA needs expression order."""
    path = tmp_path / "geno.eigenvec"
    path.write_text(
        "#FID\tIID\tPC1\tPC2\n"
        "NA12878\tNA12878\t0.5\t0.6\n"
        "HG00096\tHG00096\t0.1\t0.2\n"
    )
    cohort = pd.DataFrame({"sample_id": ["HG00096", "NA12878"]})
    pcs = load_eigenvec(path, cohort)
    assert pcs.shape == (2, 2)
    assert np.allclose(pcs[0], [0.1, 0.2])
    assert np.allclose(pcs[1], [0.5, 0.6])


def test_load_eigenvec_rejects_missing_samples(tmp_path):
    path = tmp_path / "geno.eigenvec"
    path.write_text("#FID\tIID\tPC1\nHG00096\tHG00096\t0.1\n")
    cohort = pd.DataFrame({"sample_id": ["HG00096", "MISSING"]})
    with pytest.raises(KeyError, match="MISSING"):
        load_eigenvec(path, cohort)


def test_load_dosages_parses_raw_and_mean_imputes(tmp_path):
    """plink2 --export A writes 6 leading columns then one column per variant."""
    path = tmp_path / "chr1.raw"
    path.write_text(
        "FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE\t1:100:A:G_G\t1:200:C:T_T\n"
        "HG00096\tHG00096\t0\t0\t0\t-9\t0\t1\n"
        "HG00097\tHG00097\t0\t0\t0\t-9\t2\tNA\n"
    )
    ids, dosages = load_dosages(path)
    assert list(ids) == ["1:100:A:G", "1:200:C:T"]
    assert dosages.shape == (2, 2)
    assert np.allclose(dosages[:, 0], [0.0, 2.0])
    # The NA is replaced by the column mean, which is 1.0 from the single value.
    assert np.allclose(dosages[:, 1], [1.0, 1.0])


@pytest.mark.network
def test_fetch_phase3_urls_finds_every_autosome():
    urls = fetch_phase3_urls()
    for chrom in range(1, 23):
        assert f"chr{chrom}_phase3.pgen.zst" in urls
        assert f"chr{chrom}_phase3.pvar.zst" in urls
    assert "phase3_orig.psam" in urls
    assert urls["chr1_phase3.pgen.zst"].startswith("https://")
```

Register the `network` marker by creating `pytest.ini` at the repo root:

```ini
[pytest]
markers =
    network: test requires internet access
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_genotype.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pcca.genotype'`

- [ ] **Step 4: Implement genotype.py**

`playground/pcca/genotype.py`:

```python
"""plink2 driver for the genotype side of PCCA.

Everything works one chromosome at a time: the machine has ~22 GiB free on a
98%-full volume, and a whole-cohort decompressed callset does not fit.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from pcca import paths

PLINK2 = str(Path("~/bin/plink2").expanduser())
AUTOSOMES = list(range(1, 23))


def fetch_phase3_urls() -> dict[str, str]:
    """Scrape the plink2 resources page for build-37 Phase 3 download links.

    The files live on Dropbox behind opaque hashes that change between releases,
    so hardcoding them would rot silently.
    """
    with urllib.request.urlopen(paths.PLINK2_RESOURCES_URL, timeout=120) as response:
        html = response.read().decode("utf-8", errors="replace")
    urls = {}
    for match in re.finditer(r'https://www\.dropbox\.com/s/[^"\'>]+', html):
        url = match.group(0)
        filename = url.split("/")[-1].split("?")[0]
        urls.setdefault(filename, url)
    return urls


def write_keep_file(cohort: pd.DataFrame, path: Path) -> Path:
    """Write a plink2 --keep file (FID IID) for the cohort's 1KGP sample IDs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [f"{s}\t{s}" for s in cohort["sample_id"]]
    path.write_text("\n".join(lines) + "\n")
    return path


def _run(args: list[str]) -> None:
    subprocess.run([str(a) for a in args], check=True)


def download_phase3_chromosome(chrom: int, urls: dict[str, str], out_dir: Path) -> Path:
    """Download and decompress one chromosome's pgen/pvar, plus the psam."""
    out_dir.mkdir(parents=True, exist_ok=True)
    psam = out_dir / "phase3_orig.psam"
    if not psam.exists():
        urllib.request.urlretrieve(urls["phase3_orig.psam"], psam)

    prefix = out_dir / f"chr{chrom}_phase3"
    for suffix in ("pgen", "pvar"):
        zst = out_dir / f"chr{chrom}_phase3.{suffix}.zst"
        plain = out_dir / f"chr{chrom}_phase3.{suffix}"
        if plain.exists():
            continue
        if not zst.exists():
            urllib.request.urlretrieve(urls[f"chr{chrom}_phase3.{suffix}.zst"], zst)
        _run([PLINK2, "--zst-decompress", zst, plain])
        zst.unlink()
    shutil.copy(psam, prefix.with_suffix(".psam"))
    return prefix


def filter_chromosome(chrom: int, keep_file: Path, work_dir: Path, out_dir: Path,
                      urls: dict[str, str]) -> Path:
    """Subset one chromosome to the cohort and to common biallelic autosomal SNPs.

    Deletes the decompressed source immediately afterwards so the next chromosome
    has room. No LD pruning: the paper does not prune, and pruning would change
    the PC structure being replicated.
    """
    out_prefix = out_dir / f"geuv_chr{chrom}"
    if out_prefix.with_suffix(".pgen").exists():
        return out_prefix

    src = download_phase3_chromosome(chrom, urls, work_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    _run([
        PLINK2, "--pfile", src,
        "--keep", keep_file,
        "--snps-only", "--max-alleles", "2", "--min-alleles", "2",
        "--maf", "0.05", "--autosome",
        "--set-all-var-ids", "@:#:$r:$a",
        "--new-id-max-allele-len", "100", "truncate",
        "--make-pgen", "--out", out_prefix,
    ])
    for suffix in (".pgen", ".pvar", ".psam"):
        src.with_suffix(suffix).unlink(missing_ok=True)
    return out_prefix


def merge_and_pca(out_dir: Path, prefix_glob: str = "geuv_chr*", n_pcs: int = 5) -> Path:
    """Merge per-chromosome sets and compute genotype PCs."""
    prefixes = sorted(
        out_dir.glob(f"{prefix_glob}.pgen"),
        key=lambda p: int(re.search(r"chr(\d+)", p.name).group(1)),
    )
    merge_list = out_dir / "merge_list.txt"
    merge_list.write_text("\n".join(str(p.with_suffix("")) for p in prefixes) + "\n")

    merged = out_dir / "geuv_merged"
    # The trailing `pfile` tells plink2 the list holds pgen/pvar/psam prefixes.
    _run([PLINK2, "--pmerge-list", merge_list, "pfile", "--make-pgen", "--out", merged])
    _run([PLINK2, "--pfile", merged, "--pca", str(n_pcs), "--out", merged])
    return merged.with_suffix(".eigenvec")


def load_eigenvec(path: Path, cohort: pd.DataFrame) -> np.ndarray:
    """Read plink2 PCs, reordered to match the cohort's row order.

    plink2 emits samples in .psam order, which is not the expression matrix order.
    Misaligned rows would produce a plausible-looking but meaningless CCA.
    """
    df = pd.read_csv(path, sep="\t")
    # plink2 writes a leading '#' on the header and may emit either "#IID" or
    # "#FID<TAB>IID" depending on the input .psam, so locate IID by name.
    df.columns = [c.lstrip("#") for c in df.columns]
    id_col = "IID" if "IID" in df.columns else df.columns[0]
    df[id_col] = df[id_col].astype(str)
    indexed = df.set_index(id_col)
    missing = [s for s in cohort["sample_id"] if s not in indexed.index]
    if missing:
        raise KeyError(f"{len(missing)} cohort samples absent from {path.name}: {missing[:5]}")
    pc_cols = [c for c in indexed.columns if c.upper().startswith("PC")]
    return indexed.loc[cohort["sample_id"], pc_cols].to_numpy(dtype=float)


def export_dosages(prefix: Path, out_prefix: Path) -> Path:
    """Export additive dosages for one chromosome as a plink2 .raw file."""
    _run([PLINK2, "--pfile", prefix, "--export", "A", "--out", out_prefix])
    return out_prefix.with_suffix(".raw")


def load_dosages(raw_path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Parse a plink2 .raw file into (variant_ids, samples x variants dosages).

    plink2 appends the counted allele to each column name (`1:100:A:G_G`); it is
    stripped so IDs join back to the .pvar. Missing calls are mean-imputed, which
    is standard for eQTL scans and keeps the regression well-defined.
    """
    df = pd.read_csv(raw_path, sep="\t")
    leading = ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]
    dosage_cols = [c for c in df.columns if c not in leading]
    ids = np.array([c.rsplit("_", 1)[0] for c in dosage_cols])
    mat = df[dosage_cols].to_numpy(dtype=np.float32)
    col_means = np.nanmean(mat, axis=0)
    inds = np.where(np.isnan(mat))
    mat[inds] = np.take(col_means, inds[1])
    return ids, mat
```

- [ ] **Step 5: Run the unit tests**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_genotype.py -v -m "not network"`
Expected: PASS, 4 passed.

Then the network test: `/opt/anaconda3/bin/python -m pytest tests/test_genotype.py -v -m network`
Expected: PASS, 1 passed.

- [ ] **Step 6: Build the genotype callset**

```python
# Run as a script or in the notebook.
from pathlib import Path
from pcca import paths
from pcca.genotype import fetch_phase3_urls, filter_chromosome, merge_and_pca, write_keep_file
from pcca.metadata import load_sdrf, replication_cohort

cohort = replication_cohort(load_sdrf(paths.SDRF))
geno_dir = paths.GENOTYPE
work = geno_dir / "work"
keep = write_keep_file(cohort, geno_dir / "cohort_keep.txt")

urls = fetch_phase3_urls()
for chrom in range(1, 23):
    print(f"chr{chrom}…", flush=True)
    filter_chromosome(chrom, keep, work, geno_dir / "phase3", urls)

eigenvec = merge_and_pca(geno_dir / "phase3")
print(eigenvec)
```

Monitor free space between chromosomes with `df -h .`; it should stay roughly flat because each decompressed source is deleted after filtering.

- [ ] **Step 7: Sanity-check the genotype PCs**

```python
import numpy as np
from pcca.genotype import load_eigenvec

pcs = load_eigenvec(eigenvec, cohort)
assert pcs.shape == (343, 5)
for pc in range(2):
    yri = pcs[(cohort["pop"] == "YRI").to_numpy(), pc]
    eur = pcs[(cohort["pop"] != "YRI").to_numpy(), pc]
    print(f"PC{pc+1}: YRI mean {yri.mean():+.4f}, EUR mean {eur.mean():+.4f}")
```

Expected: one of the first two PCs separates YRI from the European populations by a wide margin — the classic African/European axis. If neither does, the sample ordering or the filtering is wrong; stop and fix before proceeding.

---

### Task 6: Figure 1 and Figure 2 — the core replication

This is the task that either reproduces the paper or tells us the preprocessing is wrong.

**Files:**
- Create: `playground/pcca/figures.py`, `tests/test_figures.py`
- Create: `playground/geuvadis_pcca.ipynb`

**Interfaces:**
- Consumes: everything from Tasks 1–5
- Produces: `scatter_by_population(coords, cohort, title, xlabel, ylabel, path) -> Path`. `POP_COLORS: dict[str, str]`.

- [ ] **Step 1: Write the failing test**

`tests/test_figures.py`:

```python
import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")

from pcca.figures import POP_COLORS, scatter_by_population


@pytest.fixture
def cohort():
    return pd.DataFrame({
        "sample_id": [f"S{i}" for i in range(8)],
        "pop": ["GBR", "FIN", "TSI", "YRI"] * 2,
    })


def test_every_geuvadis_population_has_a_colour():
    assert set(POP_COLORS) >= {"GBR", "FIN", "TSI", "YRI", "CEU"}


def test_scatter_writes_a_png(tmp_path, cohort):
    coords = np.arange(16, dtype=float).reshape(8, 2)
    out = scatter_by_population(coords, cohort, "t", "x", "y", tmp_path / "f.png")
    assert out.exists() and out.stat().st_size > 1000


def test_scatter_rejects_mismatched_lengths(tmp_path, cohort):
    coords = np.zeros((3, 2))
    with pytest.raises(ValueError, match="rows"):
        scatter_by_population(coords, cohort, "t", "x", "y", tmp_path / "f.png")
```

- [ ] **Step 2: Run to verify it fails**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_figures.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pcca.figures'`

- [ ] **Step 3: Implement the scatter helper**

`playground/pcca/figures.py`:

```python
"""Figures for the PCCA replication."""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

POP_COLORS = {
    "GBR": "#4C72B0",
    "FIN": "#DD8452",
    "TSI": "#55A868",
    "YRI": "#C44E52",
    "CEU": "#8172B3",
}
POP_ORDER = ["GBR", "FIN", "TSI", "YRI", "CEU"]


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
            coords[mask, 0], coords[mask, 1],
            c=POP_COLORS[pop], s=24, alpha=0.8, linewidths=0, label=pop,
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
```

- [ ] **Step 4: Run to verify it passes**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_figures.py -v`
Expected: PASS, 3 passed.

- [ ] **Step 5: Create the notebook and produce Figure 1**

Create `playground/geuvadis_pcca.ipynb` with a first cell that sets up imports and loads everything:

```python
import sys
from pathlib import Path

# The notebook lives in playground/, so its own directory is the import root.
sys.path.insert(0, str(Path.cwd()))

import numpy as np
from pcca import paths
from pcca.core import cca, loo_projection, orthonormalize, pca
from pcca.expression import build_expression_matrix, load_gd660
from pcca.figures import scatter_by_population
from pcca.genotype import load_eigenvec
from pcca.metadata import load_sdrf, replication_cohort

cohort = replication_cohort(load_sdrf(paths.SDRF))
ann, expr = load_gd660(paths.GD660_GENE_RPKM)
ann_f, X = build_expression_matrix(ann, expr, cohort)
geno_pcs = load_eigenvec(paths.GENOTYPE / "phase3" / "geuv_merged.eigenvec", cohort)
print(f"{X.shape[0]} individuals, {X.shape[1]} genes, {geno_pcs.shape[1]} genotype PCs")
```

Then Figure 1 — the paper's negative result:

```python
u_x, _ = pca(X, dim=30)
scatter_by_population(
    u_x[:, :2], cohort,
    "Figure 1A — PCA of expression fails to cluster by population",
    "Expression PC1", "Expression PC2",
    paths.FIGURES / "pcca_fig1a_expression_pca.png",
)
scatter_by_population(
    geno_pcs[:, :2], cohort,
    "Figure 1B — PCA of genotype clusters by population",
    "Genotype PC1", "Genotype PC2",
    paths.FIGURES / "pcca_fig1b_genotype_pca.png",
)
```

Expected: 1A shows populations intermixed; 1B shows YRI cleanly separated from the Europeans. If 1A already separates populations, the confounder regression or the input file is wrong.

- [ ] **Step 6: Run PCCA and check the canonical correlations**

```python
u_y = orthonormalize(geno_pcs)
l_x, l_y, rho = cca(u_x, u_y, dim=2)

print(f"canonical correlations: {rho[0]:.3f}, {rho[1]:.3f}")
print("paper reports:          0.963, 0.766")

scatter_by_population(
    l_x, cohort,
    f"Figure 2A — PCCA projection of expression (rho = {rho[0]:.3f}, {rho[1]:.3f})",
    "Coordinate 1", "Coordinate 2",
    paths.FIGURES / "pcca_fig2a_projection.png",
)
```

**This is the replication gate.** Expected: `rho[0]` near 0.963 and `rho[1]` near 0.766, and the scatter separating YRI from the Europeans with FIN distinguishable.

If the correlations are materially low, the mitigations, in order:

1. **Swap the quantification.** Use `GD660.TrQuantRPKM.txt.gz` summed to genes instead of `GD660.GeneQuantRPKM.txt.gz`. This reproduces the paper's transcript→gene step more closely and is the only remaining lever on the expression side.
2. **Check the confounder regression.** Confirm all 7 labs are present in the design and that the residuals are orthogonal to it (`test_build_expression_matrix_replication_cohort` already asserts this).
3. **Vary the expression PC count.** The paper chose 30 by the elbow method on its own quantification; a different quantifier can put the elbow elsewhere. Sweep 20–40 and report the curve rather than silently tuning to the target.

Note that the Phase 1 genotype route is **not** available as a mitigation — it was ruled out on download size. If the correlations remain low after the above, report the observed values with the genotype release named as the untested remaining explanation, rather than presenting a miss as a match.

- [ ] **Step 7: Leave-one-out cross-validation**

```python
# Genotype Gram matrix from the orthonormal PCs: the CCA only ever sees this
# subspace, so it is the right object to subset for holdouts.
gram_y = geno_pcs @ geno_pcs.T
loo = loo_projection(X, gram_y, dim_x=30, dim_y=5, dim_z=2)

scatter_by_population(
    loo, cohort,
    "Figure 2B — leave-one-out projection",
    "Coordinate 1", "Coordinate 2",
    paths.FIGURES / "pcca_fig2b_loo.png",
)
```

Expected: held-out individuals land near their own population's cluster. Runs in seconds thanks to the Gram trick; if it takes minutes, `loo_projection` is recomputing full PCAs and the implementation is wrong.

- [ ] **Step 8: Produce the production-config embedding for the fusion work**

The replication cohort deliberately drops CEU and the 41 individuals without Phase 1 genotypes. The follow-on fusion with Antoni's SNP clustering needs all 462, so run the same code path over the production cohort and persist the coordinates.

This requires a second genotype callset keyed to all 462 individuals:

```python
from pcca.metadata import production_cohort

prod = production_cohort(load_sdrf(paths.SDRF))
prod_dir = paths.GENOTYPE / "phase3_all462"
keep_all = write_keep_file(prod, paths.GENOTYPE / "cohort_keep_all462.txt")
for chrom in range(1, 23):
    filter_chromosome(chrom, keep_all, paths.GENOTYPE / "work", prod_dir, urls)
eigenvec_all = merge_and_pca(prod_dir)
```

Then the same PCCA, and write the result out:

```python
import pandas as pd

ann_p, X_p = build_expression_matrix(ann, expr, prod)
geno_pcs_p = load_eigenvec(eigenvec_all, prod)
u_x_p, _ = pca(X_p, dim=30)
l_x_p, l_y_p, rho_p = cca(u_x_p, orthonormalize(geno_pcs_p), dim=2)
print(f"production canonical correlations: {rho_p[0]:.3f}, {rho_p[1]:.3f}")

scatter_by_population(
    l_x_p, prod,
    f"PCCA projection, all 462 individuals (rho = {rho_p[0]:.3f}, {rho_p[1]:.3f})",
    "Coordinate 1", "Coordinate 2",
    paths.FIGURES / "pcca_production_462.png",
)

out = pd.DataFrame({
    "sample_id": prod["sample_id"],
    "pop": prod["pop"],
    "pcca_c1": l_x_p[:, 0],
    "pcca_c2": l_x_p[:, 1],
})
for i in range(5):
    out[f"geno_pc{i + 1}"] = geno_pcs_p[:, i]
out.to_csv(paths.METADATA / "geuvadis_pcca_coords.tsv", sep="\t", index=False)
print(f"wrote {paths.METADATA / 'geuvadis_pcca_coords.tsv'}")
```

Expected: 462 rows with both transcriptional coordinates and genotype PCs — the transcriptional vector the fusion work consumes. Including CEU may lower the canonical correlations slightly, since CEU is the population the paper excluded for biased expression; note the value rather than treating a drop as a failure.

---

### Task 7: cis-eQTL scan

**Files:**
- Create: `playground/pcca/eqtl.py`, `tests/test_eqtl.py`

**Interfaces:**
- Consumes: `pcca.genotype.load_dosages`, `pcca.expression`
- Produces: `residualize(mat, C) -> np.ndarray`. `association_pvalues(y, G, dof) -> (r, p)`. `cis_scan(ann, X, chrom_dosages, covariates, window=1_000_000) -> pd.DataFrame` with columns `gene, best_p, n_snps`.

- [ ] **Step 1: Write the failing tests**

`tests/test_eqtl.py`:

```python
import numpy as np
import pandas as pd
import pytest
from scipy import stats

from pcca.eqtl import association_pvalues, cis_scan, residualize


def test_association_pvalues_match_scipy_ols_without_covariates():
    """Ground truth: with only an intercept, this reduces to Pearson correlation."""
    rng = np.random.default_rng(0)
    n = 200
    G = rng.integers(0, 3, size=(n, 5)).astype(float)
    y = G[:, 0] * 0.4 + rng.normal(size=n)

    y_c = y - y.mean()
    G_c = G - G.mean(axis=0)
    _, p = association_pvalues(y_c, G_c, dof=n - 2)

    for j in range(5):
        expected = stats.pearsonr(y, G[:, j])[1]
        assert p[j] == pytest.approx(expected, rel=1e-6)


def test_association_pvalues_detect_a_planted_effect():
    rng = np.random.default_rng(1)
    n = 300
    G = rng.integers(0, 3, size=(n, 10)).astype(float)
    y = 1.5 * G[:, 3] + rng.normal(size=n)
    y_c = y - y.mean()
    G_c = G - G.mean(axis=0)
    _, p = association_pvalues(y_c, G_c, dof=n - 2)
    assert np.argmin(p) == 3
    assert p[3] < 1e-20


def test_residualize_matches_regress_out_semantics():
    rng = np.random.default_rng(2)
    n = 50
    C = np.column_stack([np.ones(n), rng.normal(size=n)])
    M = rng.normal(size=(n, 8))
    resid = residualize(M, C)
    assert np.allclose(C.T @ resid, 0.0, atol=1e-9)


def test_cis_scan_only_tests_snps_inside_the_window():
    """A SNP 5 Mb away must not be tested even if it is perfectly correlated."""
    rng = np.random.default_rng(3)
    n = 100
    ann = pd.DataFrame({"TargetID": ["G1"], "Chr": ["1"], "Coord": [5_000_000]})
    driver = rng.normal(size=n)
    X = (driver + 0.01 * rng.normal(size=n)).reshape(n, 1)
    # Column 0 is in-window; column 1 is far away but identical to the driver.
    G = np.column_stack([rng.normal(size=n), driver])
    positions = np.array([5_100_000, 20_000_000])
    ids = np.array(["1:5100000:A:G", "1:20000000:A:G"])

    C = np.ones((n, 1))
    result = cis_scan(ann, X, {"1": (ids, positions, G)}, C, window=1_000_000)
    assert len(result) == 1
    assert result.loc[0, "n_snps"] == 1


def test_cis_scan_reports_the_best_p_per_gene():
    rng = np.random.default_rng(4)
    n = 200
    ann = pd.DataFrame({"TargetID": ["G1"], "Chr": ["1"], "Coord": [1_000_000]})
    G = rng.integers(0, 3, size=(n, 4)).astype(float)
    X = (2.0 * G[:, 2] + rng.normal(size=n)).reshape(n, 1)
    ids = np.array([f"1:{1_000_000 + i}:A:G" for i in range(4)])
    positions = np.array([1_000_000 + i for i in range(4)])

    C = np.ones((n, 1))
    result = cis_scan(ann, X, {"1": (ids, positions, G)}, C, window=1_000_000)
    assert result.loc[0, "n_snps"] == 4
    assert result.loc[0, "best_p"] < 1e-10
```

- [ ] **Step 2: Run to verify they fail**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_eqtl.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pcca.eqtl'`

- [ ] **Step 3: Implement eqtl.py**

`playground/pcca/eqtl.py`:

```python
"""Vectorized cis-eQTL scan.

The comparison the paper draws is between two covariate strategies, so the scan
runs twice over identical genotypes and phenotypes with only the covariates
swapped. A naive per-(gene, SNP) loop is far too slow; instead the genotypes are
residualized against the covariates once per chromosome, after which each test is
a dot product.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

CIS_WINDOW = 1_000_000


def residualize(mat: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Residuals of `mat` (samples x k) after regressing on covariates C."""
    beta, *_ = np.linalg.lstsq(C, mat, rcond=None)
    return mat - C @ beta


def association_pvalues(
    y: np.ndarray, G: np.ndarray, dof: int
) -> tuple[np.ndarray, np.ndarray]:
    """Two-sided p-values for the correlation of y with every column of G.

    Both inputs must already be covariate-residualized (and therefore centered).
    For simple linear regression the t-statistic depends only on the correlation,
    so this is exact rather than an approximation.
    """
    y_norm = np.linalg.norm(y)
    g_norms = np.linalg.norm(G, axis=0)
    g_norms = np.where(g_norms == 0, np.inf, g_norms)
    r = (G.T @ y) / (g_norms * y_norm)
    r = np.clip(r, -1 + 1e-12, 1 - 1e-12)
    t = r * np.sqrt(dof / (1 - r**2))
    p = 2 * stats.t.sf(np.abs(t), dof)
    return r, p


def cis_scan(
    ann: pd.DataFrame,
    X: np.ndarray,
    chrom_dosages: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]],
    covariates: np.ndarray,
    window: int = CIS_WINDOW,
) -> pd.DataFrame:
    """Best cis association per gene.

    Args:
        ann: gene annotation with TargetID, Chr, Coord (GRCh37).
        X: samples x genes expression, columns aligned to `ann` rows.
        chrom_dosages: chromosome -> (variant_ids, positions, samples x variants).
        covariates: samples x k design matrix, must include an intercept.
        window: cis distance either side of the gene's coordinate.

    Returns:
        DataFrame with gene, best_p, n_snps — one row per gene that had at least
        one cis SNP.
    """
    n_samples = X.shape[0]
    # Rank, not column count: the PCCA covariate set mixes canonical components
    # with expression PCs drawn from the same 30-PC subspace, so the design can be
    # rank-deficient. Using the column count there would overstate the degrees of
    # freedom and quietly inflate every p-value in the comparison.
    dof = n_samples - int(np.linalg.matrix_rank(covariates)) - 1
    X_resid = residualize(X, covariates)

    records = []
    for chrom, (ids, positions, G) in chrom_dosages.items():
        gene_rows = np.where(ann["Chr"].to_numpy() == chrom)[0]
        if len(gene_rows) == 0:
            continue
        G_resid = residualize(G.astype(float), covariates)

        for row in gene_rows:
            coord = int(ann["Coord"].iloc[row])
            in_window = np.abs(positions - coord) <= window
            n_snps = int(in_window.sum())
            if n_snps == 0:
                continue
            _, p = association_pvalues(X_resid[:, row], G_resid[:, in_window], dof)
            records.append(
                {
                    "gene": ann["TargetID"].iloc[row],
                    "chrom": chrom,
                    "best_p": float(np.min(p)),
                    "n_snps": n_snps,
                }
            )
    return pd.DataFrame.from_records(records).reset_index(drop=True)


def count_egenes(scan: pd.DataFrame, alpha: float) -> int:
    """Genes with at least one cis SNP below alpha."""
    return int((scan["best_p"] < alpha).sum())
```

- [ ] **Step 4: Run to verify they pass**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_eqtl.py -v`
Expected: PASS, 5 passed.

---

### Task 8: Panels B, C and D

**Files:**
- Modify: `playground/pcca/figures.py`, `tests/test_figures.py`
- Modify: `playground/geuvadis_pcca.ipynb`

**Interfaces:**
- Produces: `qq_plot(pvalue_sets, path)`, `egene_curve(scans, alphas, path)`, `egene_venn(sets_by_method, alpha, path)`

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_figures.py`:

```python
from pcca.figures import egene_curve, egene_venn, qq_plot


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
```

- [ ] **Step 2: Run to verify they fail**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_figures.py -v -k "qq or curve or venn"`
Expected: FAIL — `ImportError: cannot import name 'qq_plot'`

- [ ] **Step 3: Implement the three panels**

Append to `playground/pcca/figures.py`:

```python
METHOD_COLORS = {"pca": "#4C72B0", "pcca": "#DD8452"}


def qq_plot(pvalue_sets: dict[str, np.ndarray], path: Path) -> Path:
    """Observed vs expected -log10(p) under the uniform null.

    The paper's claim is that the PCCA covariates reduce inflation at the high
    end, so the interesting region is the top right against the diagonal.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5, 4.5))
    limit = 0.0
    for name, pvals in pvalue_sets.items():
        pvals = np.asarray(pvals)
        pvals = pvals[np.isfinite(pvals) & (pvals > 0)]
        observed = -np.log10(np.sort(pvals))[::-1]
        expected = -np.log10(
            (np.arange(len(observed), 0, -1) - 0.5) / len(observed)
        )
        ax.scatter(expected, observed, s=6, alpha=0.6, linewidths=0,
                   label=name, c=METHOD_COLORS.get(name))
        limit = max(limit, expected.max(), observed.max())
    ax.plot([0, limit], [0, limit], color="#666666", linewidth=1, linestyle="--")
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
        ax.plot(alphas, counts, marker="o", markersize=4, linewidth=1.5,
                label=name, color=METHOD_COLORS.get(name))
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
```

- [ ] **Step 4: Run to verify they pass**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_figures.py -v`
Expected: PASS, 6 passed.

- [ ] **Step 5: Run both eQTL scans in the notebook**

```python
from pcca.eqtl import cis_scan, count_egenes
from pcca.genotype import export_dosages, load_dosages
from pcca.figures import egene_curve, egene_venn, qq_plot

# PEER substitute: the top 10 expression PCs (see the spec's Deviations table).
peer_surrogate = pca(X, dim=10)[0]

# The paper uses 5 canonical components for the eQTL covariates, not the 2 used
# for the figure, and takes them from both projections.
l_x5, l_y5, _ = cca(u_x, u_y, dim=5)

intercept = np.ones((X.shape[0], 1))
cov_pca = np.column_stack([intercept, geno_pcs, peer_surrogate])
cov_pcca = np.column_stack([intercept, l_x5, l_y5, peer_surrogate])
print(f"covariate counts — pca: {cov_pca.shape[1]}, pcca: {cov_pcca.shape[1]}")

# Load dosages one chromosome at a time to stay inside memory.
geno_dir = paths.GENOTYPE / "phase3"
scans = {"pca": [], "pcca": []}
for chrom in range(1, 23):
    raw = export_dosages(geno_dir / f"geuv_chr{chrom}", geno_dir / f"dose_chr{chrom}")
    ids, G = load_dosages(raw)
    positions = np.array([int(v.split(":")[1]) for v in ids])
    payload = {str(chrom): (ids, positions, G)}
    for name, cov in (("pca", cov_pca), ("pcca", cov_pcca)):
        scans[name].append(cis_scan(ann_f, X, payload, cov))
    raw.unlink()
    del G
    print(f"chr{chrom} done", flush=True)

import pandas as pd
scans = {k: pd.concat(v, ignore_index=True) for k, v in scans.items()}
for name, scan in scans.items():
    print(f"{name}: {count_egenes(scan, 1e-6)} eGenes at 1e-6 ({len(scan)} genes tested)")
```

Expected: eGene counts in the low thousands. The paper reports 2,818 (PCA) and 2,732 (PCCA); exact agreement is not expected given the quantification and genotype-release differences, but PCCA should find slightly fewer.

- [ ] **Step 6: Draw panels B, C and D**

```python
qq_plot({k: v["best_p"].to_numpy() for k, v in scans.items()},
        paths.FIGURES / "pcca_panelB_qq.png")

egene_curve(scans, np.logspace(-9, -3, 13), paths.FIGURES / "pcca_panelC_curve.png")

egene_sets = {k: set(v.loc[v["best_p"] < 1e-6, "gene"]) for k, v in scans.items()}
_, counts = egene_venn(egene_sets, 1e-6, paths.FIGURES / "pcca_panelD_venn.png")
print(counts)
```

---

### Task 9: Correct the `tbooo map geuvadis` stage

**Files:**
- Modify: `src/tbooo/pipeline/geuvadis.py`
- Create: `tests/test_pipeline_geuvadis.py`

**Interfaces:**
- Consumes: the GD660 file downloaded in Task 2
- Produces: unchanged outputs — `geuvadis_expression_pcs.tsv`, `geuvadis_pca_variance.tsv`, `geuvadis_pc1..20` columns in `participant.parquet`

- [ ] **Step 1: Write the failing test**

`tests/test_pipeline_geuvadis.py`:

```python
import gzip

import pandas as pd
import pytest

from tbooo.pipeline.geuvadis import _is_sample_id, _load_rpkm


def test_is_sample_id_accepts_1kgp_ids():
    assert _is_sample_id("HG00096")
    assert _is_sample_id("NA12878")


def test_is_sample_id_rejects_annotation_columns():
    assert not _is_sample_id("TargetID")
    assert not _is_sample_id("Chr")
    assert not _is_sample_id("Coord")


def test_load_rpkm_collapses_gd660_run_columns_to_sample_ids(tmp_path):
    """GD660 headers are run IDs like HG00096.1.M_111124_6, not bare sample IDs."""
    path = tmp_path / "gd660.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(
            "TargetID\tGene_Symbol\tChr\tCoord\t"
            "HG00096.1.M_111124_6\tNA12878.2.M_120209_7\n"
        )
        handle.write("ENSG1\tENSG1\t1\t1000\t5.0\t7.0\n")
    frame = _load_rpkm(path)
    assert list(frame.columns) == ["HG00096", "NA12878"]
    assert frame.shape == (1, 2)
```

- [ ] **Step 2: Run to verify it fails**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_pipeline_geuvadis.py -v`
Expected: FAIL on the third test — the current `_load_rpkm` keeps only columns whose names are exactly 7 characters, so every GD660 run column is discarded and it raises `ValueError: No 1KGP-style sample ID columns found`.

- [ ] **Step 3: Update the stage**

In `src/tbooo/pipeline/geuvadis.py`:

Change the filename constant:

```python
_RPKM_FILENAME = "GD660.GeneQuantRPKM.txt.gz"
```

Replace `_load_rpkm` so it maps run columns to sample IDs and averages the technical replicates:

```python
def _load_rpkm(path) -> pd.DataFrame:
    """Return DataFrame with genes as rows, 1KGP sample IDs as columns.

    GD660 column headers are run identifiers (`HG00096.1.M_111124_6`) covering
    462 individuals across 660 runs. Runs are averaged per individual so the
    result is one column per person.
    """
    df = pd.read_csv(path, sep="\t", compression="gzip", index_col=0)

    run_cols = [c for c in df.columns if _is_sample_id(c.split(".")[0])]
    if not run_cols:
        raise ValueError(
            f"No 1KGP-style sample columns found in {path}. "
            "Expected run identifiers like 'HG00096.1.M_111124_6'."
        )
    sample_ids = [c.split(".")[0] for c in run_cols]
    return df[run_cols].groupby(sample_ids, axis=1).mean()
```

Update the module docstring's step 1 to read:

```
    1. Read GD660 raw RPKM matrix (genes × runs), averaging technical replicates
       per individual. Raw rather than PEER-residualized: PEER factors absorb
       population structure, which is exactly the signal downstream work needs.
```

- [ ] **Step 4: Run to verify the tests pass**

Run: `/opt/anaconda3/bin/python -m pytest tests/test_pipeline_geuvadis.py -v`
Expected: PASS, 3 passed.

If pandas 1.5.3 warns that `groupby(axis=1)` is deprecated, switch to `df[run_cols].T.groupby(sample_ids).mean().T`, which is equivalent and forward-compatible.

- [ ] **Step 5: Re-run the stage end to end**

The existing `table_ok` guard skips recomputation when the old outputs are intact, so delete them first:

```bash
cd /Users/nk2/Desktop/biobank/tbooo
rm -f data/metadata/geuvadis_expression_pcs.tsv data/metadata/geuvadis_pca_variance.tsv
/opt/anaconda3/bin/tbooo map geuvadis
head -2 data/metadata/geuvadis_expression_pcs.tsv
```

Expected: 462 samples matched to EIDs, and a variance table whose PC1 share is reported in the log.

- [ ] **Step 6: Regenerate the EDA figures**

```bash
/opt/anaconda3/bin/tbooo run --target eda 2>/dev/null || /opt/anaconda3/bin/tbooo map eda
```

If neither target exists, call `run_eda` directly:

```bash
/opt/anaconda3/bin/python -c "
from tbooo.config import Config
from tbooo.pipeline.eda import run_eda
run_eda(Config.load('config.yaml'))
"
```

Expected: `data/metadata/eda/rna_pca_grid.png` regenerated from raw RPKM. It should still show populations intermixed — that is the paper's Figure 1A result, now honestly derived.

---

### Task 10: Three-paper methodology document

**Files:**
- Create: `docs/4_pcca_methods.md`
- Modify: `README.md` (documentation table)

- [ ] **Step 1: Write the methodology document**

Create `docs/4_pcca_methods.md` covering, for each of the three papers: the question it asks, its inputs, its preprocessing, its algorithm, its headline results, and what it means for this project. Required content:

1. **Lappalainen et al., Nature 2013** — the GEUVADIS resource. 462 LCL RNA-seq samples, 5 populations (CEU 91, FIN 95, GBR 94, TSI 93, YRI 89 by the SDRF), 7 sequencing labs, 660 runs including technical replicates. Quantification files and what the name parts mean: `GD660` = all QC-passed runs, `GD462` = nonredundant individuals, `50FN` = quantile normalized over genes expressed in >50% of samples, `resk10` = residuals after removing 10 PEER factors. State plainly that PEER residualization removes population structure, which is why this project uses GD660.

2. **Brown, Bray & Pachter, PLOS Genetics 2018** — PCCA. 343 individuals, CEU excluded. kallisto/GENCODE v27, 14,070 autosomal non-MHC genes, Phase 1 genotypes with MAF ≥ 5% and no LD pruning, 7 lab indicators + sex regressed out. 30 expression PCs × 5 genotype PCs → 2 canonical components; canonical correlations 0.963 and 0.766. Leave-one-out validation. eQTL comparison: 2,818 vs 2,732 eGenes at α = 1e-6, replicating in GTEx at 49.92% vs 51.02%.

3. **Fachrul et al., Communications Biology 2023** — RGStraP. Genetic PCs called directly from RNA-seq reads: STAR two-pass to GRCh38, GATK4, MAPQ ≥ 20, autosomal SNPs, HLA removed, MAF > 0.05, LD r² < 0.05, intersected with HapMap3. Nepal cohort of 376 (362 post-QC), validated on GEUVADIS. Genotype concordance 0.925 mean; GEUVADIS canonical correlations 0.994 / 0.942 / 0.752. Differential-expression inflation reduced from m = 0.935 to array-PC-equivalent levels.

Include a comparison table with rows for input data, whether genotypes are required, dimensionality reduction, what population structure is recovered from, and computational cost.

Then a **replicated vs deviated** table copied from the spec's Deviations section, plus the actual canonical correlations obtained, so the meeting audience can see exactly what was and was not reproduced.

Close with a section on what this enables: the transcriptional vector is now defined, Antoni's SNP-based clustering supplies the genetic vector, and the two can be fused into a single latent representation.

- [ ] **Step 2: Add it to the README documentation table**

In `README.md`, add a row to the table under `## Documentation`:

```markdown
| [docs/4_pcca_methods.md](docs/4_pcca_methods.md) | PCCA replication - methodology of the three population-structure papers and what this repo reproduces |
```

- [ ] **Step 3: Render the PDF**

Use the `pdf` skill to convert `docs/4_pcca_methods.md` into `docs/4_pcca_methods.pdf`, embedding the figures from `data/metadata/eda/pcca_*.png`.

- [ ] **Step 4: Final verification**

```bash
cd /Users/nk2/Desktop/biobank/tbooo
/opt/anaconda3/bin/python -m pytest tests/ -v
ls -la data/metadata/eda/pcca_*.png
ls -la data/metadata/geuvadis_pcca_coords.tsv
ls -la docs/4_pcca_methods.pdf
```

Expected: all tests pass; eight PCCA figures present (`fig1a`, `fig1b`, `fig2a`, `fig2b`, `production_462`, `panelB`, `panelC`, `panelD`); `geuvadis_pcca_coords.tsv` has 462 rows plus a header; the PDF exists.

**Do not commit.** Report what was produced and let the user decide.

---

## Notes for the implementer

- **The single most important number is the pair of canonical correlations.** Everything else can be approximately right and still be useful; if 0.963 / 0.766 are badly missed, the replication has not succeeded and the cause must be identified before the meeting.
- **Orthonormality of the genotype PCs is the most likely silent failure.** `cca` raises if its inputs are not orthonormal, which is deliberate — do not weaken that check.
- **Row alignment between expression and genotype is the second most likely silent failure.** `load_eigenvec` reorders to the cohort; do not bypass it.
- **Watch disk space.** `df -h .` between chromosomes. If free space drops below ~3 GB, stop and clean `data/raw/geuvadis/genotype/work/`.
