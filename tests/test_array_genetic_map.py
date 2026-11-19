"""Tests for the genetic map lookup of the array stage."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from tbooo.pipeline.array import _find_genetic_map, _inject_cm, _read_genetic_map


class _Cfg:
    def __init__(self, reference_dir: Path):
        self.reference_dir = reference_dir


def _omni(maps: Path, chrom: str, rows: list[tuple[str, int, float]]) -> Path:
    d = maps / "interpolated_OMNI"
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"chr{chrom}.OMNI.interpolated_genetic_map.gz"
    with gzip.open(p, "wt") as f:
        for rs, pos, cm in rows:
            f.write(f"{rs} {pos} {cm}\n")
    return p


def test_the_map_of_chr1_is_not_the_map_of_chr10(tmp_path):
    maps = tmp_path / "genetic_maps"
    _omni(maps, "10", [("rs1", 100, 0.0)])
    one = _omni(maps, "1", [("rs2", 100, 0.0)])
    assert _find_genetic_map(maps, "1") == one
    assert _find_genetic_map(maps, "2") is None


def test_interpolated_maps_are_read_without_a_header(tmp_path):
    p = _omni(tmp_path / "genetic_maps", "22", [("rs9", 300, 2.0), ("rs8", 100, 0.0), ("rs7", 200, 1.0)])
    gmap = _read_genetic_map(p)
    assert list(gmap["pos"]) == [100, 200, 300]
    assert list(gmap["cm"]) == [0.0, 1.0, 2.0]


def test_phase3_style_text_maps_are_read_with_their_header(tmp_path):
    p = tmp_path / "genetic_map_GRCh37_chr22.txt"
    p.write_text("position COMBINED_rate(cM/Mb) Genetic_Map(cM)\n100 1.5 0.0\n200 1.5 0.5\n")
    gmap = _read_genetic_map(p)
    assert list(gmap["pos"]) == [100, 200] and list(gmap["cm"]) == [0.0, 0.5]


def test_inject_cm_interpolates_between_map_points(tmp_path):
    _omni(tmp_path / "genetic_maps", "22", [("rs8", 100, 0.0), ("rs9", 300, 2.0)])
    bim = tmp_path / "x.bim"
    bim.write_text("22\tsnpA\t0\t200\tA\tG\n22\tsnpB\t0\t300\tC\tT\n")
    _inject_cm(_Cfg(tmp_path), "22", bim)
    got = pd.read_csv(bim, sep="\t", header=None)
    assert list(got[2]) == [1.0, 2.0]
    assert list(got[1]) == ["snpA", "snpB"]
