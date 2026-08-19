"""Generic, registry-driven ingestion of externally-held cohorts.

This package contains NO cohort-specific knowledge. Everything about a given
cohort — where its data lives, how samples are named, how metadata is parsed —
comes from a `cohorts.yaml` registry that is never committed. See
`cohorts.example.yaml` for the schema.
"""

from __future__ import annotations

from tbooo.cohorts.registry import Cohort, RegistryError, Sample, SourceFile, load_registry

__all__ = ["Cohort", "Sample", "SourceFile", "RegistryError", "load_registry"]
