"""
Transcript-confidence filters.

The transcript set selected before classification is the single most consequential
parameter: the count of intron-containing 3'UTR genes and of Class III genes both
depend strongly on transcript-support stringency (TSL / APPRIS / GENCODE 'basic').
Each level below is defined purely on GENCODE attributes so results are exactly
reproducible. ``protein_coding`` (the default) keeps transcripts annotated
``transcript_type == "protein_coding"`` and is a good general starting point;
raise stringency with ``basic`` / ``tsl12`` / ``tsl1`` / ``strict``.
"""

from __future__ import annotations

from typing import Callable, Dict

from .gtf import Transcript

# Each filter takes a Transcript and returns True to keep it.
FILTERS: Dict[str, Callable[[Transcript], bool]] = {
    # everything with an annotated CDS (GENCODE comprehensive coding set)
    "all_coding": lambda t: t.is_coding,
    # transcript_type == protein_coding: excludes NMD, non_stop_decay,
    # retained_intron and other flagged biotypes.
    "protein_coding": lambda t: t.transcript_type == "protein_coding",
    # + GENCODE 'basic' tag (representative transcript set)
    "basic": lambda t: t.transcript_type == "protein_coding" and "basic" in t.tags,
    # + high transcript support level (well-supported splice junctions)
    "tsl12": lambda t: t.transcript_type == "protein_coding" and t.tsl in ("1", "2"),
    "tsl1": lambda t: t.transcript_type == "protein_coding" and t.tsl == "1",
    # highest stringency: basic + TSL1
    "strict": lambda t: (t.transcript_type == "protein_coding"
                         and "basic" in t.tags and t.tsl == "1"),
}

# Default transcript set.
DEFAULT_FILTER = "protein_coding"

# Order for the sensitivity sweep (loose -> strict).
SWEEP_ORDER = ["all_coding", "protein_coding", "basic", "tsl12", "tsl1", "strict"]


def apply_filter(transcripts: Dict[str, Transcript], level: str) -> Dict[str, Transcript]:
    keep = FILTERS[level]
    return {tid: t for tid, t in transcripts.items() if keep(t)}
