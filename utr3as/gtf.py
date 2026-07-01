"""
GENCODE GTF parser for the AS-3'UTR classification pipeline.

We parse the GENCODE annotation into a per-transcript representation that keeps
exactly the features needed to reconstruct 3'UTR architecture:

    * exon intervals (with GENCODE ``exon_number``)
    * CDS intervals
    * stop_codon interval(s)
    * UTR intervals (GENCODE labels 5' and 3' UTR jointly as ``UTR``; we split
      them by position relative to the CDS in :mod:`utr3as.classify`)

The parser is intentionally dependency-free (standard library only) so the
pipeline is trivially reproducible on any machine with Python >= 3.9.

Coordinate convention
---------------------
All coordinates are 1-based, fully-closed genomic coordinates, exactly as they
appear in the GTF. Strand is retained so downstream code can reason in
transcription order.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ---- attribute parsing -------------------------------------------------------

# Only a handful of attribute keys are needed. A targeted regex is far faster
# than fully tokenising the attribute column for every one of ~1.9M lines.
_ATTR_RE = re.compile(r'(\w+) "?([^";]+)"?;')
_TAG_RE = re.compile(r'tag "([^"]+)"')


def _parse_attributes(attr_field: str) -> Dict[str, str]:
    return {m.group(1): m.group(2) for m in _ATTR_RE.finditer(attr_field)}


Interval = Tuple[int, int]  # (start, end), 1-based inclusive


@dataclass
class Transcript:
    """All features of a single transcript needed for 3'UTR classification."""

    transcript_id: str
    gene_id: str
    gene_name: str
    gene_type: str
    transcript_type: str
    chrom: str
    strand: str
    tags: set = field(default_factory=set)
    tsl: str = ""  # transcript_support_level ("1".."5", "NA", or "")

    # feature intervals
    exons: List[Interval] = field(default_factory=list)          # (start, end)
    exon_numbers: List[int] = field(default_factory=list)         # parallel to exons
    cds: List[Interval] = field(default_factory=list)
    stop_codon: List[Interval] = field(default_factory=list)
    start_codon: List[Interval] = field(default_factory=list)
    utr: List[Interval] = field(default_factory=list)

    # ---- convenience -----------------------------------------------------
    @property
    def is_coding(self) -> bool:
        """Protein-coding = has at least one annotated CDS feature."""
        return len(self.cds) > 0

    @property
    def n_exons(self) -> int:
        return len(self.exons)

    def sorted_exons(self) -> List[Tuple[int, Interval]]:
        """Exons as ``(exon_number, (start, end))`` sorted by exon_number."""
        return sorted(zip(self.exon_numbers, self.exons), key=lambda x: x[0])

    def cds_signature(self) -> Tuple[Interval, ...]:
        """Sorted tuple of CDS intervals.

        Two transcripts with an identical CDS signature encode an identical
        protein and terminate at the same stop codon; this is the key used to
        detect isoforms that share coding sequence but differ in the 3'UTR
        (Class III).
        """
        return tuple(sorted(self.cds))


def parse_gencode(gtf_path: str,
                  coding_only: bool = False,
                  main_chroms_only: bool = True) -> Dict[str, Transcript]:
    """Parse a (optionally gzipped) GENCODE GTF into ``{transcript_id: Transcript}``.

    Parameters
    ----------
    gtf_path : str
        Path to ``gencode.vNN.annotation.gtf`` or ``.gtf.gz``.
    coding_only : bool
        If True, drop transcripts with no CDS after parsing.
    main_chroms_only : bool
        If True, keep only chr1-22, chrX, chrY, chrM (drop scaffolds/patches).
    """
    opener = gzip.open if gtf_path.endswith(".gz") else open
    main_chroms = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y", "M"]}

    transcripts: Dict[str, Transcript] = {}

    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _source, feature, start, end, _score, strand, _frame, attrs = f
            if feature == "gene":
                continue
            if main_chroms_only and chrom not in main_chroms:
                continue

            a = _parse_attributes(attrs)
            tid = a.get("transcript_id")
            if tid is None:
                continue

            t = transcripts.get(tid)
            if t is None:
                t = Transcript(
                    transcript_id=tid,
                    gene_id=a.get("gene_id", ""),
                    gene_name=a.get("gene_name", ""),
                    gene_type=a.get("gene_type", ""),
                    transcript_type=a.get("transcript_type", ""),
                    chrom=chrom,
                    strand=strand,
                )
                transcripts[tid] = t

            s, e = int(start), int(end)
            if feature == "exon":
                t.exons.append((s, e))
                t.exon_numbers.append(int(a.get("exon_number", 0)))
            elif feature == "CDS":
                t.cds.append((s, e))
            elif feature == "stop_codon":
                t.stop_codon.append((s, e))
            elif feature == "start_codon":
                t.start_codon.append((s, e))
            elif feature == "UTR":
                t.utr.append((s, e))
            elif feature == "transcript":
                # a transcript line can carry several `tag "..."` attributes
                for m in _TAG_RE.finditer(attrs):
                    t.tags.add(m.group(1))
                t.tsl = a.get("transcript_support_level", "")

    if coding_only:
        transcripts = {k: v for k, v in transcripts.items() if v.is_coding}

    return transcripts


def group_by_gene(transcripts: Dict[str, Transcript]) -> Dict[str, List[Transcript]]:
    """Group transcripts by ``gene_id``."""
    genes: Dict[str, List[Transcript]] = {}
    for t in transcripts.values():
        genes.setdefault(t.gene_id, []).append(t)
    return genes


if __name__ == "__main__":
    import sys
    import time

    path = sys.argv[1] if len(sys.argv) > 1 else \
        "data/raw/gencode.v38.annotation.gtf.gz"
    t0 = time.time()
    tx = parse_gencode(path)
    print(f"parsed {len(tx)} transcripts in {time.time()-t0:.1f}s")
    coding = {k: v for k, v in tx.items() if v.is_coding}
    genes = group_by_gene(coding)
    print(f"coding transcripts: {len(coding)}; coding genes: {len(genes)}")
