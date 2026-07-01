"""
3'UTR architecture and Class I/II/III classification.

Part 1 - per-transcript 3'UTR architecture
------------------------------------------
GENCODE labels 5' and 3' UTR jointly as ``UTR``; they are split here by position
relative to the CDS (transcription-aware). A transcript carries an *internal
3'UTR intron* when its 3'UTR spans >= 2 exons -- i.e. both the terminal exon and
an upstream exon carry 3'UTR sequence. Splicing out an intron located in the
3'UTR is what makes the 3'UTR span two exons; retaining that intron leaves a
single, longer 3'UTR exon.

Part 2 - gene-level Class I/II/III classification
-------------------------------------------------
Every gene that harbours >= 1 intron-containing-3'UTR transcript is partitioned
into exactly one class along two axes -- the protein C-terminus (stop codon) and
the mRNA 3' end (poly(A) site):

  * Class I  - Alternative Stop Codon (ASC):
        distinct C-terminus AND distinct 3' end. Operationally, the gene uses
        >= 2 distinct stop codons and an intron-containing-3'UTR isoform's 3'UTR
        intron overlaps the CDS of a different-stop isoform (the "3'UTR intron"
        is coding sequence in another isoform).
  * Class II - Alternative Polyadenylation (APA):
        same C-terminus, distinct 3' end -- isoform diversity driven by
        alternative 3' ends rather than 3'UTR splicing or alternative stop codons.
  * Class III - Alternative Splicing within the 3'UTR (AS-3'UTR):
        same C-terminus AND same 3' end, but internally 3'UTR-spliced: isoforms
        sharing both the stop codon and the 3' end include >= 2 distinct 3'UTR
        intron signatures with >= 1 intron-contained isoform (intron retention,
        cassette/MXE exon, or alternative 5'/3' splice site).

Assignment priority is III > I > II. Redundant transcripts are collapsed to
distinct 3'UTR-intron signatures before counting, so counts are not inflated by
the many identical-structure transcript records GENCODE lists per gene. The
grouping unit is configurable via ``classify_gene(class3_group=...)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from .gtf import Interval, Transcript

# ===========================================================================
# Part 1: per-transcript 3'UTR architecture
# ===========================================================================


def _overlaps(a: Interval, b: Interval) -> bool:
    return a[0] <= b[1] and b[0] <= a[1]


def _within(inner: Interval, outer: Interval) -> bool:
    return outer[0] <= inner[0] and inner[1] <= outer[1]


def cds_bounds(t: Transcript) -> Optional[Tuple[int, int]]:
    """Genomic (min_start, max_end) over all CDS intervals, or None."""
    if not t.cds:
        return None
    return min(s for s, _ in t.cds), max(e for _, e in t.cds)


def stop_codon_pos(t: Transcript) -> Optional[int]:
    """3'-most coding coordinate (transcription-aware) marking the protein C-terminus.

    Uses the annotated stop_codon when present, else the 3' end of the CDS.
    Identical values across two isoforms => identical C-terminus.
    """
    b = cds_bounds(t)
    if b is None:
        return None
    cmin, cmax = b
    if t.stop_codon:
        smin = min(s for s, _ in t.stop_codon)
        smax = max(e for _, e in t.stop_codon)
        return smax if t.strand == "+" else smin
    return cmax if t.strand == "+" else cmin


def three_prime_utr(t: Transcript) -> List[Interval]:
    """UTR intervals lying 3' of the CDS (the 3'UTR), transcription-aware."""
    b = cds_bounds(t)
    if b is None:
        return []
    cmin, cmax = b
    if t.strand == "+":
        return sorted(u for u in t.utr if u[0] > cmax)
    return sorted(u for u in t.utr if u[1] < cmin)


def five_prime_utr(t: Transcript) -> List[Interval]:
    """UTR intervals lying 5' of the CDS (the 5'UTR)."""
    b = cds_bounds(t)
    if b is None:
        return []
    cmin, cmax = b
    if t.strand == "+":
        return sorted(u for u in t.utr if u[1] < cmin)
    return sorted(u for u in t.utr if u[0] > cmax)


def utr3_introns(t: Transcript) -> List[Interval]:
    """Introns internal to the 3'UTR, i.e. between consecutive 3'UTR exons.

    An intron counts as a 3'UTR intron only when it is flanked by 3'UTR-bearing
    exons on both sides. Selecting introns merely by position 3' of the stop
    codon would wrongly include the intron that *separates* the coding exon from
    a downstream 3'UTR-only exon (which is not internal to the 3'UTR); this
    exon-adjacency definition avoids that.
    """
    tputr = three_prime_utr(t)
    if not tputr:
        return []
    ex = sorted(e for e in t.exons if any(_overlaps(e, u) for u in tputr))
    return [(ex[i][1] + 1, ex[i + 1][0] - 1)
            for i in range(len(ex) - 1)
            if ex[i + 1][0] - 1 >= ex[i][1] + 1]


def utr3_exons(t: Transcript) -> List[Interval]:
    """Whole exons that carry 3'UTR sequence (3' of the stop codon), strand-aware."""
    sp = stop_codon_pos(t)
    b = cds_bounds(t)
    if b is None:
        return []
    cmin, cmax = b
    ref = sp if sp is not None else (cmax if t.strand == "+" else cmin)
    out = []
    for ex in sorted(t.exons):
        if t.strand == "+" and ex[1] > ref:
            out.append(ex)
        elif t.strand == "-" and ex[0] < ref:
            out.append(ex)
    return out


def n_three_prime_utr_exons(t: Transcript) -> int:
    """Number of exons that overlap the 3'UTR."""
    tputr = three_prime_utr(t)
    if not tputr:
        return 0
    return sum(1 for ex in t.exons if any(_overlaps(ex, u) for u in tputr))


def three_prime_utr_length(t: Transcript) -> int:
    """Summed exonic length of the 3'UTR (nt)."""
    return sum(e - s + 1 for s, e in three_prime_utr(t))


def has_3utr_intron(t: Transcript) -> bool:
    """True iff the 3'UTR spans >= 2 exons (an internal 3'UTR intron exists)."""
    return t.is_coding and n_three_prime_utr_exons(t) >= 2


def utr3_intron_signature(t: Transcript) -> Tuple[Interval, ...]:
    """Canonical key for the 3'UTR splicing pattern (its set of 3'UTR introns).

    Two isoforms with the same signature have identically-spliced 3'UTRs; the
    empty tuple denotes an unspliced (single-exon / intron-retained) 3'UTR.
    Poly(A)-site differences alone do not change this signature, so collapsing
    by it removes GENCODE's redundant transcript records and isolates *splicing*
    differences from APA.
    """
    return tuple(utr3_introns(t))


@dataclass
class TxArchitecture:
    """Cached 3'UTR architecture summary for one transcript."""

    transcript_id: str
    gene_id: str
    gene_name: str
    strand: str
    is_coding: bool
    transcript_type: str
    n_3utr_exons: int
    utr3_len: int
    has_3utr_intron: bool
    stop_pos: Optional[int]
    three_end: Optional[int]           # 3'-most coordinate of the 3'UTR (mRNA 3' end)
    cds_sig: Tuple[Interval, ...]
    intron_sig: Tuple[Interval, ...]
    three_utr: Tuple[Interval, ...]


def three_prime_end(t: Transcript) -> Optional[int]:
    """3'-most coordinate of the 3'UTR (the mRNA 3' end / poly(A) site), strand-aware."""
    tputr = three_prime_utr(t)
    if not tputr:
        return None
    return max(e for _, e in tputr) if t.strand == "+" else min(s for s, _ in tputr)


def architecture(t: Transcript) -> TxArchitecture:
    return TxArchitecture(
        transcript_id=t.transcript_id,
        gene_id=t.gene_id,
        gene_name=t.gene_name,
        strand=t.strand,
        is_coding=t.is_coding,
        transcript_type=t.transcript_type,
        n_3utr_exons=n_three_prime_utr_exons(t),
        utr3_len=three_prime_utr_length(t),
        has_3utr_intron=has_3utr_intron(t),
        stop_pos=stop_codon_pos(t),
        three_end=three_prime_end(t),
        cds_sig=t.cds_signature(),
        intron_sig=utr3_intron_signature(t),
        three_utr=tuple(three_prime_utr(t)),
    )


# ===========================================================================
# Part 2: gene-level Class I / II / III classification
# ===========================================================================

CLASS_III = "III_alternative_splicing"
CLASS_I = "I_alternative_stop_codon"
CLASS_II = "II_alternative_polyadenylation"


@dataclass
class GeneClass:
    gene_id: str
    gene_name: str
    strand: str
    utr_class: str
    n_coding_tx: int
    n_intron_3utr_tx: int              # transcripts with a spliced (multi-exon) 3'UTR
    n_distinct_alt_utr: int            # distinct spliced 3'UTR structures (intron sigs != ())
    n_class3_pairs: int                # distinct (reference, alternative) 3'UTR pairs
    class3_pairs: List[Tuple[str, str]]  # representative transcript IDs per pair
    n_stop_codons: int


def classify_gene(arches: List[TxArchitecture],
                  class3_group: str = "stop_end") -> Optional[GeneClass]:
    """Classify one gene from its coding-transcript architectures.

    Returns None if the gene harbours no intron-containing-3'UTR transcript
    (i.e. it is not part of the AS-3'UTR gene set).

    ``class3_group`` selects the Class III grouping unit. Class III is defined on
    two axes -- protein C-terminus (stop codon) and mRNA 3' end -- so the default
    is ``"stop_end"``:

      * ``"stop_end"`` (default) - group isoforms sharing BOTH the same stop codon
        AND the same mRNA 3' end; Class III = internal 3'UTR splicing within such
        a group. This isolates alternative splicing from alternative polyadenylation.
      * ``"stop"`` - group by stop codon only (same C-terminus, any 3' end).
        Looser; conflates APA with splicing.
      * ``"cds"`` - group by identical CDS. Stricter; a same-CDS isoform pair is
        required, which can miss genes whose 3'UTR-spliced isoforms differ
        elsewhere in coding while sharing a C-terminus.
    """
    coding = [a for a in arches if a.is_coding]
    if not coding:
        return None
    intron_tx = [a for a in coding if a.has_3utr_intron]
    if not intron_tx:
        return None

    gene_id = coding[0].gene_id
    gene_name = coding[0].gene_name
    strand = coding[0].strand
    n_stops = len({a.stop_pos for a in coding if a.stop_pos is not None})

    # ---- Group isoforms by the chosen unit --------------------------------
    def group_key(a):
        if class3_group == "stop_end":
            return (a.stop_pos, a.three_end)   # same C-terminus AND same 3' end
        if class3_group == "stop":
            return a.stop_pos
        return a.cds_sig

    by_stop: Dict[object, List[TxArchitecture]] = {}
    for a in coding:
        by_stop.setdefault(group_key(a), []).append(a)

    # ---- Class III: same C-terminus & same 3' end, but internally 3'UTR- ----
    #      spliced: >= 2 distinct 3'UTR intron signatures with >= 1 intron-
    #      contained isoform. Requiring the shared 3' end excludes APA (Class II)
    #      and isolates pure alternative splicing within the 3'UTR.
    class3_pairs: List[Tuple[str, str]] = []
    alt_sigs: set = set()
    for stop, group in by_stop.items():
        rep: Dict[Tuple[Interval, ...], TxArchitecture] = {}
        for a in group:
            rep.setdefault(a.intron_sig, a)
        if len(rep) < 2 or not any(a.has_3utr_intron for a in rep.values()):
            continue
        # reference = the isoform with the fewest 3'UTR introns (most canonical)
        sigs = sorted(rep.keys(), key=lambda s: (len(s), s))
        ref = rep[sigs[0]]
        for s in sigs[1:]:
            class3_pairs.append((ref.transcript_id, rep[s].transcript_id))
            alt_sigs.add((stop, s))

    if class3_pairs:
        return GeneClass(
            gene_id=gene_id, gene_name=gene_name, strand=strand,
            utr_class=CLASS_III, n_coding_tx=len(coding),
            n_intron_3utr_tx=len(intron_tx),
            n_distinct_alt_utr=len(alt_sigs),
            n_class3_pairs=len(class3_pairs), class3_pairs=class3_pairs,
            n_stop_codons=n_stops,
        )

    # ---- Class I: alternative stop codon linking 3'UTR intron to coding ----
    class1 = False
    if n_stops >= 2:
        for it in intron_tx:
            for other in coding:
                if other.stop_pos == it.stop_pos or not other.cds_sig:
                    continue
                ocds_min = min(s for s, _ in other.cds_sig)
                ocds_max = max(e for _, e in other.cds_sig)
                if any(us <= ocds_max and ocds_min <= ue for us, ue in it.three_utr):
                    class1 = True
                    break
            if class1:
                break
    if class1:
        return GeneClass(
            gene_id=gene_id, gene_name=gene_name, strand=strand,
            utr_class=CLASS_I, n_coding_tx=len(coding),
            n_intron_3utr_tx=len(intron_tx), n_distinct_alt_utr=0,
            n_class3_pairs=0, class3_pairs=[], n_stop_codons=n_stops,
        )

    # ---- Class II: alternative polyadenylation / remaining -----------------
    return GeneClass(
        gene_id=gene_id, gene_name=gene_name, strand=strand,
        utr_class=CLASS_II, n_coding_tx=len(coding),
        n_intron_3utr_tx=len(intron_tx), n_distinct_alt_utr=0,
        n_class3_pairs=0, class3_pairs=[], n_stop_codons=n_stops,
    )
