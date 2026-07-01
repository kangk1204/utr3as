"""
Classification of the alternative-splicing event type for Class III 3'UTR pairs.

Given two coding isoforms that share an identical CDS but differ by a splicing
event inside the 3'UTR (a Class III pair), we assign the canonical AS event type
by comparing their intron chains restricted to the 3'UTR:

    * IR    - intron retention: one isoform retains (as a single exon) an intron
              that the other splices out.
    * ES    - cassette / mutually-exclusive exon: one isoform includes an
              internal 3'UTR exon lying entirely within an intron of the other.
    * A3SS  - alternative 3' splice site: two introns share the 5' donor but use
              different 3' acceptors.
    * A5SS  - alternative 5' splice site: two introns share the 3' acceptor but
              use different 5' donors.
    * complex - none of the clean single-event patterns above (e.g. combined
              retention + cassette, or multi-intron rearrangements).

All calls are transcription-strand aware.
"""

from __future__ import annotations

from collections import Counter
from typing import List

from .gtf import Interval, Transcript
from .classify import utr3_introns, utr3_exons, _within

IR = "IR_intron_retention"
ES = "ES_cassette_or_MXE"
A3SS = "A3SS_alt_3prime_ss"
A5SS = "A5SS_alt_5prime_ss"
COMPLEX = "complex"

ALL_TYPES = [IR, ES, A3SS, A5SS, COMPLEX]


def classify_pair(a: Transcript, b: Transcript) -> str:
    """Return the AS event type distinguishing two Class III isoforms."""
    Ia, Ib = utr3_introns(a), utr3_introns(b)
    sIa, sIb = set(Ia), set(Ib)
    only_a = sorted(sIa - sIb)
    only_b = sorted(sIb - sIa)

    if not only_a and not only_b:
        return COMPLEX

    strand = a.strand

    # ---- Intron retention -------------------------------------------------
    def is_retention(host_introns_only, retainer: Transcript) -> bool:
        rex = utr3_exons(retainer)
        return any(_within(intr, ex) for intr in host_introns_only for ex in rex)

    a_retains = bool(only_b) and not only_a and is_retention(only_b, a)
    b_retains = bool(only_a) and not only_b and is_retention(only_a, b)
    if a_retains or b_retains:
        return IR

    # ---- Cassette / MXE exon ---------------------------------------------
    def has_cassette(exon_donor: Transcript, intron_host: Transcript) -> bool:
        host_introns = utr3_introns(intron_host)
        return any(_within(ex, intr)
                   for ex in utr3_exons(exon_donor) for intr in host_introns)

    if has_cassette(a, b) or has_cassette(b, a):
        return ES

    # ---- Alt 5'/3' splice site -------------------------------------------
    for ia in only_a:
        for ib in only_b:
            same_start = ia[0] == ib[0]   # shared donor on + strand
            same_end = ia[1] == ib[1]     # shared acceptor on + strand
            if same_start and not same_end:
                return A3SS if strand == "+" else A5SS
            if same_end and not same_start:
                return A5SS if strand == "+" else A3SS

    return COMPLEX


def primary_gene_as_type(pair_types: List[str]) -> str:
    """Collapse the AS types observed across a gene's Class III pairs to one call.

    Most frequent type wins; ties broken by specificity (ES/A3SS/A5SS over IR
    over complex).
    """
    if not pair_types:
        return COMPLEX
    specificity = {ES: 4, A3SS: 3, A5SS: 3, IR: 2, COMPLEX: 1}
    c = Counter(pair_types)
    return max(c.items(), key=lambda kv: (kv[1], specificity.get(kv[0], 0)))[0]
