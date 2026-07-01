"""Tests for the AS-3'UTR classifier.

Synthetic unit tests run with no data. The ground-truth tests parse a GENCODE
GTF if one is present under ``data/`` (download it with ``run.py`` or
``download_gencode.sh``); otherwise they skip.

Run:  pytest -q
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from utr3as.gtf import Transcript, parse_gencode  # noqa: E402
from utr3as import classify as C  # noqa: E402
from utr3as import as_types as A  # noqa: E402

GTF = os.path.join(os.path.dirname(__file__), "..", "data",
                   "gencode.v38.annotation.gtf.gz")


def _tx(strand, exons, cds, utr, stop=None):
    return Transcript(
        transcript_id="t", gene_id="g", gene_name="G", gene_type="protein_coding",
        transcript_type="protein_coding", chrom="chr1", strand=strand,
        exons=list(exons), exon_numbers=list(range(1, len(exons) + 1)),
        cds=list(cds), utr=list(utr), stop_codon=[stop] if stop else [])


# ---- 3'UTR architecture ----------------------------------------------------
def test_single_exon_3utr_is_canonical():
    t = _tx("+", [(100, 500)], [(100, 300)], [(304, 500)], stop=(301, 303))
    assert C.n_three_prime_utr_exons(t) == 1
    assert not C.has_3utr_intron(t)


def test_multi_exon_3utr_has_intron():
    t = _tx("+", [(100, 310), (400, 500)], [(100, 300)],
            [(304, 310), (400, 500)], stop=(301, 303))
    assert C.n_three_prime_utr_exons(t) == 2
    assert C.has_3utr_intron(t)
    assert C.utr3_intron_signature(t) == ((311, 399),)


def test_separating_intron_not_counted_as_3utr_intron():
    # 3'UTR entirely in a single downstream exon, separated from the coding exon
    # by an intron -> that separating intron is NOT internal to the 3'UTR.
    t = _tx("+", [(100, 303), (400, 500)], [(100, 300)],
            [(400, 500)], stop=(301, 303))
    assert C.n_three_prime_utr_exons(t) == 1
    assert C.has_3utr_intron(t) is False
    assert C.utr3_introns(t) == []
    assert C.has_3utr_intron(t) == (len(C.utr3_introns(t)) > 0)


def test_minus_strand_3utr_split():
    t = _tx("-", [(100, 200), (300, 460)], [(400, 460)],
            [(100, 200), (300, 396)], stop=(397, 399))
    assert C.n_three_prime_utr_exons(t) == 2
    assert C.has_3utr_intron(t)


# ---- AS event typing -------------------------------------------------------
def test_intron_retention_call():
    spliced = _tx("+", [(100, 310), (400, 500)], [(100, 300)],
                  [(304, 310), (400, 500)], stop=(301, 303))
    retained = _tx("+", [(100, 500)], [(100, 300)], [(304, 500)], stop=(301, 303))
    assert A.classify_pair(spliced, retained) == A.IR


def test_a3ss_call():
    a = _tx("+", [(100, 310), (400, 500)], [(100, 300)],
            [(304, 310), (400, 500)], stop=(301, 303))
    b = _tx("+", [(100, 310), (420, 500)], [(100, 300)],
            [(304, 310), (420, 500)], stop=(301, 303))
    assert A.classify_pair(a, b) == A.A3SS


# ---- ground truth on real annotation (optional) ----------------------------
def _classify_named(coding, name):
    txs = [t for t in coding.values() if t.gene_name == name
           and t.transcript_type == "protein_coding"]
    return C.classify_gene([C.architecture(t) for t in txs])


def test_known_example_genes_classify_correctly():
    if not os.path.exists(GTF):
        import pytest
        pytest.skip("no GENCODE GTF under data/; run `python run.py` first")
    tx = parse_gencode(GTF)
    coding = {k: v for k, v in tx.items() if v.is_coding}
    # known example genes, one per class
    assert _classify_named(coding, "KRAS").utr_class == C.CLASS_I
    assert _classify_named(coding, "CHEK1").utr_class == C.CLASS_II
    assert _classify_named(coding, "CAMK2B").utr_class == C.CLASS_III
