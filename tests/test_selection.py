"""Gene-level selection, interval arithmetic, and CLI output contracts."""

import csv
import json
import os
from pathlib import Path
import random
import subprocess
import sys

import pytest

from run import as_type_calls, classify_set
from utr3as.classify import architecture, classify_gene
from utr3as.gtf import Transcript, parse_gencode
from utr3as.selection import select_class3_genes, symmetric_difference_nt


@pytest.mark.parametrize("first,second,expected", [
    ([], [], 0),
    ([(1, 1)], [], 1),
    ([(1, 20)], [(1, 15)], 5),
    ([(1, 20)], [(1, 14)], 6),
    ([(1, 5)], [(1, 5)], 0),
    ([(1, 5)], [(5, 9)], 8),
    ([(1, 4)], [(10, 14)], 9),  # total 4 + 5, not net length difference 1
    ([(1, 5), (3, 8), (1, 5)], [(1, 8)], 0),
    ([(1, 3), (4, 6)], [(1, 6)], 0),
])
def test_symmetric_difference_boundaries(first, second, expected):
    assert symmetric_difference_nt(first, second) == expected
    assert symmetric_difference_nt(second, first) == expected


def test_interval_sweep_matches_independent_base_set_oracle():
    rng = random.Random(731)
    for _ in range(200):
        inputs = [[tuple(sorted(rng.sample(range(1, 50), 2)))
                   for _ in range(rng.randrange(8))] for _ in range(2)]
        bases = [{base for start, end in intervals for base in range(start, end + 1)}
                 for intervals in inputs]
        assert symmetric_difference_nt(*inputs) == len(bases[0] ^ bases[1])


def test_reversed_interval_rejected():
    with pytest.raises(ValueError, match="start"):
        symmetric_difference_nt([(5, 4)], [])


def transcript(tid, acceptor=300, gene="g", strand="+", end=500):
    t = Transcript(
        transcript_id=tid, gene_id=gene, gene_name=gene.upper(),
        gene_type="protein_coding", transcript_type="protein_coding",
        chrom="chr1", strand=strand, exons=[(100, 250), (acceptor, end)],
        exon_numbers=[1, 2], cds=[(100, 200)], stop_codon=[(201, 203)],
        utr=[(204, 250), (acceptor, end)])
    if strand == "-":
        for feature in ("exons", "cds", "stop_codon", "utr"):
            setattr(t, feature, [(601 - e, 601 - s) for s, e in getattr(t, feature)])
    return t


def select(transcripts, group="stop_end"):
    tx = {t.transcript_id: t for t in transcripts}
    arch = {tid: architecture(t) for tid, t in tx.items()}
    gc = classify_gene(list(arch.values()), group)
    return gc, select_class3_genes([gc], arch, tx)


@pytest.mark.parametrize("difference,excluded", [(1, True), (3, True), (5, True), (6, False)])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_inclusive_gene_exclusion_threshold(difference, excluded, strand):
    _, (retained, removed, rows) = select([
        transcript("ref", strand=strand),
        transcript("alt", acceptor=300 + difference, strand=strand)])
    assert len(removed) == int(excluded)
    assert len(retained) == int(not excluded)
    assert rows[0]["total_3utr_difference_nt"] == difference
    assert rows[0]["include_in_subsequent_analyses"] == int(not excluded)


def test_one_short_pair_excludes_entire_gene_including_long_pairs():
    gc, (retained, removed, rows) = select([
        transcript("ref"), transcript("short", 303), transcript("long", 348)])
    assert gc.n_class3_pairs == 2  # original classification remains intact
    assert retained == [] and removed == [gc]
    assert {r["total_3utr_difference_nt"] for r in rows} == {3, 48}
    assert all(r["include_in_subsequent_analyses"] == 0 for r in rows)


def test_only_classifier_reference_pairs_determine_exclusion():
    # The two alternatives differ by 3 nt, but each differs >5 nt from reference.
    gc, (retained, removed, rows) = select([
        transcript("ref"), transcript("alt1", 340), transcript("alt2", 343)])
    assert gc.class3_pairs == [("ref", "alt1"), ("ref", "alt2")]
    assert {r["total_3utr_difference_nt"] for r in rows} == {40, 43}
    assert retained == [gc] and not removed


def test_duplicate_signature_preserves_existing_first_representative():
    gc, (_, _, rows) = select([
        transcript("first"), transcript("duplicate"), transcript("alt", 310)])
    assert gc.class3_pairs == [("first", "alt")]
    assert [(r["ref_transcript"], r["alt_transcript"]) for r in rows] == gc.class3_pairs


def test_short_pair_in_another_stop_end_group_excludes_whole_gene():
    gc, (retained, removed, rows) = select([
        transcript("r1"), transcript("a1", 320),
        transcript("r2", end=600), transcript("a2", 303, end=600)])
    assert gc.n_class3_pairs == 2
    assert not retained and removed == [gc]
    assert all(r["include_in_subsequent_analyses"] == 0 for r in rows)


@pytest.mark.parametrize("group", ["stop", "cds"])
def test_nondefault_group_uses_its_classifier_pairs(group):
    gc, (_, removed, rows) = select([
        transcript("ref"), transcript("alt", 303, end=501)], group=group)
    assert gc.n_class3_pairs == 1
    assert rows[0]["total_3utr_difference_nt"] == 4
    assert removed == [gc]


def write_gtf(path, transcripts):
    with path.open("w") as f:
        for t in transcripts:
            attrs = (f'gene_id "{t.gene_id}"; gene_name "{t.gene_name}"; '
                     f'gene_type "protein_coding"; transcript_id "{t.transcript_id}"; '
                     'transcript_type "protein_coding";')
            for feature, intervals in (("exon", t.exons), ("CDS", t.cds),
                                       ("stop_codon", t.stop_codon), ("UTR", t.utr)):
                for start, end in intervals:
                    f.write(f"chr1\ttest\t{feature}\t{start}\t{end}\t.\t{t.strand}\t.\t{attrs}\n")


def read_tsv(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_cli_writes_both_sets_and_replaces_stale_selected_outputs(tmp_path):
    gtf, out = tmp_path / "test.gtf", tmp_path / "results"
    short = [transcript("s0", gene="short"), transcript("s1", 303, gene="short")]
    long = [transcript("l0", gene="long"), transcript("l1", 320, gene="long")]
    write_gtf(gtf, short + long)
    script = Path(__file__).resolve().parents[1] / "run.py"
    command = [sys.executable, str(script), "--gtf", str(gtf), "--outdir", str(out)]
    result = subprocess.run(command, check=True, capture_output=True, text=True)
    summary = json.loads((out / "summary.json").read_text())
    assert summary["class_III"] == 2 and summary["class3_pairs"] == 2
    assert summary["class3_selection"]["retained_genes"] == 1
    assert summary["class3_selection"]["excluded_gene_symbols"] == ["SHORT"]
    assert "retained for subsequent analyses: 1" in result.stdout
    assert len(read_tsv(out / "class3_genes.tsv")) == 2
    assert (out / "class3_retained_gene_list.txt").read_text() == "LONG\n"
    assert [x["gene_name"] for x in read_tsv(out / "class3_excluded_genes.tsv")] == ["SHORT"]
    rows = read_tsv(out / "class3_reference_pairs.tsv")
    assert sum(int(x["include_in_subsequent_analyses"]) for x in rows) == 1
    # Reusing an output directory must not leave the old retained list behind.
    write_gtf(gtf, short)
    subprocess.run(command, check=True, capture_output=True, text=True)
    assert (out / "class3_retained_gene_list.txt").read_text() == ""
    assert read_tsv(out / "class3_retained_genes.tsv") == []
    # A run with no Class III genes still produces every table with its header.
    write_gtf(gtf, [transcript("single")])
    subprocess.run(command, check=True, capture_output=True, text=True)
    assert read_tsv(out / "class3_reference_pairs.tsv") == []
    assert read_tsv(out / "class3_excluded_genes.tsv") == []


def test_gencode_v38_selected_membership_and_reference_pairs():
    gtf = Path(os.environ.get("UTR3AS_TEST_GTF", str(
        Path(__file__).resolve().parents[1] / "data/gencode.v38.annotation.gtf.gz")))
    if not gtf.exists():
        pytest.skip("no GENCODE v38 GTF; set UTR3AS_TEST_GTF or download into data/")
    tx = parse_gencode(str(gtf))
    sub, by_gene, arch, genes = classify_set(
        {tid: t for tid, t in tx.items() if t.is_coding}, "protein_coding", "stop_end")
    class3 = [g for g in genes.values() if g.class3_pairs]
    retained, excluded, rows = select_class3_genes(class3, arch, sub)
    expected = read_tsv(Path(__file__).parent / "data/gencode_v38_retained_genes.tsv")
    assert {(g.gene_id, g.gene_name) for g in retained} == {
        (r["gene_id"], r["gene_name"]) for r in expected}
    assert len(class3) == 127 and len(retained) == 121 and len(rows) == 155
    assert {g.gene_name for g in excluded} == {
        "GNAS", "KCNH4", "KLHDC4", "KRTCAP3", "PNCK", "TEX51"}
    assert sum(g.n_class3_pairs for g in retained) == 147
    assert sum(r["include_in_subsequent_analyses"] for r in rows) == 147
    assert sum(r["total_3utr_difference_nt"] <= 5 for r in rows) == 6
    assert {r["total_3utr_difference_nt"] for r in rows if r["gene_name"] == "GNAS"} == {3, 45, 48}
    assert all(r["include_in_subsequent_analyses"] for r in rows if r["gene_name"] == "CACTIN")
    all_pairs, primary, _ = as_type_calls(class3, by_gene, arch, sub)
    assert len(all_pairs) == 172
    from collections import Counter
    assert Counter(primary[g.gene_id] for g in retained) == {
        "IR_intron_retention": 79, "A3SS_alt_3prime_ss": 26,
        "ES_cassette_or_MXE": 11, "A5SS_alt_5prime_ss": 5}
