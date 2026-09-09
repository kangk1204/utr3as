"""Select Class III genes for subsequent analyses after initial classification."""

from collections import defaultdict

from .as_types import classify_pair

MAX_SHORT_DIFFERENCE_NT = 5


def symmetric_difference_nt(first, second):
    """Count bases present in exactly one union of 1-based inclusive intervals.

    Coverage changes at start and end + 1. Between consecutive boundaries,
    count the span only when one transcript, but not both, covers it. Duplicate
    or overlapping intervals within a transcript do not count bases twice.
    """
    changes = defaultdict(lambda: [0, 0])
    for axis, intervals in enumerate((first, second)):
        for start, end in intervals:
            if start > end:
                raise ValueError("interval start must not exceed end")
            changes[start][axis] += 1
            changes[end + 1][axis] -= 1
    coverage = [0, 0]
    previous = None
    total = 0
    for position in sorted(changes):
        if previous is not None and bool(coverage[0]) != bool(coverage[1]):
            total += position - previous
        coverage[0] += changes[position][0]
        coverage[1] += changes[position][1]
        previous = position
    return total


def select_class3_genes(class3, arch, transcripts):
    """Return retained genes, excluded genes, and reference-pair evidence.

    Reuse the classifier's representative pairs, including its reference and
    duplicate-signature choices. Exclude an entire gene if any such pair has
    a total 3'UTR sequence difference of <=5 nt. Initial classes are unchanged.
    """
    pair_rows = []
    excluded_ids = set()
    for gene in sorted(class3, key=lambda g: g.gene_name):
        for ref_id, alt_id in gene.class3_pairs:
            ref, alt = arch[ref_id], arch[alt_id]
            difference = symmetric_difference_nt(ref.three_utr, alt.three_utr)
            if difference <= MAX_SHORT_DIFFERENCE_NT:
                excluded_ids.add(gene.gene_id)
            pair_rows.append({
                "gene_id": gene.gene_id,
                "gene_name": gene.gene_name,
                "ref_transcript": ref_id,
                "alt_transcript": alt_id,
                "as_event_type": classify_pair(transcripts[ref_id], transcripts[alt_id]),
                "ref_utr3_nt": ref.utr3_len,
                "alt_utr3_nt": alt.utr3_len,
                "total_3utr_difference_nt": difference,
            })
    for row in pair_rows:
        row["include_in_subsequent_analyses"] = int(row["gene_id"] not in excluded_ids)
    retained = [g for g in class3 if g.gene_id not in excluded_ids]
    excluded = [g for g in class3 if g.gene_id in excluded_ids]
    return retained, excluded, pair_rows
