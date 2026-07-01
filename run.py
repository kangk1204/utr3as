"""Command-line driver: download a GENCODE annotation, classify AS-3'UTRs, write tables.

Examples
--------
    python run.py                                  # GENCODE v38, protein_coding -> results/
    python run.py --gencode-version 44 --level tsl12
    python run.py --gtf /path/to/annotation.gtf.gz --sweep

Outputs (in --outdir, default ``results/``):
    transcript_architecture.tsv   per-transcript 3'UTR exon count, length, intron flag
    gene_classification.tsv       per-gene Class I/II/III + supporting counts
    class3_genes.tsv              Class III genes with their primary AS event type
    class3_pairs.tsv              every Class III isoform pair with its AS event type
    class3_gene_list.txt          plain list of Class III gene symbols
    summary.json                  headline counts (also printed to stdout)
"""
import argparse
import csv
import json
import os
import sys
import urllib.request
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utr3as.gtf import parse_gencode  # noqa: E402
from utr3as.classify import (  # noqa: E402
    architecture, classify_gene, CLASS_I, CLASS_II, CLASS_III,
)
from utr3as.as_types import classify_pair, primary_gene_as_type  # noqa: E402
from utr3as.filters import apply_filter, DEFAULT_FILTER, SWEEP_ORDER  # noqa: E402

GENCODE_URL = ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
               "release_{v}/gencode.v{v}.annotation.gtf.gz")


def ensure_gencode(version, datadir):
    """Return a local path to the GENCODE GTF, downloading it if necessary."""
    os.makedirs(datadir, exist_ok=True)
    path = os.path.join(datadir, f"gencode.v{version}.annotation.gtf.gz")
    if os.path.exists(path):
        print(f"[data] using cached {path}")
        return path
    url = GENCODE_URL.format(v=version)
    print(f"[data] downloading {url}")
    tmp = path + ".part"
    urllib.request.urlretrieve(url, tmp)
    os.replace(tmp, path)
    print(f"[data] saved -> {path} ({os.path.getsize(path)/1e6:.0f} MB)")
    return path


def classify_set(coding, level, group):
    sub = apply_filter(coding, level)
    by_gene = defaultdict(list)
    for t in sub.values():
        by_gene[t.gene_id].append(t)
    arch = {tid: architecture(t) for tid, t in sub.items()}
    gene_classes = {}
    for gid, txs in by_gene.items():
        gc = classify_gene([arch[t.transcript_id] for t in txs], group)
        if gc is not None:
            gene_classes[gid] = gc
    return sub, by_gene, arch, gene_classes


def as_type_calls(class3, by_gene, arch, sub):
    """Assign AS event types to Class III pairs (grouped by stop codon + 3' end)."""
    pair_rows, gene_primary = [], {}
    clean_by_gene = Counter()
    for gc in class3:
        by_key = defaultdict(dict)  # (stop_pos, three_end) -> {intron_sig: transcript}
        for t in by_gene[gc.gene_id]:
            a = arch[t.transcript_id]
            by_key[(a.stop_pos, a.three_end)].setdefault(a.intron_sig, t)
        types = []
        for _key, rep in by_key.items():
            reps = list(rep.values())
            for i in range(len(reps)):
                for j in range(i + 1, len(reps)):
                    et = classify_pair(reps[i], reps[j])
                    types.append(et)
                    pair_rows.append([gc.gene_id, gc.gene_name,
                                      reps[i].transcript_id, reps[j].transcript_id, et])
        gene_primary[gc.gene_id] = primary_gene_as_type(types)
        clean = [t for t in types if t != "complex"]
        if clean:
            clean_by_gene[primary_gene_as_type(clean)] += 1
    return pair_rows, gene_primary, clean_by_gene


def write_tables(outdir, arch, gene_classes, class3, gene_primary, pair_rows):
    os.makedirs(outdir, exist_ok=True)

    with open(os.path.join(outdir, "transcript_architecture.tsv"), "w",
              newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["transcript_id", "gene_id", "gene_name", "strand",
                    "transcript_type", "n_3utr_exons", "utr3_len",
                    "has_3utr_intron", "n_3utr_introns"])
        for a in arch.values():
            w.writerow([a.transcript_id, a.gene_id, a.gene_name, a.strand,
                        a.transcript_type, a.n_3utr_exons, a.utr3_len,
                        int(a.has_3utr_intron), len(a.intron_sig)])

    with open(os.path.join(outdir, "gene_classification.tsv"), "w",
              newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "gene_name", "strand", "utr_class", "n_coding_tx",
                    "n_intron_3utr_tx", "n_distinct_alt_utr", "n_class3_pairs",
                    "n_stop_codons", "primary_as_type"])
        for gc in sorted(gene_classes.values(), key=lambda g: g.gene_name):
            w.writerow([gc.gene_id, gc.gene_name, gc.strand, gc.utr_class,
                        gc.n_coding_tx, gc.n_intron_3utr_tx, gc.n_distinct_alt_utr,
                        gc.n_class3_pairs, gc.n_stop_codons,
                        gene_primary.get(gc.gene_id, "")])

    with open(os.path.join(outdir, "class3_genes.tsv"), "w",
              newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "gene_name", "strand", "n_class3_pairs",
                    "n_distinct_alt_utr", "primary_as_type"])
        for gc in sorted(class3, key=lambda g: g.gene_name):
            w.writerow([gc.gene_id, gc.gene_name, gc.strand, gc.n_class3_pairs,
                        gc.n_distinct_alt_utr, gene_primary[gc.gene_id]])

    with open(os.path.join(outdir, "class3_pairs.tsv"), "w",
              newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "gene_name", "ref_transcript", "alt_transcript",
                    "as_event_type"])
        w.writerows(pair_rows)

    syms = sorted({gc.gene_name for gc in class3 if gc.gene_name})
    with open(os.path.join(outdir, "class3_gene_list.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(syms) + "\n")


def main():
    ap = argparse.ArgumentParser(
        description="Classify alternative splicing within 3'UTRs from a GENCODE annotation.")
    ap.add_argument("--gencode-version", default="38",
                    help="GENCODE human release number (default: 38)")
    ap.add_argument("--gtf", default=None,
                    help="path to a GENCODE GTF(.gz); skips the download")
    ap.add_argument("--level", default=DEFAULT_FILTER, choices=SWEEP_ORDER,
                    help="transcript-confidence filter (default: protein_coding)")
    ap.add_argument("--group", default="stop_end", choices=["stop_end", "stop", "cds"],
                    help="Class III grouping unit (default: stop_end)")
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--datadir", default="data")
    ap.add_argument("--sweep", action="store_true",
                    help="also print a transcript-stringency sensitivity sweep")
    args = ap.parse_args()

    gtf = args.gtf or ensure_gencode(args.gencode_version, args.datadir)
    print(f"[parse] {gtf}")
    tx = parse_gencode(gtf)
    coding = {k: v for k, v in tx.items() if v.is_coding}
    pc_genes = len({t.gene_id for t in coding.values()
                    if t.transcript_type == "protein_coding"})
    print(f"[parse] {len(coding):,} coding transcripts; {pc_genes:,} protein-coding genes")

    if args.sweep:
        print("\n[sweep] class counts by transcript-confidence level:")
        print(f"  {'level':14s}{'intron3UTR':>12}{'ClassI':>8}{'ClassII':>9}{'ClassIII':>10}")
        for lvl in SWEEP_ORDER:
            _, _, _, gcs = classify_set(coding, lvl, args.group)
            cc = Counter(g.utr_class for g in gcs.values())
            print(f"  {lvl:14s}{len(gcs):>12,}{cc[CLASS_I]:>8,}"
                  f"{cc[CLASS_II]:>9,}{cc[CLASS_III]:>10,}")

    sub, by_gene, arch, gene_classes = classify_set(coding, args.level, args.group)
    class3 = [gc for gc in gene_classes.values() if gc.utr_class == CLASS_III]
    pair_rows, gene_primary, clean_by_gene = as_type_calls(class3, by_gene, arch, sub)
    write_tables(args.outdir, arch, gene_classes, class3, gene_primary, pair_rows)

    cc = Counter(g.utr_class for g in gene_classes.values())
    n_intron = len(gene_classes)
    summary = {
        "gencode_version": args.gencode_version,
        "level": args.level,
        "class3_group": args.group,
        "protein_coding_genes": pc_genes,
        "intron_3utr_genes": n_intron,
        "pct_intron_3utr_genes": round(100 * n_intron / pc_genes, 2) if pc_genes else 0,
        "class_I": cc[CLASS_I],
        "class_II": cc[CLASS_II],
        "class_III": cc[CLASS_III],
        "class3_pairs": sum(gc.n_class3_pairs for gc in class3),
        "class3_as_type_clean": dict(clean_by_gene),
    }
    with open(os.path.join(args.outdir, "summary.json"), "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)

    print("\n================ SUMMARY ================")
    print(f"  intron-containing 3'UTR genes : {n_intron:,} "
          f"({summary['pct_intron_3utr_genes']}% of {pc_genes:,} protein-coding genes)")
    print(f"  Class I   (alt stop codon)    : {cc[CLASS_I]:,}")
    print(f"  Class II  (alt polyadenylation): {cc[CLASS_II]:,}")
    print(f"  Class III (alt splicing)      : {cc[CLASS_III]:,}")
    print(f"  Class III AS event types      : {dict(clean_by_gene)}")
    print(f"\n  wrote tables + summary.json to {args.outdir}/")


if __name__ == "__main__":
    main()
