# utr3as

Genome-wide classification of alternative splicing within 3′UTRs (AS-3′UTR) from a GENCODE annotation.

`utr3as` identifies genes with an intron-containing 3′UTR in a GENCODE annotation
and assigns them to one of three classes. The default settings are
`--level protein_coding --group stop_end`.

The classifier evaluates the classes in the order III, I, then II. With the
default grouping, it uses these rules:

| Class | Assignment rule |
|-------|-----------------|
| III (alternative splicing) | At least one group of transcripts with the same stop-codon position and mRNA 3′ end contains at least two distinct 3′UTR-intron signatures, including an intron-containing isoform. |
| I (alternative stop codon) | Among genes not assigned to III, at least two stop-codon positions are present, and a 3′UTR segment of an intron-containing isoform overlaps the genomic span of the CDS of another isoform with a different stop. |
| II (alternative polyadenylation) | All remaining genes with an intron-containing 3′UTR. The code does not separately require a shared stop position and distinct 3′ ends for this class. |

Class I does not require distinct mRNA 3′ ends. For Class III, the default
`stop_end` rule requires matching stop and 3′-end coordinates, but does not
require identical complete CDS structures or transcript 5′ ends. These are
structural classifications based on the annotation.

Class III isoform comparisons are typed as intron retention (IR), cassette /
mutually-exclusive exon (ES/MXE), alternative 5′/3′ splice site (A5SS/A3SS), or
`complex` when none of these rules applies.

The parser and classifier require Python 3.9 or newer and use only the standard
library. To reproduce a run, use the same GTF file, code revision, and options.
The GTF checksum identifies the input file independently of its filename.

## Install

```bash
git clone https://github.com/kangk1204/utr3as.git && cd utr3as
python -m pip install -r requirements.txt   # only needed to run the tests
```

`run.py` works with a plain Python installation
(use `python3` instead of `python` on systems where that is the interpreter name).

## Quick start

```bash
python run.py                 # downloads GENCODE v38, classifies, writes results/
```

`run.py` fetches the annotation automatically (into `data/`) if it is not already
present, then writes the tables below and prints a summary. Each run parses and
classifies the GTF again, reusing the downloaded file when available. Download
size and runtime vary by release and machine. To use a different release or a
local file:

```bash
python run.py --gencode-version 44          # any GENCODE human release
python run.py --gtf /path/to/annotation.gtf.gz
python run.py --sweep                        # also print a stringency sweep
```

## Outputs (`results/`)

| file | contents |
|------|----------|
| `gene_classification.tsv` | one row per intron-3′UTR gene ID with its Class (I/II/III) and supporting counts |
| `class3_genes.tsv` | Class III genes with their primary AS event type |
| `class3_gene_list.txt` | plain list of Class III gene symbols |
| `class3_pairs.tsv` | all unordered pairs of representative 3′UTR-intron signatures within each stop/3′-end group of Class III genes, with AS event types |
| `transcript_architecture.tsv` | per-transcript 3′UTR exon count, length, intron flag |
| `summary.json` | headline counts (also printed to stdout) |

The `n_class3_pairs` fields and `summary.json`'s `class3_pairs` count use one
reference-to-alternative pair for each non-reference signature in a qualifying
Class III group. The reference has the fewest 3′UTR introns, with ties resolved
by intron coordinates. `class3_pairs.tsv` instead contains all pairwise
comparisons of representatives within stop/3′-end groups. With the GENCODE v38
defaults, these give 155 reference-to-alternative pairs and 172 pair-table rows.
The `ref_transcript` and `alt_transcript` columns follow enumeration order;
either transcript can carry the retained intron.

## Options

| flag | default | meaning |
|------|---------|---------|
| `--gencode-version` | `38` | GENCODE human release to download |
| `--gtf` | none | use a local GTF(.gz) instead of downloading |
| `--level` | `protein_coding` | transcript-confidence filter: `all_coding`, `protein_coding`, `basic`, `tsl12`, `tsl1`, `strict` |
| `--group` | `stop_end` | Class III gene-assignment grouping: `stop_end` (stop position + 3′ end), `stop` (stop position only), or `cds` (identical sorted CDS intervals) |
| `--sweep` | off | print class counts across all stringency levels |
| `--outdir` / `--datadir` | `results` / `data` | output and download locations |

The transcript filter and grouping affect the gene counts. All levels operate
on transcripts with at least one CDS feature. `protein_coding` additionally
requires `transcript_type == "protein_coding"`; `all_coding` includes other
CDS-bearing transcript biotypes. `basic` requires the GENCODE basic tag;
`tsl12` and `tsl1` require TSL 1/2 or TSL 1, respectively; `strict` requires
both basic and TSL 1. These four filters also require the protein-coding
transcript type. `basic` and the TSL filters are separate criteria, not
successive subsets of one another.

`--group cds` requires identical CDS intervals and allows different mRNA 3′
ends. `stop_end` requires matching stop and 3′-end coordinates and allows
differences elsewhere in the CDS. The driver performs AS event typing within
`stop_end` groups for all three settings. With `--group stop` or `--group cds`,
some pairs supporting a Class III assignment may therefore be absent from the
event-type outputs.

## Method

1. Parse GENCODE on chr1-22, chrX, chrY and chrM. Keep transcripts with at least
   one CDS feature, then apply `--level`. Coordinates are 1-based and inclusive.
2. Derive the 3′UTR from annotated UTR intervals downstream of the CDS in the
   direction of transcription. The stop position is the transcriptionally
   terminal stop-codon coordinate, or the terminal CDS coordinate when no stop
   codon is annotated. The mRNA 3′ end is the terminal 3′UTR coordinate.
3. An internal 3′UTR intron is a gap between consecutive 3′UTR-bearing exons.
   A transcript is marked as intron-containing when its annotated 3′UTR overlaps
   at least two exons. The separating intron upstream of a single downstream
   3′UTR exon is not counted as internal to the 3′UTR.
4. Group transcripts within each gene using `--group`. Within each group,
   collapse duplicate 3′UTR-intron signatures to representative transcripts.
   An empty signature represents a 3′UTR without an internal intron. Apply the
   Class III, I, and II rules above in that order. The Class I overlap uses the
   minimum-to-maximum genomic span of the other isoform's CDS, including gaps
   between CDS intervals; it does not test overlap with individual CDS exons.
5. For Class III genes, compare representative intron chains within stop/3′-end
   groups and assign strand-aware AS event types. The most frequent pair type
   gives the primary gene type, with specificity used to resolve ties.
   `summary.json`'s `class3_as_type_clean` first excludes `complex` pair calls
   and counts only genes with at least one remaining call; table primary types
   use all pair calls.

See [`utr3as/classify.py`](utr3as/classify.py) for classification,
[`utr3as/filters.py`](utr3as/filters.py) for transcript filters, and
[`run.py`](run.py) for pair enumeration and output summaries.

## GENCODE v38 reference counts

With GENCODE v38 (GRCh38), `--level protein_coding --group stop_end` gives:

| Quantity | Count |
|----------|------:|
| Protein-coding gene IDs | 19,681 |
| Intron-containing 3′UTR gene IDs | 2,714 |
| Class I | 1,611 |
| Class II | 976 |
| Class III | 127 |
| Class III reference-to-alternative pairs | 155 |

Gene classification counts use GENCODE gene IDs. In this release, the 2,714
intron-containing gene records correspond to 2,711 unique gene symbols because
CD99, CRLF2 and CSF2RA each have separate X and PAR_Y records. Gene-symbol
comparisons require a separately deduplicated symbol set.

The clean primary AS-type counts are IR 79, A3SS 30, ES/MXE 12, and A5SS 6
(127 genes in total). The CLI reports the complete Class III set without an
exon-length filter. If a downstream analysis uses a subset, report its gene
list and filtering criteria separately.

## Analysis scope

The input is an existing GENCODE annotation. RNA-seq transcript assembly and
comparisons with the published Riley or Chan sets are outside the CLI. Those
comparisons need separate reference files and gene-identifier mapping scripts.
Gene-symbol comparisons measure gene overlap; event-level agreement requires
transcript or splice-junction matching.

## Tests

```bash
pytest -q
```

Synthetic unit tests run with no data; ground-truth tests run automatically once a
GENCODE GTF has been downloaded.

## License

MIT. See [LICENSE](LICENSE).
