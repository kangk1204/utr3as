# utr3as

**Genome-wide classification of alternative splicing within 3′UTRs (AS-3′UTR) from a GENCODE annotation.**

Many protein-coding genes carry an **intron inside their 3′UTR**. `utr3as` scans a
GENCODE annotation, finds every protein-coding gene with an intron-containing
3′UTR, and classifies each into one of three mechanistic classes based on two
axes — the protein C-terminus (stop codon) and the mRNA 3′ end (poly(A) site):

| Class | protein C-terminus | mRNA 3′ end | mechanism |
|-------|--------------------|-------------|-----------|
| **I**   | distinct | distinct | alternative stop codon |
| **II**  | same     | distinct | alternative polyadenylation |
| **III** | **same** | **same** | **alternative splicing within the 3′UTR** |

Class III — two isoforms that start and end at the same points but differ by an
internal 3′UTR splicing event — is the case of interest, and each Class III event
is further typed as intron retention (IR), cassette / mutually-exclusive exon
(ES/MXE), or alternative 5′/3′ splice site (A5SS/A3SS).

The parser and classifier use only the Python standard library, so results are
exactly reproducible on any machine with Python ≥ 3.9.

## Install

```bash
git clone <your-repo-url> utr3as && cd utr3as
python -m pip install -r requirements.txt   # only needed to run the tests
```

No compiled dependencies; `run.py` works against a plain Python install.

## Quick start

```bash
python run.py                 # downloads GENCODE v38, classifies, writes results/
```

`run.py` fetches the annotation automatically (into `data/`) if it is not already
present, then writes the tables below and prints a summary. To use a different
release or a local file:

```bash
python run.py --gencode-version 44          # any GENCODE human release
python run.py --gtf /path/to/annotation.gtf.gz
python run.py --sweep                        # also print a stringency sweep
```

## Outputs (`results/`)

| file | contents |
|------|----------|
| `gene_classification.tsv` | every intron-3′UTR gene with its Class (I/II/III) and supporting counts |
| `class3_genes.tsv` | Class III genes with their primary AS event type |
| `class3_gene_list.txt` | plain list of Class III gene symbols |
| `class3_pairs.tsv` | each Class III isoform pair with its AS event type |
| `transcript_architecture.tsv` | per-transcript 3′UTR exon count, length, intron flag |
| `summary.json` | headline counts (also printed to stdout) |

## Options

| flag | default | meaning |
|------|---------|---------|
| `--gencode-version` | `38` | GENCODE human release to download |
| `--gtf` | – | use a local GTF(.gz) instead of downloading |
| `--level` | `protein_coding` | transcript-confidence filter: `all_coding`, `protein_coding`, `basic`, `tsl12`, `tsl1`, `strict` |
| `--group` | `stop_end` | Class III grouping unit: `stop_end` (stop codon + 3′ end), `stop`, or `cds` |
| `--sweep` | off | print class counts across all stringency levels |
| `--outdir` / `--datadir` | `results` / `data` | output and download locations |

The transcript-confidence filter (`--level`) is the single most consequential
parameter — the number of intron-3′UTR and Class III genes depends strongly on it.
It is an explicit, named option so any result can be reproduced exactly.

## Method

1. Parse GENCODE; keep transcripts with ≥ 1 CDS. Derive the **3′UTR** as the UTR
   segments 3′ of the CDS (strand-aware; GENCODE labels 5′/3′ UTR jointly).
2. A transcript has an **internal 3′UTR intron** when its 3′UTR spans ≥ 2 exons.
3. Collapse redundant transcripts to distinct **3′UTR-intron signatures** so counts
   are not inflated by GENCODE's many identical-structure transcript records.
4. Group a gene's transcripts by shared **stop codon + 3′ end**; assign **Class III**
   when such a group shows ≥ 2 distinct 3′UTR-intron structures with ≥ 1
   intron-contained isoform. Otherwise assign **Class I** (distinct stop codon) or
   **Class II** (same stop, distinct 3′ end = APA).
5. Type each Class III event (IR / ES-MXE / A3SS / A5SS) by comparing 3′UTR intron
   chains.

## Tests

```bash
pytest -q
```

Synthetic unit tests run with no data; ground-truth tests run automatically once a
GENCODE GTF has been downloaded.

## License

MIT — see [LICENSE](LICENSE).
