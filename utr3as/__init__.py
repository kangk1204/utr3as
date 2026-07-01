"""utr3as - genome-wide classification of alternative splicing within 3'UTRs.

Identifies protein-coding genes whose 3'UTR contains an internal intron from a
GENCODE annotation and classifies them into three mechanistic classes:

    Class I   - alternative stop codon        (distinct C-terminus, distinct 3' end)
    Class II  - alternative polyadenylation    (same C-terminus, distinct 3' end)
    Class III - alternative splicing in 3'UTR  (same C-terminus, same 3' end)

Dependency-light (standard library only for parsing and classification).
"""

__version__ = "1.0.0"
