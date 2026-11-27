# Supported nomenclatures

Variation Normalizer accepts HGVS, free-text, and gnomAD-style VCF input for
GRCh37 and GRCh38. The supported change types depend on the input format and
coordinate type:

| Input | Coordinate type | Supported change types |
| --- | --- | --- |
| HGVS or free text | Protein (p.) | Substitution, stop gain, reference agree, deletion, insertion, deletion-insertion |
| HGVS or free text | Coding DNA (c.) | Substitution, reference agree, deletion, insertion, deletion-insertion |
| HGVS or free text | Genomic (g.) | Substitution, reference agree, deletion, ambiguous deletion, insertion, deletion-insertion, duplication, ambiguous duplication |
| gnomAD-style VCF | Genomic (g.) | Substitution, reference agree, deletion, insertion, deletion-insertion |
| Free text | Categorical | Amplification |

## Input formats

| Format | Description | Example |
| --- | --- | --- |
| HGVS | An HGVS expression with an accession and coordinate type | `NP_004324.2:p.Val600Glu` |
| Free text | A gene name with a protein, coding DNA, or genomic change, or a categorical variant | `BRAF V600E`, `BRAF c.1799T>A`, `BRAF amplification` |
| gnomAD-style VCF | Chromosome, position, reference allele, and alternate allele, separated by hyphens | `7-140753336-A-T` |

## Examples

### Protein (p.)

| Variant type | Examples |
| --- | --- |
| Substitution | `BRAF V600E`, `NP_004324.2:p.Val600Glu` |
| Stop gain | `NP_000542.1:p.Tyr185Ter`, `NP_000542.1:p.Tyr185*` |
| Reference agree | `NP_000542.1:p.Pro61=` |
| Deletion | `EGFR L747_T751del`, `NP_004439.2:p.Leu755_Thr759del` |
| Insertion | `EGFR D770_N771insGL`, `NP_005219.2:p.Asp770_Asn771insGlyLeu` |
| Deletion-insertion | `EGFR L747_T751delinsP`, `NP_001333827.1:p.Leu747_Thr751delinsPro` |

### Coding DNA (c.)

| Variant type | Examples |
| --- | --- |
| Substitution | `BRAF c.1799T>A`, `BRAF V600E (c.1799T>A)`, `NM_004333.4:c.1799T>A`, `ENST00000288602.10:c.1799T>A` |
| Reference agree | `BRAF c.1799=`, `BRAF V600E c.1799=`, `NM_004333.4:c.1799=` |
| Deletion | `ERBB2 c.2264_2278delTGAGGGAAAACACAT`, `ERBB2 Leu755_Thr759del c.2264_2278del`, `NM_004448.3:c.2264_2278del` |
| Insertion | `VHL C77_N78insL (c.230_231insTCT)`, `NM_000551.3:c.230_231insTCT` |
| Deletion-insertion | `VHL I206fs (c.615delinsAA)`, `NM_000551.3:c.615delinsAA` |

### Genomic (g.)

| Variant type | Examples |
| --- | --- |
| Substitution | `BRAF g.140453136A>T`, `BRAF V600E (g.140453136A>T)`, `NC_000007.13:g.140453136A>T`, `7-140753336-A-T` |
| Reference agree | `BRAF g.140453136=`, `NC_000007.13:g.140453136=`, `7-140753336-A-A` |
| Deletion | `NC_000003.12:g.10146527_10146528del`, `3-10146526-TCT-T` |
| Ambiguous deletion | `NC_000023.11:g.(31060227_31100351)_(33274278_33417151)del`, `NC_000023.11:g.(?_31120496)_(33339477_?)del`, `NC_000023.11:g.(?_18575354)_18653629del`, `NC_000006.12:g.133462764_(133464858_?)del` |
| Insertion | `ERBB2 g.37880993_37880994insGCTTACGTGATG`, `NC_000017.10:g.37880993_37880994insGCTTACGTGATG`, `17-37880993-G-GGCTTACGTGATG` |
| Deletion-insertion | `BRAF g.140453135_140453136delinsAT`, `NC_000007.13:g.140453135_140453136delinsAT`, `3-10149938-C-AA` |
| Duplication | `NC_000003.12:g.49531262dup`, `NC_000016.10:g.2087938_2087948dup` |
| Ambiguous duplication | `NC_000023.11:g.(31060227_31100351)_(33274278_33417151)dup`, `NC_000020.11:g.(?_30417576)_(31394018_?)dup`, `NC_000023.11:g.(?_154021812)_154092209dup`, `NC_000023.11:g.154021812_(154092209_?)dup` |

### Categorical

| Variant type | Examples |
| --- | --- |
| Amplification | `BRAF amplification` |

See [HGVS deletion and duplication modes](hgvs_dup_del_mode.md) for copy number interpretation.
