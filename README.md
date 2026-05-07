# TRBSeqMock — Immunorepertoire Sequence Simulation Platform

A Python-based simulation tool for generating **T-cell receptor beta (TRB)** and
**immunoglobulin heavy chain (IGH)** rearranged sequences from germline gene
reference databases.  The tool models the complete V(D)J recombination process
including gene selection, exonuclease deletion, non-templated nucleotide
insertion, sequencing error, and (for IGH) somatic hypermutation (SHM).

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Architecture](#architecture)
3. [Simulation Pipeline](#simulation-pipeline)
4. [SHM Model (IGH only)](#shm-model-igh-only)
5. [Gene Selection & Pool Construction](#gene-selection--pool-construction)
6. [Usage](#usage)
7. [Output Files](#output-files)
8. [Dependencies](#dependencies)
9. [References](#references)

---

## Project Overview

**Purpose**: Generate realistic, fully annotated, rearranged immune receptor
sequences for evaluating germline V/J gene inference tools (IMPre, TIgGER,
IgDiscover, partis) and other AIRR-seq bioinformatics pipelines.

**Scope**:

| Chain | Recombination | SHM | Status |
|-------|--------------|-----|--------|
| TRB   | V-D-J        | —   | Stable |
| IGH   | V-D-J        | ✓   | Stable |

---

## Architecture

```
TRBSeqMock/
├── TRBSeqMock.py              # Unified CLI entry point & module dispatcher
├── Source_TRB/                # TRB-specific simulation module
│   ├── TRBSeqMock.py
│   ├── data_loader/           # Reference gene loading & distribution parsing
│   │   ├── config.py          # File paths for TRB reference data
│   │   ├── read_fasta.py      # IMGT FASTA parser & gene pool builder
│   │   ├── read_distributions.py  # Deletion/insertion distribution reader
│   │   └── __init__.py
│   └── seq_sim/               # Core recombination simulation engine
│       ├── main.py            # Top-level simulation orchestrator
│       ├── gene_select.py     # Random gene selection from pool
│       ├── delete_process.py  # Exonucleolytic deletion simulation
│       ├── create_insertion.py# N-nucleotide addition simulation
│       ├── gene_assemble.py   # Fragment concatenation
│       └── __init__.py
├── Source_IGH/                # IGH-specific simulation module
│   ├── IGHSeqMock.py          # Entry script with SHM support
│   ├── data_loader/           # IGH-specific config & gene loading
│   └── seq_sim/               # IGH recombination + SHM integration
├── Source_SHM/                # ★ Somatic Hypermutation module
│   ├── __init__.py
│   ├── shm_constants.py       # WRC/GYW motifs, IMGT CDR/FR boundaries,
│   │                          #   substitution weight matrix, default config
│   └── shm_simulator.py       # SHMConfig dataclass, scan_hotspots(),
│                              #   apply_shm(), compute_hotspot_enrichment()
├── Source_error/              # Sequencing error simulator
│   ├── __init__.py
│   └── error_simulator.py     # Per-base substitution (p=0.001), numba-accelerated
│                              #   with pure-Python fallback
├── Source_write/              # Output writer module
│   ├── __init__.py
│   └── write.py               # FASTA, TSV, SHM stats, ground truth JSON
├── Data/
│   ├── Data_Ref/              # IMGT germline reference FASTA files
│   └── Data_Distribute/       # V(D)J deletion/insertion length distributions
└── Ref/                       # Additional reference data
```

---

## Simulation Pipeline

### Step-by-step process (for each simulated sequence)

```
 ┌──────────────┐
1│ Select V/D/J │  Random draw from gene pool (uniform probability)
 │ genes        │
 └──────┬───────┘
        ▼
 ┌──────────────┐
2│ V 3' deletion│  Random draw from deletion length distribution (V3)
 │              │  → truncate V gene 3' end
 └──────┬───────┘
        ▼
 ┌──────────────┐
3│ D 5' deletion│  Random draw → truncate D gene 5' end
 │              │  → produces intermediate d_gene_mid
 └──────┬───────┘
        ▼
 ┌──────────────┐
4│ D 3' deletion│  Random draw → truncate D gene 3' end from d_gene_mid
 │              │  → produces final d_gene
 └──────┬───────┘
        ▼
 ┌──────────────┐
5│ J 5' deletion│  Random draw → truncate J gene 5' end
 └──────┬───────┘
        ▼
 ┌──────────────┐
6│ VD insertion │  Random draw from insertion length distribution (DJ)
 │              │  → generates random nucleotide sequence (A/C/G/T)
 └──────┬───────┘
        ▼
 ┌──────────────┐
7│ DJ insertion │  Same process → generates DJ junction insert
 └──────┬───────┘
        ▼
 ┌──────────────┐
8│   Assemble   │  V_processed + VD_insert + D_processed + DJ_insert + J_processed
 └──────┬───────┘
        ▼
 ┌──────────────┐
9│  Sequencing  │  Per-base substitution error (p=0.001 default)
 │  error       │  Numba JIT acceleration (graceful fallback if unavailable)
 └──────┬───────┘
        ▼
 ┌──────────────┐
10│  SHM (IGH   │  ★ Position-specific probability model applied to IGHV
 │   only)      │  region only (FR1-CDR1-FR2-CDR2-FR3, IMGT 1–104)
 └──────────────┘
```

---

## SHM Model (IGH only)

### Design principles

The SHM module implements a **position-specific probability model** grounded in
the biochemistry of Activation-Induced Deaminase (AID) and empirical
observations from AIRR-seq data.  SHM is applied **only to the IGHV region**
(FR1–CDR1–FR2–CDR2–FR3, IMGT positions 1–104); CDR3, D, J, and N-insertion
regions are left unmutated.

### Three-component mutation rate

For each nucleotide position *i* in the IGHV region:

<p align="center">
<b>P<sub>mut</sub>(i) = base_rate × f<sub>hotspot</sub>(i) × f<sub>CDR</sub>(i)</b>
</p>

| Factor | Symbol | Description | Default Value | Reference |
|--------|--------|-------------|---------------|-----------|
| Global SHM rate | *μ* | Per-base mutation probability | 0.05 (5%) | Zheng et al., 2005; Kleinstein et al., 2003 |
| Hotspot multiplier | *f<sub>hotspot</sub>* | ×4.0 if position is in WRC/GYW motif; else background rate (0.002) | 4.0 | Rogozin & Diaz, 2004 |
| CDR bias | *f<sub>CDR</sub>* | ×1.5 if position is in CDR1 or CDR2 (IMGT) | 1.5 | Sheng et al., 2017 |

### Hotspot motifs

The model scans germline IGHV sequences for four classes of AID-targeting motifs:

| Motif | Pattern | Hotspot Base | Reference |
|-------|---------|-------------|-----------|
| WRC | [AT][AG]C | C (position 2) | Rogozin & Diaz, 2004 |
| GYW | G[CT][AT] | G (position 0) | ibid. |
| DGYW | [AGT]G[CT][AT] | G (position 1) | ibid. (extended) |
| WRCH | [AT][AG]C[ACT] | C (position 2) | ibid. (extended) |

### Substitution bias

When a mutation event occurs, the replacement base is selected with
transition/transversion weights reflecting AID biochemistry:

| From | To A | To C | To G | To T | Dominant pathway |
|------|------|------|------|------|------------------|
| A | — | 0.30 | 0.35 | 0.35 | A to G (transition) |
| C | 0.15 | — | 0.10 | **0.75** | C to T (AID deamination) |
| G | **0.55** | 0.20 | — | 0.25 | G to A (complementary C to T) |
| T | 0.30 | 0.35 | 0.35 | — | balanced |

*Reference*: Di Noia & Neuberger, 2007; Sheng et al., 2017.

### CDR/FR boundaries (IMGT numbering)

| Region | IMGT Positions | Approx. bp |
|--------|---------------|------------|
| FR1 | 1–26 | 78 |
| CDR1 | 27–38 | 36 |
| FR2 | 39–55 | 51 |
| CDR2 | 56–65 | 30 |
| FR3 | 66–104 | 117 |
| **CDR3** | >104 | **excluded** |

*References*: Lefranc et al., 2003 (IMGT numbering); Ohlin et al., 2019 (IARC standard).

---

## Gene Selection & Pool Construction

### V and J gene pool (filtered)

1. Parse IMGT-standard FASTA reference files (header format `>...|gene|...|F|...`)
2. Retain only **functional** (`|F|`) and **non-partial** genes
3. Include **all `*01` alleles**
4. Randomly sample **20% of `*02` alleles** to reflect natural population variation
5. Uniform random selection from final pool

### D gene pool

All D genes from the reference FASTA are included without filtering.

---

## Usage

### Basic usage

```bash
# TRB simulation (no SHM)
python3 TRBSeqMock.py -n <N> -t TRB -o <output_dir>

# IGH simulation (with default SHM at 5%)
python3 TRBSeqMock.py -n <N> -t IGH -o <output_dir>

# IGH simulation with custom SHM parameters
python3 TRBSeqMock.py -n 50000 -t IGH -o output/ \
    --shm-rate 0.05 \
    --shm-hotspot-rr 4.0 \
    --shm-cdr-bias 1.5 \
    --shm-cold-spot-rate 0.002 \
    --shm-seed 42
```

### SHM CLI parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--shm-rate` | float | 0.05 | Global per-base SHM rate (0.0 = naive, 0.02 = low, 0.05 = medium, 0.10 = high) |
| `--shm-hotspot-rr` | float | 4.0 | Relative risk multiplier for WRC/GYW hotspot positions |
| `--shm-cdr-bias` | float | 1.5 | CDR mutation rate multiplier vs FR |
| `--shm-cold-spot-rate` | float | 0.002 | Background mutation rate outside hotspot positions |
| `--shm-seed` | int | None | Random seed for SHM reproducibility |

---

## Output Files

Each run produces the following files in the output directory:

| File | Format | Description |
|------|--------|-------------|
| `result.fasta` | FASTA | Simulated sequences. Header: `>Vgene\|Dgene\|Jgene\|length\|shm_count\|v_identity\|shm_rate` |
| `detail_result.tsv` | TSV | Per-sequence annotation with 24 columns including V/D/J genes, deletion/insertion lengths & sequences, and full SHM report |
| `shm_statistics.tsv` | TSV | Per-V-gene SHM summary: avg mutations, CDR/FR ratio, hotspot enrichment (IGH only) |
| `ground_truth.json` | JSONL | One JSON record per sequence with complete germline truth: V/D/J names, full germline sequences, SHM parameters, rearrangement details |
| `insertion_distribution.tsv` | TSV | Output insertion length distribution (VD and DJ) |
| `deletion_distribution.tsv` | TSV | Output deletion length distribution (V3, D5, D3, J5) |

### detail_result.tsv columns

```
V_gene  D_gene  J_gene  V3_del  D5_del  D3_del  J5_del
VD_ins_len  DJ_ins_len  VD_ins_seq  DJ_ins_seq
SHM_count  SHM_v_identity  SHM_cdr_count  SHM_fr_count
SHM_hotspot_count  SHM_transitions  SHM_transversions
SHM_rate_used  SHM_positions  SHM_details
Germline_V_seq  Germline_J_seq
```

### ground_truth.json example

```json
{
  "seq_id": 0,
  "germline_v_gene": "IGHV4-61*01",
  "germline_d_gene": "IGHD3-22*01",
  "germline_j_gene": "IGHJ4*01",
  "germline_v_seq": "CAGGTGCAG...",
  "germline_j_seq": "ACTACTTTGAC...",
  "shm_applied": true,
  "shm_rate_used": 0.05,
  "shm_count": 3,
  "shm_v_identity": 98.99,
  "rearrangement": {
    "v3_del": 1, "d5_del": 1, "d3_del": 1, "j5_del": 0,
    "vd_insert_len": 0, "dj_insert_len": 13,
    "vd_insert_seq": "", "dj_insert_seq": "ACTAGAATCCTAC"
  }
}
```

---

## Dependencies

| Package | Required? | Purpose |
|---------|-----------|---------|
| Python 3.8+ | ✓ | Runtime |
| numpy | ✓ | Numerical computation for sequencing error simulation probability arrays |
| numba | optional | JIT acceleration for `Source_error`; gracefully falls back to pure Python if not installed |
| biopython | optional | Used only in `scripts/sequencing_error_final.py` (standalone error tool); `data_loader/` uses pure Python FASTA parsing and does NOT require biopython |

---

## References

### SHM biology & mechanisms

1. **Di Noia JM, Neuberger MS.** Molecular mechanisms of antibody somatic hypermutation. *Annu Rev Biochem*. 2007;76:1–22. doi:10.1146/annurev.biochem.76.061705.090740
2. **Rogozin IB, Diaz M.** Cutting edge: DGYW/WRCH is a better predictor of mutability at G:C bases in Ig hypermutation than the widely accepted RGYW/WRCY motif. *J Immunol*. 2004;172(6):3382–4. doi:10.4049/jimmunol.172.6.3382
3. **Kleinstein SH, Louzoun Y, Shlomchik MJ.** Estimating hypermutation rates from clonal tree data. *J Immunol*. 2003;171(9):4639–49. doi:10.4049/jimmunol.171.9.4639
4. **Zheng NY, Wilson K, Jared M, Wilson PC.** Intricate targeting of immunoglobulin somatic hypermutation maximizes the efficiency of affinity maturation. *J Exp Med*. 2005;201(9):1467–78. doi:10.1084/jem.20042483

### SHM substitution profiles & targeting

5. **Sheng Z, Schramm CA, Kong R, NISC Comparative Sequencing Program, Mullikin JC, Mascola JR, et al.** Gene-specific substitution profiles describe the types and frequencies of amino acid changes during antibody somatic hypermutation. *Front Immunol*. 2017;8:537. doi:10.3389/fimmu.2017.00537
6. **Kirik U, Persson H, Levander F, Greiff L, Ohlin M.** Antibody heavy chain variable domains of different germline gene origins diversify through different paths. *Front Immunol*. 2017;8:1433. doi:10.3389/fimmu.2017.01433
7. **Schramm CA, Douek DC.** Beyond hot spots: biases in antibody somatic hypermutation and implications for vaccine design. *Front Immunol*. 2018;9:1876. doi:10.3389/fimmu.2018.01876

### Germline gene inference & evaluation (project motivation)

8. **Ohlin M, Scheepers C, Corcoran M, Lees WD, Busse CE, Bagnara D, et al.** Inferred allelic variants of immunoglobulin receptor genes: a system for their evaluation, documentation, and naming. *Front Immunol*. 2019;10:435. doi:10.3389/fimmu.2019.00435
9. **Gadala-Maria D, Yaari G, Uduman M, Kleinstein SH.** Automated analysis of high-throughput B-cell sequencing data reveals a high frequency of novel immunoglobulin V gene segment alleles. *Proc Natl Acad Sci USA*. 2015;112(8):E862–70. doi:10.1073/pnas.1417683112
10. **Zhang W, Wang IM, Wang C, Lin L, Chai X, Wu J, et al.** IMPre: an accurate and efficient software for prediction of T- and B-cell receptor germline genes and alleles from rearranged repertoire data. *Front Immunol*. 2016;7:457. doi:10.3389/fimmu.2016.00457

### AIRR-seq methodology & tools

11. **Yaari G, Kleinstein SH.** Practical guidelines for B-cell receptor repertoire sequencing analysis. *Genome Med*. 2015;7:121. doi:10.1186/s13073-015-0243-2
12. **Gupta NT, Vander Heiden JA, Uduman M, Gadala-Maria D, Yaari G, Kleinstein SH.** Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data. *Bioinformatics*. 2015;31(20):3356–8. doi:10.1093/bioinformatics/btv359
13. **Corcoran MM, Phad GE, Vázquez Bernat N, Stahl-Hennig C, Sumida N, Persson MA, et al.** Production of individualized V gene databases reveals high levels of immunoglobulin genetic diversity. *Nat Commun*. 2016;7:13642. doi:10.1038/ncomms13642

### IMGT nomenclature & numbering

14. **Lefranc M-P, Giudicelli V, Duroux P, Jabado-Michaloud J, Folch G, Aouinti S, et al.** IMGT®, the international ImMunoGeneTics information system® 25 years on. *Nucleic Acids Res*. 2015;43(D1):D413–22. doi:10.1093/nar/gku1056
15. **Brochet X, Lefranc M-P, Giudicelli V.** IMGT/V-QUEST: the highly customized and integrated system for IG and TR standardized V-J and V-D-J sequence analysis. *Nucleic Acids Res*. 2008;36(Web Server):W503–8. doi:10.1093/nar/gkn316

---

*For questions, contact: 2210240103@csu.edu.cn*
*GitHub: https://github.com/NitroMint/TRBSeqMock*
