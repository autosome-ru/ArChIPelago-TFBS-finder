# ArChIPelago-TFBS-finder

**Command-line tool for scanning DNA sequences for transcription factor binding sites (TFBS) using pre-trained ArChIPelago Random Forest models.**

ArChIPelago-TFBS-finder addresses the false-positive problem inherent in
repurposing a classifier as a scanner by constructing an **empirical null
distribution** from dinucleotide-shuffled sequences and applying
**Benjamini–Hochberg FDR control**.

> **Paper:** *ArChIPelago — an Automated Pipeline for Comprehensive ChIP-seq
> Data Analysis* (Kravchenko et al.)

---

## Table of Contents

1. [Features](#features)
2. [Supported Transcription Factors](#supported-transcription-factors)
3. [Installation](#installation)
4. [Downloading Data from Zenodo](#downloading-data-from-zenodo)
5. [Directory Layout](#directory-layout)
6. [Usage](#usage)
7. [Output Files](#output-files)
8. [Examples](#examples)
9. [Testing](#testing)
10. [How It Works](#how-it-works)
11. [Citation](#citation)
12. [License](#license)

---

## Features

| Feature | Description |
|---|---|
| **36 pre-trained RF models** | One model per TF, trained on ChIP-seq data from HOCOMOCO |
| **Mono- & di-nucleotide PWMs** | Scans with SARUS using both PWM types |
| **Sliding-window scanning** | Configurable frame size (default 301 bp) and step (default 150 bp) |
| **Empirical null calibration** | Dinucleotide-shuffled null distribution to estimate p-values |
| **BH-FDR control** | Benjamini–Hochberg correction; default FDR ≤ 0.10 |
| **Multiple output formats** | TSV (full + significant) and BED |
| **Single-TF or batch mode** | Scan for one TF or all 36 in one run |

---

## Supported Transcription Factors

```
ANDR   AP2A   CEBPB  COE1   CTCF   E2F4   ERG    ESR1
FLI1   GATA1  GATA2  GATA3  GCR    HNF4A  IRF1   IRF4
JUND   MAFK   MAX    MYC    P53    PPARG  PRGR   REST
RUNX1  RXRA   SOX2   SPI1   SRF    STA5A  STAT1  STAT3
TAL1   TF65   TFE2   USF2
```

Run `python scanning_tool.py --list_tfs` to see the full list.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/ArChIPelago-TFBS-finder.git
cd ArChIPelago-TFBS-finder
```

### 2. Create a Python environment

Python **≥ 3.8** is required. We recommend **conda**:

```bash
conda create -n archipelago python=3.8
conda activate archipelago
pip install -r requirements.txt
```

<details>
<summary>requirements.txt contents</summary>

```
numpy>=1.19.0
pandas>=1.2.0
scikit-learn>=1.3.0
joblib>=1.0.0
biopython>=1.79
pytest>=7.0.0
```

</details>

### 3. Install Java (for SARUS)

The tool calls [SARUS](https://github.com/autosome-ru/sarus) to score sequences
against PWMs. You need **Java ≥ 8**:

```bash
java -version          # check if already installed
# macOS:  brew install openjdk
# Ubuntu: sudo apt install default-jre
```

### 4. Download the SARUS jar

```bash
mkdir -p sarus/releases
# Download sarus-2.0.1.jar from:
# https://github.com/autosome-ru/sarus/releases
# and place it in sarus/releases/sarus-2.0.1.jar
```

The tool auto-detects the jar at `sarus/releases/sarus-2.0.1.jar` relative to
the script location. You can also pass `--sarus_jar /path/to/sarus.jar`.

---

## Downloading Data from Zenodo

**Pre-trained models and PWM matrices are NOT included in this repository.**
Download them from Zenodo:

> **DOI:** [10.5281/zenodo.14927304](https://doi.org/10.5281/zenodo.14927304)

After downloading and extracting the archive, place three directories in the
repository root:

```bash
# From the Zenodo archive, copy these three directories:
ArChIPelago-TFBS-finder/
├── Models_sklearn13/      # Pre-trained Random Forest models (.sav files)
├── PWMs_mono_HUMAN/       # Mono-nucleotide PWMs (.pwm files)
└── PWMs_di_HUMAN/         # Di-nucleotide PWMs (.dpwm files)
```

Each directory contains one sub-folder per TF (e.g. `Models_sklearn13/CTCF/`).

Verify the data is in place:

```bash
python scanning_tool.py --list_models --tf CTCF
```

---

## Directory Layout

After setup, your working directory should look like this:

```
ArChIPelago-TFBS-finder/
├── scanning_tool.py           # Main CLI tool
├── requirements.txt           # Python dependencies
├── pytest.ini                 # Test configuration
├── README.md                  # This file
├── .gitignore
│
├── Models_sklearn13/          # ⬇ from Zenodo
│   ├── CTCF/
│   │   ├── model_CTCF_m1.sav
│   │   └── ...
│   └── .../
│
├── PWMs_mono_HUMAN/           # ⬇ from Zenodo
│   ├── CTCF/
│   │   └── *.pwm
│   └── .../
│
├── PWMs_di_HUMAN/             # ⬇ from Zenodo
│   ├── CTCF/
│   │   └── *.dpwm
│   └── .../
│
├── sarus/
│   └── releases/
│       └── sarus-2.0.1.jar   # SARUS scanner (download separately)
│
├── synthetic_CTCF_demo.fasta  # Demo input (shipped with repo)
├── realdata_CTCF_test.fasta   # Real CTCF test set (shipped with repo)
│
└── tests/                     # Unit & integration tests
    ├── test_scanning_tool.py
    └── test_integration.py
```

---

## Usage

```
python scanning_tool.py -f <FASTA> --tf <TF_NAME> [OPTIONS]
```

### Required arguments

| Argument | Description |
|---|---|
| `-f`, `--fasta` | Path to input FASTA file |
| `--tf` | Transcription factor name (e.g. `CTCF`). Use `all` to scan all 36 TFs |

### Optional arguments

| Argument | Default | Description |
|---|---|---|
| `-o`, `--output` | `Results/` | Output directory |
| `--pwm_type` | `mono_di` | PWM type: `mono`, `di`, or `mono_di` |
| `--frame` | `301` | Sliding window size (bp) |
| `--step` | `150` | Sliding window step (bp) |
| `--fdr` | `0.1` | FDR threshold for significant hits |
| `--prob` | `0.7` | Probability threshold (used when null is disabled) |
| `--n_null` | `50` | Number of dinucleotide-shuffled sequences for null |
| `--no-null` | off | Skip null calibration (faster, no FDR) |
| `-v`, `--verbose` | off | Verbose logging |
| `--models_dir` | auto | Path to models directory |
| `--pwm_mono_dir` | auto | Path to mono-nucleotide PWM directory |
| `--pwm_di_dir` | auto | Path to di-nucleotide PWM directory |
| `--sarus_jar` | auto | Path to SARUS jar file |

### Info-only arguments (no scanning)

| Argument | Description |
|---|---|
| `--list_tfs` | Print all supported TF names and exit |
| `--list_models` | Print available model files for the specified TF and exit |

---

## Output Files

All outputs are written to the directory specified by `-o` (default `Results/`).

| File | Description |
|---|---|
| `<TF>_predictions_full.tsv` | All scanned windows with probabilities, p-values, and q-values |
| `<TF>_significant_FDR0.1.tsv` | Windows passing FDR threshold |
| `<TF>_significant.bed` | BED file of significant hits (for genome browsers) |

### TSV columns

| Column | Description |
|---|---|
| `window_id` | `<sequence_id>@window_<N>` |
| `position` | Start position of the window (0-based) |
| `predicted_probability` | RF model probability of being a true binding site |
| `empirical_pvalue` | p-value from empirical null distribution |
| `qvalue` | Benjamini–Hochberg adjusted p-value |

---

## Examples

### Quick scan with demo data

```bash
# Scan synthetic sequences for CTCF
python scanning_tool.py -f synthetic_CTCF_demo.fasta --tf CTCF

# Scan real data with stricter FDR
python scanning_tool.py -f realdata_CTCF_test.fasta --tf CTCF --fdr 0.05

# Fast scan without null calibration
python scanning_tool.py -f synthetic_CTCF_demo.fasta --tf CTCF --no-null

# Scan with only mono-nucleotide PWMs, verbose output
python scanning_tool.py -f synthetic_CTCF_demo.fasta --tf CTCF --pwm_type mono -v

# Batch: scan all 36 TFs
python scanning_tool.py -f sequences.fasta --tf all -o Results/batch_scan/

# Custom window size and step
python scanning_tool.py -f sequences.fasta --tf CTCF --frame 501 --step 250
```

### Inspecting available resources

```bash
# List supported TFs
python scanning_tool.py --list_tfs

# Show model files for a specific TF
python scanning_tool.py --list_models --tf CTCF
```

---

## Testing

The test suite includes 122 tests (84 unit + 38 integration).

```bash
# Run all tests
pytest

# Run only unit tests
pytest tests/test_scanning_tool.py

# Run only integration tests (requires models + PWMs)
pytest tests/test_integration.py

# Verbose output
pytest -v
```

---

## How It Works

1. **Parse input FASTA** — reads sequences, filters ambiguous nucleotides.
2. **Sliding-window decomposition** — each sequence is split into overlapping
   windows of `--frame` bp with `--step` bp stride.
3. **PWM scoring with SARUS** — each window is scored against mono- and/or
   di-nucleotide PWMs using the SARUS Java tool.
4. **Feature matrix construction** — PWM scores form the feature vector for each
   window.
5. **Random Forest prediction** — pre-trained RF model outputs the probability
   that the window contains a true binding site.
6. **Empirical null calibration** *(unless `--no-null`)* — dinucleotide-shuffled
   versions of each input sequence are scanned identically; the resulting
   probability distribution forms the null.
7. **Empirical p-value calculation** — for each real window, the fraction of
   null probabilities ≥ the observed probability gives the p-value.
8. **Benjamini–Hochberg FDR correction** — q-values are computed; windows with
   q-value ≤ `--fdr` are reported as significant.
9. **Export** — full results (TSV), significant hits (TSV), and BED file.

---

## Citation

If you use ArChIPelago-TFBS-finder in your research, please cite:

```bibtex
@article{kravchenko2025archipelago,
  title   = {ArChIPelago — an Automated Pipeline for Comprehensive ChIP-seq Data Analysis},
  author  = {Kravchenko, Pavel and others},
  year    = {2025},
}
```

**Data:**
> Kravchenko, P. (2025). ArChIPelago pre-trained models and PWMs [Data set].
> Zenodo. https://doi.org/10.5281/zenodo.14927304

---

## License

MIT License. See the source header of `scanning_tool.py` for details.
