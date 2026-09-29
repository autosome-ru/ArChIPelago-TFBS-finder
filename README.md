# ArChIPelago-TFBS-finder

Command-line tool that scans DNA sequences for transcription factor binding sites with the pre-trained
ArChIPelago Random Forest models: each model combines the scores of all human mono- and dinucleotide PWMs of a
transcription factor (TF). Hits are called against an empirical null distribution built from
dinucleotide-shuffled sequences, with Benjamini-Hochberg FDR control.

> Kravchenko P., Vorontsov I.E., Grosse I., Makeev V.J., Kulakovskiy I.V., and Penzar D.D. (2026).
> *Classic machine learning on top of multiple position weight matrices improves genomic prediction of
> transcription factor binding sites.* Pipeline and manuscript analyses:
> [autosome-ru/ArChIPelago](https://github.com/autosome-ru/ArChIPelago).

---

## Contents

1. [Models](#models)
2. [Supported transcription factors](#supported-transcription-factors)
3. [Installation](#installation)
4. [Models and PWMs](#models-and-pwms)
5. [Directory layout](#directory-layout)
6. [Usage](#usage)
7. [Output files](#output-files)
8. [Examples](#examples)
9. [Testing](#testing)
10. [How it works](#how-it-works)
11. [Citation](#citation)
12. [License](#license)

---

## Models

Three models per TF, one for each PWM set: monoPWMs (`mono`), diPWMs (`di`) and monoPWMs + diPWMs
(`mono_di`, the default).

| | |
|---|---|
| Algorithm | `RandomForestClassifier(max_depth=6, max_samples=0.8, n_estimators=100, random_state=0)`, scikit-learn 1.3 |
| Training data | human ChIP-Seq peak regions of the TF (positives) and GC-matched peak regions of unrelated TF families (negatives), human training chromosomes 2-7, 9, 10, 13-20 |
| Features | best-hit log-odds score of every PWM of the TF in the 300-bp sequence (SPRY-SARUS), standardised with the mean and standard deviation of that feature on the training set |
| Files | `Models/<TF>/ArChIPelago_<TF>_<mono\|di\|mono_di>.sav` (joblib) and `.json` (feature specification) |

The `.json` file lists the features in the order the model expects (`mono_<k>` = `PWMs_mono_HUMAN/<TF>/<k>.pwm`,
`di_<k>` = `PWMs_di_HUMAN/<TF>/<k>.dpwm`), the training mean and standard deviation of each feature
(`scaler_mean`, `scaler_scale`), the size of the training set and the auROC / auPRC of the model on the human test
set (chromosomes 1, 8, 21). Over the 36 TFs the `mono_di` models reach a median auROC of 0.8915 and a median auPRC of
0.3026 on the human test set (`mono`: 0.8796 / 0.2734; `di`: 0.8847 / 0.2911).

The tool stops with an error if a model, its `.json`, a PWM file of one of its features or a SARUS scan is missing.

---

## Supported transcription factors

```
ANDR   AP2A   CEBPB  COE1   CTCF   E2F4   ERG    ESR1
FLI1   GATA1  GATA2  GATA3  GCR    HNF4A  IRF1   IRF4
JUND   MAFK   MAX    MYC    P53    PPARG  PRGR   REST
RUNX1  RXRA   SOX2   SPI1   SRF    STA5A  STAT1  STAT3
TAL1   TF65   TFE2   USF2
```

`python scanning_tool.py --list_tfs` prints the list.

---

## Installation

### 1. Clone the repository

```bash
git clone --recurse-submodules https://github.com/autosome-ru/ArChIPelago-TFBS-finder.git
cd ArChIPelago-TFBS-finder
```

`--recurse-submodules` also fetches SPRY-SARUS (`sarus/`, pinned at release 2.2.3). In a clone made without it,
run `git submodule update --init --recursive`.

### 2. Python environment

Python 3.8 or later; the models need scikit-learn 1.3 or later.

```bash
conda create -n archipelago python=3.8
conda activate archipelago
pip install -r requirements.txt
```

### 3. Java (for SPRY-SARUS)

The PWM scores are computed with [SPRY-SARUS](https://github.com/autosome-ru/sarus), which needs Java 8 or later:

```bash
java -version
# macOS:  brew install openjdk
# Ubuntu: sudo apt install default-jre
```

### 4. SARUS jar

The `sarus` submodule provides `sarus/releases/sarus-2.2.3.jar`. The tool uses the first jar it finds among `../sarus/releases/sarus-2.2.3.jar` (the `sarus`
submodule when the tool is used inside the ArChIPelago repository), `sarus/releases/sarus-2.2.3.jar`,
`sarus/sarus-2.2.3.jar`, `../sarus/releases/sarus-2.0.1.jar` and `sarus/releases/sarus-2.0.1.jar`;
`--sarus_jar /path/to/sarus.jar` sets it explicitly.

---

## Models and PWMs

The models (`Models/`, 36 TFs x 3 models) and the human PWMs (`PWMs_mono_HUMAN/`, 1,495 monoPWMs;
`PWMs_di_HUMAN/`, 780 diPWMs) are part of this repository; the same files are in the ArChIPelago Zenodo record
([10.5281/zenodo.14927303](https://doi.org/10.5281/zenodo.14927303): `Models.tar.gz`, `PWMs_mono_HUMAN.tar.gz`,
`PWMs_di_HUMAN.tar.gz`).

Check the models of a TF:

```bash
python scanning_tool.py --list_models --tf CTCF
```

---

## Directory layout

```
ArChIPelago-TFBS-finder/
├── scanning_tool.py           # command-line tool
├── requirements.txt
├── pytest.ini
├── README.md
│
├── Models/                    # Random Forest models (.sav) and feature specifications (.json)
│   ├── CTCF/
│   │   ├── ArChIPelago_CTCF_mono.sav
│   │   ├── ArChIPelago_CTCF_mono.json
│   │   ├── ArChIPelago_CTCF_di.sav
│   │   ├── ArChIPelago_CTCF_di.json
│   │   ├── ArChIPelago_CTCF_mono_di.sav
│   │   └── ArChIPelago_CTCF_mono_di.json
│   └── .../
│
├── PWMs_mono_HUMAN/           # monoPWMs
│   ├── CTCF/
│   │   └── <k>.pwm
│   └── .../
│
├── PWMs_di_HUMAN/             # diPWMs
│   ├── CTCF/
│   │   └── <k>.dpwm
│   └── .../
│
├── sarus/
│   └── releases/
│       └── sarus-2.2.3.jar    # SPRY-SARUS (git submodule)
│
├── synthetic_CTCF_demo.fasta  # demo input
├── realdata_CTCF_test.fasta   # CTCF test input
│
└── tests/
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
| `-f`, `--fasta` | input FASTA file |
| `--tf` | TF name (e.g. `CTCF`), or `all` for all 36 TFs |

### Optional arguments

| Argument | Default | Description |
|---|---|---|
| `-o`, `--output` | `Results/` | output directory |
| `--pwm_type` | `mono_di` | model: `mono`, `di` or `mono_di` |
| `--frame` | `300` | window size in bp (the length of the training sequences) |
| `--step` | `150` | window step in bp |
| `--fdr` | `0.1` | FDR threshold for significant windows |
| `--prob` | `0.7` | probability threshold, used with `--no-null` |
| `--n_null` | `50` | number of dinucleotide-shuffled sequences for the null distribution |
| `--no-null` | off | skip the null calibration (no p- and q-values) |
| `-v`, `--verbose` | off | verbose logging |
| `--models_dir` | `Models/` | models directory |
| `--pwm_mono_dir` | `PWMs_mono_HUMAN/` | monoPWM directory |
| `--pwm_di_dir` | `PWMs_di_HUMAN/` | diPWM directory |
| `--sarus_jar` | see [SARUS jar](#4-sarus-jar) | SARUS jar file |

Default directories are relative to `scanning_tool.py`.

### Information only

| Argument | Description |
|---|---|
| `--list_tfs` | print the supported TFs and exit |
| `--list_models` | print the models available for `--tf` and exit |

---

## Output files

Written to the directory given by `-o` (default `Results/`):

| File | Content |
|---|---|
| `<TF>_predictions_full.tsv` | every window with its probability (and p- and q-value with the null calibration) |
| `<TF>_significant_FDR<fdr>.tsv` | windows with q-value <= `--fdr`, sorted by q-value |
| `<TF>_significant_prob<prob>.tsv` | with `--no-null`: windows with probability >= `--prob`, sorted by probability |
| `<TF>_significant.bed` | the significant windows as BED (sequence id, start, start + frame, `<TF>_hit`, 1000 x probability, `.`) |

### TSV columns

| Column | Content |
|---|---|
| `window_id` | `<sequence_id>@<window index>` (index 0, 1, 2, ... along the sequence) |
| `position` | 0-based start of the window in the sequence |
| `predicted_probability` | Random Forest probability of a binding site in the window |
| `empirical_pvalue` | p-value against the null distribution: (number of null probabilities >= observed + 1) / (null size + 1) |
| `qvalue` | Benjamini-Hochberg q-value |

---

## Examples

```bash
# CTCF on the demo sequence
python scanning_tool.py -f synthetic_CTCF_demo.fasta --tf CTCF

# stricter FDR
python scanning_tool.py -f realdata_CTCF_test.fasta --tf CTCF --fdr 0.05

# without null calibration
python scanning_tool.py -f synthetic_CTCF_demo.fasta --tf CTCF --no-null

# monoPWM model, verbose
python scanning_tool.py -f synthetic_CTCF_demo.fasta --tf CTCF --pwm_type mono -v

# all 36 TFs
python scanning_tool.py -f sequences.fasta --tf all -o Results/batch_scan/

# other window size and step
python scanning_tool.py -f sequences.fasta --tf CTCF --frame 500 --step 250

# supported TFs and the models of one TF
python scanning_tool.py --list_tfs
python scanning_tool.py --list_models --tf CTCF
```

---

## Testing

The test suite has 115 tests: 81 unit tests (`tests/test_scanning_tool.py`, no data needed) and 34 integration
tests (`tests/test_integration.py`) that use `Models/`, `PWMs_mono_HUMAN/`, `PWMs_di_HUMAN/`, the SARUS jar and
the FASTA files of the repository; they need scikit-learn 1.3 or later to load the models.

```bash
pytest                                # all tests
pytest tests/test_scanning_tool.py    # unit tests
pytest tests/test_integration.py      # integration tests
```

---

## How it works

1. **Input.** The FASTA sequences are read and upper-cased.
2. **Windows.** Each sequence is cut into windows of `--frame` bp with a step of `--step` bp; a sequence shorter
   than the frame is one window. Windows with IUPAC ambiguity codes (including N) are skipped.
3. **PWM scores.** Every window is scanned with each PWM of the model by SPRY-SARUS
   (`--skipn --show-non-matching --output-scoring-mode score besthit`), giving the best-hit log-odds score.
4. **Features.** The scores form one column per PWM, in the order listed in the model's `.json`, and are
   standardised with the training mean and standard deviation of each column.
5. **Prediction.** The Random Forest gives the probability of a binding site for every window.
6. **Null distribution** (unless `--no-null`). All input sequences are concatenated and shuffled `--n_null` times
   with the dinucleotide counts preserved exactly (Altschul-Erickson; seed 42); the shuffled sequences are scanned
   together in one run in the same way.
7. **p-values.** Each window is compared with the pooled null probabilities.
8. **FDR.** Benjamini-Hochberg q-values; windows with q-value <= `--fdr` are significant.
9. **Export.** Full table, significant windows and BED file.

---

## Citation

Kravchenko P., Vorontsov I.E., Grosse I., Makeev V.J., Kulakovskiy I.V., and Penzar D.D. (2026). Classic machine
learning on top of multiple position weight matrices improves genomic prediction of transcription factor binding
sites.

Models, PWMs and data: Zenodo [10.5281/zenodo.14927303](https://doi.org/10.5281/zenodo.14927303).

---

## License

MIT License (see the header of `scanning_tool.py`).
