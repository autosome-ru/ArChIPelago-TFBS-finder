# ArChIPelago-TFBS-finder

[![tests](https://github.com/autosome-ru/ArChIPelago-TFBS-finder/actions/workflows/tests.yml/badge.svg)](https://github.com/autosome-ru/ArChIPelago-TFBS-finder/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14927303.svg)](https://doi.org/10.5281/zenodo.14927303)
[![python](https://img.shields.io/badge/python-3.8--3.12-blue.svg)](requirements.txt)
[![scikit-learn](https://img.shields.io/badge/scikit--learn-1.3-orange.svg)](requirements.txt)
[![License: MIT](https://img.shields.io/badge/license-MIT-lightgrey.svg)](scanning_tool.py)

Command-line tool that scans DNA sequences for transcription factor binding sites with the pre-trained
ArChIPelago Random Forest models: each model combines the scores of all mono- and dinucleotide PWMs of a
transcription factor (TF). Models trained on human data (human PWMs) and on mouse data (mouse PWMs) are
available for 35 TFs; `--species` selects them. Hits are called against an empirical null distribution built from
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

Six models per TF: for each species (human, mouse) three models, one for each PWM set: monoPWMs (`mono`), diPWMs
(`di`) and monoPWMs + diPWMs (`mono_di`, the default). Human models use the human PWMs and are trained on human
peaks; mouse models use the mouse PWMs and are trained on mouse peaks. The tool uses the human models unless
`--species mouse` is given.

| | |
|---|---|
| Algorithm | `RandomForestClassifier(max_depth=6, max_samples=0.8, n_estimators=100, random_state=0)`, scikit-learn 1.3 |
| Positives | per peak set the top 1,000 GTRD peaks (tags >= 10, ranked by -10log10 P), +-500 bp windows, repeat filter, at most 10,000 positives per TF and species, 300-bp sequences |
| Negatives (v3) | peaks of TFs of other TFClass classes; windows overlapping the TF's own positives are removed; the rest are GC-matched to the positives with BiasAway |
| Training chromosomes | human chr2-7, 9-20, 22; mouse chr2-7, 9-18 |
| Test chromosomes | human chr1, 8, 21; mouse chr1, 8, 19 |
| Features | best-hit log-odds score of every PWM of the TF in the 300-bp sequence (SPRY-SARUS), standardised with the mean and standard deviation of that feature on the training set |
| Files | `Models/<TF>/ArChIPelago_<TF>_<mono\|di\|mono_di>.sav` (human) and `ArChIPelago_<TF>_MOUSE_<mono\|di\|mono_di>.sav` (mouse), joblib, each with a `.json` (feature specification) |

The `.json` file lists the features in the order the model expects (`mono_<k>` = `PWMs_mono_<SPECIES>/<TF>/<k>.pwm`,
`di_<k>` = `PWMs_di_<SPECIES>/<TF>/<k>.dpwm`, with `<SPECIES>` = `HUMAN` or `MOUSE`), the training mean and standard
deviation of each feature (`scaler_mean`, `scaler_scale`), the sizes of the training and test sets and the auROC and
auPRC of the model on the test chromosomes of its own species and of the other species. In the mouse `.json` files the
`pwm_files` entry names the PWM folders `PWMs_MOUSE_mono` and `PWMs_MOUSE_di`; in this repository they are
`PWMs_mono_MOUSE/` and `PWMs_di_MOUSE/`.

Median test auROC / auPRC over the 35 TFs on the test chromosomes of the training species (from the `.json` files):

| model | human `mono` | human `di` | human `mono_di` | mouse `mono` | mouse `di` | mouse `mono_di` |
|---|---|---|---|---|---|---|
| auROC | 0.8906 | 0.9021 | 0.9021 | 0.8948 | 0.8968 | 0.8973 |
| auPRC | 0.3142 | 0.3251 | 0.3360 | 0.3274 | 0.3441 | 0.3563 |

The tool stops with an error if a model, its `.json`, a PWM file of one of its features or a SARUS scan is missing.

---

## Supported transcription factors

```
ANDR   AP2A   CEBPB  COE1   CTCF   E2F4   ERG    ESR1
FLI1   GATA1  GATA2  GATA3  GCR    HNF4A  IRF1   IRF4
JUND   MAFK   MAX    MYC    PPARG  PRGR   REST   RUNX1
RXRA   SOX2   SPI1   SRF    STA5A  STAT1  STAT3  TAL1
TF65   TFE2   USF2
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

Python 3.8-3.12 and scikit-learn 1.3 (the models were saved with 1.3.0; scikit-learn 1.4 and later cannot use them, and the tool stops with an error).

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

The models (`Models/`, 35 TFs x 2 species x 3 models), the human PWMs (`PWMs_mono_HUMAN/`, 1,447 monoPWMs;
`PWMs_di_HUMAN/`, 754 diPWMs) and the mouse PWMs (`PWMs_mono_MOUSE/`, 1,178 monoPWMs; `PWMs_di_MOUSE/`, 617 diPWMs)
are part of this repository. The same files, with the training and test sets and the other data of the revision
(R2), are in the ArChIPelago Zenodo record
([10.5281/zenodo.14927303](https://doi.org/10.5281/zenodo.14927303), concept DOI of all versions): `Models.tar.gz`,
`PWMs_mono_HUMAN.tar.gz`, `PWMs_di_HUMAN.tar.gz`, `PWMs_mono_MOUSE.tar.gz`, `PWMs_di_MOUSE.tar.gz`.

Check the models of a TF:

```bash
python scanning_tool.py --list_models --tf CTCF
python scanning_tool.py --list_models --tf CTCF --species mouse
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
│   │   ├── ArChIPelago_CTCF_mono_di.json
│   │   ├── ArChIPelago_CTCF_MOUSE_mono.sav       # mouse-trained models
│   │   ├── ArChIPelago_CTCF_MOUSE_mono.json
│   │   └── ...                                    # MOUSE_di, MOUSE_mono_di
│   └── .../
│
├── PWMs_mono_HUMAN/           # human monoPWMs
│   ├── CTCF/
│   │   └── <k>.pwm
│   └── .../
│
├── PWMs_di_HUMAN/             # human diPWMs
│   ├── CTCF/
│   │   └── <k>.dpwm
│   └── .../
│
├── PWMs_mono_MOUSE/           # mouse monoPWMs (same layout)
├── PWMs_di_MOUSE/             # mouse diPWMs (same layout)
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
| `--tf` | TF name (e.g. `CTCF`), or `all` for all 35 TFs |

### Optional arguments

| Argument | Default | Description |
|---|---|---|
| `-o`, `--output` | `Results/` | output directory |
| `--species` | `human` | species of the training set and the PWMs: `human` or `mouse` |
| `--pwm_type` | `mono_di` | model: `mono`, `di` or `mono_di` |
| `--frame` | `300` | window size in bp (the length of the training sequences) |
| `--step` | `150` | window step in bp |
| `--fdr` | `0.1` | FDR threshold for significant windows |
| `--prob` | `0.7` | probability threshold, used with `--no-null` |
| `--n_null` | `50` | number of dinucleotide-shuffled sequences for the null distribution |
| `--no-null` | off | skip the null calibration (no p- and q-values) |
| `-v`, `--verbose` | off | verbose logging |
| `--models_dir` | `Models/` | models directory |
| `--pwm_mono_dir` | `PWMs_mono_HUMAN/` (`PWMs_mono_MOUSE/` with `--species mouse`) | monoPWM directory |
| `--pwm_di_dir` | `PWMs_di_HUMAN/` (`PWMs_di_MOUSE/` with `--species mouse`) | diPWM directory |
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

With `--species mouse` the file names read `<TF>_MOUSE_predictions_full.tsv`, `<TF>_MOUSE_significant_FDR<fdr>.tsv`,
`<TF>_MOUSE_significant_prob<prob>.tsv` and `<TF>_MOUSE_significant.bed`.

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

# mouse-trained model with the mouse PWMs
python scanning_tool.py -f synthetic_CTCF_demo.fasta --tf CTCF --species mouse

# all 35 TFs
python scanning_tool.py -f sequences.fasta --tf all -o Results/batch_scan/

# other window size and step
python scanning_tool.py -f sequences.fasta --tf CTCF --frame 500 --step 250

# supported TFs and the models of one TF
python scanning_tool.py --list_tfs
python scanning_tool.py --list_models --tf CTCF
```

---

## Testing

The test suite has 140 tests: 99 unit tests (`tests/test_scanning_tool.py`, no data needed) and 41 integration
tests (`tests/test_integration.py`) that use `Models/`, the human and mouse PWM directories, the SARUS jar and
the FASTA files of the repository; they need scikit-learn 1.3 to load the models.

```bash
pytest                                # all tests
pytest tests/test_scanning_tool.py    # unit tests
pytest tests/test_integration.py      # integration tests
```

GitHub Actions runs the tests on Python 3.9 and 3.12 and a CTCF scan with null calibration on every push
(`.github/workflows/tests.yml`).

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

Models, PWMs and data: Zenodo [10.5281/zenodo.14927303](https://doi.org/10.5281/zenodo.14927303) (all versions); release v1.0.0: [10.5281/zenodo.23036475](https://doi.org/10.5281/zenodo.23036475).

---

## License

MIT License (see the header of `scanning_tool.py`).
