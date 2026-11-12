# ArChIPelago-TFBS-finder

[![tests](https://github.com/autosome-ru/ArChIPelago-TFBS-finder/actions/workflows/tests.yml/badge.svg)](https://github.com/autosome-ru/ArChIPelago-TFBS-finder/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14927303.svg)](https://doi.org/10.5281/zenodo.14927303)
[![python](https://img.shields.io/badge/python-3.10%20%7C%203.12-blue.svg)](.github/workflows/tests.yml)
[![scikit-learn](https://img.shields.io/badge/scikit--learn-1.3-orange.svg)](requirements.txt)
[![License: MIT](https://img.shields.io/badge/license-MIT-lightgrey.svg)](scanning_tool.py)

`scanning_tool.py` scans DNA sequences for binding sites of 35 transcription factors (TFs). It cuts each sequence
into windows, scores every window with all mononucleotide and dinucleotide position weight matrices (monoPWMs and
diPWMs) of the TF using SPRY-SARUS, and passes the scores to a pre-trained Random Forest model, which gives the
probability of a binding site. The tool then compares these probabilities with those of windows from
dinucleotide-shuffled copies of the input and reports the windows that pass a Benjamini-Hochberg false discovery rate
(FDR) threshold. Models trained on human data and on mouse data are included; `--species` selects them.

## Contents

1. [Installation](#installation)
2. [Models](#models)
3. [PWMs](#pwms)
4. [Directory layout](#directory-layout)
5. [Usage](#usage)
6. [Output files](#output-files)
7. [Examples](#examples)
8. [Testing](#testing)
9. [How it works](#how-it-works)
10. [Citation](#citation)
11. [Licence](#licence)

## Installation

### 1. Clone with the SARUS submodule

```bash
git clone --recurse-submodules https://github.com/autosome-ru/ArChIPelago-TFBS-finder.git
cd ArChIPelago-TFBS-finder
```

The `sarus` submodule is SPRY-SARUS ([autosome-ru/sarus](https://github.com/autosome-ru/sarus)) at tag v2.2.3 and
provides `sarus/releases/sarus-2.2.3.jar`. In a clone made without `--recurse-submodules`, run
`git submodule update --init --recursive`.

### 2. Python packages

```bash
python3 -m venv archipelago-env
source archipelago-env/bin/activate
pip install -r requirements.txt
```

`requirements.txt` asks for numpy >= 1.19.0, pandas >= 1.2.0, scikit-learn >= 1.3.0 and < 1.4, joblib >= 1.0.0,
biopython >= 1.79 and pytest >= 7.0.0. The models were saved with scikit-learn 1.3. With another minor version of
scikit-learn the tool stops with an error when it loads a model, and with a version older than 1.3 it fails on import.
The test workflow runs on Python 3.10 and 3.12 (`.github/workflows/tests.yml`).

### 3. Java

SPRY-SARUS is a Java program, so `java` must be on the `PATH` (check with `java -version`). The test workflow uses
Java 17. Each SARUS call is started with `-Xmx2G`.

### 4. SARUS jar

Without `--sarus_jar`, the tool uses the first of these files that exists, with paths relative to the folder of
`scanning_tool.py`:

1. `../sarus/releases/sarus-2.2.3.jar`
2. `sarus/releases/sarus-2.2.3.jar`
3. `sarus/sarus-2.2.3.jar`
4. `../sarus/releases/sarus-2.0.1.jar`
5. `sarus/releases/sarus-2.0.1.jar`

In a clone made as in step 1, this is `sarus/releases/sarus-2.2.3.jar`.

## Models

`Models/` holds 210 models: 35 TFs, two training species (human, mouse) and three feature sets per species.

| feature set | `--pwm_type` |
|---|---|
| monoPWMs of the TF | `mono` |
| diPWMs of the TF | `di` |
| monoPWMs and diPWMs of the TF | `mono_di` (default) |

Human models are trained on human data and use the human PWMs; mouse models are trained on mouse data and use the
mouse PWMs. The tool uses the human models unless `--species mouse` is given. Each model is a scikit-learn
`RandomForestClassifier(max_depth=6, max_samples=0.8, n_estimators=100, random_state=0)` saved with joblib, with a
`.json` file of the same name next to it:

```
Models/<TF>/ArChIPelago_<TF>_<mono|di|mono_di>.sav          human
Models/<TF>/ArChIPelago_<TF>_MOUSE_<mono|di|mono_di>.sav    mouse
```

The `.json` file lists the features in the order the model expects (`mono_<k>` is
`PWMs_mono_<SPECIES>/<TF>/<k>.pwm`, `di_<k>` is `PWMs_di_<SPECIES>/<TF>/<k>.dpwm`, `<SPECIES>` is `HUMAN` or
`MOUSE`) and the training mean and standard deviation of each feature (`scaler_mean`, `scaler_scale`). It also gives
the training and test chromosomes (human: training chr2-7, 9-20, 22, test chr1, 8, 21; mouse: training chr2-7, 9-18,
test chr1, 8, 19), the sizes of the training and test sets, and the auROC and auPRC of the model on the test
chromosomes of both species. In the mouse `.json` files the `pwm_files` entry names the folders `PWMs_MOUSE_mono` and
`PWMs_MOUSE_di`; in this repository they are `PWMs_mono_MOUSE/` and `PWMs_di_MOUSE/`.

Median auROC and auPRC of the 35 models of each kind on the test chromosomes of their training species (field
`test_set_same_species` of the `.json` files):

| | human `mono` | human `di` | human `mono_di` | mouse `mono` | mouse `di` | mouse `mono_di` |
|---|---|---|---|---|---|---|
| auROC | 0.8906 | 0.9021 | 0.9021 | 0.8948 | 0.8968 | 0.8973 |
| auPRC | 0.3142 | 0.3251 | 0.3360 | 0.3274 | 0.3441 | 0.3563 |

The tool stops with an error if a model, its `.json` file or the PWM file of one of its features is missing, or if a
SARUS run fails.

### Transcription factors

```
ANDR   AP2A   CEBPB  COE1   CTCF   E2F4   ERG    ESR1
FLI1   GATA1  GATA2  GATA3  GCR    HNF4A  IRF1   IRF4
JUND   MAFK   MAX    MYC    PPARG  PRGR   REST   RUNX1
RXRA   SOX2   SPI1   SRF    STA5A  STAT1  STAT3  TAL1
TF65   TFE2   USF2
```

`python scanning_tool.py --list_tfs` prints this list. `python scanning_tool.py --list_models --tf CTCF` prints the
feature sets available for one TF (add `--species mouse` for the mouse models).

## PWMs

Four folders, each with one subfolder per TF (35 TFs):

| folder | content | files |
|---|---|---|
| `PWMs_mono_HUMAN/` | human monoPWMs, `<TF>/<k>.pwm` | 1,447 |
| `PWMs_di_HUMAN/` | human diPWMs, `<TF>/<k>.dpwm` | 754 |
| `PWMs_mono_MOUSE/` | mouse monoPWMs, `<TF>/<k>.pwm` | 1,178 |
| `PWMs_di_MOUSE/` | mouse diPWMs, `<TF>/<k>.dpwm` | 617 |

The models and PWMs are also archived on Zenodo, concept DOI
[10.5281/zenodo.14927303](https://doi.org/10.5281/zenodo.14927303) (all versions).

## Directory layout

```
ArChIPelago-TFBS-finder/
├── scanning_tool.py            command-line tool
├── requirements.txt
├── pytest.ini
├── Models/<TF>/                Random Forest models (.sav) and feature specifications (.json)
├── PWMs_mono_HUMAN/<TF>/       human monoPWMs (.pwm)
├── PWMs_di_HUMAN/<TF>/         human diPWMs (.dpwm)
├── PWMs_mono_MOUSE/<TF>/       mouse monoPWMs (.pwm)
├── PWMs_di_MOUSE/<TF>/         mouse diPWMs (.dpwm)
├── sarus/                      SPRY-SARUS (git submodule)
├── synthetic_CTCF_demo.fasta   demo input, one sequence of 10,000 bp
├── realdata_CTCF_test.fasta    demo input, one sequence of 12,020 bp
├── tests/
│   ├── test_scanning_tool.py   unit tests
│   └── test_integration.py     integration tests
└── .github/workflows/tests.yml
```

## Usage

```
python scanning_tool.py -f <FASTA> --tf <TF> [options]
```

### Required

| option | meaning |
|---|---|
| `-f`, `--fasta` | input FASTA file |
| `--tf`, `--tf_name` | TF name as in the list above (e.g. `CTCF`), or `all` for all 35 TFs one after another |

### Optional

| option | default | meaning |
|---|---|---|
| `-o`, `--output` | `Results/` | output directory |
| `--species` | `human` | `human` or `mouse`: species of the model and of the default PWM folders |
| `--pwm_type` | `mono_di` | feature set of the model: `mono`, `di` or `mono_di` |
| `--frame` | `300` | window length in bp (300 bp is the length of the training sequences) |
| `--step` | `150` | window step in bp |
| `--fdr` | `0.1` | q-value threshold for significant windows |
| `--prob` | `0.7` | probability threshold for significant windows, used with `--no-null` |
| `--n_null` | `50` | number of dinucleotide-shuffled copies for the null distribution |
| `--no-null` | off | skip the null distribution (no p-values or q-values) |
| `--models_dir` | `Models/` | models folder |
| `--pwm_mono_dir` | `PWMs_mono_HUMAN/`, or `PWMs_mono_MOUSE/` with `--species mouse` | monoPWM folder |
| `--pwm_di_dir` | `PWMs_di_HUMAN/`, or `PWMs_di_MOUSE/` with `--species mouse` | diPWM folder |
| `--sarus_jar` | see [SARUS jar](#4-sarus-jar) | SARUS jar file |
| `-v`, `--verbose` | off | debug logging |

Default folders, including `Results/`, are relative to the folder of `scanning_tool.py`, not to the working directory.

### Information only

| option | meaning |
|---|---|
| `--list_tfs` | print the 35 TFs and exit |
| `--list_models` | print the feature sets available for `--tf` (and `--species`) and exit |

`python scanning_tool.py --help` prints all options.

## Output files

Files written to the output directory (`<TF>` becomes `<TF>_MOUSE` with `--species mouse`):

| file | content |
|---|---|
| `<TF>_predictions_full.tsv` | all windows |
| `<TF>_significant_FDR<fdr>.tsv` | windows with q-value <= `--fdr`, sorted by q-value (null distribution used) |
| `<TF>_significant_prob<prob>.tsv` | windows with probability >= `--prob`, sorted by probability (with `--no-null`) |
| `<TF>_significant.bed` | the significant windows in BED format |

`<fdr>` and `<prob>` are the threshold values as Python prints them, for example `CTCF_significant_FDR0.1.tsv`. If
the null distribution yields no windows, the tool logs a warning and selects windows by `--prob` instead.

Columns of the TSV files:

| column | content |
|---|---|
| `window_id` | `<sequence id>@<index>`; the index counts the scanned windows of the sequence (0, 1, 2, ...) |
| `position` | 0-based start of the window in the sequence |
| `predicted_probability` | Random Forest probability of a binding site in the window |
| `empirical_pvalue` | (number of null probabilities >= the observed one + 1) / (number of null windows + 1); not written with `--no-null` |
| `qvalue` | Benjamini-Hochberg q-value; not written with `--no-null` |

The BED file has six columns: sequence id, window start, window start + `--frame`, `<TF>_hit` (without `_MOUSE`), the
probability x 1000 truncated to an integer, and `.`. A short summary (number of windows, probability mean, median and
maximum, null pool size and its 95th percentile, number of significant windows) is printed to standard output.

## Examples

Run from the repository folder:

```bash
# CTCF, human model, null distribution with 50 shuffles, FDR 0.1
python scanning_tool.py -f synthetic_CTCF_demo.fasta --tf CTCF

# FDR 0.05, output in another folder
python scanning_tool.py -f realdata_CTCF_test.fasta --tf CTCF --fdr 0.05 -o Results/realdata

# no null distribution, probability threshold 0.7
python scanning_tool.py -f synthetic_CTCF_demo.fasta --tf CTCF --no-null

# monoPWM model, debug logging
python scanning_tool.py -f synthetic_CTCF_demo.fasta --tf CTCF --pwm_type mono -v

# mouse model with the mouse PWMs
python scanning_tool.py -f synthetic_CTCF_demo.fasta --tf CTCF --species mouse

# window 500 bp, step 250 bp
python scanning_tool.py -f realdata_CTCF_test.fasta --tf CTCF --frame 500 --step 250

# all 35 TFs (35 scans, each with its own null distribution)
python scanning_tool.py -f realdata_CTCF_test.fasta --tf all -o Results/all_tfs
```

With `--frame 300 --step 150`, `synthetic_CTCF_demo.fasta` gives 65 windows and `realdata_CTCF_test.fasta` gives 79.
The test workflow checks the second number.

## Testing

```bash
pytest                                # all tests
pytest tests/test_scanning_tool.py    # unit tests
pytest tests/test_integration.py      # integration tests
```

`tests/test_scanning_tool.py` has 99 unit tests, which build their own small input files. `tests/test_integration.py`
has 41 integration tests, which use `Models/`, the four PWM folders, the two FASTA files and
`sarus/releases/sarus-2.2.3.jar` (or `../sarus/releases/sarus-2.2.3.jar`); the tests that need the jar fail if it is
missing. All tests need the packages of `requirements.txt`, and the integration tests also need Java. On every push,
`.github/workflows/tests.yml` runs all tests on Python 3.10 and 3.12 and, on Python 3.12, a CTCF scan of
`realdata_CTCF_test.fasta` with the null distribution.

## How it works

1. **Input.** The FASTA sequences are read and converted to upper case.
2. **Windows.** Each sequence is cut into windows of `--frame` bp with a step of `--step` bp. A sequence shorter than
   `--frame` is one window. Windows that contain IUPAC ambiguity codes (B, D, H, K, M, N, R, S, V, W, Y) are skipped.
3. **PWM scores.** SPRY-SARUS scans every window with each PWM of the model
   (`--skipn --show-non-matching --output-scoring-mode score besthit`) and gives the log-odds score of the best hit.
4. **Features.** The scores form one column per PWM, in the order of the model's `.json` file, and each column is
   standardised with the training mean and standard deviation from the `.json` file.
5. **Prediction.** The Random Forest gives the probability of a binding site for each window.
6. **Null distribution** (skipped with `--no-null`). All input sequences are joined into one sequence, which is
   shuffled `--n_null` times with its dinucleotide counts kept exactly (Altschul-Erickson algorithm, seed 42). The
   shuffled copies are cut into windows and scored in the same way.
7. **p-values.** Each window's probability is compared with the pooled probabilities of the null windows.
8. **FDR.** p-values are converted to Benjamini-Hochberg q-values; windows with q-value <= `--fdr` are significant.
9. **Export.** The full table, the significant windows and the BED file are written.

## Citation

Kravchenko P., Vorontsov I.E., Grosse I., Makeev V.J., Kulakovskiy I.V., and Penzar D.D. (2026). Classic machine
learning on top of multiple position weight matrices improves genomic prediction of transcription factor binding
sites.

Zenodo: [10.5281/zenodo.14927303](https://doi.org/10.5281/zenodo.14927303) (all versions); release v1.0.0:
[10.5281/zenodo.23036475](https://doi.org/10.5281/zenodo.23036475).

## Licence

MIT licence (see the header of `scanning_tool.py`).
