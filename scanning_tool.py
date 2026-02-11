#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ArChIPelago-TFBS-finder: Command-Line Scanning Tool for Transcription Factor Binding Sites

Scans DNA sequences for transcription factor binding sites using pre-trained
ArChIPelago Random Forest models that aggregate multiple PWM (Position Weight
Matrix) scores.

The tool addresses the false-positive problem inherent in repurposing a classifier
as a scanner by constructing an **empirical null distribution** from
dinucleotide-shuffled sequences and applying Benjamini-Hochberg FDR control.

Author: Pavel Kravchenko
License: MIT

Usage examples:
    # Scan for CTCF binding sites with default FDR threshold
    python scanning_tool.py -f sequences.fasta --tf CTCF

    # Scan with specific FDR and PWM type
    python scanning_tool.py -f sequences.fasta --tf CTCF --fdr 0.05 --pwm_type mono_di

    # List available transcription factors
    python scanning_tool.py --list_tfs

    # Scan with custom window size and disable null calibration
    python scanning_tool.py -f sequences.fasta --tf CTCF --frame 301 --step 50 --no-null
"""

import os
import sys
import re
import argparse
import subprocess
import tempfile
import logging
import random
from copy import copy
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
from collections import defaultdict

import numpy as np
import pandas as pd
import joblib
from Bio import SeqIO


# =============================================================================
# Constants
# =============================================================================

AVAILABLE_TFS = [
    'ANDR', 'AP2A', 'CEBPB', 'COE1', 'CTCF', 'E2F4', 'ERG', 'ESR1',
    'FLI1', 'GATA1', 'GATA2', 'GATA3', 'GCR', 'HNF4A', 'IRF1', 'IRF4',
    'JUND', 'MAFK', 'MAX', 'MYC', 'P53', 'PPARG', 'PRGR', 'REST',
    'RUNX1', 'RXRA', 'SOX2', 'SPI1', 'SRF', 'STA5A', 'STAT1', 'STAT3',
    'TAL1', 'TF65', 'TFE2', 'USF2',
]

AMBIGUOUS_NUCLEOTIDES = set('BDHKMSVWYNR')

PWM_TYPES = ['mono', 'di', 'mono_di']


# =============================================================================
# Logging
# =============================================================================

def setup_logging(verbose: bool = False) -> logging.Logger:
    """Configure logging for the scanning tool."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    return logging.getLogger('archipelago')


# =============================================================================
# Data Classes
# =============================================================================

@dataclass
class ScanConfig:
    """All parameters for a single scanning run."""
    fasta_file: Path
    tf_name: str
    output_dir: Path
    models_dir: Path
    pwm_mono_dir: Path
    pwm_di_dir: Path
    sarus_jar: Path
    pwm_type: str = 'mono_di'
    frame: int = 301
    step: int = 150
    fdr_threshold: float = 0.1
    n_null_shuffles: int = 50
    prob_threshold: float = 0.7
    skip_null: bool = False


@dataclass
class ScanResult:
    """Container for scanning output."""
    tf_name: str
    predictions_df: pd.DataFrame
    n_windows: int
    n_positive: int
    model_is_real: bool = True


# =============================================================================
# Sequence Processing
# =============================================================================

def parse_fasta(fasta_file: Path) -> List[Tuple[str, str]]:
    """Parse FASTA file; returns list of (id, sequence) tuples."""
    sequences = []
    for record in SeqIO.parse(str(fasta_file), "fasta"):
        sequences.append((str(record.id), str(record.seq).upper()))
    return sequences


def has_ambiguous(seq: str) -> bool:
    """True if sequence contains IUPAC ambiguity codes."""
    return bool(set(seq.upper()) & AMBIGUOUS_NUCLEOTIDES)


def generate_windows(sequences, frame, step):
    """
    Sliding-window generator.

    Returns list of (window_id, sequence, position) tuples.
    """
    windows = []
    for seq_id, seq in sequences:
        if len(seq) >= frame:
            idx = 0
            for pos in range(0, len(seq) - frame + 1, step):
                w = seq[pos:pos + frame]
                if not has_ambiguous(w):
                    windows.append((f"{seq_id}@{idx}", w, pos))
                    idx += 1
        else:
            if not has_ambiguous(seq):
                windows.append((f"{seq_id}@0", seq, 0))
    return windows


def windows_to_fasta(windows, output_file):
    """Write windows to a temporary FASTA file."""
    with open(output_file, 'w') as f:
        for wid, seq, _ in windows:
            f.write(f">{wid}\n{seq}\n")
    return len(windows)


def reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(comp.get(b, 'N') for b in reversed(seq))


# =============================================================================
# PWM Handling and SARUS Scoring
# =============================================================================

def _extract_pwm_number(path: Path) -> int:
    """Numerical sort key for PWM filenames (0.pwm, 1.pwm, ..., 10.pwm)."""
    m = re.match(r'^(\d+)', path.stem)
    return int(m.group(1)) if m else 999999


def get_pwm_files(pwm_dir: Path, tf_name: str, ext: str = '.pwm') -> List[Path]:
    """Get numerically sorted list of PWM files for a TF."""
    d = pwm_dir / tf_name
    if not d.exists():
        return []
    return sorted(d.glob(f"*{ext}"), key=_extract_pwm_number)


def find_model_file(models_dir: Path, tf_name: str,
                    pwm_type: str) -> Optional[Path]:
    """Locate the appropriate .sav model file."""
    d = models_dir / tf_name
    if not d.exists():
        return None
    suffix_map = {
        'mono':    'mono_RF_on_all_PWMs',
        'di':      'di_RF_on_all_PWMs',
        'mono_di': 'mono_di_full_MODEL_all_features',
    }
    suffix = suffix_map.get(pwm_type)
    if not suffix:
        return None
    # Prefer base model (no SLIM/ChIPMunk extensions)
    for f in d.glob("*.sav"):
        if suffix in f.name and not any(
                x in f.name for x in ['SLIM', 'ChIPMunk']):
            return f
    for f in d.glob("*.sav"):
        if suffix in f.name:
            return f
    return None


def get_available_models(models_dir: Path, tf_name: str) -> Dict[str, Path]:
    """Get available model files for a transcription factor."""
    d = models_dir / tf_name
    if not d.exists():
        return {}
    models = {}
    for f in d.glob("*.sav"):
        key = f.stem.replace(f'model_{tf_name}_HUMAN_', '').replace(
            f'finalized_model_{tf_name}_HUMAN_', '')
        models[key] = f
    return models


def validate_paths(config: ScanConfig, logger: logging.Logger) -> bool:
    """Validate that all required paths exist."""
    errors = []
    if not Path(config.fasta_file).exists():
        errors.append(f"FASTA file not found: {config.fasta_file}")
    if not Path(config.models_dir).exists():
        errors.append(f"Models directory not found: {config.models_dir}")
    if not Path(config.sarus_jar).exists():
        errors.append(f"SARUS JAR not found: {config.sarus_jar}")
    if config.pwm_type in ('mono', 'mono_di'):
        if not Path(config.pwm_mono_dir).exists():
            errors.append(
                f"Mono PWM directory not found: {config.pwm_mono_dir}")
    if config.pwm_type in ('di', 'mono_di'):
        if not Path(config.pwm_di_dir).exists():
            errors.append(
                f"Di PWM directory not found: {config.pwm_di_dir}")
    for error in errors:
        logger.error(error)
    return len(errors) == 0


def run_sarus(sarus_jar, fasta_file, pwm_file, output_file, is_di=False):
    """Run SARUS PWM scanner; returns True on success."""
    cls = "ru.autosome.di.SARUS" if is_di else "ru.autosome.SARUS"
    cmd = (f'java -Xmx2G -cp {sarus_jar} {cls} {fasta_file} {pwm_file} '
           f'--skipn --show-non-matching --output-scoring-mode score besthit '
           f'2>/dev/null | grep -v ">" > {output_file}')
    r = subprocess.run(cmd, shell=True, capture_output=True)
    return r.returncode == 0 and Path(output_file).exists()


def scan_all_pwms(sarus_jar, fasta_file, pwm_files, out_dir, is_di=False):
    """Scan with all PWMs; returns {pwm_stem: score_file_path}."""
    out_dir.mkdir(parents=True, exist_ok=True)
    results = {}
    tag = 'di' if is_di else 'mono'
    for pwm in pwm_files:
        out_f = out_dir / f"{pwm.stem}_{tag}.tab"
        if run_sarus(sarus_jar, fasta_file, pwm, out_f, is_di):
            results[pwm.stem] = out_f
    return results


# =============================================================================
# Feature Matrix Construction and Model Interface
# =============================================================================

def build_feature_matrix(windows, mono_scores, di_scores,
                         mono_pwms, di_pwms):
    """
    Build the feature matrix in the exact column order used during training:
    [mono_0, mono_1, ..., mono_N, di_0, di_1, ..., di_M].

    Returns (DataFrame, mono_col_names, di_col_names).
    """
    base = pd.DataFrame({
        'window_id': [w[0] for w in windows],
        'sequence':  [w[1] for w in windows],
        'position':  [w[2] for w in windows],
    })
    mono_cols, di_cols = [], []
    score_series = {}
    for pwm in mono_pwms:
        if pwm.stem in mono_scores:
            col = f"mono_{pwm.stem}"
            try:
                scores = pd.read_csv(
                    mono_scores[pwm.stem], header=None, sep='\t')[0].values
                if len(scores) == len(windows):
                    score_series[col] = scores
                    mono_cols.append(col)
            except Exception:
                pass
    for pwm in di_pwms:
        if pwm.stem in di_scores:
            col = f"di_{pwm.stem}"
            try:
                scores = pd.read_csv(
                    di_scores[pwm.stem], header=None, sep='\t')[0].values
                if len(scores) == len(windows):
                    score_series[col] = scores
                    di_cols.append(col)
            except Exception:
                pass
    if score_series:
        df = pd.concat([base, pd.DataFrame(score_series)], axis=1)
    else:
        df = base
    return df, mono_cols, di_cols


def map_model_features(model, mono_cols, di_cols, pwm_type):
    """
    Map model feature names to the scanning feature columns.

    Training data layout: columns 0-3 = metadata, column 4+ = PWM scores
    in order [mono_0, ..., mono_N, di_0, ..., di_M].
    """
    if not hasattr(model, 'feature_names_in_'):
        all_cols = mono_cols + di_cols
        n_expected = getattr(model, 'n_features_in_', len(all_cols))
        return all_cols[:n_expected]

    numeric_feats = [n for n in model.feature_names_in_
                     if not n.startswith('pred')]
    if pwm_type == 'mono_di':
        n_mono_train = (len(numeric_feats) + 1) // 2
    elif pwm_type == 'mono':
        n_mono_train = len(numeric_feats)
    else:
        n_mono_train = 0

    selected = []
    for name in model.feature_names_in_:
        if name.startswith('pred'):
            continue
        try:
            idx = int(name) - 4
        except ValueError:
            continue
        if idx < 0:
            continue
        if n_mono_train > 0 and idx < n_mono_train:
            if idx < len(mono_cols):
                selected.append(mono_cols[idx])
        else:
            di_idx = idx - n_mono_train
            if di_idx < len(di_cols):
                selected.append(di_cols[di_idx])

    n_expected = getattr(model, 'n_features_in_', 0)
    return selected if selected else (mono_cols + di_cols)[:n_expected]


def compute_pwm_summary(df, mono_cols, di_cols, model=None, feat_cols=None):
    """
    Compute per-window PWM summary statistics.

    Adds columns: pwm_mean_mono, pwm_mean_di, pwm_best_mono, pwm_best_di,
    pwm_consensus, pwm_n_high_mono, pwm_n_high_di.
    """
    if mono_cols:
        mono_vals = df[mono_cols].values
        df['pwm_mean_mono'] = mono_vals.mean(axis=1)
        df['pwm_best_mono'] = mono_vals.max(axis=1)
        mono_z = (mono_vals - mono_vals.mean(axis=0)) / (
            mono_vals.std(axis=0) + 1e-10)
        df['pwm_n_high_mono'] = (mono_z > 2).sum(axis=1)
    else:
        df['pwm_mean_mono'] = 0.0
        df['pwm_best_mono'] = 0.0
        df['pwm_n_high_mono'] = 0

    if di_cols:
        di_vals = df[di_cols].values
        df['pwm_mean_di'] = di_vals.mean(axis=1)
        df['pwm_best_di'] = di_vals.max(axis=1)
        di_z = (di_vals - di_vals.mean(axis=0)) / (
            di_vals.std(axis=0) + 1e-10)
        df['pwm_n_high_di'] = (di_z > 2).sum(axis=1)
    else:
        df['pwm_mean_di'] = 0.0
        df['pwm_best_di'] = 0.0
        df['pwm_n_high_di'] = 0

    if (model is not None and feat_cols is not None
            and hasattr(model, 'feature_importances_')):
        imp = model.feature_importances_
        available = [c for c in feat_cols if c in df.columns]
        imp_map = dict(zip(feat_cols, imp))
        weights = np.array([imp_map.get(c, 0) for c in available])
        weights = weights / (weights.sum() + 1e-10)
        df['pwm_consensus'] = df[available].values @ weights
    else:
        all_pwm = mono_cols + di_cols
        df['pwm_consensus'] = df[all_pwm].mean(axis=1) if all_pwm else 0.0

    return df


# =============================================================================
# Complete Scanning Pipeline
# =============================================================================

def scan_sequence(config: ScanConfig, verbose=True,
                  logger=None) -> Optional[ScanResult]:
    """
    Complete pipeline: FASTA -> windows -> SARUS scoring -> RF prediction.

    Returns all windows with their predicted probabilities and PWM summaries.
    """
    if logger is None:
        logger = logging.getLogger('archipelago')

    # 1. Parse and window
    sequences = parse_fasta(config.fasta_file)
    windows = generate_windows(sequences, config.frame, config.step)
    if not windows:
        logger.warning("No valid windows generated from input sequences.")
        return None

    # 2. Get PWMs
    mono_pwms = (get_pwm_files(config.pwm_mono_dir, config.tf_name, '.pwm')
                 if config.pwm_type in ('mono', 'mono_di') else [])
    di_pwms = (get_pwm_files(config.pwm_di_dir, config.tf_name, '.dpwm')
               if config.pwm_type in ('di', 'mono_di') else [])
    if not mono_pwms and not di_pwms:
        logger.error(f"No PWM files found for {config.tf_name}.")
        return None

    # 3. Load model
    model_path = find_model_file(config.models_dir, config.tf_name,
                                 config.pwm_type)
    if model_path is None:
        logger.error(f"No model file found for {config.tf_name} "
                     f"({config.pwm_type}).")
        return None
    model = joblib.load(str(model_path))

    # 4. Scan with SARUS
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        win_fasta = tmp / "windows.fasta"
        windows_to_fasta(windows, win_fasta)
        feat_dir = tmp / "features"

        mono_scores = (scan_all_pwms(
            config.sarus_jar, win_fasta, mono_pwms, feat_dir, is_di=False
        ) if mono_pwms else {})
        di_scores = (scan_all_pwms(
            config.sarus_jar, win_fasta, di_pwms, feat_dir, is_di=True
        ) if di_pwms else {})

        feature_df, mono_cols, di_cols = build_feature_matrix(
            windows, mono_scores, di_scores, mono_pwms, di_pwms)

    # 5. Select features matching model expectations
    feat_cols = map_model_features(model, mono_cols, di_cols, config.pwm_type)
    if verbose:
        logger.info(f"  Windows: {len(windows)}, Features: {len(feat_cols)}")

    # 6. Predict
    X = feature_df[feat_cols].values
    probs = model.predict_proba(X)[:, 1]
    feature_df['predicted_probability'] = probs

    # 7. Compute per-window PWM summary scores
    feature_df = compute_pwm_summary(feature_df, mono_cols, di_cols,
                                     model, feat_cols)

    n_pos = int(np.sum(probs >= config.prob_threshold))
    return ScanResult(
        tf_name=config.tf_name,
        predictions_df=feature_df,
        n_windows=len(windows),
        n_positive=n_pos,
        model_is_real=True,
    )


# =============================================================================
# Empirical Null Calibration
# =============================================================================

def dinucleotide_shuffle(seq: str, rng=None) -> str:
    """
    Shuffle a DNA sequence while preserving dinucleotide frequencies.

    Uses the Altschul-Erickson algorithm: build an Eulerian path over the
    dinucleotide graph, which preserves exact edge (dinucleotide) counts.
    Falls back to mononucleotide shuffle if the graph is disconnected.
    """
    if rng is None:
        rng = np.random.default_rng()
    seq = seq.upper()
    n = len(seq)
    if n <= 2:
        return seq

    # Build adjacency lists
    edges = defaultdict(list)
    for i in range(n - 1):
        edges[seq[i]].append(seq[i + 1])

    # Shuffle each adjacency list, reserving last edge for completion
    last_edge = {}
    for nuc in edges:
        lst = edges[nuc]
        rng.shuffle(lst)
        last_edge[nuc] = lst.pop()

    # Build Eulerian path
    result = [seq[0]]
    current = seq[0]
    for _ in range(n - 2):
        if edges[current]:
            nxt = edges[current].pop()
        else:
            nxt = last_edge.pop(current, None)
            if nxt is None:
                # Fallback: mononucleotide shuffle
                lst = list(seq)
                rng.shuffle(lst)
                return ''.join(lst)
        result.append(nxt)
        current = nxt

    if current in last_edge:
        result.append(last_edge[current])
    else:
        result.append(rng.choice(list('ACGT')))

    return ''.join(result)


def generate_null_sequences(seq: str, n_shuffles: int = 50,
                            seed: int = 42) -> List[str]:
    """Generate *n* dinucleotide-shuffled null sequences."""
    rng = np.random.default_rng(seed)
    return [dinucleotide_shuffle(seq, rng) for _ in range(n_shuffles)]


def compute_empirical_pvalues(observed: np.ndarray,
                              null_pool: np.ndarray) -> np.ndarray:
    """
    Empirical p-value for each observed score against the null pool.

    Uses the (rank + 1) / (N + 1) formulation to avoid zero p-values.
    """
    n_null = len(null_pool)
    sorted_null = np.sort(null_pool)
    counts_ge = n_null - np.searchsorted(sorted_null, observed, side='left')
    pvals = (counts_ge + 1) / (n_null + 1)
    return pvals


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values with monotonicity enforcement."""
    n = len(pvalues)
    if n == 0:
        return np.array([])
    order = np.argsort(pvalues)
    sorted_p = pvalues[order]
    q = np.minimum(1.0, sorted_p * n / np.arange(1, n + 1, dtype=float))
    for i in range(n - 2, -1, -1):
        q[i] = min(q[i], q[i + 1])
    qvalues = np.empty_like(pvalues)
    qvalues[order] = q
    return qvalues


def build_null_distribution(config: ScanConfig, input_sequence: str,
                            logger: logging.Logger) -> np.ndarray:
    """
    Build empirical null distribution by scanning dinucleotide-shuffled
    versions of the input sequence with the same model.

    Returns array of null predicted probabilities.
    """
    null_probs_all = []
    rng = np.random.default_rng(42)

    for i in range(config.n_null_shuffles):
        null_seq = dinucleotide_shuffle(input_sequence, rng)
        with tempfile.NamedTemporaryFile(
                mode='w', suffix='.fasta', delete=False) as tmp_f:
            tmp_f.write(f">null_{i}\n{null_seq}\n")
            null_fasta = Path(tmp_f.name)

        cfg_null = copy(config)
        cfg_null.fasta_file = null_fasta

        null_result = scan_sequence(cfg_null, verbose=False, logger=logger)
        if null_result is not None:
            null_probs_all.extend(
                null_result.predictions_df['predicted_probability'].values)

        null_fasta.unlink(missing_ok=True)

        if (i + 1) % 10 == 0:
            logger.info(
                f"  Null calibration: {i + 1}/{config.n_null_shuffles}")

    return np.array(null_probs_all)


# =============================================================================
# Export Functions
# =============================================================================

def export_results(preds: pd.DataFrame, config: ScanConfig,
                   logger: logging.Logger) -> Dict[str, Path]:
    """
    Export prediction results to TSV and BED files.

    Returns dict of output file paths.
    """
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}

    base_cols = ['window_id', 'position', 'predicted_probability']
    extra_cols = ['empirical_pvalue', 'qvalue']
    export_cols = base_cols + [c for c in extra_cols if c in preds.columns]

    # Full predictions
    out_full = output_dir / f"{config.tf_name}_predictions_full.tsv"
    preds[export_cols].to_csv(out_full, sep='\t', index=False,
                              float_format='%.6g')
    logger.info(f"Full predictions: {out_full}  ({len(preds)} windows)")
    outputs['full'] = out_full

    # Significant hits (by q-value if available, else by prob threshold)
    if 'qvalue' in preds.columns:
        sig = preds[preds['qvalue'] <= config.fdr_threshold].sort_values(
            'qvalue')
        tag = f"FDR{config.fdr_threshold}"
    else:
        sig = preds[
            preds['predicted_probability'] >= config.prob_threshold
        ].sort_values('predicted_probability', ascending=False)
        tag = f"prob{config.prob_threshold}"

    out_sig = output_dir / f"{config.tf_name}_significant_{tag}.tsv"
    sig[export_cols].to_csv(out_sig, sep='\t', index=False,
                            float_format='%.6g')
    logger.info(f"Significant hits: {out_sig}  ({len(sig)} windows)")
    outputs['significant'] = out_sig

    # BED format
    out_bed = output_dir / f"{config.tf_name}_significant.bed"
    with open(out_bed, 'w') as f:
        for _, row in sig.iterrows():
            seq_id = row['window_id'].split('@')[0]
            start = int(row['position'])
            end = start + config.frame
            score = int(row['predicted_probability'] * 1000)
            f.write(f"{seq_id}\t{start}\t{end}\t"
                    f"{config.tf_name}_hit\t{score}\t.\n")
    logger.info(f"BED file: {out_bed}")
    outputs['bed'] = out_bed

    return outputs


def print_summary(preds: pd.DataFrame, config: ScanConfig,
                  null_pool: Optional[np.ndarray] = None):
    """Print a human-readable summary to stdout."""
    probs = preds['predicted_probability'].values
    print("=" * 64)
    print("  ArChIPelago-TFBS-finder: Scan Summary")
    print("=" * 64)
    print(f"  Transcription factor:       {config.tf_name}")
    print(f"  PWM type:                   {config.pwm_type}")
    print(f"  Window / step:              {config.frame} / {config.step} bp")
    print(f"  Total windows scanned:      {len(preds)}")
    print(f"  Prob mean:                  {probs.mean():.4f}")
    print(f"  Prob median:                {np.median(probs):.4f}")
    print(f"  Prob max:                   {probs.max():.4f}")

    if null_pool is not None and len(null_pool) > 0:
        print(f"  Null shuffles:              {config.n_null_shuffles}")
        print(f"  Null pool size:             {len(null_pool)}")
        print(f"  Null 95th percentile:       "
              f"{np.percentile(null_pool, 95):.4f}")

    if 'qvalue' in preds.columns:
        n_sig = int((preds['qvalue'] <= config.fdr_threshold).sum())
        print(f"  FDR threshold:              {config.fdr_threshold}")
        print(f"  FDR-significant windows:    {n_sig}")
    else:
        n_above = int((probs >= config.prob_threshold).sum())
        print(f"  Prob threshold:             {config.prob_threshold}")
        print(f"  Windows above threshold:    {n_above}")
    print("=" * 64)


# =============================================================================
# CLI
# =============================================================================

def parse_arguments(argv=None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        prog='scanning_tool',
        description=(
            'ArChIPelago-TFBS-finder: Scan DNA sequences for '
            'transcription factor binding sites using pre-trained '
            'Random Forest models with empirical FDR calibration.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -f sequences.fasta --tf CTCF
  %(prog)s -f sequences.fasta --tf CTCF --fdr 0.05 --pwm_type mono_di
  %(prog)s -f sequences.fasta --tf CTCF --no-null --prob 0.7
  %(prog)s --list_tfs
  %(prog)s --list_models --tf CTCF
        """,
    )

    # Required
    parser.add_argument('-f', '--fasta', type=str,
                        help='Path to input FASTA file')
    parser.add_argument('--tf', '--tf_name', dest='tf_name', type=str,
                        help='Transcription factor name (or "all")')

    # Output
    parser.add_argument('-o', '--output', default=None,
                        help='Output directory (default: ./Results)')

    # Paths
    parser.add_argument('--models_dir', default=None,
                        help='Models directory (default: auto-detect)')
    parser.add_argument('--pwm_mono_dir', default=None,
                        help='Mononucleotide PWM directory')
    parser.add_argument('--pwm_di_dir', default=None,
                        help='Dinucleotide PWM directory')
    parser.add_argument('--sarus_jar', default=None,
                        help='Path to SARUS JAR file')

    # Scanning parameters
    parser.add_argument('--pwm_type', choices=PWM_TYPES, default='mono_di',
                        help='PWM type (default: mono_di)')
    parser.add_argument('--frame', type=int, default=301,
                        help='Sliding window size in bp (default: 301)')
    parser.add_argument('--step', type=int, default=150,
                        help='Sliding window step in bp (default: 150)')

    # FDR / thresholds
    parser.add_argument('--fdr', type=float, default=0.1,
                        help='FDR threshold for significance (default: 0.1)')
    parser.add_argument('--prob', type=float, default=0.7,
                        help='Probability threshold when --no-null '
                             '(default: 0.7)')
    parser.add_argument('--n_null', type=int, default=50,
                        help='Number of null shuffles (default: 50)')
    parser.add_argument('--no-null', dest='skip_null', action='store_true',
                        help='Skip empirical null calibration (faster, '
                             'less accurate)')

    # Misc
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')
    parser.add_argument('--list_tfs', action='store_true',
                        help='List available transcription factors and exit')
    parser.add_argument('--list_models', action='store_true',
                        help='List available models for --tf and exit')

    return parser.parse_args(argv)


def _resolve_default_paths(script_dir: Path) -> dict:
    """Auto-detect resource directories relative to the script."""
    candidates_models = [
        script_dir / 'Models_sklearn13',
        script_dir / 'Models',
    ]
    candidates_sarus = [
        script_dir.parent / 'sarus' / 'releases' / 'sarus-2.0.1.jar',
        script_dir / 'sarus' / 'releases' / 'sarus-2.0.1.jar',
        script_dir / 'sarus' / 'sarus-2.0.1.jar',
    ]
    defaults = {
        'models_dir': next(
            (p for p in candidates_models if p.exists()), candidates_models[0]),
        'pwm_mono_dir': script_dir / 'PWMs_mono_HUMAN',
        'pwm_di_dir': script_dir / 'PWMs_di_HUMAN',
        'sarus_jar': next(
            (p for p in candidates_sarus if p.exists()), candidates_sarus[0]),
        'output_dir': script_dir / 'Results',
    }
    return defaults


def main(argv=None) -> int:
    """Main entry point."""
    args = parse_arguments(argv)
    logger = setup_logging(args.verbose)

    # --list_tfs
    if args.list_tfs:
        print("\nAvailable Transcription Factors (36 from HOCOMOCO v11):")
        print("-" * 52)
        for i, tf in enumerate(AVAILABLE_TFS):
            print(f"  {tf:10}", end='' if (i + 1) % 4 else '\n')
        if len(AVAILABLE_TFS) % 4:
            print()
        print()
        return 0

    # Resolve default paths
    script_dir = Path(__file__).resolve().parent
    defaults = _resolve_default_paths(script_dir)

    models_dir = (Path(args.models_dir) if args.models_dir
                  else defaults['models_dir'])

    # --list_models
    if args.list_models:
        if not args.tf_name:
            print("Error: --tf is required with --list_models")
            return 1
        models = get_available_models(models_dir, args.tf_name)
        if not models:
            print(f"No models found for {args.tf_name} in {models_dir}")
            return 1
        print(f"\nAvailable models for {args.tf_name}:")
        print("-" * 60)
        for key in sorted(models.keys()):
            print(f"  {key}")
        print()
        return 0

    # Validate required arguments
    if not args.fasta or not args.tf_name:
        logger.error("Missing required arguments: --fasta and --tf")
        logger.error("Run with --help for usage information.")
        return 1

    # Build config
    config = ScanConfig(
        fasta_file=Path(args.fasta).resolve(),
        tf_name=args.tf_name,
        output_dir=(Path(args.output).resolve() if args.output
                    else defaults['output_dir']),
        models_dir=models_dir,
        pwm_mono_dir=(Path(args.pwm_mono_dir) if args.pwm_mono_dir
                      else defaults['pwm_mono_dir']),
        pwm_di_dir=(Path(args.pwm_di_dir) if args.pwm_di_dir
                    else defaults['pwm_di_dir']),
        sarus_jar=(Path(args.sarus_jar) if args.sarus_jar
                   else defaults['sarus_jar']),
        pwm_type=args.pwm_type,
        frame=args.frame,
        step=args.step,
        fdr_threshold=args.fdr,
        n_null_shuffles=args.n_null,
        prob_threshold=args.prob,
        skip_null=args.skip_null,
    )

    if not validate_paths(config, logger):
        return 1

    # Determine TFs to scan
    if args.tf_name.lower() == 'all':
        tf_list = AVAILABLE_TFS
    else:
        if args.tf_name not in AVAILABLE_TFS:
            logger.error(
                f"Unknown TF: {args.tf_name}. "
                f"Use --list_tfs to see options.")
            return 1
        tf_list = [args.tf_name]

    # Scan each TF
    for tf in tf_list:
        tf_config = copy(config)
        tf_config.tf_name = tf

        logger.info(f"{'=' * 60}")
        logger.info(f"Scanning for {tf}")
        logger.info(f"{'=' * 60}")

        # Step 1: scan real sequence
        result = scan_sequence(tf_config, verbose=args.verbose, logger=logger)
        if result is None:
            logger.error(f"Scanning failed for {tf}. Skipping.")
            continue

        preds = result.predictions_df

        # Step 2: empirical null calibration (unless --no-null)
        null_pool = None
        if not config.skip_null:
            logger.info("Building empirical null distribution...")
            sequences = parse_fasta(tf_config.fasta_file)
            full_seq = ''.join(s for _, s in sequences)
            null_pool = build_null_distribution(tf_config, full_seq, logger)

            if len(null_pool) > 0:
                probs = preds['predicted_probability'].values
                emp_pvals = compute_empirical_pvalues(probs, null_pool)
                qvals = benjamini_hochberg(emp_pvals)
                preds['empirical_pvalue'] = emp_pvals
                preds['qvalue'] = qvals
            else:
                logger.warning(
                    "Null calibration produced no scores; "
                    "falling back to probability threshold.")

        # Step 3: export
        export_results(preds, tf_config, logger)

        # Step 4: print summary
        print_summary(preds, tf_config, null_pool)

    return 0


if __name__ == '__main__':
    sys.exit(main())
