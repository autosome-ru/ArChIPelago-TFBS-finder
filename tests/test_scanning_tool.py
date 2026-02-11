#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for ArChIPelago-TFBS-finder scanning tool.

Tests for:
- Constants and configuration
- Data classes (ScanConfig, ScanResult)
- Sequence processing (parse_fasta, has_ambiguous, generate_windows, etc.)
- PWM file handling (get_pwm_files, find_model_file, etc.)
- Feature matrix construction
- Empirical null calibration (dinucleotide_shuffle, compute_empirical_pvalues, BH)
- Export functions
- CLI argument parsing
- Model predictions (mocked)
"""

import os
import sys
import pytest
import tempfile
import numpy as np
import pandas as pd
from pathlib import Path
from unittest.mock import patch, MagicMock
from collections import Counter

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scanning_tool import (
    # Constants
    AVAILABLE_TFS,
    AMBIGUOUS_NUCLEOTIDES,
    PWM_TYPES,
    # Data classes
    ScanConfig,
    ScanResult,
    # Logging
    setup_logging,
    # Sequence processing
    parse_fasta,
    has_ambiguous,
    generate_windows,
    windows_to_fasta,
    reverse_complement,
    # PWM handling
    _extract_pwm_number,
    get_pwm_files,
    find_model_file,
    get_available_models,
    validate_paths,
    # Feature matrix
    build_feature_matrix,
    map_model_features,
    compute_pwm_summary,
    # Null calibration
    dinucleotide_shuffle,
    generate_null_sequences,
    compute_empirical_pvalues,
    benjamini_hochberg,
    # Export
    export_results,
    print_summary,
    # CLI
    parse_arguments,
    main,
)


# =============================================================================
# Fixtures
# =============================================================================

@pytest.fixture
def sample_fasta(tmp_path):
    """Create a sample FASTA file with known content."""
    content = (
        ">seq1\n"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
        "ATCG\n"
        ">seq2\n"
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n"
        ">seq3_ambiguous\n"
        "ATCGATCGATCGNNNNATCGATCGATCG\n"
    )
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(content)
    return fasta_file


@pytest.fixture
def sample_pwm_dir(tmp_path):
    """Create a mock PWM directory structure."""
    ctcf_dir = tmp_path / "PWMs_mono" / "CTCF"
    ctcf_dir.mkdir(parents=True)
    for i in range(5):
        pwm = ctcf_dir / f"{i}.pwm"
        pwm.write_text(
            f">PWM_{i}\n-0.5\t0.8\t-0.3\t0.1\n0.6\t-1.2\t0.9\t-0.4\n")
    return tmp_path / "PWMs_mono"


@pytest.fixture
def sample_dpwm_dir(tmp_path):
    """Create a mock dinucleotide PWM directory."""
    ctcf_dir = tmp_path / "PWMs_di" / "CTCF"
    ctcf_dir.mkdir(parents=True)
    for i in range(3):
        dpwm = ctcf_dir / f"{i}.dpwm"
        cols = "\t".join(["0.1"] * 16)
        dpwm.write_text(f">DPWM_{i}\n{cols}\n{cols}\n")
    return tmp_path / "PWMs_di"


@pytest.fixture
def mock_model_dir(tmp_path):
    """Create a mock model directory with dummy .sav files."""
    ctcf_dir = tmp_path / "Models" / "CTCF"
    ctcf_dir.mkdir(parents=True)
    files = [
        "model_CTCF_HUMAN_RandomForestClassifier_mono_RF_on_all_PWMs.sav",
        "model_CTCF_HUMAN_RandomForestClassifier_di_RF_on_all_PWMs.sav",
        "finalized_model_CTCF_HUMAN_RandomForestClassifier_"
        "mono_di_full_MODEL_all_features.sav",
    ]
    for fn in files:
        (ctcf_dir / fn).touch()
    return tmp_path / "Models"


@pytest.fixture
def mock_sarus_jar(tmp_path):
    """Create a dummy SARUS jar path."""
    jar = tmp_path / "sarus.jar"
    jar.touch()
    return jar


@pytest.fixture
def sample_config(tmp_path, sample_fasta, mock_model_dir, sample_pwm_dir,
                  sample_dpwm_dir, mock_sarus_jar):
    """Create a ScanConfig with all paths pointing at test fixtures."""
    return ScanConfig(
        fasta_file=sample_fasta,
        tf_name="CTCF",
        output_dir=tmp_path / "output",
        models_dir=mock_model_dir,
        pwm_mono_dir=sample_pwm_dir,
        pwm_di_dir=sample_dpwm_dir,
        sarus_jar=mock_sarus_jar,
        pwm_type="mono_di",
        frame=200,
        step=50,
    )


# =============================================================================
# Test Constants
# =============================================================================

class TestConstants:
    """Test module-level constants."""

    def test_available_tfs_count(self):
        assert len(AVAILABLE_TFS) == 36

    def test_available_tfs_contains_common_factors(self):
        for tf in ['CTCF', 'P53', 'MYC', 'MAX', 'STAT1', 'GATA1']:
            assert tf in AVAILABLE_TFS

    def test_ambiguous_nucleotides_set(self):
        expected = set('BDHKMSVWYNR')
        assert AMBIGUOUS_NUCLEOTIDES == expected

    def test_pwm_types(self):
        assert PWM_TYPES == ['mono', 'di', 'mono_di']


# =============================================================================
# Test Data Classes
# =============================================================================

class TestScanConfig:
    def test_creation_defaults(self, sample_config):
        assert sample_config.tf_name == "CTCF"
        assert sample_config.pwm_type == "mono_di"
        assert sample_config.frame == 200
        assert sample_config.step == 50
        assert sample_config.fdr_threshold == 0.1
        assert sample_config.skip_null is False

    def test_custom_values(self, tmp_path):
        cfg = ScanConfig(
            fasta_file=tmp_path / "x.fa",
            tf_name="P53",
            output_dir=tmp_path,
            models_dir=tmp_path,
            pwm_mono_dir=tmp_path,
            pwm_di_dir=tmp_path,
            sarus_jar=tmp_path / "s.jar",
            pwm_type="mono",
            frame=500,
            step=10,
            fdr_threshold=0.05,
            n_null_shuffles=20,
            prob_threshold=0.8,
            skip_null=True,
        )
        assert cfg.tf_name == "P53"
        assert cfg.frame == 500
        assert cfg.skip_null is True


class TestScanResult:
    def test_creation(self):
        df = pd.DataFrame({'a': [1, 2]})
        r = ScanResult(tf_name="CTCF", predictions_df=df,
                       n_windows=10, n_positive=3)
        assert r.tf_name == "CTCF"
        assert r.n_windows == 10
        assert r.n_positive == 3
        assert r.model_is_real is True


# =============================================================================
# Test Sequence Processing
# =============================================================================

class TestParseFasta:
    def test_parse_basic(self, sample_fasta):
        seqs = parse_fasta(sample_fasta)
        assert len(seqs) == 3
        assert seqs[0][0] == "seq1"
        assert seqs[1][0] == "seq2"
        assert all(s.upper() == s for _, s in seqs)

    def test_empty_file(self, tmp_path):
        fa = tmp_path / "empty.fa"
        fa.write_text("")
        seqs = parse_fasta(fa)
        assert seqs == []


class TestHasAmbiguous:
    @pytest.mark.parametrize("seq,expected", [
        ("ATCGATCG", False),
        ("ACGT", False),
        ("ATCNGATC", True),
        ("YATCGATCG", True),
        ("RATCG", True),
        ("BATCG", True),
        ("atcg", False),
        ("atcNg", True),
    ])
    def test_various(self, seq, expected):
        assert has_ambiguous(seq) == expected


class TestGenerateWindows:
    def test_basic_windowing(self):
        sequences = [("s1", "A" * 500)]
        wins = generate_windows(sequences, frame=200, step=100)
        # positions: 0, 100, 200, 300 -> 4 windows
        assert len(wins) == 4
        assert all(len(w[1]) == 200 for w in wins)

    def test_short_sequence_returned_as_single(self):
        sequences = [("s1", "ACGT")]
        wins = generate_windows(sequences, frame=200, step=50)
        assert len(wins) == 1
        assert wins[0][1] == "ACGT"
        assert wins[0][2] == 0

    def test_ambiguous_windows_skipped(self):
        seq = "A" * 100 + "N" * 50 + "A" * 200
        sequences = [("s1", seq)]
        wins = generate_windows(sequences, frame=200, step=50)
        for wid, wseq, pos in wins:
            assert 'N' not in wseq

    def test_window_ids_sequential(self):
        sequences = [("seq1", "A" * 1000)]
        wins = generate_windows(sequences, frame=200, step=100)
        ids = [w[0] for w in wins]
        for i, wid in enumerate(ids):
            assert wid == f"seq1@{i}"

    def test_exact_frame_length(self):
        sequences = [("s1", "A" * 200)]
        wins = generate_windows(sequences, frame=200, step=50)
        assert len(wins) == 1
        assert wins[0][2] == 0

    def test_multiple_sequences(self):
        sequences = [("s1", "A" * 400), ("s2", "C" * 400)]
        wins = generate_windows(sequences, frame=200, step=200)
        s1_wins = [w for w in wins if w[0].startswith("s1")]
        s2_wins = [w for w in wins if w[0].startswith("s2")]
        assert len(s1_wins) == 2
        assert len(s2_wins) == 2


class TestWindowsToFasta:
    def test_writes_correct_format(self, tmp_path):
        windows = [
            ("s1@0", "ACGTACGT", 0),
            ("s1@1", "TGCATGCA", 100),
        ]
        out = tmp_path / "windows.fa"
        n = windows_to_fasta(windows, out)
        assert n == 2
        content = out.read_text()
        assert ">s1@0\nACGTACGT\n" in content
        assert ">s1@1\nTGCATGCA\n" in content


class TestReverseComplement:
    @pytest.mark.parametrize("seq,expected", [
        ("ATCG", "CGAT"),
        ("AAAA", "TTTT"),
        ("ACGT", "ACGT"),
        ("N", "N"),
        ("", ""),
    ])
    def test_rc(self, seq, expected):
        assert reverse_complement(seq) == expected


# =============================================================================
# Test PWM Handling
# =============================================================================

class TestExtractPwmNumber:
    def test_numeric(self):
        assert _extract_pwm_number(Path("5.pwm")) == 5
        assert _extract_pwm_number(Path("10.pwm")) == 10
        assert _extract_pwm_number(Path("0.pwm")) == 0

    def test_non_numeric(self):
        assert _extract_pwm_number(Path("abc.pwm")) == 999999


class TestGetPwmFiles:
    def test_gets_sorted_files(self, sample_pwm_dir):
        files = get_pwm_files(sample_pwm_dir, "CTCF", ".pwm")
        assert len(files) == 5
        stems = [f.stem for f in files]
        assert stems == ['0', '1', '2', '3', '4']

    def test_missing_tf(self, sample_pwm_dir):
        files = get_pwm_files(sample_pwm_dir, "NONEXISTENT", ".pwm")
        assert files == []


class TestFindModelFile:
    def test_find_mono(self, mock_model_dir):
        p = find_model_file(mock_model_dir, "CTCF", "mono")
        assert p is not None
        assert "mono_RF_on_all_PWMs" in p.name

    def test_find_di(self, mock_model_dir):
        p = find_model_file(mock_model_dir, "CTCF", "di")
        assert p is not None
        assert "di_RF_on_all_PWMs" in p.name

    def test_find_mono_di(self, mock_model_dir):
        p = find_model_file(mock_model_dir, "CTCF", "mono_di")
        assert p is not None
        assert "mono_di_full_MODEL_all_features" in p.name

    def test_invalid_tf(self, mock_model_dir):
        assert find_model_file(mock_model_dir, "INVALID", "mono") is None

    def test_invalid_pwm_type(self, mock_model_dir):
        assert find_model_file(mock_model_dir, "CTCF", "bad") is None

    def test_prefers_base_model(self, tmp_path):
        d = tmp_path / "M" / "CTCF"
        d.mkdir(parents=True)
        (d / "model_CTCF_HUMAN_RandomForestClassifier_"
             "mono_RF_on_all_PWMs_SLIM_m1.sav").touch()
        (d / "model_CTCF_HUMAN_RandomForestClassifier_"
             "mono_RF_on_all_PWMs.sav").touch()
        p = find_model_file(tmp_path / "M", "CTCF", "mono")
        assert "SLIM" not in p.name


class TestGetAvailableModels:
    def test_returns_models(self, mock_model_dir):
        models = get_available_models(mock_model_dir, "CTCF")
        assert len(models) == 3

    def test_missing_tf(self, mock_model_dir):
        models = get_available_models(mock_model_dir, "NONEXISTENT")
        assert models == {}


class TestValidatePaths:
    def test_all_valid(self, sample_config):
        logger = setup_logging(False)
        assert validate_paths(sample_config, logger) is True

    def test_missing_fasta(self, sample_config):
        logger = setup_logging(False)
        sample_config.fasta_file = Path("/nonexistent/file.fa")
        assert validate_paths(sample_config, logger) is False

    def test_missing_models_dir(self, sample_config):
        logger = setup_logging(False)
        sample_config.models_dir = Path("/nonexistent/models")
        assert validate_paths(sample_config, logger) is False

    def test_missing_sarus_jar(self, sample_config):
        logger = setup_logging(False)
        sample_config.sarus_jar = Path("/nonexistent/sarus.jar")
        assert validate_paths(sample_config, logger) is False


# =============================================================================
# Test Feature Matrix
# =============================================================================

class TestBuildFeatureMatrix:
    def test_builds_correct_shape(self, tmp_path):
        windows = [("s1@0", "ACGT", 0), ("s1@1", "TGCA", 100)]
        mono_dir = tmp_path / "scores"
        mono_dir.mkdir()
        score_f = mono_dir / "0_mono.tab"
        score_f.write_text("1.5\n2.3\n")

        mono_pwms = [Path("0.pwm")]
        mono_scores = {"0": score_f}
        df, mcols, dcols = build_feature_matrix(
            windows, mono_scores, {}, mono_pwms, [])
        assert len(df) == 2
        assert "mono_0" in mcols
        assert len(dcols) == 0
        assert df["mono_0"].tolist() == [1.5, 2.3]

    def test_length_mismatch_skipped(self, tmp_path):
        windows = [("s1@0", "ACGT", 0)]
        score_f = tmp_path / "0_mono.tab"
        score_f.write_text("1.5\n2.3\n")  # 2 scores but 1 window
        mono_pwms = [Path("0.pwm")]
        df, mcols, dcols = build_feature_matrix(
            windows, {"0": score_f}, {}, mono_pwms, [])
        assert "mono_0" not in mcols


class TestMapModelFeatures:
    def test_no_feature_names(self):
        model = MagicMock(spec=[])
        model.n_features_in_ = 3
        del model.feature_names_in_
        cols = map_model_features(
            model, ["mono_0", "mono_1"], ["di_0", "di_1"], "mono_di")
        assert cols == ["mono_0", "mono_1", "di_0"]

    def test_with_feature_names(self):
        model = MagicMock()
        model.feature_names_in_ = np.array(['4', '5', '6', '7'])
        model.n_features_in_ = 4
        cols = map_model_features(
            model,
            ["mono_0", "mono_1", "mono_2"],
            ["di_0", "di_1"],
            "mono_di",
        )
        assert len(cols) > 0


class TestComputePwmSummary:
    def test_with_mono_cols(self):
        df = pd.DataFrame({
            "mono_0": [1.0, 2.0, 3.0],
            "mono_1": [0.5, 1.5, 2.5],
        })
        result = compute_pwm_summary(df, ["mono_0", "mono_1"], [])
        assert "pwm_mean_mono" in result.columns
        assert "pwm_best_mono" in result.columns
        assert "pwm_n_high_mono" in result.columns
        assert "pwm_consensus" in result.columns

    def test_empty_cols(self):
        df = pd.DataFrame({"x": [1, 2]})
        result = compute_pwm_summary(df, [], [])
        assert result["pwm_mean_mono"].tolist() == [0.0, 0.0]
        assert result["pwm_mean_di"].tolist() == [0.0, 0.0]


# =============================================================================
# Test Empirical Null Calibration
# =============================================================================

class TestDinucleotideShuffle:
    def test_preserves_length(self):
        seq = "ACGTACGTACGTACGT"
        shuffled = dinucleotide_shuffle(seq)
        assert len(shuffled) == len(seq)

    def test_preserves_composition(self):
        seq = "AACCGGTTAACCGGTT"
        shuffled = dinucleotide_shuffle(seq)
        assert Counter(shuffled) == Counter(seq)

    def test_preserves_dinucleotide_frequencies(self):
        seq = "ACGTACGTACGTACGTACGT" * 10
        rng = np.random.default_rng(42)
        shuffled = dinucleotide_shuffle(seq, rng)
        orig_di = Counter(seq[i:i+2] for i in range(len(seq) - 1))
        shuf_di = Counter(
            shuffled[i:i+2] for i in range(len(shuffled) - 1))
        assert orig_di == shuf_di

    def test_short_seq_returned_as_is(self):
        assert dinucleotide_shuffle("AC") == "AC"
        assert dinucleotide_shuffle("A") == "A"

    def test_deterministic_with_seed(self):
        seq = "ACGTACGTACGTACGT"
        r1 = dinucleotide_shuffle(seq, np.random.default_rng(0))
        r2 = dinucleotide_shuffle(seq, np.random.default_rng(0))
        assert r1 == r2


class TestGenerateNullSequences:
    def test_correct_count(self):
        seqs = generate_null_sequences("ACGTACGTACGT", n_shuffles=5)
        assert len(seqs) == 5

    def test_all_same_length(self):
        seq = "AACCGGTTAACCGGTT"
        seqs = generate_null_sequences(seq, n_shuffles=10)
        assert all(len(s) == len(seq) for s in seqs)


class TestComputeEmpiricalPvalues:
    def test_basic(self):
        observed = np.array([0.9, 0.5, 0.1])
        null_pool = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        pvals = compute_empirical_pvalues(observed, null_pool)
        assert len(pvals) == 3
        assert pvals[0] < pvals[1] < pvals[2]

    def test_no_zero_pvalues(self):
        observed = np.array([100.0])
        null_pool = np.array([0.1, 0.2, 0.3])
        pvals = compute_empirical_pvalues(observed, null_pool)
        assert pvals[0] > 0

    def test_all_below_null(self):
        observed = np.array([0.01])
        null_pool = np.array([0.5, 0.6, 0.7, 0.8, 0.9])
        pvals = compute_empirical_pvalues(observed, null_pool)
        assert pvals[0] == pytest.approx(1.0)

    def test_range(self):
        observed = np.random.rand(100)
        null_pool = np.random.rand(1000)
        pvals = compute_empirical_pvalues(observed, null_pool)
        assert np.all(pvals > 0)
        assert np.all(pvals <= 1)


class TestBenjaminiHochberg:
    def test_empty(self):
        q = benjamini_hochberg(np.array([]))
        assert len(q) == 0

    def test_single(self):
        q = benjamini_hochberg(np.array([0.05]))
        assert q[0] == pytest.approx(0.05)

    def test_monotonicity(self):
        pvals = np.array([0.01, 0.04, 0.03, 0.20])
        q = benjamini_hochberg(pvals)
        sorted_q = q[np.argsort(pvals)]
        for i in range(len(sorted_q) - 1):
            assert sorted_q[i] <= sorted_q[i + 1]

    def test_qvalues_bounded(self):
        pvals = np.random.rand(50)
        q = benjamini_hochberg(pvals)
        assert np.all(q >= 0)
        assert np.all(q <= 1)

    def test_known_values(self):
        pvals = np.array([0.01, 0.04, 0.03, 0.20])
        q = benjamini_hochberg(pvals)
        assert np.argmin(q) == np.argmin(pvals)


# =============================================================================
# Test Export Functions
# =============================================================================

class TestExportResults:
    def test_creates_output_files(self, tmp_path):
        preds = pd.DataFrame({
            'window_id': ['s1@0', 's1@1'],
            'position': [0, 100],
            'predicted_probability': [0.9, 0.3],
            'empirical_pvalue': [0.001, 0.5],
            'qvalue': [0.01, 0.8],
        })
        config = ScanConfig(
            fasta_file=tmp_path / "x.fa",
            tf_name="CTCF",
            output_dir=tmp_path / "out",
            models_dir=tmp_path,
            pwm_mono_dir=tmp_path,
            pwm_di_dir=tmp_path,
            sarus_jar=tmp_path / "s.jar",
            frame=200,
        )
        logger = setup_logging(False)
        outputs = export_results(preds, config, logger)

        assert outputs['full'].exists()
        assert outputs['significant'].exists()
        assert outputs['bed'].exists()

        full = pd.read_csv(outputs['full'], sep='\t')
        assert len(full) == 2
        assert 'predicted_probability' in full.columns

        sig = pd.read_csv(outputs['significant'], sep='\t')
        assert len(sig) == 1
        assert sig.iloc[0]['window_id'] == 's1@0'

    def test_bed_format(self, tmp_path):
        preds = pd.DataFrame({
            'window_id': ['chr1@0'],
            'position': [500],
            'predicted_probability': [0.95],
            'qvalue': [0.001],
        })
        config = ScanConfig(
            fasta_file=tmp_path / "x.fa",
            tf_name="CTCF",
            output_dir=tmp_path / "out",
            models_dir=tmp_path,
            pwm_mono_dir=tmp_path,
            pwm_di_dir=tmp_path,
            sarus_jar=tmp_path / "s.jar",
            frame=301,
        )
        logger = setup_logging(False)
        outputs = export_results(preds, config, logger)
        bed_lines = outputs['bed'].read_text().strip().split('\n')
        assert len(bed_lines) == 1
        fields = bed_lines[0].split('\t')
        assert fields[0] == 'chr1'
        assert fields[1] == '500'
        assert fields[2] == '801'  # 500 + 301
        assert fields[3] == 'CTCF_hit'

    def test_no_null_uses_prob_threshold(self, tmp_path):
        preds = pd.DataFrame({
            'window_id': ['s1@0', 's1@1'],
            'position': [0, 100],
            'predicted_probability': [0.9, 0.3],
        })
        config = ScanConfig(
            fasta_file=tmp_path / "x.fa",
            tf_name="CTCF",
            output_dir=tmp_path / "out",
            models_dir=tmp_path,
            pwm_mono_dir=tmp_path,
            pwm_di_dir=tmp_path,
            sarus_jar=tmp_path / "s.jar",
            prob_threshold=0.5,
        )
        logger = setup_logging(False)
        outputs = export_results(preds, config, logger)
        sig = pd.read_csv(outputs['significant'], sep='\t')
        assert len(sig) == 1


class TestPrintSummary:
    def test_runs_without_error(self, capsys):
        preds = pd.DataFrame({
            'predicted_probability': [0.9, 0.5, 0.3],
            'qvalue': [0.01, 0.5, 0.9],
        })
        config = ScanConfig(
            fasta_file=Path("x.fa"),
            tf_name="CTCF",
            output_dir=Path("out"),
            models_dir=Path("m"),
            pwm_mono_dir=Path("p"),
            pwm_di_dir=Path("p"),
            sarus_jar=Path("s.jar"),
        )
        print_summary(preds, config, np.array([0.1, 0.2, 0.3]))
        captured = capsys.readouterr()
        assert "CTCF" in captured.out
        assert "FDR-significant" in captured.out


# =============================================================================
# Test CLI Argument Parsing
# =============================================================================

class TestParseArguments:
    def test_basic(self):
        args = parse_arguments(['-f', 'in.fa', '--tf', 'CTCF'])
        assert args.fasta == 'in.fa'
        assert args.tf_name == 'CTCF'
        assert args.pwm_type == 'mono_di'
        assert args.frame == 301
        assert args.step == 150
        assert args.fdr == 0.1
        assert args.skip_null is False

    def test_all_options(self):
        args = parse_arguments([
            '-f', 'in.fa', '--tf', 'P53',
            '-o', '/tmp/out', '--pwm_type', 'mono',
            '--frame', '500', '--step', '100',
            '--fdr', '0.05', '--prob', '0.8',
            '--n_null', '20', '--no-null', '-v',
        ])
        assert args.tf_name == 'P53'
        assert args.pwm_type == 'mono'
        assert args.frame == 500
        assert args.step == 100
        assert args.fdr == 0.05
        assert args.prob == 0.8
        assert args.n_null == 20
        assert args.skip_null is True
        assert args.verbose is True

    def test_list_tfs(self):
        args = parse_arguments(['--list_tfs'])
        assert args.list_tfs is True

    def test_list_models(self):
        args = parse_arguments(['--list_models', '--tf', 'CTCF'])
        assert args.list_models is True
        assert args.tf_name == 'CTCF'


class TestMain:
    def test_list_tfs(self, capsys):
        rc = main(['--list_tfs'])
        assert rc == 0
        captured = capsys.readouterr()
        assert "CTCF" in captured.out
        assert "P53" in captured.out

    def test_missing_args(self):
        rc = main([])
        assert rc == 1

    def test_list_models(self, mock_model_dir, capsys):
        rc = main(['--list_models', '--tf', 'CTCF',
                    '--models_dir', str(mock_model_dir)])
        assert rc == 0
        captured = capsys.readouterr()
        assert "CTCF" in captured.out

    def test_list_models_no_tf(self, capsys):
        rc = main(['--list_models'])
        assert rc == 1


# =============================================================================
# Test Model Predictions (mocked)
# =============================================================================

class TestModelPredictions:
    def test_prediction_shape_and_range(self):
        from sklearn.ensemble import RandomForestClassifier
        model = RandomForestClassifier(n_estimators=5, random_state=42)
        X = np.random.randn(50, 10)
        y = np.random.randint(0, 2, 50)
        model.fit(X, y)
        X_test = np.random.randn(20, 10)
        probs = model.predict_proba(X_test)[:, 1]
        assert probs.shape == (20,)
        assert np.all(probs >= 0) and np.all(probs <= 1)

    def test_model_save_load(self, tmp_path):
        from sklearn.ensemble import RandomForestClassifier
        import joblib as jl
        model = RandomForestClassifier(n_estimators=3, random_state=0)
        X = np.random.randn(30, 5)
        y = np.random.randint(0, 2, 30)
        model.fit(X, y)
        path = tmp_path / "model.sav"
        jl.dump(model, path)
        loaded = jl.load(path)
        assert hasattr(loaded, 'predict_proba')
        p1 = model.predict_proba(X)
        p2 = loaded.predict_proba(X)
        np.testing.assert_array_equal(p1, p2)


# =============================================================================
# Test Performance
# =============================================================================

class TestPerformance:
    def test_large_feature_matrix(self):
        import time
        t0 = time.time()
        data = np.random.randn(10000, 200)
        df = pd.DataFrame(data)
        assert time.time() - t0 < 5.0
        assert df.shape == (10000, 200)

    def test_bh_performance(self):
        import time
        pvals = np.random.rand(100000)
        t0 = time.time()
        q = benjamini_hochberg(pvals)
        assert time.time() - t0 < 2.0
        assert len(q) == 100000

    def test_empirical_pvalue_performance(self):
        import time
        obs = np.random.rand(1000)
        null = np.random.rand(50000)
        t0 = time.time()
        p = compute_empirical_pvalues(obs, null)
        assert time.time() - t0 < 1.0
        assert len(p) == 1000


# =============================================================================
# Main
# =============================================================================

if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
