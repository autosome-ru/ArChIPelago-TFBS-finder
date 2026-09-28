#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integration tests for ArChIPelago-TFBS-finder scanning tool.

These tests exercise real PWM files, real pre-trained models, and real
FASTA data that ship with the repository.  They are skipped automatically
when the required data directories are absent.

Run only integration tests:
    pytest -m integration -v

Run everything *except* integration tests:
    pytest -m "not integration" -v
"""

import os
import sys
import pytest
import tempfile
import numpy as np
import pandas as pd
from pathlib import Path
from collections import Counter

# Ensure the package root is importable
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    import joblib
    from Bio import SeqIO
    from sklearn.ensemble import RandomForestClassifier
    IMPORTS_AVAILABLE = True
except ImportError:
    IMPORTS_AVAILABLE = False

from scanning_tool import (
    AVAILABLE_TFS,
    PWM_TYPES,
    ScanConfig,
    ScanResult,
    parse_fasta,
    has_ambiguous,
    generate_windows,
    windows_to_fasta,
    reverse_complement,
    find_model_file,
    load_model,
    model_pwm_files,
    get_available_models,
    validate_paths,
    build_feature_matrix,
    dinucleotide_shuffle,
    generate_null_sequences,
    compute_empirical_pvalues,
    benjamini_hochberg,
    export_results,
    print_summary,
    setup_logging,
    main,
)


# =============================================================================
# Fixtures – real data paths
# =============================================================================

@pytest.fixture
def project_root():
    """Return project root (parent of tests/)."""
    return Path(__file__).resolve().parent.parent


@pytest.fixture
def real_mono_pwm_dir(project_root):
    d = project_root / "PWMs_mono_HUMAN"
    if not d.exists():
        pytest.skip("PWMs_mono_HUMAN directory not found")
    return d


@pytest.fixture
def real_di_pwm_dir(project_root):
    d = project_root / "PWMs_di_HUMAN"
    if not d.exists():
        pytest.skip("PWMs_di_HUMAN directory not found")
    return d


@pytest.fixture
def real_models_dir(project_root):
    d = project_root / "Models"
    if not d.exists():
        pytest.skip("Models directory not found")
    return d


@pytest.fixture
def real_fasta_file(project_root):
    """Use the synthetic demo FASTA shipped with the repo."""
    f = project_root / "synthetic_CTCF_demo.fasta"
    if not f.exists():
        # Fallback to the real test FASTA
        f = project_root / "realdata_CTCF_test.fasta"
    if not f.exists():
        f = project_root / "CTCF_data" / "test_pos_id_HUMAN.fasta"
    if not f.exists():
        pytest.skip("No suitable FASTA file found")
    return f


@pytest.fixture
def sarus_jar(project_root):
    jar = project_root.parent / "sarus" / "releases" / "sarus-2.2.3.jar"
    if not jar.exists():
        pytest.skip("SARUS jar not found")
    return jar


# =============================================================================
# PWM File Validation
# =============================================================================

@pytest.mark.integration
@pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Required imports not available")
class TestRealPWMFiles:
    """Validate the on-disk PWM files expected by ArChIPelago."""

    def test_mono_pwm_structure(self, real_mono_pwm_dir):
        """Mono PWM files must have a >header and rows of 4 tab-separated floats."""
        ctcf_dir = real_mono_pwm_dir / "CTCF"
        assert ctcf_dir.exists(), "CTCF directory not found"
        pwm_files = list(ctcf_dir.glob("*.pwm"))
        assert len(pwm_files) > 0

        with open(pwm_files[0]) as f:
            lines = f.readlines()
        assert lines[0].startswith(">")
        for line in lines[1:]:
            if line.strip():
                vals = line.strip().split("\t")
                assert len(vals) == 4, f"Expected 4 cols, got {len(vals)}"
                for v in vals:
                    float(v)  # should not raise

    def test_di_pwm_structure(self, real_di_pwm_dir):
        """Di-PWM files must have a >header and rows of 16 tab-separated floats."""
        ctcf_dir = real_di_pwm_dir / "CTCF"
        assert ctcf_dir.exists()
        dpwm_files = list(ctcf_dir.glob("*.dpwm"))
        assert len(dpwm_files) > 0

        with open(dpwm_files[0]) as f:
            lines = f.readlines()
        assert lines[0].startswith(">")
        for line in lines[1:]:
            if line.strip():
                vals = line.strip().split("\t")
                assert len(vals) == 16, f"Expected 16 cols, got {len(vals)}"

    def test_all_tfs_have_mono_pwms(self, real_mono_pwm_dir):
        """Every TF in AVAILABLE_TFS should have at least one mono PWM."""
        present = {d.name for d in real_mono_pwm_dir.iterdir() if d.is_dir()}
        for tf in AVAILABLE_TFS:
            if tf not in present:
                continue  # some TFs may not ship in the test data
            pwm_files = list((real_mono_pwm_dir / tf).glob("*.pwm"))
            assert len(pwm_files) > 0, f"No .pwm files for {tf}"

    def test_ctcf_has_many_pwms(self, real_mono_pwm_dir):
        """CTCF should have a large number of PWMs (from HOCOMOCO)."""
        ctcf = real_mono_pwm_dir / "CTCF"
        assert len(list(ctcf.glob("*.pwm"))) >= 50



# =============================================================================
# Model File Validation
# =============================================================================

@pytest.mark.integration
@pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Required imports not available")
class TestRealModelFiles:
    """Validate pre-trained model files: one model per TF and PWM type."""

    def test_every_tf_has_three_models(self, real_models_dir):
        for tf in AVAILABLE_TFS:
            assert sorted(get_available_models(real_models_dir, tf)) == sorted(PWM_TYPES), tf

    @pytest.mark.parametrize("pwm_type", ["mono", "di", "mono_di"])
    def test_ctcf_model_and_spec(self, real_models_dir, real_mono_pwm_dir,
                                 real_di_pwm_dir, pwm_type):
        model, spec = load_model(real_models_dir, "CTCF", pwm_type)
        assert hasattr(model, "predict_proba")
        files = model_pwm_files(spec, real_mono_pwm_dir, real_di_pwm_dir, "CTCF")
        assert len(files) == model.n_features_in_
        kinds = {f.split("_")[0] for f in spec["features"]}
        assert kinds == ({"mono", "di"} if pwm_type == "mono_di" else {pwm_type})

    def test_every_model_matches_its_pwms(self, real_models_dir, real_mono_pwm_dir,
                                          real_di_pwm_dir):
        """Each model uses every PWM of its type(s) of the TF, each exactly once."""
        for tf in AVAILABLE_TFS:
            n_mono = len(list((real_mono_pwm_dir / tf).glob("*.pwm")))
            n_di = len(list((real_di_pwm_dir / tf).glob("*.dpwm")))
            for pwm_type, n in (("mono", n_mono), ("di", n_di), ("mono_di", n_mono + n_di)):
                model, spec = load_model(real_models_dir, tf, pwm_type)
                assert len(set(spec["features"])) == model.n_features_in_ == n, (tf, pwm_type)


# =============================================================================
# FASTA Processing with Real Data
# =============================================================================

@pytest.mark.integration
@pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Required imports not available")
class TestFASTAProcessing:
    """Test FASTA handling on real files."""

    def test_parse_fasta_nonempty(self, real_fasta_file):
        seqs = parse_fasta(real_fasta_file)
        assert len(seqs) > 0

    def test_sequences_are_valid_dna(self, real_fasta_file):
        valid = set("ACGTNBDHKMSVWYR")
        for sid, seq in parse_fasta(real_fasta_file):
            chars = set(seq)
            bad = chars - valid
            assert len(bad) == 0, f"Invalid chars in {sid}: {bad}"

    def test_sliding_window_on_real_data(self, real_fasta_file):
        seqs = parse_fasta(real_fasta_file)
        frame, step = 301, 150
        windows = generate_windows(seqs, frame, step)
        assert len(windows) > 0, "No windows generated"
        for wid, wseq, pos in windows:
            # Windows longer than frame or with ambiguity should not appear
            if len(wseq) == frame:
                assert not has_ambiguous(wseq)

    def test_windows_to_fasta_roundtrip(self, real_fasta_file, tmp_path):
        seqs = parse_fasta(real_fasta_file)
        windows = generate_windows(seqs, frame=301, step=150)
        out = tmp_path / "windows.fa"
        n = windows_to_fasta(windows, out)
        assert n == len(windows)
        # Read back
        back = list(SeqIO.parse(str(out), "fasta"))
        assert len(back) == n


# =============================================================================
# Data Consistency Across Directories
# =============================================================================

@pytest.mark.integration
@pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Required imports not available")
class TestDataConsistency:
    """Cross-directory consistency checks."""

    def test_mono_di_tf_sets_match(self, real_mono_pwm_dir, real_di_pwm_dir):
        mono_tfs = {d.name for d in real_mono_pwm_dir.iterdir() if d.is_dir()}
        di_tfs = {d.name for d in real_di_pwm_dir.iterdir() if d.is_dir()}
        assert mono_tfs == di_tfs, (
            f"Mismatch: only-mono={mono_tfs - di_tfs}, "
            f"only-di={di_tfs - mono_tfs}")

    def test_tfs_have_pwms_and_models(self, real_mono_pwm_dir, real_models_dir):
        mono_tfs = {d.name for d in real_mono_pwm_dir.iterdir() if d.is_dir()}
        model_tfs = {d.name for d in real_models_dir.iterdir() if d.is_dir()}
        common = mono_tfs & model_tfs
        assert len(common) > 0, "No TFs with both PWMs and models"


# =============================================================================
# Validate Config & Paths with Real Directories
# =============================================================================

@pytest.mark.integration
@pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Required imports not available")
class TestRealValidatePaths:

    def test_validate_paths_all_present(
        self, real_fasta_file, real_models_dir, real_mono_pwm_dir,
        real_di_pwm_dir, sarus_jar, tmp_path,
    ):
        cfg = ScanConfig(
            fasta_file=real_fasta_file,
            tf_name="CTCF",
            output_dir=tmp_path / "out",
            models_dir=real_models_dir,
            pwm_mono_dir=real_mono_pwm_dir,
            pwm_di_dir=real_di_pwm_dir,
            sarus_jar=sarus_jar,
        )
        logger = setup_logging(False)
        assert validate_paths(cfg, logger) is True


# =============================================================================
# Pipeline Pieces with Real Data
# =============================================================================

@pytest.mark.integration
@pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Required imports not available")
class TestPipelineIntegration:
    """Test pipeline pieces using real models."""

    def test_mono_model_prediction_with_real_features(
        self, real_models_dir, real_mono_pwm_dir
    ):
        """Load the mono model, predict on random standardised features."""
        model, spec = load_model(real_models_dir, "CTCF", "mono")
        X = np.random.randn(20, len(spec["features"]))
        probs = model.predict_proba(X)[:, 1]
        assert probs.shape == (20,)
        assert np.all((probs >= 0) & (probs <= 1))

    def test_prediction_aggregation_by_sequence(self):
        """Verify that window-level predictions can be aggregated to per-sequence."""
        preds = pd.DataFrame({
            "window_id": ["s1@0", "s1@1", "s1@2", "s2@0", "s2@1"],
            "predicted_probability": [0.3, 0.8, 0.5, 0.6, 0.4],
        })
        preds["orig_seq"] = preds["window_id"].str.rsplit("@", n=1).str[0]
        agg = preds.groupby("orig_seq")["predicted_probability"].agg(
            ["max", "mean", "count"]).reset_index()
        assert len(agg) == 2
        row_s1 = agg[agg["orig_seq"] == "s1"].iloc[0]
        assert row_s1["max"] == 0.8
        assert row_s1["count"] == 3

    def test_output_tsv_roundtrip(self, tmp_path):
        """Write predictions to TSV and re-read unchanged."""
        df = pd.DataFrame({
            "window_id": ["s1@0", "s2@0"],
            "position": [0, 150],
            "predicted_probability": [0.85, 0.72],
        })
        p = tmp_path / "preds.tsv"
        df.to_csv(p, sep="\t", index=False, float_format="%.6g")
        back = pd.read_csv(p, sep="\t")
        pd.testing.assert_frame_equal(df, back)


# =============================================================================
# Null Calibration with Realistic Sequences
# =============================================================================

@pytest.mark.integration
@pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Required imports not available")
class TestNullCalibrationIntegration:
    """Longer / heavier null-calibration tests."""

    def test_dinucleotide_shuffle_on_real_sequence(self, real_fasta_file):
        seqs = parse_fasta(real_fasta_file)
        if not seqs:
            pytest.skip("No sequences")
        # Take first sequence
        _, seq = seqs[0]
        if has_ambiguous(seq):
            seq = seq.replace("N", "")
        shuffled = dinucleotide_shuffle(seq)
        assert len(shuffled) == len(seq)
        assert Counter(shuffled) == Counter(seq)

    def test_null_sequence_generation(self, real_fasta_file):
        seqs = parse_fasta(real_fasta_file)
        if not seqs:
            pytest.skip("No sequences")
        _, seq = seqs[0]
        seq = seq.replace("N", "")  # remove ambiguous
        nulls = generate_null_sequences(seq[:500], n_shuffles=10, seed=42)
        assert len(nulls) == 10
        assert all(len(n) == len(seq[:500]) for n in nulls)

    def test_empirical_pvalues_monotone(self):
        """Higher observed → lower p-value."""
        null = np.random.rand(5000)
        obs = np.array([0.05, 0.50, 0.95])
        pv = compute_empirical_pvalues(obs, null)
        assert pv[0] > pv[1] > pv[2]

    def test_bh_fdr_on_mixed_pvalues(self):
        """Mix of significant and non-significant; q-values should be valid."""
        pvals = np.concatenate([
            np.random.uniform(0.0, 0.01, 20),   # significant
            np.random.uniform(0.2, 1.0, 80),     # non-significant
        ])
        q = benjamini_hochberg(pvals)
        assert np.all(q >= 0) and np.all(q <= 1)
        # Significant p-values should generally yield lower q-values
        assert q[:20].mean() < q[20:].mean()


# =============================================================================
# Export Functions with Real-ish Data
# =============================================================================

@pytest.mark.integration
@pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Required imports not available")
class TestExportIntegration:

    def test_full_export_pipeline(self, tmp_path):
        preds = pd.DataFrame({
            "window_id": [f"chr1@{i}" for i in range(50)],
            "position": list(range(0, 7500, 150)),
            "predicted_probability": np.random.rand(50),
            "empirical_pvalue": np.random.rand(50),
            "qvalue": np.concatenate([
                np.random.uniform(0.0, 0.05, 10),
                np.random.uniform(0.2, 1.0, 40),
            ]),
        })
        cfg = ScanConfig(
            fasta_file=tmp_path / "in.fa",
            tf_name="CTCF",
            output_dir=tmp_path / "out",
            models_dir=tmp_path,
            pwm_mono_dir=tmp_path,
            pwm_di_dir=tmp_path,
            sarus_jar=tmp_path / "s.jar",
            frame=301,
            fdr_threshold=0.1,
        )
        logger = setup_logging(False)
        outputs = export_results(preds, cfg, logger)

        # Check files exist
        assert outputs["full"].exists()
        assert outputs["significant"].exists()
        assert outputs["bed"].exists()

        # Full file should have all rows
        full = pd.read_csv(outputs["full"], sep="\t")
        assert len(full) == 50

        # Significant file should only have q <= 0.1
        sig = pd.read_csv(outputs["significant"], sep="\t")
        assert all(sig["qvalue"] <= 0.1)

        # BED format sanity
        bed_text = outputs["bed"].read_text().strip()
        if bed_text:
            for line in bed_text.split("\n"):
                fields = line.split("\t")
                assert len(fields) == 6
                start, end = int(fields[1]), int(fields[2])
                assert end - start == 301  # frame

    def test_export_without_null(self, tmp_path):
        """When no qvalue column, export should use prob_threshold."""
        preds = pd.DataFrame({
            "window_id": ["s1@0", "s1@1", "s1@2"],
            "position": [0, 150, 300],
            "predicted_probability": [0.9, 0.6, 0.3],
        })
        cfg = ScanConfig(
            fasta_file=tmp_path / "in.fa",
            tf_name="CTCF",
            output_dir=tmp_path / "out",
            models_dir=tmp_path,
            pwm_mono_dir=tmp_path,
            pwm_di_dir=tmp_path,
            sarus_jar=tmp_path / "s.jar",
            prob_threshold=0.5,
        )
        logger = setup_logging(False)
        outputs = export_results(preds, cfg, logger)
        sig = pd.read_csv(outputs["significant"], sep="\t")
        assert len(sig) == 2  # 0.9 and 0.6 are ≥ 0.5


# =============================================================================
# CLI Integration
# =============================================================================

@pytest.mark.integration
@pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Required imports not available")
class TestCLIIntegration:
    """Test the CLI entry-point with real resource paths."""

    def test_list_tfs_via_main(self, capsys):
        rc = main(["--list_tfs"])
        assert rc == 0
        out = capsys.readouterr().out
        for tf in ["CTCF", "P53", "MYC"]:
            assert tf in out

    def test_list_models_via_main(self, real_models_dir, capsys):
        rc = main(["--list_models", "--tf", "CTCF",
                    "--models_dir", str(real_models_dir)])
        assert rc == 0
        out = capsys.readouterr().out
        assert "CTCF" in out

    def test_missing_fasta_returns_error(self, real_models_dir):
        rc = main(["-f", "/nonexistent/file.fasta", "--tf", "CTCF",
                    "--models_dir", str(real_models_dir)])
        assert rc == 1

    def test_unknown_tf_returns_error(self, real_fasta_file, real_models_dir):
        rc = main(["-f", str(real_fasta_file), "--tf", "NOTREAL",
                    "--models_dir", str(real_models_dir)])
        assert rc == 1


# =============================================================================
# Edge Cases
# =============================================================================

@pytest.mark.integration
@pytest.mark.skipif(not IMPORTS_AVAILABLE, reason="Required imports not available")
class TestEdgeCases:

    def test_empty_fasta(self, tmp_path):
        fa = tmp_path / "empty.fa"
        fa.write_text("")
        seqs = parse_fasta(fa)
        assert seqs == []

    def test_single_short_sequence(self, tmp_path):
        fa = tmp_path / "short.fa"
        fa.write_text(">tiny\nACGT\n")
        seqs = parse_fasta(fa)
        wins = generate_windows(seqs, frame=301, step=150)
        # Too short for a full window, should still produce one entry
        assert len(wins) == 1
        assert wins[0][1] == "ACGT"

    def test_all_ambiguous_sequence(self, tmp_path):
        fa = tmp_path / "ambig.fa"
        fa.write_text(">ambig\n" + "N" * 500 + "\n")
        seqs = parse_fasta(fa)
        wins = generate_windows(seqs, frame=301, step=150)
        # All windows should be filtered out
        assert len(wins) == 0

    def test_threshold_boundaries(self):
        probs = np.array([0.0, 0.5, 0.500001, 1.0])
        above_half = probs >= 0.5
        assert sum(above_half) == 3
        above_zero = probs >= 0.0
        assert sum(above_zero) == 4
        above_one = probs >= 1.0
        assert sum(above_one) == 1

    def test_reverse_complement_palindrome(self):
        """A sequence that is its own reverse complement."""
        seq = "ACGT"
        assert reverse_complement(seq) == seq


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
