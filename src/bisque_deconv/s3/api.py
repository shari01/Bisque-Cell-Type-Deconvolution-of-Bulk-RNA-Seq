#!/usr/bin/env python3
"""
S3 (analysis) — Public wrapper
s3_enhanced_analysis
"""

from __future__ import annotations

from typing import Any, Dict, Optional

from .analyzer import DeconvolutionAnalyzer

__all__ = ["s3_enhanced_analysis"]


def s3_enhanced_analysis(
    input_tsv: str,
    metadata_path: Optional[str] = None,
    sample_col: str = "sample_id",
    condition_col: str = "condition",
    control_label: str = "Control",
    patient_label: str = "Disease",
    out_dir: str = "deconvolution_analysis",
    plots: bool = True,
) -> Dict[str, Any]:
    """
    High-level convenience API for S3 analysis.

    Parameters
    ----------
    input_tsv : str
        Path to `bisque_bulk_proportions.tsv` (or CSV/TSV with cell types in the first column).
    metadata_path : str | None
        Optional table with cohort labels. Must include columns specified by
        `sample_col` and `condition_col`.
    sample_col : str
        Column in metadata that holds sample IDs (matching columns of `input_tsv`).
    condition_col : str
        Column in metadata with cohort labels (e.g., Control vs Disease).
    control_label : str
        Label in `condition_col` that denotes controls.
    patient_label : str
        Label in `condition_col` that denotes patients/cases.
    out_dir : str
        Output directory root for plots/, data/, reports/.
    plots : bool
        If True, attempts to render figures via `s3.plots` (optional module).

    Returns
    -------
    Dict[str, Any]
        Result dictionary from `DeconvolutionAnalyzer.analyze(...)`.
    """
    analyzer = DeconvolutionAnalyzer(output_dir=out_dir, plots_enabled=plots)
    return analyzer.analyze(
        proportions_file=input_tsv,
        metadata_file=metadata_path,
        sample_col=sample_col,
        condition_col=condition_col,
        control_label=control_label,
        patient_label=patient_label,
    )
