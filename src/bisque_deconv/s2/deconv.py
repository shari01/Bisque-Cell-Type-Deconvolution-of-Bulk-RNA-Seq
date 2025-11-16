# src/bisque_deconv/s2/deconv.py
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path
import pandas as pd


def s2_bisque_deconv(
    *,
    ref_dir: str,
    out_dir: str,
    drop_mito: bool = True,
    allow_duplicate_genes: bool = True,
    warn_min_overlap: int = 200,
    max_prop_deviation: float = 0.05,
    top_markers_per_ct: int = 25,
    top_genes_per_sample: int = 40,
    auto_install: bool = True,
) -> str:
    """
    Stage 2: Run BisqueRNA (via R) and normalize outputs for Stage 3.

    Returns
    -------
    str
        Path to bisque_bulk_proportions.tsv
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    print(f"[S2] Using reference dir: {ref_dir}")
    print(f"[S2] Writing to:         {out_dir}")

    # --- Locate the R script ---
    r_script = Path(__file__).with_name("r_code.R")
    if not r_script.exists():
        raise FileNotFoundError(f"[S2] Missing R script: {r_script}")

    # --- Environment for R ---
    env = os.environ.copy()
    env.update({
        "REF_DIR_PY": os.path.abspath(ref_dir),
        "OUT_DIR_PY": os.path.abspath(out_dir),
        "DROP_MITO_PY": "true" if drop_mito else "false",
        "ALLOW_DUP_PY": "true" if allow_duplicate_genes else "false",
        "WARN_MIN_OVERLAP_PY": str(int(warn_min_overlap)),
        "MAX_PROP_DEV_PY": str(float(max_prop_deviation)),
        "TOP_MARKERS_PER_CT_PY": str(int(top_markers_per_ct)),
        "TOP_GENES_PER_SAMPLE_PY": str(int(top_genes_per_sample)),
        "AUTO_INSTALL_PY": "true" if auto_install else "false",
    })

    # If overlap TSV exists, pass it explicitly (optional for R; we read it inside R anyway)
    bulk_overlap = Path(ref_dir) / "bulk_counts_overlap.tsv"
    if bulk_overlap.exists():
        env["BULK_FILE_PY"] = str(bulk_overlap.resolve())

    # --- Log files ---
    log_file = out / "R_console.log"
    err_file = out / "R_stderr.log"

    # --- Run R script ---
    cmd = ["Rscript", str(r_script)]
    print(f"[S2] Running: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            check=True,
            env=env,
            capture_output=True,
            text=True,
        )
        log_file.write_text(result.stdout or "")
        err_file.write_text(result.stderr or "")
        print("[S2] ✓ R completed successfully.")
    except FileNotFoundError:
        raise RuntimeError(
            "Rscript not found on PATH. Install R and ensure 'Rscript' is available."
        )
    except subprocess.CalledProcessError as e:
        # Always save R output
        try:
            log_file.write_text(e.stdout or "")
            err_file.write_text(e.stderr or "")
        except Exception:
            pass
        # Also echo last stderr lines to console for faster debugging
        tail = (e.stderr or "").splitlines()[-20:]
        print("\n[S2][R stderr tail]\n" + "\n".join(tail))
        raise RuntimeError(
            f"S2 R step failed (exit {e.returncode}). "
            f"See logs:\n  {log_file}\n  {err_file}"
        )

    # --- Normalize outputs from R ---
    # Prefer the TSV (canonical), fallback to CSV if present
    tsv = out / "bisque_bulk_proportions.tsv"
    csv = out / "Bisque_proportions.csv"

    if tsv.exists():
        df = pd.read_csv(tsv, sep="\t", index_col=0)
    elif csv.exists():
        df = pd.read_csv(csv, index_col=0)
    else:
        raise FileNotFoundError(
            "[S2] Expected R output not found: neither "
            f"{tsv} nor {csv} exist."
        )

    # Heuristic: ensure rows are cell types, columns are samples
    # If looks like samples x celltypes, transpose
    if df.shape[0] < df.shape[1]:
        df = df.T

    bulk_props_tsv = out / "bisque_bulk_proportions.tsv"
    df.to_csv(bulk_props_tsv, sep="\t")

    # Write placeholders for downstream (if not already created by R)
    (out / "signature_matrix.tsv").touch(exist_ok=True)
    (out / "pseudobulk_donor_level.tsv").touch(exist_ok=True)

    # --- Write a JSON summary ---
    summary = {
        "ref_dir": str(Path(ref_dir).resolve()),
        "out_dir": str(out.resolve()),
        "files": {
            "bisque_bulk_proportions": str(bulk_props_tsv),
            "bisque_result_rds": str(out / "Bisque_result.rds"),
            "bisque_summary_json": str(out / "bisque_summary.json"),
            "r_console_log": str(log_file),
            "r_stderr_log": str(err_file),
        },
        "shape_bisque_bulk_proportions": list(map(int, df.shape)),
        "index_is_cell_type": True,
    }

    with open(out / "bisque_python_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[S2] Wrote {bulk_props_tsv}")
    return str(bulk_props_tsv)
