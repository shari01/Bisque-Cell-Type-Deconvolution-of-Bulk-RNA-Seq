# bisque-deconv — Quick Windows Guide (Python 3.12 + R 4.5 / Bioc 3.21)

This guide shows **exact steps** to install and run the pipeline on Windows with your paths and settings. It includes a Bisque installation recipe for **R 4.5** (Bioconductor **3.21**).

---

## 1) Install (Python)

```powershell
# Create and activate a virtual environment
python -m venv .venv
. .\.venv\Scripts\Activate.ps1

# Install the package (editable)
pip install -e .

python -m bisque_deconv.cli `
  --ref_h5ad "enter_your_Dataset_location\thymus.h5ad" `
  --bulk     "enter_your_Dataset_location" `
  --bulk_gene_col "Gene" `
  --metadata "enter_your_Dataset_location" `
  --outdir   "output_location\Bisque-Preprocessing" `
  --s2_outdir "output_location\Bisque_deconvolution_results" `
  --enh_outdir "output_location\Enhanced-Deconv-Reports" `
  --species  "human"
```
```

> Ensure Microsoft **Visual C++** Redistributable / Build Tools are present (they’re required by several dependencies).

---

## 2) Run (CLI — copy/paste)

Use **PowerShell** (backtick ` for line continuation). Adjust paths if needed.

```powershell
python -m bisque_deconv.cli `
  --ref_h5ad "enter_your_Dataset_location\thymus.h5ad" `
  --bulk     "enter_your_Dataset_location" `
  --bulk_gene_col "Gene" `
  --metadata "enter_your_Dataset_location" `
  --outdir   "enter_your_Dataset_location\thymus-test\Bisque-Preprocessing" `
  --s2_outdir "enter_your_Dataset_location\thymus-test\Bisque_deconvolution_results" `
  --enh_outdir "enter_your_Dataset_location\thymus-test\Enhanced-Deconv-Reports" `
  --species  "human"
```

**Notes**
- Provide **either** `--ref_h5ad` **or** `--root`. The above uses `--ref_h5ad`.
- `--bulk_gene_col` must match the column in your bulk CSV (here: **Gene**).

---

## 3) Results (expected outputs)

- **S1 → `--outdir`**: `sc_counts.tsv`, `sc_metadata.tsv`, `reference_by_celltype.tsv`, `markers.tsv`, `bulk_counts_overlap.tsv`, `bulk_counts_unmapped.tsv`, `Bisque_ready.h5ad`
- **S2 → `--s2_outdir`**: `bisque_bulk_proportions.tsv`, `bisque_sc_proportions.tsv`, `signature_matrix.tsv`, `pseudobulk_donor_level.tsv`, plus QC plots
- **S3 → `--enh_outdir`**: summary/plots/report

---

## 4) Status on data sizes

- **Tested on a smaller dataset:** Completed end‑to‑end successfully.
- **Larger datasets:** Memory allocation issues may appear on Windows. Recommendations:
  - Subset the reference `.h5ad` to fewer cells / HVGs (e.g., top 2k genes).
  - Filter to tissues of interest via `--tissue_include` and relevant `--tissue_keys`.
  - Run stages separately with `--skip_s1/--skip_s2/--skip_s3` to avoid recomputation.
  - Increase system paging file (virtual memory) if possible.

---

## 5) Bisque (R) installation — R 4.5 / Bioc 3.21

Run the following **in R**. It installs BisqueRNA and its common dependencies in your **user library** (no admin).

```r
# ---- 0) Use your user library (no admin needed) ----
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

# ---- 1) Build tools check (R 4.5 needs Rtools45) ----
if (!requireNamespace("pkgbuild", quietly = TRUE)) install.packages("pkgbuild")
pkgbuild::check_build_tools(debug = TRUE)  # should say TRUE
print(Sys.which("make"))                   # should show a path, not ""

# ---- 2) Get Bioconductor deps that BisqueRNA uses ----
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# Core deps commonly needed by BisqueRNA (safe to install)
BiocManager::install(c("Biobase", "BiocGenerics", "preprocessCore"), ask = FALSE, update = FALSE)

# ---- 3) Install from GitHub (official repo) ----
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("cozygene/bisque", upgrade = "never", build_vignettes = FALSE)

# ---- 4) Smoke test ----
library(BisqueRNA)
print(packageVersion("BisqueRNA"))
```

**If `rpy2` can’t find R**: set `R_HOME` (e.g., `C:\Program Files\R\R-4.5.1`) in Windows Environment Variables and reopen your terminal.

---

## 6) Optional: stage-only runs

If you’ve already produced S1 outputs:
```powershell
python -m bisque_deconv.cli `
  --ref_h5ad "D:\...\thymus.h5ad" `
  --bulk "C:\...\count.csv" `
  --bulk_gene_col "Gene" `
  --meta 'metadata-location'
  --outdir "D:\...\Bisque-Preprocessing" `
  --s2_outdir "D:\...\Bisque_deconvolution_results" `
  --enh_outdir "D:\...\Enhanced-Deconv-Reports" `
  --species human `
  --skip_s1
```
Use `--skip_s2` or `--skip_s3` similarly when re-running specific parts.

---

## 7) Versions used

- Python 3.12 (venv)  
- R 4.5.1 / Bioconductor 3.21  
- Key Python libs: `anndata>=0.12`, `scanpy>=1.9`, `numpy>=1.24`, `pandas>=2.0`, `rpy2>=3.5`, `scikit-learn`, `scrublet` (optional on Windows/Py≥3.13)

---

**That’s it.** Install with `pip install -e .`, then run the provided CLI block.
