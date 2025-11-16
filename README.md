<p align="center">
  <a href="https://hub.docker.com/r/sheryar09/bisque-deconv">
    <img src="https://img.shields.io/docker/pulls/sheryar09/bisque-deconv.svg" alt="Docker Pulls" />
  </a>
  <a href="https://hub.docker.com/r/sheryar09/bisque-deconv">
    <img src="https://img.shields.io/docker/v/sheryar09/bisque-deconv/latest?label=docker%20image" alt="Docker Image Version" />
  </a>
  <a href="https://github.com/shari01/Bisque-Cell-Type-Deconvolution-of-Bulk-RNA-Seq/actions">
    <img src="https://github.com/shari01/Bisque-Cell-Type-Deconvolution-of-Bulk-RNA-Seq/workflows/Docker%20Image%20CI/badge.svg" alt="CI Status" />
  </a>
</p>

<h1 align="center">bisque-deconv — Bulk RNA-Seq Cell-Type Deconvolution</h1>
<h3 align="center">Python 3.12 + R 4.5 + Bioconductor 3.21</h3>

<p align="center">
  A complete S1 → S2 → S3 Bisque RNA deconvolution workflow using Python (Scanpy, AnnData)
  and R (BisqueRNA via rpy2). Generates cell-type proportions, signature matrices, QC metrics,
  and enhanced tumor microenvironment (TME) reports from bulk RNA-Seq data.
</p>

<hr />

<h2>1. Python Installation</h2>

<p><strong>Recommended: use a virtual environment</strong></p>

<pre><code>python -m venv .venv
. .venv/Scripts/activate   # Windows PowerShell
pip install -e .
</code></pre>

<hr />

<h2>2. Run the Pipeline (CLI)</h2>

<p>Core entry point: <code>python -m bisque_deconv.cli</code></p>

<pre><code>python -m bisque_deconv.cli `
  --ref_h5ad "path/to/reference.h5ad" `
  --bulk "path/to/bulk.csv" `
  --bulk_gene_col "Gene" `
  --metadata "path/to/metadata.csv" `
  --outdir "outputs/S1" `
  --s2_outdir "outputs/S2" `
  --enh_outdir "outputs/S3" `
  --species "human"
</code></pre>

<ul>
  <li><code>--bulk_gene_col</code> must match the gene column in your bulk expression file (e.g. <code>Gene</code>).</li>
  <li>Use <strong>either</strong> <code>--ref_h5ad</code> or <code>--root</code> (reference scRNA-seq input).</li>
</ul>

<hr />

<h2>3. Outputs (S1 → S2 → S3)</h2>

<h3>S1 — Preprocessing (path: <code>--outdir</code>)</h3>
<ul>
  <li><code>sc_counts.tsv</code></li>
  <li><code>sc_metadata.tsv</code></li>
  <li><code>reference_by_celltype.tsv</code></li>
  <li><code>markers.tsv</code></li>
  <li><code>bulk_counts_overlap.tsv</code></li>
  <li><code>bulk_counts_unmapped.tsv</code></li>
  <li><code>Bisque_ready.h5ad</code></li>
</ul>

<h3>S2 — BisqueRNA Deconvolution (path: <code>--s2_outdir</code>)</h3>
<ul>
  <li><code>bisque_bulk_proportions.tsv</code> (bulk-level cell-type fractions)</li>
  <li><code>bisque_sc_proportions.tsv</code> (single-cell level proportions)</li>
  <li><code>signature_matrix.tsv</code> (cell-type-specific expression signatures)</li>
  <li><code>pseudobulk_donor_level.tsv</code></li>
  <li>QC plots and diagnostic summaries</li>
</ul>

<h3>S3 — Enhanced Reporting (path: <code>--enh_outdir</code>)</h3>
<ul>
  <li>Cell-type variability and composition plots</li>
  <li>Summary reports of deconvolution results</li>
  <li>Final graphics to support TME interpretation</li>
</ul>

<hr />

<h2>4. Conceptual Pipeline Overview</h2>

<p>
  The workflow can be summarized as:
</p>

<ol>
  <li><strong>S1 (Preprocessing)</strong>: Harmonize bulk RNA-Seq and reference scRNA-Seq, generate <code>Bisque_ready.h5ad</code>.</li>
  <li><strong>S2 (Deconvolution)</strong>: Run BisqueRNA to estimate cell-type proportions and build a signature matrix.</li>
  <li><strong>S3 (Reporting)</strong>: Produce summarized tables and plots for interpretation of the tumor microenvironment.</li>
</ol>

<p>
  This design allows you to re-run individual stages (e.g., S2 only) without recomputing the entire pipeline.
</p>

<hr />

<h2>5. Docker Usage</h2>

<h3>5.1 Build the image</h3>

<pre><code>docker build -t bisque-deconv:latest .
</code></pre>

<h3>5.2 Run the pipeline using Docker</h3>

<p>
  Mount your current working directory into <code>/workspace</code> inside the container:
</p>

<pre><code>docker run --rm -it ^
  -v "%cd%":/workspace ^
  bisque-deconv:latest ^
  python -m bisque_deconv.cli \
    --ref_h5ad "/workspace/reference.h5ad" \
    --bulk "/workspace/bulk.csv" \
    --bulk_gene_col "Gene" \
    --metadata "/workspace/metadata.csv" \
    --outdir "/workspace/out/S1" \
    --s2_outdir "/workspace/out/S2" \
    --enh_outdir "/workspace/out/S3" \
    --species "human"
</code></pre>

<p>
  Inside the container, all paths are Linux-style (e.g. <code>/workspace/...</code>), even when running on Windows.
</p>

<h3>5.3 Save the Docker image to a local file</h3>

<pre><code>docker save -o bisque-deconv.tar bisque-deconv:latest
</code></pre>

<h3>5.4 Load the image from a local file</h3>

<pre><code>docker load -i bisque-deconv.tar
</code></pre>

<hr />

<h2>6. GitHub Actions — Auto-Build Docker Image</h2>

<p>
  To enable automatic Docker builds and pushes to Docker Hub, create the file:
  <code>.github/workflows/docker-publish.yml</code>
</p>

<pre><code>name: Docker Image CI

on:
  push:
    branches: ["main", "master"]
    tags: ["v*"]
  workflow_dispatch:

jobs:
  build-and-push:
    runs-on: ubuntu-latest

    steps:
      - uses: actions/checkout@v4

      - uses: docker/login-action@v3
        with:
          username: ${{ secrets.DOCKERHUB_USERNAME }}
          password: ${{ secrets.DOCKERHUB_TOKEN }}

      - uses: docker/setup-buildx-action@v3

      - uses: docker/build-push-action@v6
        with:
          context: .
          push: true
          tags: |
            sheryar09/bisque-deconv:latest
            sheryar09/bisque-deconv:${{ github.sha }}
</code></pre>

<p>
  In your GitHub repository settings, define the following secrets:
</p>

<ul>
  <li><code>DOCKERHUB_USERNAME</code></li>
  <li><code>DOCKERHUB_TOKEN</code> (a Docker Hub access token)</li>
</ul>

<hr />

<h2>7. R &amp; BisqueRNA Installation (Standalone)</h2>

<p>
  If you want to run BisqueRNA natively (outside Docker), you can use the following R setup:
</p>

<pre><code>dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages("pkgbuild")
pkgbuild::check_build_tools(debug = TRUE)

install.packages("BiocManager")
BiocManager::install(c("Biobase", "BiocGenerics", "preprocessCore"))

install.packages("remotes")
remotes::install_github("cozygene/bisque", upgrade = "never")

library(BisqueRNA)
packageVersion("BisqueRNA")
</code></pre>

<hr />

<h2>8. Stage-Only Execution</h2>

<p>
  The CLI allows skipping earlier stages if outputs already exist:
</p>

<ul>
  <li><strong>Skip S1 (preprocessing)</strong>, run only S2 + S3:</li>
</ul>

<pre><code>python -m bisque_deconv.cli --skip_s1 \
  --ref_h5ad "reference.h5ad" \
  --bulk "bulk.csv" \
  --bulk_gene_col "Gene" \
  --s2_outdir "out/S2" \
  --enh_outdir "out/S3" \
  --species "human"
</code></pre>

<ul>
  <li><strong>Skip S1 and S2</strong>, run only S3 (assuming S1 and S2 outputs already exist):</li>
</ul>

<pre><code>python -m bisque_deconv.cli --skip_s1 --skip_s2 \
  --enh_outdir "out/S3" \
  --species "human"
</code></pre>

<hr />

<h2>9. Requirements Summary</h2>

<ul>
  <li><strong>Python</strong>: 3.12</li>
  <li><strong>R</strong>: 4.5.1</li>
  <li><strong>Bioconductor</strong>: 3.21</li>
  <li><strong>Core Python dependencies</strong>: <code>numpy</code>, <code>pandas</code>, <code>scipy</code>, <code>scanpy</code>, <code>anndata</code>, <code>scikit-learn</code>, <code>rpy2</code>, <code>matplotlib</code>, <code>seaborn</code>, <code>tqdm</code>, <code>statsmodels</code></li>
</ul>

<hr />

<h2>10. Summary</h2>

<p>
  <strong>bisque-deconv</strong> provides a robust, reproducible, and container-friendly pipeline for
  cell-type deconvolution of bulk RNA-Seq data using BisqueRNA. It integrates modern Python tooling
  (Scanpy, AnnData) with R-based deconvolution, producing publication-ready outputs for TME analysis
  and cell-type composition profiling.
</p>
