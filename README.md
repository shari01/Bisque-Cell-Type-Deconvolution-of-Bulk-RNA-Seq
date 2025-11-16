<p align="center">
  <a href="https://hub.docker.com/r/sheryar09/bisque-deconv" target="_blank">
    <img src="https://img.shields.io/docker/pulls/sheryar09/bisque-deconv.svg" alt="Docker Pulls">
  </a>
  <a href="https://hub.docker.com/r/sheryar09/bisque-deconv" target="_blank">
    <img src="https://img.shields.io/docker/v/sheryar09/bisque-deconv/latest?label=docker%20image" alt="Docker Image Version">
  </a>
  <a href="https://github.com/YOUR_GITHUB_USERNAME/bisque-deconv/actions" target="_blank">
    <img src="https://github.com/YOUR_GITHUB_USERNAME/bisque-deconv/workflows/Docker%20Image%20CI/badge.svg" alt="CI Status">
  </a>
</p>

<h1 align="center">bisque-deconv — Bulk RNA-Seq Cell-Type Deconvolution<br>(Python 3.12 + R 4.5 + Bioconductor 3.21)</h1>

<p align="center">
A full S1→S2→S3 Bisque RNA deconvolution pipeline with Python orchestration and R integration (via rpy2).  
Provides high-resolution cellular composition estimates, signature matrices, QC diagnostics, and enhanced reporting.
</p>

<hr/>

<h2>📦 1) Python Installation</h2>

<h3>Virtual environment (recommended)</h3>

<pre>
python -m venv .venv
. .venv/Scripts/activate      # Windows PowerShell
pip install -e .
</pre>

<hr/>

<h2>▶ 2) Run the Pipeline (CLI)</h2>

<pre>
python -m bisque_deconv.cli `
    --ref_h5ad "path/to/reference.h5ad" `
    --bulk "path/to/bulk.csv" `
    --bulk_gene_col "Gene" `
    --metadata "path/to/metadata.csv" `
    --outdir "outputs/S1" `
    --s2_outdir "outputs/S2" `
    --enh_outdir "outputs/S3" `
    --species "human"
</pre>

<ul>
<li>The first column in the bulk dataset must match <code>--bulk_gene_col</code>.</li>
<li>Use <code>--ref_h5ad</code> <strong>or</strong> <code>--root</code>.</li>
</ul>

<hr/>

<h2>📁 3) Outputs (Expected)</h2>

<h3>S1 — Preprocessing ( <code>--outdir</code> )</h3>
<ul>
  <li>sc_counts.tsv</li>
  <li>sc_metadata.tsv</li>
  <li>reference_by_celltype.tsv</li>
  <li>markers.tsv</li>
  <li>bulk_counts_overlap.tsv</li>
  <li>bulk_counts_unmapped.tsv</li>
  <li>Bisque_ready.h5ad</li>
</ul>

<h3>S2 — BisqueRNA Deconvolution ( <code>--s2_outdir</code> )</h3>
<ul>
  <li>bisque_bulk_proportions.tsv</li>
  <li>bisque_sc_proportions.tsv</li>
  <li>signature_matrix.tsv</li>
  <li>pseudobulk_donor_level.tsv</li>
  <li>QC plots</li>
</ul>

<h3>S3 — Enhanced Reporting ( <code>--enh_outdir</code> )</h3>
<ul>
  <li>Cell-type variability plots</li>
  <li>High-level summaries</li>
  <li>TME composition graphics</li>
  <li>Final report</li>
</ul>

<hr/>

<h2>🧬 4) Pipeline Diagram</h2>

```mermaid
flowchart LR
    A[Bulk RNA-seq + Metadata + scRNA-seq Reference] --> B[S1: Preprocessing<br/>(QC, Alignment, Bisque_ready.h5ad)]
    B --> C[S2: BisqueRNA Deconvolution<br/>(Proportions, Signature Matrix)]
    C --> D[S3: Enhanced Reporting<br/>(Plots, Summaries, QC)]
    D --> E[TME Interpretation / Biological Insights]
<hr/> <h2>🐋 5) Docker Usage</h2> <h3>Build Docker image</h3> <pre> docker build -t bisque-deconv:latest . </pre> <h3>Run the pipeline using Docker</h3>
Mount your project directory into <code>/workspace</code>:

<pre> docker run --rm -it ^ -v "%cd%":/workspace ^ bisque-deconv:latest ^ python -m bisque_deconv.cli \ --ref_h5ad "/workspace/reference.h5ad" \ --bulk "/workspace/bulk.csv" \ --bulk_gene_col "Gene" \ --metadata "/workspace/metadata.csv" \ --outdir "/workspace/out/S1" \ --s2_outdir "/workspace/out/S2" \ --enh_outdir "/workspace/out/S3" \ --species "human" </pre> <h3>Save the image to a local .tar file</h3> <pre> docker save -o bisque-deconv.tar bisque-deconv:latest </pre> <h3>Load the image from a .tar file</h3> <pre> docker load -i bisque-deconv.tar </pre> <hr/> <h2>🚀 6) Docker Hub Auto-Build (GitHub Actions)</h2> <p>Create the file: <code>.github/workflows/docker-publish.yml</code></p> <pre> name: Docker Image CI on: push: branches: [ "main", "master" ] tags: [ "v*" ] workflow_dispatch: jobs: build-and-push: runs-on: ubuntu-latest steps: - uses: actions/checkout@v4 - name: Log in to Docker Hub uses: docker/login-action@v3 with: username: ${{ secrets.DOCKERHUB_USERNAME }} password: ${{ secrets.DOCKERHUB_TOKEN }} - uses: docker/setup-buildx-action@v3 - name: Build and push Docker image uses: docker/build-push-action@v6 with: context: . push: true tags: | sheryar09/bisque-deconv:latest sheryar09/bisque-deconv:${{ github.sha }} </pre> <p><b>Note:</b> Add secrets in GitHub → Settings → Secrets → Actions:</p> <ul> <li><code>DOCKERHUB_USERNAME</code></li> <li><code>DOCKERHUB_TOKEN</code></li> </ul> <hr/> <h2>📙 7) Install BisqueRNA (R 4.5 / Bioc 3.21)</h2> <pre> dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE) .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths())) options(repos = c(CRAN = "https://cloud.r-project.org")) if (!requireNamespace("pkgbuild", quietly = TRUE)) install.packages("pkgbuild") pkgbuild::check_build_tools(debug = TRUE) if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager") BiocManager::install(c("Biobase", "BiocGenerics", "preprocessCore"), ask = FALSE) if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes") remotes::install_github("cozygene/bisque", upgrade = "never", build_vignettes = FALSE) library(BisqueRNA) print(packageVersion("BisqueRNA")) </pre> <hr/> <h2>🧩 8) Stage-Only Execution</h2>
Run only S2:

<pre> python -m bisque_deconv.cli --skip_s1 \ --ref_h5ad "reference.h5ad" \ --bulk "bulk.csv" \ --bulk_gene_col "Gene" \ --s2_outdir "out/S2" \ --species human </pre>
Run only S3:

<pre> python -m bisque_deconv.cli --skip_s1 --skip_s2 ... </pre> <hr/> <h2>📌 9) Requirements</h2> <ul> <li>Python 3.12</li> <li>R 4.5.1</li> <li>Bioconductor 3.21</li> <li>Dependencies: anndata, scanpy, numpy, pandas, rpy2, scikit-learn</li> </ul> <hr/> <h2>✔ Summary</h2> <p> <b>bisque-deconv</b> provides a robust, reproducible, and fully containerized workflow for cell-type deconvolution of bulk RNA-Seq data, integrating modern Python tooling with the BisqueRNA R package to produce high-quality TME profiles and biological insights. </p> ```
