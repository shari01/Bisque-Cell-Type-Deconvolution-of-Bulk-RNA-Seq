# Python 3.12 base image (Debian-based)
FROM python:3.12-slim

ENV PYTHONUNBUFFERED=1 \
    DEBIAN_FRONTEND=noninteractive

# ---------------------------
# System deps (R + build tools)
# ---------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base r-base-dev \
    build-essential \
    gfortran \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    git \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------
# R: install Bisque and deps
# ---------------------------
# This script installs:
# - BiocManager + Biobase/BiocGenerics/preprocessCore
# - remotes
# - BisqueRNA from GitHub (cozygene/bisque)
RUN echo 'dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)' >> /tmp/install_bisque.R && \
    echo '.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))' >> /tmp/install_bisque.R && \
    echo 'options(repos = c(CRAN = "https://cloud.r-project.org"))' >> /tmp/install_bisque.R && \
    echo 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")' >> /tmp/install_bisque.R && \
    echo 'BiocManager::install(c("Biobase","BiocGenerics","preprocessCore"), ask = FALSE, update = FALSE)' >> /tmp/install_bisque.R && \
    echo 'if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")' >> /tmp/install_bisque.R && \
    echo 'remotes::install_github("cozygene/bisque", upgrade = "never", build_vignettes = FALSE)' >> /tmp/install_bisque.R && \
    Rscript /tmp/install_bisque.R && \
    rm /tmp/install_bisque.R

# ---------------------------
# Python: install bisque-deconv package
# ---------------------------
WORKDIR /app

# Copy project metadata and source
COPY pyproject.toml README.md ./
COPY src/ ./src/

# Install your package (and its Python deps) into the image
RUN pip install --upgrade pip && \
    pip install . && \
    pip cache purge

# Default working dir for user
WORKDIR /workspace

# ---------------------------
# Default entrypoint
# ---------------------------
# You can also use "bisque-deconv" since it's defined as a console_script,
# but python -m is very explicit.
ENTRYPOINT ["python", "-m", "bisque_deconv.cli"]
