# Chromnitron

This repository contains the official implementation of **Chromnitron** and its associated data-processing pipeline, **Chrom2vec**. For early testing, only  the inference pipeline for Chromnitron and ATAC-seq processing component for Chrom2Vec are included. Please stay tuned for upcoming full release.

Read our pre-print on [biorxiv](https://www.biorxiv.org/content/10.1101/2025.08.17.670761v1.article-metrics).

## Quick Start
For those looking to get started with Chromnitron quickly, we recommend the following resources:
* **Google Colab:** Explore Chromnitron's features in an interactive environment [Colab](https://colab.research.google.com/drive/1NEgKCjJ4ipGMrtq8MuTmWy4LcahGswjv?usp=sharing).
* **Local Inference:** For local setup, refer to the detailed instructions in the [Chromnitron README](chromnitron).
* **Custom Analysis:** If you want to use Chromnitron with your own ATAC-seq data, please refer to the [Chrom2Vec README](chrom2vec).

## Repository Overview

This repository is organized into two sub-repositories:

* **`chromnitron/`**: This directory holds the core Chromnitron model, designed for predicting chromatin-associated protein binding.
* **`chrom2vec/`**: This directory contains the Chrom2vec pipeline, essential for processing raw sequencing data into a format suitable for Chromnitron.

For more detailed information on each component, please see the respective README files under each subrepo:

* `chromnitron/README.md`
* `chrom2vec/README.md`
