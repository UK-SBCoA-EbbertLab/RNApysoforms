# RNApysoforms <img src="./assets/RNA-pysoforms-logo.svg" align="right" height="80" />


<!-- badges: start -->
[![Run Tests](https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/actions/workflows/main.yml/badge.svg)](https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/actions/workflows/main.yml)
[![Codecov test coverage](https://codecov.io/gh/UK-SBCoA-EbbertLab/RNApysoforms/branch/main/graph/badge.svg)](https://app.codecov.io/gh/UK-SBCoA-EbbertLab/RNApysoforms?branch=main)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![GitHub issues](https://img.shields.io/github/issues/UK-SBCoA-EbbertLab/RNApysoforms)](https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/issues)
[![GitHub pulls](https://img.shields.io/github/issues-pr/UK-SBCoA-EbbertLab/RNApysoforms)](https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/pulls)
[![Documentation Status](https://readthedocs.org/projects/rna-pysoforms/badge/?version=latest)](https://rna-pysoforms.readthedocs.io/en/latest/?badge=latest)

<!-- badges: end -->



`RNApysoforms` is a Python package that provides Plotly-based implementations of functionality similar to the R package `ggtranscript`. It's designed to help visualize genomic transcript structures using Plotly.

## Features

- `calculate_exon_number()`: Assigns exon numbers to exons, CDS, and introns within a genomic annotation dataset based on transcript structure and strand direction.

- `gene_filtering()`: Filters genomic annotations and optionally an expression matrix for a specific gene, with options to order and select top expressed transcripts.

- `make_plot()`: Creates a multi-panel Plotly figure combining transcript structure plots and expression data plots.

- `make_traces()`: Generates Plotly traces for visualizing transcript structures and expression data.

- `read_expression_matrix()`: Loads and processes an expression matrix, optionally merging with metadata, performing CPM normalization, and calculating relative transcript abundance.

- `read_ensembl_gtf()`: Reads a GTF (Gene Transfer Format) file and returns the data as a Polars DataFrame.

- `shorten_gaps()`: Shortens intron and transcript start gaps between exons in genomic annotations to enhance visualization.

- `to_intron()`: Converts exon coordinates into corresponding intron coordinates within a genomic annotation dataset.

## Installation

You can install `RNApysoforms` using pip:

```bash
pip install RNApysoforms
```

## Usage

Here's a basic example of how to use `RNApysoforms`:

```python

```

## Contributing

Contributions to `RNApysoforms` are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License.
