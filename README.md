# RNApysoforms <img src="./assets/RNA-pysoforms-logo.png" align="right" height="80" />


<!-- badges: start -->
[![Run Tests](https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/actions/workflows/main.yml/badge.svg)](https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/actions/workflows/main.yml)
[![Codecov test coverage](https://codecov.io/gh/UK-SBCoA-EbbertLab/RNApysoforms/branch/main/graph/badge.svg)](https://app.codecov.io/gh/UK-SBCoA-EbbertLab/RNApysoforms?branch=main)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![GitHub issues](https://img.shields.io/github/issues/UK-SBCoA-EbbertLab/RNApysoforms)](https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/issues)
[![GitHub pulls](https://img.shields.io/github/issues-pr/UK-SBCoA-EbbertLab/RNApysoforms)](https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/pulls)
[![Documentation Status](https://readthedocs.org/projects/rna-pysoforms/badge/?version=latest)](https://rna-pysoforms.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/863499343.svg)](https://doi.org/10.5281/zenodo.14052788)

<!-- badges: end -->


`RNApysoforms` is a Python package designed for visualizing RNA isoform structures and expression levels. Leveraging Plotly for interactive plotting and Polars for efficient data manipulation, it enables the creation of fast-rendering, interactive plots suitable for both local and web applications. Inspired by the R package [ggtranscript](https://github.com/dzhang32/ggtranscript), RNApysoforms brings similar RNA visualization capabilities to the Python ecosystem, facilitating effective exploration and presentation of RNA sequencing data.

## Pre-print for citation

https://www.biorxiv.org/content/10.1101/2024.11.06.622357v1


## Important

### RNApysoforms expects feature start and end coordinates in GTF format, where coordinates are 1-indexed and inclusive on both ends.


## Installation

You can install `RNApysoforms` using pip:

```bash
pip install RNApysoforms
```


## Quick Start

[Basic usage (quick start)](https://rna-pysoforms.readthedocs.io/en/latest/examples/basic_usage.html)


##  More vignettes (usage examples)

[Rescaling introns for a prettier RNA isoform structure plot](https://rna-pysoforms.readthedocs.io/en/latest/examples/rescaled_introns.html)

[Plotting RNA isoform structure and expression](https://rna-pysoforms.readthedocs.io/en/latest/examples/expression_plot.html)

[Plotting RNA isoform structure and normalized expression](https://rna-pysoforms.readthedocs.io/en/latest/examples/expression_plot_filtered_and_ordered.html)



## Test data and documentation

[Download small test dataset](https://zenodo.org/records/13961009/files/RNApysoforms_test_data.zip?download=1)

[Function documentation and vignettes](https://rna-pysoforms.readthedocs.io/en/latest/index.html)


## Issues

Please go through the [documentation and vignettes](https://rna-pysoforms.readthedocs.io/en/latest/index.html) before submitting an issue.


## Contributing

Contributions to `RNApysoforms` are welcome! Please feel free to submit a Pull Request.

The function implementations are under the `src/RNApysoforms` directory.


## Functions

- [calculate_exon_number()](https://rna-pysoforms.readthedocs.io/en/latest/_autosummary/RNApysoforms.calculate_exon_number.html): Assigns exon numbers to exons, CDS, and introns within a genomic annotation dataset based on transcript structure and strand direction.

- [gene_filtering()](https://rna-pysoforms.readthedocs.io/en/latest/_autosummary/RNApysoforms.gene_filtering.html): Filters genomic annotations and optionally an expression matrix for a specific gene, with options to order and select top expressed transcripts.

- [make_plot()](https://rna-pysoforms.readthedocs.io/en/latest/_autosummary/RNApysoforms.make_plot.html): Creates a Plotly figure panel for transcript structure plots and/or expression data plots.

- [make_traces()](https://rna-pysoforms.readthedocs.io/en/latest/_autosummary/RNApysoforms.make_traces.html): Generates Plotly traces for visualizing transcript structures and expression data.

- [read_expression_matrix()](https://rna-pysoforms.readthedocs.io/en/latest/_autosummary/RNApysoforms.read_expression_matrix.html): Loads and processes an expression matrix, optionally merging with metadata, performing CPM normalization, and calculating relative transcript abundance.

- [read_ensembl_gtf()](https://rna-pysoforms.readthedocs.io/en/latest/_autosummary/RNApysoforms.read_ensembl_gtf.html): Reads an ENSEMBL GTF (Gene Transfer Format) file and returns the data as a Polars DataFrame.

- [shorten_gaps()](https://rna-pysoforms.readthedocs.io/en/latest/_autosummary/RNApysoforms.shorten_gaps.html): Shortens intron and transcript start gaps between exons in genomic annotations to enhance visualization.

- [to_intron()](https://rna-pysoforms.readthedocs.io/en/latest/_autosummary/RNApysoforms.to_intron.html): Converts exon coordinates into corresponding intron coordinates within a genomic annotation dataset.


## License

This project is licensed under the MIT License.
