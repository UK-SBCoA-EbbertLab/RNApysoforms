# RNApysoforms <img src="./assets/RNA-pysoforms-logo.svg" align="right" height="100" />


<!-- badges: start -->
[![Run Tests](https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/actions/workflows/main.yml/badge.svg)](https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/actions/workflows/main.yml)
[![GitHub issues](https://img.shields.io/github/issues/UK-SBCoA-EbbertLab/RNApysoforms)](https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/issues)
[![GitHub pulls](https://img.shields.io/github/issues-pr/UK-SBCoA-EbbertLab/RNApysoforms)](https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/pulls)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check-bioc](https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/actions/workflows/R-CMD-check-bioc.yml/badge.svg)](https://github.com/UK-SBCoA-EbbertLab/RNApysoforms/actions)
[![Codecov test coverage](https://codecov.io/gh/UK-SBCoA-EbbertLab/RNApysoforms/branch/main/graph/badge.svg)](https://app.codecov.io/gh/UK-SBCoA-EbbertLab/RNApysoforms?branch=main)

<!-- badges: end -->



`RNApysoforms` is a Python package that provides Plotly-based implementations of functionality similar to the R package `ggtranscript`. It's designed to help visualize genomic transcript structures using Plotly.

## Features

- `geom_range`: Visualize exons or other range-based genomic features
- `geom_intron`: Draw introns with strand arrows
- `shorten_gaps`: Improve transcript structure visualization by shortening gaps
- `to_intron`: Convert exon coordinates to intron coordinates


## Installation

You can install `RNApysoforms` directly from the source:

```bash
pip install -e /path/to/RNApysoforms
```

## Usage

Here's a basic example of how to use `RNApysoforms`:

```python
import pandas as pd
from RNApysoforms import geom_range, geom_intron, to_intron

# Prepare your data
exons = pd.DataFrame({
    'start': [100, 300, 600],
    'end': [200, 400, 700],
    'transcript_name': ['tx1', 'tx1', 'tx1'],
    'strand': ['+', '+', '+']
})

# Convert exons to introns
introns = to_intron(exons, group_var='transcript_name')

# Create a plot
fig = go.Figure()

# Add exons
exon_traces = geom_range(exons, x_start='start', x_end='end', y='transcript_name')
for trace in exon_traces:
    fig.add_trace(trace)

# Add introns
intron_traces = geom_intron(introns, x_start='start', x_end='end', y='transcript_name', strand='strand')
for trace in intron_traces:
    fig.add_trace(trace)

# Show the plot
fig.show()
```

## Contributing

Contributions to `RNApysoforms` are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License.
