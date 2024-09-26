# plotly_ggtranscript

`plotly_ggtranscript` is a Python package that provides Plotly-based implementations of functionality similar to the R package `ggtranscript`. It's designed to help visualize genomic transcript structures using Plotly.

## Features

- `geom_range`: Visualize exons or other range-based genomic features
- `geom_intron`: Draw introns with strand arrows
- `shorten_gaps`: Improve transcript structure visualization by shortening gaps
- `to_intron`: Convert exon coordinates to intron coordinates

## Installation

You can install `plotly_ggtranscript` directly from the source:

```bash
pip install -e /path/to/plotly_ggtranscript
```

## Usage

Here's a basic example of how to use `plotly_ggtranscript`:

```python
import pandas as pd
from plotly_ggtranscript import geom_range, geom_intron, to_intron

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

Contributions to `plotly_ggtranscript` are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License.