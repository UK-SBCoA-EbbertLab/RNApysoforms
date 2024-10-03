# Import necessary functions from local modules
from .shorten_gaps import shorten_gaps  # Function to shorten gaps in data
from .to_intron import to_intron        # Function to convert exons to introns
from .set_axis import set_axis          # Function to set axis properties for plot
from .read_gtf import read_gtf          # Function to read and parse GTF (Gene Transfer Format) files
from .make_traces import make_traces    # Function to generate traces for plotting
from .load_counts_matrix import load_counts_matrix # Function to load counts matrix
from .gene_filtering import gene_filtering # Function to filter by gene_name 

# Define the public API of this module by specifying which functions to expose when imported
__all__ = ['shorten_gaps', 'to_intron', 'set_axis', "read_gtf", "make_traces"
           "load_counts_matrix", "gene_filtering"]

