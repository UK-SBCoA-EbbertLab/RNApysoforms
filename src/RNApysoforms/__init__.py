# Import necessary functions from local modules
from .shorten_gaps import shorten_gaps  # Function to shorten gaps in data
from .to_intron import to_intron        # Function to convert exons to introns
from .read_ensembl_gtf import read_ensembl_gtf          # Function to read and parse GTF (Gene Transfer Format) files
from .read_expression_matrix import read_expression_matrix # Function to load counts matrix
from .read_expression_matrix import process_expression_matrix # Function to process already-loaded counts matrix
from .gene_filtering import gene_filtering # Function to filter by gene_name 
from .calculate_exon_number import calculate_exon_number ## Function to calculate exon number if missing
from .make_plot import make_plot
from .make_traces import make_traces

# Define the public API of this module by specifying which functions to expose when imported
__all__ = ['shorten_gaps', 'to_intron', "read_ensembl_gtf", "make_traces",
           "read_expression_matrix", "process_expression_matrix", "gene_filtering", "calculate_exon_number", "make_plot"]

__version__ = "1.2.11"

