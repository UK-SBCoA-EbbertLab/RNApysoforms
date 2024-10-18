RNApysoforms documentation
===========================

Welcome to the documentation for the RNApysoforms package. Below are the key functions included in this library.


Functions Overview
====================

.. autosummary::
   :toctree: _autosummary
   :recursive:

   RNApysoforms.calculate_exon_number
   RNApysoforms.gene_filtering
   RNApysoforms.make_plot
   RNApysoforms.make_traces
   RNApysoforms.read_expression_matrix
   RNApysoforms.read_ensembl_gtf
   RNApysoforms.shorten_gaps
   RNApysoforms.to_intron


Example Vignettes
================

.. toctree::
   :maxdepth: 2
   :caption: Examples:
   :recursive:

   0.basic_usage.html
   01.rescaled_introns.html
   02.expression_plot.html
   03.expression_plot_with_metadata.html
   04.expression_plot_filtered_and_ordered.html
   05.plot_specific_transcripts.html
   06.custom_color_palettes.html
   07.custom_color_maps.html
   08.separate_CDS_interactivity.html
   09.autoscale_plots.html
   10.dealing_with_different_gtf_files.html
   11.making_expression_plot_only.html
   12.making_changes_to_figure_after_rendering.html