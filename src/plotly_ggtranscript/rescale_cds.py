import pandas as pd

def rescale_cds(cds_exon_diff, gene_rescaled_exons):
    """
    Rescales CDS regions based on exon positions and the calculated differences between CDS and exon positions.
    This function aligns the CDS regions to the rescaled exon positions and adjusts their start and end points
    accordingly.

    Parameters
    ----------
    cds_exon_diff : pd.DataFrame
        DataFrame containing the differences between the start and end positions of exons and CDS regions.
        Expected columns:
        - 'diff_start': Absolute difference between exon start and CDS start positions.
        - 'diff_end': Absolute difference between exon end and CDS end positions.
        - Any columns necessary for joining (e.g., 'transcript_id', 'gene_id').
    gene_rescaled_exons : pd.DataFrame
        DataFrame containing the rescaled exon positions.
        Expected columns:
        - 'start': Rescaled start position of the exon.
        - 'end': Rescaled end position of the exon.
        - Any columns necessary for joining.

    Returns
    -------
    gene_rescaled_cds : pd.DataFrame
        DataFrame containing the rescaled CDS positions, with adjusted 'start' and 'end' positions,
        and 'transcript_id' ordered according to 'factor_order'.
    """

    # Step 1: Prepare the CDS DataFrame
    # - Assign a new column 'type' with the value "CDS"
    # - Drop unnecessary columns: 'c_start', 'c_end', 'e_start', 'e_end' if they exist

    cds_prepared = (
        cds_exon_diff
        .assign(type="CDS")
        .drop(columns=['cds_start', 'cds_end', 'exon_start', 'exon_end'], errors='ignore')
    )

    # Step 2: Prepare the Exon DataFrame
    # - Rename 'start' to 'exon_start' and 'end' to 'exon_end' to avoid column name conflicts
    # - Drop the 'type' column if it exists

    exons_prepared = (
        gene_rescaled_exons
        .rename(columns={'start': 'exon_start', 'end': 'exon_end'})
        .drop(columns=['type'], errors='ignore')
    )

    # Step 3: Identify common columns for joining
    # - Find columns that are present in both DataFrames to use as keys for the join

    common_columns = [col for col in cds_prepared.columns if col in exons_prepared.columns]
    if not common_columns:
        raise ValueError("No common columns to perform join on. Ensure both DataFrames have common keys.")

    # Step 4: Perform the left join on the common columns
    # - This aligns the CDS data with the corresponding rescaled exon positions

    gene_rescaled_cds = pd.merge(
        cds_prepared,
        exons_prepared,
        how='left',
        on=common_columns
    )

    # Step 5: Calculate the adjusted 'start' and 'end' positions for the CDS regions
    # - 'start' is adjusted by adding 'diff_start' to 'exon_start'
    # - 'end' is adjusted by subtracting 'diff_end' from 'exon_end'

    gene_rescaled_cds['start'] = gene_rescaled_cds['exon_start'] + gene_rescaled_cds['diff_start']
    gene_rescaled_cds['end'] = gene_rescaled_cds['exon_end'] - gene_rescaled_cds['diff_end']

    # Step 6: Drop unnecessary columns used for calculations
    # - Remove 'exon_start', 'exon_end', 'diff_start', 'diff_end' as they are no longer needed

    gene_rescaled_cds = gene_rescaled_cds.drop(
        columns=['exon_start', 'exon_end', 'diff_start', 'diff_end'],
        errors='ignore'
    )

    return gene_rescaled_cds
