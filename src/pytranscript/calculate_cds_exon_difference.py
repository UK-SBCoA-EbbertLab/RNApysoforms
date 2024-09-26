import pandas as pd

def calculate_cds_exon_difference(gene_cds_regions, gene_exons):
    """
    Calculates the absolute differences between the start and end positions of exons and CDS regions.
    This function is used to prepare data for re-scaling CDS regions based on exon positions.

    Parameters:
    ----------
    gene_cds_regions : pd.DataFrame
        DataFrame containing CDS (Coding DNA Sequence) regions with at least 'start' and 'end' columns,
        along with any common columns used for joining (e.g., 'transcript_id', 'gene_id').
    gene_exons : pd.DataFrame
        DataFrame containing exon regions with at least 'start' and 'end' columns,
        along with any common columns used for joining.

    Returns:
    -------
    cds_exon_diff : pd.DataFrame
        DataFrame resulting from a left join of the CDS and exon data, including the calculated
        absolute differences 'diff_start' and 'diff_end' between exon and CDS start and end positions.
    """

    # Step 1: Rename 'start' and 'end' columns in CDS regions to 'cds_start' and 'cds_end'
    cds_regions = gene_cds_regions.rename(columns={'start': 'cds_start', 'end': 'cds_end'})

    # Remove the 'type' column if it exists
    if 'type' in cds_regions.columns:
        cds_regions = cds_regions.drop(columns=['type'])

    # Step 2: Rename 'start' and 'end' columns in exon regions to 'exon_start' and 'exon_end'
    exons = gene_exons.rename(columns={'start': 'exon_start', 'end': 'exon_end'})

    # Remove the 'type' column if it exists
    if 'type' in exons.columns:
        exons = exons.drop(columns=['type'])

    # Step 3: Identify common columns to perform the left join
    common_columns = list(set(cds_regions.columns) & set(exons.columns))
    if not common_columns:
        raise ValueError("No common columns to perform join on. Ensure both DataFrames have common keys.")

    # Step 4: Perform the left join on the common columns
    cds_exon_diff = pd.merge(cds_regions, exons, how='left', on=common_columns)

    print(cds_exon_diff)

    # Step 5: Calculate absolute differences between exon and CDS start positions
    cds_exon_diff['diff_start'] = (cds_exon_diff['exon_start'] - cds_exon_diff['cds_start']).abs()

    # Step 6: Calculate absolute differences between exon and CDS end positions
    cds_exon_diff['diff_end'] = (cds_exon_diff['exon_end'] - cds_exon_diff['cds_end']).abs()

    return cds_exon_diff
