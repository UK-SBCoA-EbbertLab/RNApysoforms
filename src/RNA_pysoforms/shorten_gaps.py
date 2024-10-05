import plotly.graph_objects as go
from typing import List, Union
from RNA_pysoforms.to_intron import to_intron
import polars as pl
from RNA_pysoforms.utils import check_df

def shorten_gaps(
    annotation: pl.DataFrame, 
    group_var: str = "transcript_id", 
    target_gap_width: int = 100
) -> pl.DataFrame:
    
    """
    Shortens intron gaps between exons to create a more compact transcript visualization.

    This function processes genomic annotations, rescaling the width of intron gaps to a target size, while preserving 
    exon and CDS regions. The goal is to make transcript visualizations clearer by reducing the visual space occupied 
    by long intron regions, while keeping the relative transcript structure intact.

    Parameters
    ----------
    annotation : pl.DataFrame
        A Polars DataFrame containing exon data and optionally CDS and/or intron data, 
        with columns like 'start', 'end', 'seqnames', 'type',  and 'strand'.
    group_var : str, optional
        Column used to group transcripts, by default "transcript_id".
    target_gap_width : int, optional
        Maximum allowed width for intron gaps in the rescaled visualization, by default 100.

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame with shortened intron gaps and rescaled coordinates for exons, introns, and CDS regions.

    Raises
    ------
    ValueError
        If the required columns ('start', 'end', 'type', 'strand', 'seqnames', and the group_var) are missing in the input DataFrame.

    Examples
    --------
    Shorten intron gaps in a genomic annotation DataFrame:

    >>> import polars as pl
    >>> from RNA_pysoforms.plot import shorten_gaps
    >>> df = pl.DataFrame({
    ...     "transcript_id": ["tx1", "tx1", "tx1"],
    ...     "start": [100, 200, 500],
    ...     "end": [150, 250, 600],
    ...     "type": ["exon", "exon", "exon"],
    ...     "strand": ["+", "+", "+"],
    ...     "seqnames": ["chr1", "chr1", "chr1"]
    ... })
    >>> shortened_df = shorten_gaps(df, group_var="transcript_id", target_gap_width=50)
    >>> print(shortened_df.head())

    This will return a DataFrame where the intron gaps have been shortened to a maximum width of 50.

    Notes
    -----
    - The function ensures that exon and CDS regions remain unchanged, while intron gaps are shortened to a defined width.
    - The function takes an input dataframe with or without intron entries. If there are no intron entries the function generates them,
       else the function uses the intron entries already provided.
    - The input dataframe must contain columns 'start', 'end', 'type', 'strand', 'seqnames', and the group_var.
    - The function can handle gaps at the start of transcripts and rescale the entire transcript structure for better visualization.
    - Rescaling is applied to the entire transcript structure after the gaps have been shortened.
    """
    
    # Validate the input DataFrame to ensure required columns are present
    check_df(annotation, ["start", "end", "type", "strand", "seqnames", group_var])


    # Check if there are intron entries in the annotation data
    if "intron" in annotation["type"].unique().to_list():
        introns = annotation.filter(pl.col("type") == "intron") # Separate intron data
    else:
        annotation = to_intron(annotation=annotation, group_var=group_var) # Add intron annotations
        introns = annotation.filter(pl.col("type") == "intron") # Separate intron data
    

    # Check if there are CDS entries in the annotation data
    if "CDS" in annotation["type"].unique().to_list():
        cds = annotation.filter(pl.col("type") == "CDS")  # Separate CDS data
    else:
        cds = None  # No CDS entries in the data

    # Separate exons from the rest of the annotation data
    exons = annotation.filter(pl.col("type") == "exon")

    # Ensure the 'type' column in exons and introns is set correctly
    exons = _get_type(exons, "exons")  # Mark the type as 'exon'
    introns = _get_type(introns, "introns")  # Mark the type as 'intron'

    # Adjust intron positions to avoid overlap with exons
    introns = introns.with_columns([
        pl.col("start") + 1,
        pl.col("end") - 1
    ])

    # Identify gaps between exons
    gaps = _get_gaps(exons)  # Gaps between exons within the same chromosome and strand

    # Map gaps to introns
    gap_map = _get_gap_map(introns, gaps)

    # Shorten gaps based on the target gap width
    introns_shortened = _get_shortened_gaps(introns, gaps, gap_map, group_var, target_gap_width)

    # Handle gaps at the start of transcripts if a grouping variable is provided
    if group_var:
        tx_start_gaps = _get_tx_start_gaps(exons, group_var)  # Gaps at the start of transcripts
        gap_map_tx_start = _get_gap_map(tx_start_gaps, gaps)
        tx_start_gaps_shortened = _get_shortened_gaps(tx_start_gaps, gaps, gap_map_tx_start, group_var, target_gap_width)
        # Remove unnecessary columns from the transcript start gaps DataFrame
        tx_start_gaps_shortened = tx_start_gaps_shortened.drop(['start', 'end', 'strand', 'seqnames'])
    else:
        tx_start_gaps_shortened = None  # No gaps at transcript start if no group_var provided

    # Rescale the coordinates of exons and introns after shortening the gaps
    rescaled_tx = _get_rescaled_txs(exons, introns_shortened, tx_start_gaps_shortened, group_var)

    # Process CDS if available
    if isinstance(cds, pl.DataFrame):
        cds_diff = _get_cds_exon_difference(exons, cds)
        rescaled_cds = _get_rescale_cds(cds_diff, rescaled_tx.filter(pl.col("type") == "exon"))
        rescaled_cds = rescaled_cds[rescaled_tx.columns]
        # Combine the rescaled CDS data into the final DataFrame
        rescaled_tx = pl.concat([rescaled_tx, rescaled_cds])
        # Sort the DataFrame by group_var (e.g., transcript) and start position
        rescaled_tx = rescaled_tx.sort(by=[group_var, 'start', 'end'] if group_var else ['start', 'end'])

    return rescaled_tx  # Return the rescaled transcript DataFrame

def _get_type(df: pl.DataFrame, df_type: str) -> pl.DataFrame:
    """
    Ensure that the DataFrame contains the correct 'type' column.

    Parameters:
    -----------
    df : pl.DataFrame
        DataFrame containing exons or introns.
    df_type : str
        Either 'exons' or 'introns', which sets the 'type' column appropriately.

    Returns:
    --------
    pl.DataFrame
        DataFrame with the 'type' column set correctly.
    """
    # If the 'type' column doesn't exist, add it with the appropriate value
    if 'type' not in df.schema:
        return df.with_columns(
            pl.lit('exon' if df_type == 'exons' else 'intron').alias('type')
        )
    # If 'type' exists and df_type is 'introns', ensure the type is correctly set to 'intron'
    elif df_type == 'introns':
        df = df.filter(pl.col('type') == 'intron')
    return df

def _get_gaps(exons: pl.DataFrame) -> pl.DataFrame:
    """
    Identify gaps between exons in a single chromosome and strand.

    Parameters:
    -----------
    exons : pl.DataFrame
        DataFrame containing exon information with 'seqnames', 'start', 'end', and 'strand'.

    Returns:
    --------
    pl.DataFrame
        DataFrame with 'start' and 'end' positions of gaps between exons.
    """
    # Ensure all exons are from a single chromosome and strand
    seqnames_unique = exons["seqnames"].n_unique()
    strand_unique = exons["strand"].n_unique()
    if seqnames_unique != 1 or strand_unique != 1:
        raise ValueError("Exons must be from a single chromosome and strand")

    # Sort exons by start position
    exons_sorted = exons.sort('start')

    # Compute cumulative maximum of 'end' shifted by 1 to identify gaps
    exons_with_cummax = exons_sorted.with_columns([
        pl.col('end').cum_max().shift(1).fill_null(0).alias('cummax_end')
    ])

    # Determine if a new group starts (i.e., no overlap with previous exons)
    exons_with_cummax = exons_with_cummax.with_columns([
        (pl.col('start') > pl.col('cummax_end')).alias('is_new_group')
    ])

    # Compute group_id as cumulative sum of 'is_new_group'
    exons_with_cummax = exons_with_cummax.with_columns([
        pl.col('is_new_group').cast(pl.Int64).cum_sum().alias('group_id')
    ])

    # Merge exons within each group to identify continuous blocks
    merged_exons = exons_with_cummax.group_by('group_id').agg([
        pl.col('start').min().alias('start'),
        pl.col('end').max().alias('end')
    ])

    # Sort merged exons by 'start'
    merged_exons = merged_exons.sort('start')

    # Compute 'prev_end' as the shifted 'end' to identify gaps
    merged_exons = merged_exons.with_columns([
        pl.col('end').shift(1).alias('prev_end')
    ])

    # Compute gap start and end positions
    merged_exons = merged_exons.with_columns([
        (pl.col('prev_end') + 1).alias('gap_start'),
        (pl.col('start') - 1).alias('gap_end')
    ])

    # Filter valid gaps where 'gap_start' is less than or equal to 'gap_end'
    gaps = merged_exons.filter(pl.col('gap_start') <= pl.col('gap_end')).select([
        pl.col('gap_start').alias('start'),
        pl.col('gap_end').alias('end')
    ])

    return gaps

def _get_tx_start_gaps(exons: pl.DataFrame, group_var: str) -> pl.DataFrame:
    """
    Identify gaps at the start of each transcript based on the first exon.

    Parameters:
    -----------
    exons : pl.DataFrame
        DataFrame containing exon information.
    group_var : str
        Column(s) used to group transcripts. Default is 'transcript_id'.

    Returns:
    --------
    pl.DataFrame
        DataFrame containing gaps at the start of each transcript.
    """
    # Get the start of the first exon for each transcript (grouped by group_var)
    tx_starts = exons.group_by(group_var).agg(pl.col('start').min())

    # Get the overall start of the first exon across all transcripts
    overall_start = exons['start'].min()

    # Use the same chromosome and strand for all transcripts
    seqnames_value = exons['seqnames'][0]
    strand_value = exons['strand'][0]

    # Create tx_start_gaps DataFrame with gap information
    tx_start_gaps = tx_starts.with_columns([
        pl.col('start').cast(pl.Int64).alias('end'),
        pl.lit(overall_start).cast(pl.Int64).alias('start'),
        pl.lit(seqnames_value).alias('seqnames'),
        pl.lit(strand_value).alias('strand'),
    ])

    return tx_start_gaps

def _get_gap_map(df: pl.DataFrame, gaps: pl.DataFrame) -> dict:
    """
    Map gaps to the corresponding introns or exons based on their start and end positions.

    Parameters:
    -----------
    df : pl.DataFrame
        DataFrame containing exons or introns.
    gaps : pl.DataFrame
        DataFrame containing gaps between exons.

    Returns:
    --------
    dict
        Dictionary containing two mappings:
        - 'equal': Gaps that exactly match with exons or introns.
        - 'pure_within': Gaps that are fully contained within exons or introns.
    """
    
    # Add an index to each gap and exon/intron row
    gaps = gaps.with_row_count("gap_index")
    df = df.with_row_count("df_index")
    
    # Find gaps where the start and end positions exactly match those of df (exons/introns)
    equal_hits = gaps.join(df, how="inner", 
                           left_on=["start", "end"], 
                           right_on=["start", "end"]).select([
                               pl.col("gap_index"), 
                               pl.col("df_index")
                           ])
    
    # Rename columns for clarity when performing the cross join
    gaps = gaps.rename({
        "start": "gaps.start",
        "end": "gaps.end"
    })

    df = df.rename({
        "start": "df.start",
        "end": "df.end"
    })
    
    # Find gaps that are fully contained within exons/introns
    within_hits = gaps.join(df, how="cross").filter(
        (pl.col("gaps.start") >= pl.col("df.start")) & 
        (pl.col("gaps.end") <= pl.col("df.end"))
    ).select([pl.col("gap_index"), pl.col("df_index")])
    
    # Remove within_hits that also appear in equal_hits (to avoid duplication)
    pure_within_hits = within_hits.join(equal_hits, how="anti", on=["df_index", "gap_index"])

    # Sort the equal_hits by gap and df index for further processing
    equal_hits = equal_hits.sort(["gap_index", "df_index"])

    # Return both the equal and pure_within mappings as a dictionary
    return {
        'equal': equal_hits,  
        'pure_within': pure_within_hits 
    }

def _get_shortened_gaps(df: pl.DataFrame, gaps: pl.DataFrame, gap_map: dict, 
                        group_var: str, target_gap_width: int) -> pl.DataFrame:
    """
    Shorten the gaps between exons or introns based on a target gap width.

    Parameters:
    -----------
    df : pl.DataFrame
        DataFrame containing exons or introns.
    gaps : pl.DataFrame
        DataFrame containing gaps between exons.
    gap_map : dict
        A dictionary mapping gaps to their corresponding exons or introns.
    group_var : str
        Column(s) used to group transcripts. Default is 'transcript_id'.
    target_gap_width : int
        The maximum allowed width for the gaps.

    Returns:
    --------
    pl.DataFrame
        DataFrame with shortened gaps and adjusted positions.
    """
    # Calculate the width of exons/introns and initialize a 'shorten_type' column
    df =  df.with_columns(
        (pl.col('end') - pl.col('start') + 1).alias('width'),  # Calculate the width
        pl.lit('none').alias('shorten_type')  # Initialize shorten_type column
    )

    # Add an index column to the df DataFrame
    df = df.with_row_count(name="df_index")

    # Update 'shorten_type' for gaps that exactly match exons/introns
    if 'equal' in gap_map and 'df_index' in gap_map['equal'].columns:
        df = df.with_columns(
            pl.when(pl.col("df_index").is_in(gap_map["equal"]["df_index"]))
            .then(pl.lit("equal"))
            .otherwise(pl.col("shorten_type"))
            .alias("shorten_type")
        )

    # Update 'shorten_type' for gaps fully within exons/introns
    if 'pure_within' in gap_map and 'df_index' in gap_map['pure_within'].columns:
        df = df.with_columns(
            pl.when(pl.col("df_index").is_in(gap_map['pure_within']['df_index']))
            .then(pl.lit("pure_within"))
            .otherwise(pl.col("shorten_type"))
            .alias("shorten_type")
        )

    # Shorten gaps that are of type 'equal' and have a width greater than the target_gap_width
    df = df.with_columns(
        pl.when((pl.col('shorten_type') == 'equal') & (pl.col('width') > target_gap_width))
        .then(pl.lit(target_gap_width))
        .otherwise(pl.col('width'))
        .alias('shortened_width')
    )

    # Handle gaps that are 'pure_within'
    if 'pure_within' in gap_map and len(gap_map['pure_within']) > 0:
        overlapping_gap_indexes = gap_map['pure_within']['gap_index']
        gaps = gaps.with_row_count(name="gap_index")

        if len(overlapping_gap_indexes) > 0:
            sum_gap_diff = pl.DataFrame({
                'intron_indexes': gap_map['pure_within']['df_index'],
                'gap_width': gaps.join(overlapping_gap_indexes.to_frame("gap_index"), on="gap_index", how="inner")
                                .with_columns((pl.col('end') - pl.col('start') + 1).alias('gap_width'))
                                .select('gap_width')
                                .to_series()
            })

            # Shorten gap width if larger than target_gap_width
            sum_gap_diff = sum_gap_diff.with_columns(
                pl.when(pl.col('gap_width') > target_gap_width)
                .then(pl.lit(target_gap_width))
                .otherwise(pl.col('gap_width'))
                .alias('shortened_gap_width')
            )

            # Calculate the gap difference
            sum_gap_diff = sum_gap_diff.with_columns(
                (pl.col('gap_width') - pl.col('shortened_gap_width')).alias('shortened_gap_diff')
            )

            # Aggregate gap differences by intron indexes
            sum_gap_diff = sum_gap_diff.group_by('intron_indexes').agg(
                pl.sum('shortened_gap_diff').alias('sum_shortened_gap_diff')
            )

            # Join the calculated gap differences with the df DataFrame
            df = df.join(sum_gap_diff, left_on='df_index', right_on='intron_indexes', how='left')

            # Adjust the width based on gap differences
            df = df.with_columns(
                pl.when(pl.col('sum_shortened_gap_diff').is_null())
                  .then(pl.col('shortened_width'))
                  .otherwise(pl.col('width') - pl.col('sum_shortened_gap_diff'))
                  .alias('shortened_width')
            )

            # Clean up unnecessary columns
            df = df.drop(['sum_shortened_gap_diff', 'shorten_type', 'width'])
            df = df.rename({'shortened_width': 'width'})

        df = df.drop('df_index')

    return df

def _get_rescaled_txs(
    exons: pl.DataFrame,
    introns_shortened: pl.DataFrame,
    tx_start_gaps_shortened: pl.DataFrame,
    group_var: str
) -> pl.DataFrame:
    """
    Rescale transcript coordinates based on shortened gaps for exons and introns.

    Parameters:
    -----------
    exons : pl.DataFrame
        DataFrame containing exon information.
    introns_shortened : pl.DataFrame
        DataFrame containing intron information with shortened gaps.
    tx_start_gaps_shortened : pl.DataFrame
        DataFrame containing rescaled transcript start gaps.
    group_var : str
        Column(s) used to group transcripts. Default is 'transcript_id'.

    Returns:
    --------
    pl.DataFrame
        Rescaled transcript DataFrame with adjusted coordinates.
    """
    # Clone exons to avoid altering the original DataFrame
    exons = exons.clone()

    # Define columns to keep for introns, including 'width'
    column_to_keep = exons.columns + ["width"]

    # Select and reorder columns for the shortened introns
    introns_shortened = introns_shortened.select(column_to_keep)

    # Add a new 'width' column to exons representing their lengths
    exons = exons.with_columns(
        (pl.col('end') - pl.col('start') + 1).alias('width')
    )

    # Concatenate exons and shortened introns into a single DataFrame
    rescaled_tx = pl.concat([exons, introns_shortened], how='vertical')

    # Sort the DataFrame by group_var (e.g., 'transcript_id') and start position
    if group_var:
        if isinstance(group_var, list):
            sort_columns = group_var + ['start', 'end']
        else:
            sort_columns = [group_var, 'start', 'end']
    else:
        sort_columns = ['start', 'end']
    rescaled_tx = rescaled_tx.sort(sort_columns)

    # Calculate cumulative sum for rescaled start and end positions
    if group_var:
        rescaled_tx = rescaled_tx.with_columns(
            pl.col('width').cum_sum().over(group_var).alias('rescaled_end')
        )
    else:
        rescaled_tx = rescaled_tx.with_columns(
            pl.col('width').cum_sum().alias('rescaled_end')
        )
    
    # Compute the rescaled start positions based on the cumulative end positions
    rescaled_tx = rescaled_tx.with_columns(
        (pl.col('rescaled_end') - pl.col('width') + 1).alias('rescaled_start')
    )

    # If no group_var, set all transcripts to start at position 1
    if not group_var:
        rescaled_tx = rescaled_tx.with_columns(
            pl.lit(1).alias('width_tx_start')
        )
    else:
        # Join rescaled transcript start gaps to adjust start positions
        rescaled_tx = rescaled_tx.join(
            tx_start_gaps_shortened, on=group_var, how='left', suffix='_tx_start'
        )
    
    # Adjust the rescaled start and end positions based on transcript start gaps
    rescaled_tx = rescaled_tx.with_columns([
        (pl.col('rescaled_end') + pl.col('width_tx_start')).alias('rescaled_end'),
        (pl.col('rescaled_start') + pl.col('width_tx_start')).alias('rescaled_start')
    ])

    # Adjust intron start and end positions to avoid overlap with exons
    rescaled_tx = rescaled_tx.with_columns([
        pl.when(pl.col('type') == 'intron')
          .then(pl.col('start') - 1)
          .otherwise(pl.col('start'))
          .alias('start'),
        pl.when(pl.col('type') == 'intron')
          .then(pl.col('end') + 1)
          .otherwise(pl.col('end'))
          .alias('end'),
        pl.when(pl.col('type') == 'intron')
          .then(pl.col('rescaled_start') - 1)
          .otherwise(pl.col('rescaled_start'))
          .alias('rescaled_start'),
        pl.when(pl.col('type') == 'intron')
          .then(pl.col('rescaled_end') + 1)
          .otherwise(pl.col('rescaled_end'))
          .alias('rescaled_end')
    ])
    
    # Drop original 'start', 'end', and 'width' columns, renaming rescaled columns
    rescaled_tx = rescaled_tx.drop(['start', 'end', 'width']).rename({
        'rescaled_start': 'start',
        'rescaled_end': 'end'
    })
    
    # Reorder columns for consistency in the output
    columns = rescaled_tx.columns
    column_order = ['seqnames', 'start', 'end', 'strand'] + [
        col for col in columns if col not in ['seqnames', 'start', 'end', 'strand']
    ]
    rescaled_tx = rescaled_tx.select(column_order)
    
    return rescaled_tx  # Return the rescaled transcript coordinates

def _get_cds_exon_difference(gene_exons: pl.DataFrame, gene_cds_regions: pl.DataFrame) -> pl.DataFrame:
    """
    Calculate the absolute differences between the start and end positions of exons and CDS regions.

    Parameters:
    ----------
    gene_exons : pl.DataFrame
        DataFrame containing exon regions.
    gene_cds_regions : pl.DataFrame
        DataFrame containing CDS (Coding DNA Sequence) regions.

    Returns:
    -------
    pl.DataFrame
        DataFrame with the absolute differences between exon and CDS start/end positions.
    """

    # Rename 'start' and 'end' columns in CDS regions for clarity
    cds_regions = gene_cds_regions.rename({'start': 'cds_start', 'end': 'cds_end'})

    # Remove the 'type' column if it exists in CDS
    if 'type' in cds_regions.columns:
        cds_regions = cds_regions.drop('type')

    # Rename 'start' and 'end' columns in exon regions for clarity
    exons = gene_exons.rename({'start': 'exon_start', 'end': 'exon_end'})

    # Remove the 'type' column if it exists in exons
    if 'type' in exons.columns:
        exons = exons.drop('type')

    # Identify common columns to join CDS and exons on (e.g., transcript_id)
    common_columns = list(set(cds_regions.columns) & set(exons.columns))
    if not common_columns:
        raise ValueError("No common columns to join on. Ensure both DataFrames have common keys.")

    # Perform left join between CDS and exon data on the common columns
    cds_exon_diff = cds_regions.join(exons, on=common_columns, how='left')

    # Calculate absolute differences between exon and CDS start positions
    cds_exon_diff = cds_exon_diff.with_columns(
        (pl.col('exon_start') - pl.col('cds_start')).abs().alias('diff_start')
    )

    # Calculate absolute differences between exon and CDS end positions
    cds_exon_diff = cds_exon_diff.with_columns(
        (pl.col('exon_end') - pl.col('cds_end')).abs().alias('diff_end')
    )

    return cds_exon_diff

def _get_rescale_cds(cds_exon_diff: pl.DataFrame, gene_rescaled_exons: pl.DataFrame) -> pl.DataFrame:
    """
    Rescale CDS regions based on exon positions and the calculated differences between them.

    Parameters:
    ----------
    cds_exon_diff : pl.DataFrame
        DataFrame with differences between exon and CDS start/end positions.
    gene_rescaled_exons : pl.DataFrame
        DataFrame containing rescaled exon positions.

    Returns:
    -------
    pl.DataFrame
        Rescaled CDS positions based on exon positions.
    """

    # Assign a 'type' column with the value "CDS" and drop unnecessary columns
    columns_to_drop = ['cds_start', 'cds_end', 'exon_start', 'exon_end']
    cds_prepared = (
        cds_exon_diff
        .with_columns(pl.lit("CDS").alias("type"))
        .drop([col for col in columns_to_drop if col in cds_exon_diff.columns])
    )

    # Rename columns in rescaled exons for consistency
    exons_prepared = gene_rescaled_exons.rename({'start': 'exon_start', 'end': 'exon_end'})

    # Drop 'type' column if present
    if 'type' in exons_prepared.columns:
        exons_prepared = exons_prepared.drop('type')

    # Identify common columns to join CDS and rescaled exons
    common_columns = [col for col in cds_prepared.columns if col in exons_prepared.columns]
    if not common_columns:
        raise ValueError("No common columns to join on. Ensure both DataFrames have common keys.")

    # Perform left join on common columns
    gene_rescaled_cds = cds_prepared.join(exons_prepared, on=common_columns, how='left')

    # Adjust start and end positions of CDS based on exon positions
    gene_rescaled_cds = gene_rescaled_cds.with_columns([
        (pl.col('exon_start') + pl.col('diff_start')).alias('start'),
        (pl.col('exon_end') - pl.col('diff_end')).alias('end')
    ])

    # Drop unnecessary columns used for the difference calculations
    gene_rescaled_cds = gene_rescaled_cds.drop(['exon_start', 'exon_end', 'diff_start', 'diff_end'])

    return gene_rescaled_cds