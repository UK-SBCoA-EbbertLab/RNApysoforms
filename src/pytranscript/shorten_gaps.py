import pandas as pd  # Import pandas for data manipulation
import numpy as np  # Import numpy for numerical operations
import plotly.graph_objects as go  # Import Plotly for creating plots
from typing import List, Union  # Import type annotations for functions
from pytranscript.to_intron import to_intron ## Import to intron function
import polars as pl


def shorten_gaps(annotation: pd.DataFrame, 
                 group_var: Union[str, List[str]] = None, 
                 target_gap_width: int = 100) -> pd.DataFrame:
    """
    Shorten the gaps between exons and introns for more compact transcript visualization.

    Parameters:
    -----------
    annotation : pd.DataFrame
        DataFrame containing exon, intron, and possibly CDS information
    introns : pd.DataFrame
        DataFrame containing intron information.
    group_var : str or List[str], optional
        Column(s) used to group transcripts. Default is None.
    target_gap_width : int, optional
        Maximum allowed width for gaps between exons. Default is 100.

    Returns:
    --------
    pd.DataFrame
        DataFrame with shortened intron gaps and rescaled coordinates.
    """

  # Separate introns, exons, and possibly CDS
    exons = annotation.filter(pl.col("type") == "exon")
    introns = to_intron(exons=exons, group_var=group_var)

    if "CDS" in annotation["type"].unique().to_list():
        cds = annotation.filter(pl.col("type") == "CDS")
    else:
        cds = None

    # Ensure required columns are present in both the exons and introns DataFrames
    for df in [exons, introns]:
        required_columns = ['start', 'end', 'strand', 'seqnames']
        assert all(col in df.columns for col in required_columns), \
            "Both exons and introns DataFrames must have 'start', 'end', 'strand', and 'seqnames' columns"

    # Ensure each DataFrame has a 'type' column indicating 'exon' or 'intron'
    exons = _get_type(exons , "exons")  # Set 'type' to 'exon' for exons DataFrame
    introns = _get_type(introns , "introns")  # Set 'type' to 'intron' for introns DataFrame

    # Adjust the start and end positions of introns slightly to avoid overlap with exons
    introns = introns.with_columns([
        pl.col("start") + 1,
        pl.col("end") - 1
    ])

    # Identify gaps between exons
    gaps = _get_gaps(exons)  # Find gaps between exons on the same chromosome and strand

    # Map the identified gaps to the corresponding introns
    gap_map = _get_gap_map(introns, gaps)  # Map gaps to introns (or exons)

    # Shorten the gaps based on the target gap width
    introns_shortened = _get_shortened_gaps(introns, gaps, gap_map,
                                                 group_var, target_gap_width)


    # If a group variable (e.g., transcript name) is provided, handle gaps at transcript starts
    if group_var:
        tx_start_gaps = _get_tx_start_gaps(exons, group_var)  # Identify gaps at transcript starts
        gap_map_tx_start = _get_gap_map(tx_start_gaps, gaps)  # Map gaps at transcript starts
        tx_start_gaps_shortened = _get_shortened_gaps(tx_start_gaps, gaps, gap_map_tx_start, group_var, target_gap_width)
        # Remove unnecessary columns from transcript start gaps
        tx_start_gaps_shortened = tx_start_gaps_shortened.drop(['start', 'end', 'strand', 'seqnames'])
    else:
        tx_start_gaps_shortened = None  # No transcript start gaps if no group_var is provided

    # Rescale the coordinates of exons and introns after shortening the gaps
    rescaled_tx= _get_rescaled_txs(exons, introns_shortened,
                                             tx_start_gaps_shortened, group_var)

    ## Process CDS if there at all
    if isinstance(cds, pl.DataFrame):
        cds_diff = _get_cds_exon_difference(cds, exons)
        rescaled_cds = _get_rescale_cds(cds_diff, rescaled_tx.filter(pl.col("type") == "exon"))
        rescaled_cds = rescaled_cds[rescaled_tx.columns]
        # Combine the CDS into final dataframe
        rescaled_tx = pl.concat([rescaled_tx, rescaled_cds])
        # Sort the DataFrame by group_var (e.g., transcript) and start position
        rescaled_tx = rescaled_tx.sort(by=[group_var, 'start', 'end'] if group_var else ['start', 'end'])

    
    return rescaled_tx  # Return the DataFrame with rescaled transcript coordinates

def _get_type(df: pl.DataFrame, df_type: str) -> pl.DataFrame:
    # If the 'type' column doesn't exist, add it with the appropriate value
    if 'type' not in df.schema:
        return df.with_columns(
            pl.lit('exon' if df_type == 'exons' else 'intron').alias('type')
        )
    # If 'type' exists and df_type is 'introns', assert all values are 'intron'
    elif df_type == 'introns':
        df = df.filter(pl.col('type') == 'intron')
    return df


def _get_gaps(exons: pl.DataFrame) -> pl.DataFrame:
    # Ensure exons are from a single chromosome and strand
    exons_single_chrom_strand = exons.filter(
        (exons["seqnames"].n_unique() == 1) & (exons["strand"].n_unique() == 1)
    )

    # Sort exons by start position
    exons_sorted = exons_single_chrom_strand.sort("start")

    # Create lagged 'end' column to find gaps
    exons_with_lag = exons_sorted.with_columns([
        pl.col("end").shift(1).alias("prev_end")
    ])

    # Compute gaps between consecutive exons
    gaps = exons_with_lag.filter(pl.col("start") > (pl.col("prev_end") + 1)).select([
        (pl.col("prev_end") + 1).alias("start"),
        (pl.col("start") - 1).alias("end")
    ])

    return gaps

def _get_tx_start_gaps(exons: pl.DataFrame, group_var: Union[str, List[str]]) -> pl.DataFrame:
    """
    Identify gaps at the start of each transcript based on the first exon.

    Parameters:
    -----------
    exons : pl.DataFrame
        DataFrame containing exon information.
    group_var : str or List[str]
        Column name(s) used to group transcripts.

    Returns:
    --------
    pl.DataFrame
        DataFrame containing gaps at the start of each transcript.
    """
    # Get the start of the first exon for each transcript (grouped by group_var)
    tx_starts = exons.group_by(group_var).agg(pl.col('start').min())

    # Get the overall start of the first exon
    overall_start = exons['start'].min()

    # Use the same chromosome and strand for all transcripts
    seqnames_value = exons['seqnames'][0]
    strand_value = exons['strand'][0]

    # Create tx_start_gaps DataFrame
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
        DataFrame containing gaps.

    Returns:
    --------
    dict
        Dictionary with 'equal' gaps (exact matches) and 'pure_within' gaps (fully within an intron or exon).
    """
    
    # Add index as a new column
    gaps = gaps.with_row_count("gap_index")
    df = df.with_row_count("df_index")

    # Exact matches: gaps where start and end positions match with df (exons/introns)
    equal_hits = gaps.join(df, how="inner", 
                           left_on=["start", "end"], 
                           right_on=["start", "end"]) \
                     .select([pl.col("gap_index"), 
                              pl.col("df_index")])
    
    
    # Rename columns for clarity (excluding index)
    gaps = gaps.rename({
        "start": "gaps.start",
        "end": "gaps.end"
    })

    df = df.rename({
        "start": "df.start",
        "end": "df.end"
    })
    
    
    # Gaps that are fully within exons/introns
    within_hits = gaps.join(df, how="cross") \
                      .filter((pl.col("gaps.start") >= pl.col("df.start")) & 
                              (pl.col("gaps.end") <= pl.col("df.end"))).select([pl.col("gap_index"), 
                               pl.col("df_index")])

    # Remove within_hits that also appear in equal_hits
    pure_within_hits = within_hits.join(equal_hits, how="anti", on=["df_index", "gap_index"])

    ## Sort by gap index and then df_index
    equal_hits = equal_hits.sort(["gap_index", "df_index"])

    # Return the results as a dictionary
    return {
        'equal': equal_hits,  
        'pure_within': pure_within_hits 
    }

def _get_shortened_gaps(df: pl.DataFrame, gaps: pl.DataFrame, gap_map: dict, 
                        group_var: Union[str, List[str]], target_gap_width: int) -> pl.DataFrame:

    # Calculate 'width' and initialize 'shorten_type'
    df =  df.with_columns(
        (pl.col('end') - pl.col('start') + 1).alias('width'),  # Calculate width
        pl.lit('none').alias('shorten_type')  # Initialize shorten_type column
    )

    ## Add df_index column
    df = df.with_row_count(name="df_index")

    # Conditionally update 'shorten_type' for 'equal'
    if 'equal' in gap_map and 'df_index' in gap_map['equal'].columns:

        # Update the 'shorten_type' column in df based on whether the index is in gap_map["equal"]["df_index"]
        df = df.with_columns(
            pl.when(pl.col("df_index").is_in(gap_map["equal"]["df_index"]))
            .then(pl.lit("equal"))
            .otherwise(pl.col("shorten_type"))
            .alias("shorten_type")
        )

    # Conditionally update 'shorten_type' for 'pure_within'
    if 'pure_within' in gap_map and 'df_index' in gap_map['pure_within']:
        df = df.with_columns(
            pl.when(pl.col("df_index").is_in(gap_map['pure_within']['df_index']))
            .then(pl.lit("pure_within"))
            .otherwise(pl.col("shorten_type"))
            .alias("shorten_type")
        )

    df = df.with_columns(
    pl.when((pl.col('shorten_type') == 'equal') & (pl.col('width') > target_gap_width))
      .then(pl.lit(target_gap_width))
      .otherwise(pl.col('width'))
      .alias('shortened_width')
    )

    if 'pure_within' in gap_map and len(gap_map['pure_within']) > 0:
        
        overlapping_gap_indexes = gap_map['pure_within']['gap_index']

    
        ## Add df_index column
        gaps = gaps.with_row_count(name="gap_index")

        if len(overlapping_gap_indexes) > 0:

            sum_gap_diff = pl.DataFrame({
                'intron_indexes': gap_map['pure_within']['df_index'],
                'gap_width': gaps.join(overlapping_gap_indexes.to_frame("gap_index"), on="gap_index", how="inner")
                                .with_columns((pl.col('end') - pl.col('start') + 1).alias('gap_width'))
                                .select('gap_width')
                                .to_series()
            })

            # Step 1: Add 'shortened_gap_width'
            sum_gap_diff = sum_gap_diff.with_columns(
                pl.when(pl.col('gap_width') > target_gap_width)
                .then(pl.lit(target_gap_width))
                .otherwise(pl.col('gap_width'))
                .alias('shortened_gap_width')
            )

            # Step 2: Add 'shortened_gap_diff' based on the newly created 'shortened_gap_width'
            sum_gap_diff = sum_gap_diff.with_columns(
                (pl.col('gap_width') - pl.col('shortened_gap_width')).alias('shortened_gap_diff')
            )

            sum_gap_diff = sum_gap_diff.group_by('intron_indexes').agg(
                pl.sum('shortened_gap_diff').alias('sum_shortened_gap_diff')
            )

            df = df.join(sum_gap_diff, left_on='df_index', right_on='intron_indexes', how='left')

            df = df.with_columns(
                pl.when(pl.col('sum_shortened_gap_diff').is_null())
                  .then(pl.col('shortened_width'))
                  .otherwise(pl.col('width') - pl.col('sum_shortened_gap_diff'))
                  .alias('shortened_width')
            )

            df = df.drop(['sum_shortened_gap_diff', 'shorten_type', 'width'])

            df = df.drop('df_index')

    return df.rename({'shortened_width': 'width'})

def _get_rescaled_txs(
    exons: pl.DataFrame,
    introns_shortened: pl.DataFrame,
    tx_start_gaps_shortened: pl.DataFrame,
    group_var: Union[str, List[str]]
) -> pl.DataFrame:
    """
    Rescale transcript coordinates based on shortened gaps.

    Parameters:
    -----------
    exons : pl.DataFrame
        DataFrame containing exon information.
    introns_shortened : pl.DataFrame
        DataFrame containing intron information with shortened gaps.
    tx_start_gaps_shortened : pl.DataFrame
        DataFrame containing rescaled transcript start gaps.
    group_var : str or List[str]
        Column(s) used to group transcripts.

    Returns:
    --------
    pl.DataFrame
        Rescaled transcript DataFrame with adjusted coordinates.
    """
    # Clone exons DataFrame to avoid modifying the original
    exons = exons.clone()

    ## Column to keep
    column_to_keep = exons.columns + ["width"]

    ## Reorder intron columns
    introns_shortened = introns_shortened.select(column_to_keep)

    # Calculate the width of each exon
    exons = exons.with_columns(
        (pl.col('end') - pl.col('start') + 1).alias('width')
    )
    
    # Combine the exons and shortened introns into one DataFrame
    rescaled_tx = pl.concat([exons, introns_shortened], how='vertical')
    
    # Sort the DataFrame by group_var (e.g., transcript) and start position
    if group_var:
        if isinstance(group_var, list):
            sort_columns = group_var + ['start', 'end']
        else:
            sort_columns = [group_var, 'start', 'end']
    else:
        sort_columns = ['start', 'end']
    rescaled_tx = rescaled_tx.sort(sort_columns)
    
    # Calculate cumulative sum to get rescaled start and end positions
    if group_var:
        rescaled_tx = rescaled_tx.with_columns(
            pl.col('width').cum_sum().over(group_var).alias('rescaled_end')
        )
    else:
        rescaled_tx = rescaled_tx.with_columns(
            pl.col('width').cum_sum().alias('rescaled_end')
        )
    rescaled_tx = rescaled_tx.with_columns(
        (pl.col('rescaled_end') - pl.col('width') + 1).alias('rescaled_start')
    )
    
    if not group_var:  # If no grouping, all transcripts start at position 1
        rescaled_tx = rescaled_tx.with_columns(
            pl.lit(1).alias('width_tx_start')
        )
    else:
        # Merge rescaled transcript start gaps to adjust start positions
        rescaled_tx = rescaled_tx.join(
            tx_start_gaps_shortened, on=group_var, how='left', suffix='_tx_start'
        )
    
    # Adjust the rescaled start and end positions based on transcript start gaps
    rescaled_tx = rescaled_tx.with_columns([
        (pl.col('rescaled_end') + pl.col('width_tx_start')).alias('rescaled_end'),
        (pl.col('rescaled_start') + pl.col('width_tx_start')).alias('rescaled_start')
    ])
    
    # Adjust start and end positions for introns to avoid exon overlap
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
    
    # Drop original start, end, and width columns, and rename rescaled columns to 'start' and 'end'
    rescaled_tx = rescaled_tx.drop(['start', 'end', 'width']).rename({
        'rescaled_start': 'start',
        'rescaled_end': 'end'
    })
    
    # Reorder columns to ensure consistency in output
    columns = rescaled_tx.columns
    column_order = ['seqnames', 'start', 'end', 'strand'] + [
        col for col in columns if col not in ['seqnames', 'start', 'end', 'strand']
    ]
    rescaled_tx = rescaled_tx.select(column_order)
    
    return rescaled_tx  # Return the DataFrame with rescaled transcript coordinates



def _get_cds_exon_difference(gene_exons: pl.DataFrame, gene_cds_regions: pl.DataFrame) -> pl.DataFrame:
    """
    Calculates the absolute differences between the start and end positions of exons and CDS regions.
    This function is used to prepare data for re-scaling CDS regions based on exon positions.

    Parameters:
    ----------
    gene_cds_regions : pl.DataFrame
        DataFrame containing CDS (Coding DNA Sequence) regions with at least 'start' and 'end' columns,
        along with any common columns used for joining (e.g., 'transcript_id', 'gene_id').
    gene_exons : pl.DataFrame
        DataFrame containing exon regions with at least 'start' and 'end' columns,
        along with any common columns used for joining.

    Returns:
    -------
    cds_exon_diff : pl.DataFrame
        DataFrame resulting from a left join of the CDS and exon data, including the calculated
        absolute differences 'diff_start' and 'diff_end' between exon and CDS start and end positions.
    """

    # Step 1: Rename 'start' and 'end' columns in CDS regions to 'cds_start' and 'cds_end'
    cds_regions = gene_cds_regions.rename({'start': 'cds_start', 'end': 'cds_end'})

    # Remove the 'type' column if it exists
    if 'type' in cds_regions.columns:
        cds_regions = cds_regions.drop('type')

    # Step 2: Rename 'start' and 'end' columns in exon regions to 'exon_start' and 'exon_end'
    exons = gene_exons.rename({'start': 'exon_start', 'end': 'exon_end'})

    # Remove the 'type' column if it exists
    if 'type' in exons.columns:
        exons = exons.drop('type')

    # Step 3: Identify common columns to perform the left join
    common_columns = list(set(cds_regions.columns) & set(exons.columns))
    if not common_columns:
        raise ValueError("No common columns to perform join on. Ensure both DataFrames have common keys.")

    # Step 4: Perform the left join on the common columns
    cds_exon_diff = cds_regions.join(exons, on=common_columns, how='left')

    # Step 5: Calculate absolute differences between exon and CDS start positions
    cds_exon_diff = cds_exon_diff.with_columns(
        (pl.col('exon_start') - pl.col('cds_start')).abs().alias('diff_start')
    )

    # Step 6: Calculate absolute differences between exon and CDS end positions
    cds_exon_diff = cds_exon_diff.with_columns(
        (pl.col('exon_end') - pl.col('cds_end')).abs().alias('diff_end')
    )

    return cds_exon_diff

def _get_rescale_cds(cds_exon_diff: pl.DataFrame, gene_rescaled_exons: pl.DataFrame) -> pl.DataFrame:
    """
    Rescales CDS regions based on exon positions and the calculated differences between CDS and exon positions.
    This function aligns the CDS regions to the rescaled exon positions and adjusts their start and end points
    accordingly.

    Parameters
    ----------
    cds_exon_diff : pl.DataFrame
        DataFrame containing the differences between the start and end positions of exons and CDS regions.
        Expected columns:
        - 'diff_start': Absolute difference between exon start and CDS start positions.
        - 'diff_end': Absolute difference between exon end and CDS end positions.
        - Any columns necessary for joining (e.g., 'transcript_id', 'gene_id').
    gene_rescaled_exons : pl.DataFrame
        DataFrame containing the rescaled exon positions.
        Expected columns:
        - 'start': Rescaled start position of the exon.
        - 'end': Rescaled end position of the exon.
        - Any columns necessary for joining.

    Returns
    -------
    gene_rescaled_cds : pl.DataFrame
        DataFrame containing the rescaled CDS positions, with adjusted 'start' and 'end' positions,
        and 'transcript_id' ordered according to 'factor_order'.
    """

    # Step 1: Prepare the CDS DataFrame
    # - Assign a new column 'type' with the value "CDS"
    # - Drop unnecessary columns: 'cds_start', 'cds_end', 'exon_start', 'exon_end' if they exist
    columns_to_drop = ['cds_start', 'cds_end', 'exon_start', 'exon_end']
    existing_columns = [col for col in columns_to_drop if col in cds_exon_diff.columns]
    cds_prepared = (
        cds_exon_diff
        .with_columns(pl.lit("CDS").alias("type"))
        .drop(existing_columns)
    )

    # Step 2: Prepare the Exon DataFrame
    # - Rename 'start' to 'exon_start' and 'end' to 'exon_end' to avoid column name conflicts
    # - Drop the 'type' column if it exists
    exons_columns_to_drop = ['type']
    exons_existing_columns = [col for col in exons_columns_to_drop if col in gene_rescaled_exons.columns]
    exons_prepared = (
        gene_rescaled_exons
        .rename({'start': 'exon_start', 'end': 'exon_end'})
        .drop(exons_existing_columns)
    )

    # Step 3: Identify common columns for joining
    # - Find columns that are present in both DataFrames to use as keys for the join
    common_columns = [col for col in cds_prepared.columns if col in exons_prepared.columns]
    if not common_columns:
        raise ValueError("No common columns to perform join on. Ensure both DataFrames have common keys.")

    # Step 4: Perform the left join on the common columns
    # - This aligns the CDS data with the corresponding rescaled exon positions
    gene_rescaled_cds = cds_prepared.join(exons_prepared, on=common_columns, how='left')

    # Step 5: Calculate the adjusted 'start' and 'end' positions for the CDS regions
    # - 'start' is adjusted by adding 'diff_start' to 'exon_start'
    # - 'end' is adjusted by subtracting 'diff_end' from 'exon_end'
    gene_rescaled_cds = gene_rescaled_cds.with_columns([
        (pl.col('exon_start') + pl.col('diff_start')).alias('start'),
        (pl.col('exon_end') - pl.col('diff_end')).alias('end')
    ])

    # Step 6: Drop unnecessary columns used for calculations
    # - Remove 'exon_start', 'exon_end', 'diff_start', 'diff_end' as they are no longer needed
    final_columns_to_drop = ['exon_start', 'exon_end', 'diff_start', 'diff_end']
    final_existing_columns = [col for col in final_columns_to_drop if col in gene_rescaled_cds.columns]
    gene_rescaled_cds = gene_rescaled_cds.drop(final_existing_columns)

    return gene_rescaled_cds