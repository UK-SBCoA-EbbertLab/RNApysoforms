import polars as pl
import plotly.graph_objects as go
from typing import List, Union
from RNApysoforms.to_intron import to_intron
from RNApysoforms.utils import check_df


def shorten_gaps(
    annotation: pl.DataFrame,
    transcript_id_column: str = "transcript_id",
    target_gap_width: int = 100
) -> pl.DataFrame:
    """
    Shortens intron and transcript start gaps between exons in genomic annotations to enhance visualization.

    This function processes genomic annotations by shortening the widths of intron gaps and gaps at the start of transcripts
    to a specified target size, while preserving exon and CDS regions. The goal is to improve the clarity of transcript
    visualizations by reducing the visual space occupied by long intron regions and aligning transcripts for consistent
    rescaling, maintaining the relative structure of the transcripts.

    Parameters
    ----------
    annotation : pl.DataFrame
        A Polars DataFrame containing genomic annotations, including exons and optionally CDS and intron data.
        Required columns include:
        - 'start': Start position of the feature.
        - 'end': End position of the feature.
        - 'type': Feature type, expected to include 'exon' and optionally 'intron' and 'CDS'.
        - 'strand': Strand information ('+' or '-').
        - 'seqnames': Chromosome or sequence name.
        - transcript_id_column: Column used to group transcripts, typically 'transcript_id'.
    transcript_id_column : str, optional
        The column used to group transcripts, by default "transcript_id". This identifies individual transcripts
        within the annotation data.
    target_gap_width : int, optional
        The maximum width for intron gaps and transcript start gaps after shortening. Gaps wider than this will be reduced
        to this size. Default is 100.

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame with shortened intron and transcript start gaps and rescaled coordinates for exons,
        introns, and CDS regions. The DataFrame includes:
        - Original columns from the input DataFrame.
        - 'rescaled_start': The rescaled start position after shortening gaps.
        - 'rescaled_end': The rescaled end position after shortening gaps.

    Raises
    ------
    TypeError
        If 'annotation' is not a Polars DataFrame.
    ValueError
        If required columns are missing in the input DataFrame.
        If exons are not from a single chromosome and strand when calculating gaps.
        If there are no common columns to join on between CDS and exons when processing CDS regions.

    Examples
    --------
    Shorten intron and transcript start gaps in a genomic annotation DataFrame:

    >>> import polars as pl
    >>> from RNApysoforms.plot import shorten_gaps
    >>> df = pl.DataFrame({
    ...     "transcript_id": ["tx1", "tx1", "tx1"],
    ...     "start": [100, 200, 500],
    ...     "end": [150, 250, 600],
    ...     "type": ["exon", "exon", "exon"],
    ...     "strand": ["+", "+", "+"],
    ...     "seqnames": ["chr1", "chr1", "chr1"]
    ... })
    >>> shortened_df = shorten_gaps(df, transcript_id_column="transcript_id", target_gap_width=50)
    >>> print(shortened_df.head())

    This will return a DataFrame where the intron and transcript start gaps have been shortened to a maximum width of 50,
    and includes rescaled coordinates for visualization.

    Notes
    -----
    - The function ensures that exon and CDS regions maintain their original lengths, while intron gaps and transcript start
      gaps are shortened.
    - If intron entries are not present in the input DataFrame, the function generates them using the 'to_intron' function.
    - The input DataFrame must contain the required columns listed above.
    - The function processes gaps at the start of transcripts to align transcripts for consistent rescaling.
    - After shortening gaps, the coordinates are rescaled to maintain the relative positions of features within and across transcripts.
    - The function returns the rescaled DataFrame with original columns plus 'rescaled_start' and 'rescaled_end'.
    """

    # Check if annotation is a Polars DataFrame
    if not isinstance(annotation, pl.DataFrame):
        raise TypeError(
            f"Expected 'annotation' to be of type pl.DataFrame, got {type(annotation)}."
            "\nYou can convert a pandas DataFrame to Polars using: polars_df = pl.from_pandas(pandas_df)"
        )

    # Validate the input DataFrame to ensure required columns are present
    check_df(annotation, ["start", "end", "type", "strand", "seqnames", transcript_id_column])

    # Check if there are intron entries in the annotation data
    if "intron" in annotation["type"].unique().to_list():
        # Separate intron data if present
        introns = annotation.filter(pl.col("type") == "intron")
    else:
        # Generate intron entries if they are not present
        annotation = to_intron(annotation=annotation, transcript_id_column=transcript_id_column)
        introns = annotation.filter(pl.col("type") == "intron")  # Separate intron data

    # Check if there are CDS entries in the annotation data
    if "CDS" in annotation["type"].unique().to_list():
        # Separate CDS data if present
        cds = annotation.filter(pl.col("type") == "CDS")
    else:
        cds = None  # No CDS entries in the data

    # Separate exons from the annotation data
    exons = annotation.filter(pl.col("type") == "exon")

    # Ensure the 'type' column in exons and introns is set correctly
    exons = _get_type(exons, "exons")  # Mark the type as 'exon'
    introns = _get_type(introns, "introns")  # Mark the type as 'intron'

    # Adjust intron positions to avoid overlap with exons
    introns = introns.with_columns([
        pl.col("start") + 1,
        pl.col("end") - 1
    ])

    # Identify gaps between exons within the same chromosome and strand
    gaps = _get_gaps(exons)

    # Map gaps to introns to identify which gaps correspond to which introns
    gap_map = _get_gap_map(introns, gaps)

    # Shorten gaps based on the target gap width
    introns_shortened = _get_shortened_gaps(
        introns, gaps, gap_map, transcript_id_column, target_gap_width
    )

    # Handle gaps at the start of transcripts to align them
    tx_start_gaps = _get_tx_start_gaps(exons, transcript_id_column)  # Gaps at the start of transcripts
    gap_map_tx_start = _get_gap_map(tx_start_gaps, gaps)
    tx_start_gaps_shortened = _get_shortened_gaps(
        tx_start_gaps, gaps, gap_map_tx_start, transcript_id_column, target_gap_width
    )
    tx_start_gaps_shortened = tx_start_gaps_shortened.drop(['start', 'end', 'strand', 'seqnames'])

    # Rescale the coordinates of exons and introns after shortening the gaps
    rescaled_tx = _get_rescaled_txs(
        exons, introns_shortened, tx_start_gaps_shortened, transcript_id_column
    )

    # Process CDS regions if available
    if isinstance(cds, pl.DataFrame):
        # Calculate differences between exons and CDS regions
        cds_diff = _get_cds_exon_difference(exons, cds)
        # Rescale CDS regions based on the rescaled exons
        rescaled_cds = _get_rescale_cds(cds_diff, rescaled_tx.filter(pl.col("type") == "exon"))
        rescaled_cds = rescaled_cds[rescaled_tx.columns]
        # Combine the rescaled CDS data into the final DataFrame
        rescaled_tx = pl.concat([rescaled_tx, rescaled_cds])

    # Sort the DataFrame by start and end positions
    rescaled_tx = rescaled_tx.sort(by=['start', 'end'])

    # Return transcripts in original order they were given
    original_order = annotation[transcript_id_column].unique(maintain_order=True).to_list()
    order_mapping = {transcript: index for index, transcript in enumerate(original_order)}
    rescaled_tx = (rescaled_tx
                   .with_columns(pl.col(transcript_id_column).replace(order_mapping).alias("order"))
                   .sort("order")
                   .drop("order"))

    # Include original columns and rescaled coordinates in the final DataFrame
    final_columns = annotation.columns + ["rescaled_start", "rescaled_end"]
    rescaled_tx = rescaled_tx[final_columns].clone()

    return rescaled_tx  # Return the rescaled transcript DataFrame


def _get_type(df: pl.DataFrame, df_type: str) -> pl.DataFrame:
    """
    Ensures that the 'type' column in the DataFrame is correctly set to 'exon' or 'intron'.

    Parameters
    ----------
    df : pl.DataFrame
        A Polars DataFrame containing genomic features.
    df_type : str
        The type to set in the 'type' column, either 'exons' or 'introns'.

    Returns
    -------
    pl.DataFrame
        The input DataFrame with the 'type' column set to 'exon' or 'intron'.

    Raises
    ------
    ValueError
        If 'df_type' is not 'exons' or 'introns'.

    Notes
    -----
    - If the 'type' column does not exist in the input DataFrame, it is added with the specified 'df_type'.
    - If 'df_type' is 'introns', the function filters the DataFrame to include only intron entries.
    """

    # Validate 'df_type' parameter
    if df_type not in ["exons", "introns"]:
        raise ValueError("df_type must be either 'exons' or 'introns'")

    # Add or set the 'type' column
    if 'type' not in df.schema:
        # If 'type' column is missing, add it with the appropriate value
        return df.with_columns(
            pl.lit('exon' if df_type == 'exons' else 'intron').alias('type')
        )
    elif df_type == 'introns':
        # If 'df_type' is 'introns', ensure only intron entries are included
        df = df.filter(pl.col('type') == 'intron')
    return df


def _get_gaps(exons: pl.DataFrame) -> pl.DataFrame:
    """
    Identifies gaps between exons within the same chromosome and strand.

    Parameters
    ----------
    exons : pl.DataFrame
        A Polars DataFrame containing exon information with 'seqnames', 'start', 'end', and 'strand'.

    Returns
    -------
    pl.DataFrame
        A DataFrame with 'start' and 'end' positions of gaps between exons.

    Raises
    ------
    ValueError
        If exons are not from a single chromosome and strand.

    Notes
    -----
    - All exons must be from the same chromosome and strand to accurately identify gaps.
    - The function merges overlapping exons and computes the gaps between them.
    """

    # Ensure all exons are from a single chromosome and strand
    seqnames_unique = exons["seqnames"].n_unique()
    strand_unique = exons["strand"].n_unique()
    if seqnames_unique != 1 or strand_unique != 1:
        raise ValueError("Exons must be from a single chromosome and strand")

    # Sort exons by start position
    exons_sorted = exons.sort('start')

    # Compute cumulative maximum of 'end' shifted by 1 to identify overlaps
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

    return gaps  # Return the DataFrame containing gap positions


def _get_tx_start_gaps(exons: pl.DataFrame, transcript_id_column: str) -> pl.DataFrame:
    """
    Identifies gaps at the start of each transcript based on the first exon.

    Parameters
    ----------
    exons : pl.DataFrame
        A Polars DataFrame containing exon information.
    transcript_id_column : str
        Column used to group transcripts (e.g., 'transcript_id').

    Returns
    -------
    pl.DataFrame
        A DataFrame containing gaps at the start of each transcript.

    Notes
    -----
    - The function calculates the gap between the overall start of the first exon across all transcripts
      and the start of each individual transcript's first exon.
    - It assumes that all exons are on the same chromosome and strand.
    """

    # Get the start of the first exon for each transcript
    tx_starts = exons.group_by(transcript_id_column).agg(pl.col('start').min())

    # Get the overall start of the first exon across all transcripts
    overall_start = exons['start'].min()

    # Use the same chromosome and strand for all transcripts
    seqnames_value = exons['seqnames'][0]
    strand_value = exons['strand'][0]

    # Create DataFrame with gaps at the start of transcripts
    tx_start_gaps = tx_starts.with_columns([
        pl.col('start').cast(pl.Int64).alias('end'),
        pl.lit(overall_start).cast(pl.Int64).alias('start'),
        pl.lit(seqnames_value).alias('seqnames'),
        pl.lit(strand_value).alias('strand'),
    ])

    return tx_start_gaps  # Return the DataFrame with transcript start gaps


def _get_gap_map(df: pl.DataFrame, gaps: pl.DataFrame) -> dict:
    """
    Maps gaps to the corresponding exons or introns based on their positions.

    Parameters
    ----------
    df : pl.DataFrame
        A DataFrame containing exons or introns, with 'start' and 'end' positions.
    gaps : pl.DataFrame
        A DataFrame containing gaps between exons, with 'start' and 'end' positions.

    Returns
    -------
    dict
        A dictionary containing mappings:
        - 'equal': DataFrame of gaps that exactly match the 'start' and 'end' of exons/introns.
        - 'pure_within': DataFrame of gaps that are fully within exons/introns but do not exactly match.

    Notes
    -----
    - The function adds row indices to both df and gaps for mapping.
    - It first identifies exact matches, then finds gaps fully within exons/introns.
    """

    # Add an index to each gap and exon/intron row
    gaps = gaps.with_row_index("gap_index")
    df = df.with_row_index("df_index")

    # Find gaps where the start and end positions exactly match those of df
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

    # Remove within_hits that also appear in equal_hits
    pure_within_hits = within_hits.join(equal_hits, how="anti", on=["df_index", "gap_index"])

    # Sort the equal_hits by gap and df index
    equal_hits = equal_hits.sort(["gap_index", "df_index"])

    # Return the mappings
    return {
        'equal': equal_hits,
        'pure_within': pure_within_hits
    }


def _get_shortened_gaps(df: pl.DataFrame, gaps: pl.DataFrame, gap_map: dict,
                        transcript_id_column: str, target_gap_width: int) -> pl.DataFrame:
    """
    Shortens the gaps between exons or introns based on a target gap width.

    Parameters
    ----------
    df : pl.DataFrame
        A DataFrame containing exons or introns.
    gaps : pl.DataFrame
        A DataFrame containing gaps between exons.
    gap_map : dict
        A dictionary mapping gaps to their corresponding exons or introns.
    transcript_id_column : str
        Column used to group transcripts (e.g., 'transcript_id').
    target_gap_width : int
        The maximum allowed width for the gaps.

    Returns
    -------
    pl.DataFrame
        A DataFrame with shortened gaps and adjusted positions.

    Notes
    -----
    - Gaps classified as 'equal' and exceeding the target width are shortened to match the target.
    - Gaps classified as 'pure_within' are adjusted based on the target width, ensuring they do not exceed the defined maximum.
    - The function updates the 'width' of each gap accordingly and removes unnecessary columns post-adjustment.
    """

    # Calculate the width of exons/introns and initialize a 'shorten_type' column
    df = df.with_columns(
        (pl.col('end') - pl.col('start') + 1).alias('width'),  # Calculate the width
        pl.lit('none').alias('shorten_type')  # Initialize shorten_type column
    )

    # Add an index column to the df DataFrame
    df = df.with_row_index(name="df_index")

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
        gaps = gaps.with_row_index(name="gap_index")

        if len(overlapping_gap_indexes) > 0:
            # Calculate the width of overlapping gaps
            overlapping_gaps = gaps.filter(pl.col("gap_index").is_in(overlapping_gap_indexes))
            overlapping_gaps = overlapping_gaps.with_columns(
                (pl.col('end') - pl.col('start') + 1).alias('gap_width')
            )

            # Shorten gap width if larger than target_gap_width
            overlapping_gaps = overlapping_gaps.with_columns(
                pl.when(pl.col('gap_width') > target_gap_width)
                .then(pl.lit(target_gap_width))
                .otherwise(pl.col('gap_width'))
                .alias('shortened_gap_width')
            )

            # Calculate the gap difference
            overlapping_gaps = overlapping_gaps.with_columns(
                (pl.col('gap_width') - pl.col('shortened_gap_width')).alias('shortened_gap_diff')
            )

            # Map the gap differences back to df
            gap_diff_df = gap_map['pure_within'].join(
                overlapping_gaps.select('gap_index', 'shortened_gap_diff'), on='gap_index', how='left'
            )

            # Aggregate gap differences by df indexes
            sum_gap_diff = gap_diff_df.group_by('df_index').agg(
                pl.sum('shortened_gap_diff').alias('sum_shortened_gap_diff')
            )

            # Join the calculated gap differences with the df DataFrame
            df = df.join(sum_gap_diff, on='df_index', how='left')

            # Adjust the width based on gap differences
            df = df.with_columns(
                pl.when(pl.col('sum_shortened_gap_diff').is_null())
                .then(pl.col('shortened_width'))
                .otherwise(pl.col('width') - pl.col('sum_shortened_gap_diff'))
                .alias('shortened_width')
            )

            # Clean up unnecessary columns
            df = df.drop('sum_shortened_gap_diff')

    df = df.drop(['shorten_type', 'width', 'df_index'])
    df = df.rename({'shortened_width': 'width'})

    return df  # Return the DataFrame with shortened gaps


def _get_rescaled_txs(
    exons: pl.DataFrame,
    introns_shortened: pl.DataFrame,
    tx_start_gaps_shortened: pl.DataFrame,
    transcript_id_column: str
) -> pl.DataFrame:
    """
    Rescales transcript coordinates based on shortened gaps for exons and introns.

    Parameters
    ----------
    exons : pl.DataFrame
        DataFrame containing exon information.
    introns_shortened : pl.DataFrame
        DataFrame containing intron information with shortened gaps.
    tx_start_gaps_shortened : pl.DataFrame
        DataFrame containing rescaled transcript start gaps.
    transcript_id_column : str
        Column used to group transcripts (e.g., 'transcript_id').

    Returns
    -------
    pl.DataFrame
        Rescaled transcript DataFrame with adjusted coordinates.

    Notes
    -----
    - The function concatenates exons and shortened introns, sorts them, and calculates rescaled start and end positions.
    - It adjusts intron positions to prevent overlap with exons.
    - Transcript start gaps are incorporated to ensure accurate rescaling across different transcripts.
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

    # Sort based on transcript_id, start, and end
    rescaled_tx = rescaled_tx.sort([transcript_id_column, 'start', 'end'])

    # Calculate cumulative sum for rescaled end positions
    rescaled_tx = rescaled_tx.with_columns(
        pl.col('width').cum_sum().over(transcript_id_column).alias('rescaled_end')
    )

    # Compute the rescaled start positions based on the cumulative end positions
    rescaled_tx = rescaled_tx.with_columns(
        (pl.col('rescaled_end') - pl.col('width') + 1).alias('rescaled_start')
    )

    # Join rescaled transcript start gaps to adjust start positions
    rescaled_tx = rescaled_tx.join(
        tx_start_gaps_shortened, on=transcript_id_column, how='left', suffix='_tx_start'
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

    # Drop 'width' column as it's no longer needed
    rescaled_tx = rescaled_tx.drop(['width'])

    # Reorder columns for consistency in the output
    columns = rescaled_tx.columns
    column_order = ['seqnames', 'start', 'end', "rescaled_start", "rescaled_end", 'strand'] + [
        col for col in columns if col not in ['seqnames', 'start', 'end', "rescaled_start", "rescaled_end", 'strand']
    ]
    rescaled_tx = rescaled_tx.select(column_order)

    return rescaled_tx  # Return the rescaled transcript coordinates


def _get_cds_exon_difference(gene_exons: pl.DataFrame, gene_cds_regions: pl.DataFrame) -> pl.DataFrame:
    """
    Calculates the absolute differences between the start and end positions of exons and CDS regions.

    Parameters
    ----------
    gene_exons : pl.DataFrame
        DataFrame containing exon regions.
    gene_cds_regions : pl.DataFrame
        DataFrame containing CDS (Coding DNA Sequence) regions.

    Returns
    -------
    pl.DataFrame
        DataFrame with the absolute differences between exon and CDS start/end positions.

    Raises
    ------
    ValueError
        If there are no common columns to join on between CDS and exons.

    Notes
    -----
    - The function joins CDS and exon DataFrames on common columns (e.g., 'transcript_id') to align corresponding regions.
    - It calculates the absolute differences between exon and CDS start and end positions to identify discrepancies.
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

    return cds_exon_diff  # Return the DataFrame with differences


def _get_rescale_cds(cds_exon_diff: pl.DataFrame, gene_rescaled_exons: pl.DataFrame) -> pl.DataFrame:
    """
    Rescales CDS regions based on exon positions and the calculated differences between them.

    Parameters
    ----------
    cds_exon_diff : pl.DataFrame
        DataFrame with differences between exon and CDS start/end positions.
    gene_rescaled_exons : pl.DataFrame
        DataFrame containing rescaled exon positions.

    Returns
    -------
    pl.DataFrame
        Rescaled CDS positions based on exon positions.

    Raises
    ------
    ValueError
        If there are no common columns to join on between CDS differences and rescaled exons.

    Notes
    -----
    - The function adjusts CDS start and end positions based on the rescaled exon positions and the previously calculated differences.
    - It ensures that CDS regions are accurately positioned relative to exons after rescaling.
    """

    # Assign a 'type' column with the value "CDS" and drop unnecessary columns
    columns_to_drop = ['exon_start', 'exon_end']
    cds_prepared = (
        cds_exon_diff
        .with_columns(pl.lit("CDS").alias("type"))
        .drop([col for col in columns_to_drop if col in cds_exon_diff.columns])
    )

    # Rename columns in rescaled exons for consistency
    exons_prepared = gene_rescaled_exons.rename({'rescaled_start': 'exon_start', 'rescaled_end': 'exon_end'})
    exons_prepared = exons_prepared.drop(["start", "end"])

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
        (pl.col('exon_start') + pl.col('diff_start')).alias('rescaled_start'),
        (pl.col('exon_end') - pl.col('diff_end')).alias('rescaled_end')
    ])

    # Drop unnecessary columns used for the difference calculations
    gene_rescaled_cds = gene_rescaled_cds.drop(['exon_start', 'exon_end', 'diff_start', 'diff_end'])

    # Rename CDS start and end to 'start' and 'end'
    gene_rescaled_cds = gene_rescaled_cds.rename({
        "cds_start": "start",
        "cds_end": "end"
    })

    return gene_rescaled_cds  # Return the rescaled CDS DataFrame
