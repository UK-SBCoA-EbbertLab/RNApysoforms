import polars as pl
from RNApysoforms.utils import check_df

def calculate_exon_number(annotation: pl.DataFrame, group_var: str = "transcript_id") -> pl.DataFrame:
    """
    Assigns exon numbers to exons, CDS, and introns within a genomic annotation dataset based on transcript structure.

    This function processes an annotation dataset, identifies exons, CDS, and introns, and assigns exon numbers to 
    these features based on the transcript's structure and strand direction. Exons are assigned numbers sequentially,
    while CDS and introns inherit exon numbers based on their overlap or position relative to exons.

    Parameters
    ----------
    annotation : pl.DataFrame
        A Polars DataFrame containing the genomic annotation data. Must include columns 'start', 'end', 'type', 
        'strand', and the grouping variable (default 'transcript_id').
    group_var : str, optional
        The column name that identifies the transcript groups (default is 'transcript_id').

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame with the original annotation data, now including a new 'exon_number' column. This column 
        assigns exon numbers to exons, CDS, and introns based on their relationships within the transcript structure.

    Raises
    ------
    ValueError
        If the input annotation is not a Polars DataFrame or if required columns are missing from the DataFrame.

    Examples
    --------
    Assign exon numbers to genomic features:

    >>> import polars as pl
    >>> from RNApysoforms.annotation import calculate_exon_number
    >>> df = pl.DataFrame({
    ...     "transcript_id": ["tx1", "tx1", "tx1", "tx2", "tx2"],
    ...     "start": [100, 200, 300, 150, 250],
    ...     "end": [150, 250, 350, 200, 300],
    ...     "type": ["exon", "intron", "CDS", "exon", "intron"],
    ...     "strand": ["+", "+", "+", "-", "-"]
    ... })
    >>> result_df = calculate_exon_number(df)
    >>> result_df.show()

    This will output a DataFrame where exons, CDS, and introns are numbered according to their order and relationships 
    within each transcript.

    Notes
    -----
    - Exons are numbered based on their relative positions within the transcript, accounting for strand orientation 
      (ascending for positive strand and descending for negative strand).
    - CDS regions are assigned the exon number of the overlapping exon.
    - Introns are assigned the exon number of the adjacent exon based on the strand direction (preceding exon for positive 
      strand and following exon for negative strand).
    
    Steps
    -----
    1. The function first checks the input DataFrame for required columns and verifies its type.
    2. Exons are extracted and assigned numbers sequentially based on their position within the transcript and strand.
    3. CDS regions are processed next, inheriting the exon number from overlapping exons.
    4. Introns are assigned exon numbers based on the nearest adjacent exon (preceding exon for '+' strand, following 
       exon for '-' strand).
    5. The processed annotations (exons, CDS, and introns) are then combined into a single DataFrame with exon numbers assigned.

    """

    # Ensure 'annotation' is a Polars DataFrame
    if not isinstance(annotation, pl.DataFrame):
        raise TypeError(f"Expected annotation to be of type pl.DataFrame, got {type(annotation)}" +
                        "\n You can use polars_df = pandas_df.from_pandas() to convert a pandas df into a polars df")
    
    # Ensure required columns are present
    check_df(annotation, ["start", "end", group_var, "type", "strand"])
    
    # Step 1: Extract exons and assign exon_numbers
    exon_annotation = annotation.filter(pl.col('type') == 'exon')
    
    exon_annotation = exon_annotation.with_columns(
        pl.when(pl.col('strand') == '+')
        .then(pl.col('start').rank('dense').over(group_var).cast(pl.Int64))
        .otherwise(pl.col('end').rank('dense', descending=True).over(group_var).cast(pl.Int64))
        .alias('exon_number')
    )
    
    # Step 2: Assign exon numbers to CDS entries based on overlapping exons
    cds_annotation = annotation.filter(pl.col('type') == 'CDS')
    
    if not cds_annotation.is_empty():
        cds_exon_annotation = cds_annotation.join(
            exon_annotation.select([group_var, 'start', 'end', 'exon_number']),
            on=group_var,
            how='left',
            suffix='_exon'
        ).filter(
            (pl.col('start') <= pl.col('end_exon')) & (pl.col('end') >= pl.col('start_exon'))
        )
        # Assign the minimum exon_number (in case of multiple overlaps)
        cds_exon_annotation = cds_exon_annotation.group_by([group_var, 'start', 'end', 'strand', 'type']).agg(
            pl.col('exon_number').min().alias('exon_number')
        )
    else:
        cds_exon_annotation = cds_annotation.with_columns(pl.lit(None).alias('exon_number'))
    
    # Step 3: Assign exon numbers to intron entries based on the exon behind them
    intron_annotation = annotation.filter(pl.col('type') == 'intron')
    intron_exon_annotation_list = []
    
    for strand in ['+', '-']:
        strand_introns = intron_annotation.filter(pl.col('strand') == strand)
        strand_exons = exon_annotation.filter(pl.col('strand') == strand)
        
        if not strand_introns.is_empty():
            intron_exon_annotation = strand_introns.join(
                strand_exons.select([group_var, 'start', 'end', 'exon_number']),
                on=group_var,
                how='left',
                suffix='_exon'
            )
            if strand == '+':
                intron_exon_annotation = intron_exon_annotation.filter(pl.col('end_exon') <= pl.col('start'))
                intron_exon_annotation = intron_exon_annotation.group_by([group_var, 'start', 'end', 'strand', 'type']).agg(
                    pl.col('exon_number').filter(pl.col('end_exon') == pl.col('end_exon').max()).first().alias('exon_number')
                )
            else:
                intron_exon_annotation = intron_exon_annotation.filter(pl.col('start_exon') >= pl.col('end'))
                intron_exon_annotation = intron_exon_annotation.group_by([group_var, 'start', 'end', 'strand', 'type']).agg(
                    pl.col('exon_number').filter(pl.col('start_exon') == pl.col('start_exon').min()).first().alias('exon_number')
                )
            intron_exon_annotation_list.append(intron_exon_annotation)
    if intron_exon_annotation_list:
        intron_exon_annotation = pl.concat(intron_exon_annotation_list)
    else:
        intron_exon_annotation = intron_annotation.with_columns(pl.lit(None).alias('exon_number'))
    
    # Combine all annotations
    result_annotation = pl.concat([
        exon_annotation.select([group_var, 'start', 'end', 'strand', 'type', 'exon_number']),
        cds_exon_annotation.select([group_var, 'start', 'end', 'strand', 'type', 'exon_number']),
        intron_exon_annotation.select([group_var, 'start', 'end', 'strand', 'type', 'exon_number'])
    ]).sort([group_var, 'start'])
    
    return result_annotation
