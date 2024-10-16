import polars as pl
from RNApysoforms.utils import check_df

def calculate_exon_number(annotation: pl.DataFrame, transcript_id_column: str = "transcript_id") -> pl.DataFrame:
    """
    Assigns exon numbers to exons, CDS, and introns within a genomic annotation dataset based on transcript structure and strand direction.

    This function processes a genomic annotation DataFrame to assign exon numbers to genomic features—specifically exons,
    CDS (Coding Sequences), and introns—based on their position within each transcript and their strand orientation.
    Exons are numbered sequentially within each transcript, accounting for the strand direction: numbering increases from
    the 5' to 3' end on the positive strand and decreases on the negative strand. CDS and introns are assigned exon numbers
    based on their overlap or adjacency to exons.

    **Required Columns in `annotation`:**
    - `start`: Start position of the feature.
    - `end`: End position of the feature.
    - `strand`: Strand direction of the feature (`"+"` or `"-"`).
    - `type`: Type of the feature (must include `"exon"`. Can also include `"CDS"`, and/or `"intron"`).
    - `transcript_id_column` (default `"transcript_id"`): Identifier for grouping features into transcripts.

    Parameters
    ----------
    annotation : pl.DataFrame
        A Polars DataFrame containing genomic annotation data. Must include columns for start and end positions,
        feature type, strand direction, and a grouping variable (default is 'transcript_id'). If a different
        grouping variable is used, specify it using the `transcript_id_column` parameter.
    transcript_id_column : str, optional
        The column name that identifies transcript groups within the DataFrame, by default "transcript_id".

    Returns
    -------
    pl.DataFrame
        A Polars DataFrame that includes all original annotation data along with a new 'exon_number' column.
        This column assigns exon numbers to exons, CDS, and introns based on their order and relationships within
        each transcript.

    Raises
    ------
    TypeError
        If the `annotation` parameter is not a Polars DataFrame.
    ValueError
        If required columns are missing from the `annotation` DataFrame based on the provided parameters.

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
    >>> print(result_df)

    This will output a DataFrame where exons, CDS, and introns are numbered according to their order and relationships
    within each transcript.

    Notes
    -----
    - **Exon Numbering:**
        - For transcripts on the positive strand (`"+"`), exons are numbered in ascending order based on their start positions.
        - For transcripts on the negative strand (`"-"`), exons are numbered in descending order based on their end positions.
    - **CDS Assignment:**
        - CDS regions inherit the exon number of the overlapping exon. CDS regions must be within the boundaries of a 
        single exon, otherwise the function might return erroneous results.
    - **Intron Assignment:**
        - Introns are assigned the exon number based on the order the introns show up in the transcript.
        the first intron will be assigned exon_number 1, so forth and so on.
    - **Data Integrity:**
        - The function ensures that all required columns are present and correctly formatted before processing.
        - If no CDS or introns are present in the data, the function handles these cases gracefully without errors.
    - **Performance:**
        - Utilizes Polars' efficient data manipulation capabilities to handle large genomic datasets effectively.

    """

    # Ensure 'annotation' is a Polars DataFrame
    if not isinstance(annotation, pl.DataFrame):
        raise TypeError(
            f"Expected 'annotation' to be of type pl.DataFrame, got {type(annotation)}. "
            "You can convert a pandas DataFrame to Polars using pl.from_pandas(pandas_df)."
        )

    # Ensure required columns are present in the DataFrame
    required_columns = [transcript_id_column, "start", "end", "type", "strand"]
    check_df(annotation, required_columns)

    # Get original column order
    column_order = annotation.columns
    column_order = column_order + ["exon_number"]

    # Step 1: Extract exons and assign exon numbers based on strand direction
    # Filter the annotation DataFrame to include only 'exon' type features
    exon_annotation = annotation.filter(pl.col('type') == 'exon')

    # Assign exon numbers within each transcript group, accounting for strand direction
    exon_annotation = exon_annotation.with_columns(
        pl.when(pl.col('strand') == '+')
        .then(
            # For positive strand, rank exons in ascending order based on 'start' position
            pl.col('start').rank('dense').over(transcript_id_column).cast(pl.Int64)
        )
        .otherwise(
            # For negative strand, rank exons in descending order based on 'end' position
            pl.col('end').rank('dense', descending=True).over(transcript_id_column).cast(pl.Int64)
        )
        .alias('exon_number')
    )

    # Step 2: Assign exon numbers to CDS entries based on overlapping exons
    # Filter the annotation DataFrame to include only 'CDS' type features
    cds_annotation = annotation.filter(pl.col('type') == 'CDS')

    if not cds_annotation.is_empty():
        # If there are CDS entries, proceed to assign exon numbers
        # Join CDS with exons on the transcript ID to find overlapping exons
        cds_exon_annotation = cds_annotation.join(
            exon_annotation.select([transcript_id_column, 'start', 'end', 'exon_number']),
            on=transcript_id_column,
            how='left',
            suffix='_exon'
        ).filter(
            # Keep only CDS and exon pairs that overlap
            (pl.col('start') <= pl.col('end_exon')) & (pl.col('end') >= pl.col('start_exon'))
        )
        # For each CDS, assign the minimum exon_number among overlapping exons
        cds_exon_annotation = cds_exon_annotation.group_by([transcript_id_column, 'start', 'end', 'strand', 'type']).agg(
            pl.col('exon_number').min().alias('exon_number')
        )
    else:
        # If no CDS entries, create an empty DataFrame with 'exon_number' as None
        cds_exon_annotation = cds_annotation.with_columns(
            pl.lit(None).cast(pl.Int64).alias('exon_number')
        )

    # Step 3: Assign exon numbers to intron entries based on adjacent exons
    # Filter the annotation DataFrame to include only 'intron' type features
    intron_annotation = annotation.filter(pl.col('type') == 'intron')

    # Initialize a list to collect intron annotations with assigned exon numbers
    intron_exon_annotation_list = []

    # Process introns separately for positive and negative strands
    for strand in ['+', '-']:
        # Filter introns and exons by strand
        strand_introns = intron_annotation.filter(pl.col('strand') == strand)
        strand_exons = exon_annotation.filter(pl.col('strand') == strand)

        if not strand_introns.is_empty():
            # Join introns with exons on the transcript ID to find adjacent exons
            intron_exon_annotation = strand_introns.join(
                strand_exons.select([transcript_id_column, 'start', 'end', 'exon_number']),
                on=transcript_id_column,
                how='left',
                suffix='_exon'
            )
            if strand == '+':
                # For positive strand, find exons where exon 'end' is less than or equal to intron 'start'
                intron_exon_annotation = intron_exon_annotation.filter(pl.col('end_exon') <= pl.col('start'))
                # For each intron, assign the exon_number of the preceding exon (max 'end_exon')
                intron_exon_annotation = intron_exon_annotation.group_by([transcript_id_column, 'start', 'end', 'strand', 'type']).agg(
                    pl.col('exon_number').filter(pl.col('end_exon') == pl.col('end_exon').max()).first().alias('exon_number')
                )
            else:
                # For negative strand, find exons where exon 'start' is greater than or equal to intron 'end'
                intron_exon_annotation = intron_exon_annotation.filter(pl.col('start_exon') >= pl.col('end'))
                # For each intron, assign the exon_number of the following exon (min 'start_exon')
                intron_exon_annotation = intron_exon_annotation.group_by([transcript_id_column, 'start', 'end', 'strand', 'type']).agg(
                    pl.col('exon_number').filter(pl.col('start_exon') == pl.col('start_exon').min()).first().alias('exon_number')
                )
            # Add the processed intron annotations to the list
            intron_exon_annotation_list.append(intron_exon_annotation)

    if not intron_annotation.is_empty():
        # Concatenate all intron annotations with assigned exon numbers
        intron_exon_annotation = pl.concat(intron_exon_annotation_list)
    else:
        # If no intron entries, create an empty DataFrame with 'exon_number' as None
        intron_exon_annotation = intron_annotation.with_columns(
            pl.lit(None).cast(pl.Int64).alias('exon_number')
        )

    # Combine exons, CDS, and introns with their assigned exon numbers
    result_annotation = pl.concat([
        exon_annotation.select([transcript_id_column, 'start', 'end', 'strand', 'type', 'exon_number']),
        cds_exon_annotation.select([transcript_id_column, 'start', 'end', 'strand', 'type', 'exon_number']),
        intron_exon_annotation.select([transcript_id_column, 'start', 'end', 'strand', 'type', 'exon_number'])
    ])

    ## Only do this part if annotation is not empty
    if not result_annotation.is_empty():

        ## Return all the original columns
        result_annotation = result_annotation.join(annotation, on=required_columns, how="inner")
        result_annotation = result_annotation[column_order]
        result_annotation = result_annotation.sort([transcript_id_column, 'start'])

    return result_annotation  # Return the final DataFrame with assigned exon numbers