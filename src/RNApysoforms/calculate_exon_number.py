import polars as pl
from RNApysoforms.utils import check_df

def calculate_exon_number(annotation: pl.DataFrame, transcript_id_column: str = "transcript_id") -> pl.DataFrame:
    """
    Assigns exon numbers to exons, CDS, and introns within a genomic annotation dataset based on transcript structure.
    
    This function processes a genomic annotation dataset to identify exons, CDS (coding sequences), and introns, 
    assigning exon numbers to these features based on their position within each transcript and the strand direction. 
    Exons are numbered sequentially, while CDS and introns inherit exon numbers based on their overlap or 
    positional relationship to exons. The numbering accounts for strand orientation, ensuring accurate representation 
    of transcript structure.
    
    **Required Columns in `annotation`:**
    - `start`: Start position of the feature.
    - `end`: End position of the feature.
    - `transcript_id_column` (default `"transcript_id"`): Identifier for grouping features into transcripts.
    - `type`: Type of the feature (`"exon"`, `"CDS"`, or `"intron"`).
    - `strand`: Strand direction of the feature (`"+"` or `"-"`).
    
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
        A Polars DataFrame that includes all original annotation data along with a new 'exon_number' column. This 
        column assigns exon numbers to exons, CDS, and introns based on their order and relationships within 
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
    >>> from RNApysoforms.utils import check_df
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
    - **Exon Numbering:**
        - For transcripts on the positive strand (`"+"`), exons are numbered in ascending order based on their start positions.
        - For transcripts on the negative strand (`"-"`), exons are numbered in descending order based on their end positions.
    - **CDS Assignment:**
        - CDS regions inherit the exon number of overlapping exons. If a CDS overlaps multiple exons, it is assigned 
          the smallest exon number among the overlapping exons.
    - **Intron Assignment:**
        - Introns are assigned the exon number of the adjacent exon based on strand direction:
            - Positive strand (`"+"`): Introns inherit the exon number of the preceding exon.
            - Negative strand (`"-"`): Introns inherit the exon number of the following exon.
    - **Data Integrity:**
        - The function ensures that all required columns are present and correctly formatted before processing.
        - If no CDS or introns are present in the data, the function handles these cases gracefully without errors.
    - **Performance:**
        - Utilizes Polars' efficient data manipulation capabilities to handle large genomic datasets effectively.
    
    Steps
    -----
    1. **Input Validation:**
        - Verifies that the input `annotation` is a Polars DataFrame.
        - Checks for the presence of required columns: 'start', 'end', the grouping variable (`transcript_id_column`), 'type', and 'strand'.
    2. **Exon Extraction and Numbering:**
        - Filters the DataFrame to extract exon entries.
        - Assigns exon numbers based on strand direction:
            - Positive strand: Exons are ranked in ascending order of their start positions.
            - Negative strand: Exons are ranked in descending order of their end positions.
    3. **CDS Numbering:**
        - Filters the DataFrame to extract CDS entries.
        - Joins CDS entries with exon annotations to identify overlapping exons.
        - Assigns the smallest overlapping exon number to each CDS.
        - If no CDS entries are present, assigns `None` to exon numbers for CDS.
    4. **Intron Numbering:**
        - Filters the DataFrame to extract intron entries.
        - For each strand (`"+"` and `"-"`), joins introns with exons to find adjacent exons.
        - Assigns exon numbers to introns based on the nearest adjacent exon:
            - Positive strand: Inherits the exon number of the preceding exon.
            - Negative strand: Inherits the exon number of the following exon.
        - If no intron entries are present, assigns `None` to exon numbers for introns.
    5. **Combining Annotations:**
        - Concatenates exon, CDS, and intron annotations with their assigned exon numbers.
        - Sorts the combined DataFrame by the grouping variable and start positions to maintain logical order.
    6. **Output:**
        - Returns the final DataFrame with exon numbers assigned to all relevant features.
    
    """
    
    # Ensure 'annotation' is a Polars DataFrame
    if not isinstance(annotation, pl.DataFrame):
        raise TypeError(
            f"Expected 'annotation' to be of type pl.DataFrame, got {type(annotation)}. "
            "You can convert a pandas DataFrame to Polars using pl.from_pandas(pandas_df)."
        )
    
    # Ensure required columns are present in the DataFrame
    required_columns = ["start", "end", transcript_id_column, "type", "strand"]
    check_df(annotation, required_columns)
    
    # Step 1: Extract exons and assign exon numbers based on strand direction
    exon_annotation = annotation.filter(pl.col('type') == 'exon')
    
    exon_annotation = exon_annotation.with_columns(
        pl.when(pl.col('strand') == '+')
        .then(pl.col('start').rank('dense').over(transcript_id_column).cast(pl.Int64))
        .otherwise(pl.col('end').rank('dense', descending=True).over(transcript_id_column).cast(pl.Int64))
        .alias('exon_number')
    )
    
    # Step 2: Assign exon numbers to CDS entries based on overlapping exons
    cds_annotation = annotation.filter(pl.col('type') == 'CDS')
    
    if not cds_annotation.is_empty():
        # Join CDS with exons to find overlapping exons
        cds_exon_annotation = cds_annotation.join(
            exon_annotation.select([transcript_id_column, 'start', 'end', 'exon_number']),
            on=transcript_id_column,
            how='left',
            suffix='_exon'
        ).filter(
            (pl.col('start') <= pl.col('end_exon')) & (pl.col('end') >= pl.col('start_exon'))
        )
        # Assign the minimum exon_number in case of multiple overlaps
        cds_exon_annotation = cds_exon_annotation.group_by([transcript_id_column, 'start', 'end', 'strand', 'type']).agg(
            pl.col('exon_number').min().alias('exon_number')
        )
    else:
        # If no CDS entries, assign None to exon_number
        cds_exon_annotation = cds_annotation.with_columns(pl.lit(None).cast(pl.Int64).alias('exon_number'))
    
    # Step 3: Assign exon numbers to intron entries based on adjacent exons
    intron_annotation = annotation.filter(pl.col('type') == 'intron')
    intron_exon_annotation_list = []
    
    for strand in ['+', '-']:
        # Filter introns and exons by strand
        strand_introns = intron_annotation.filter(pl.col('strand') == strand)
        strand_exons = exon_annotation.filter(pl.col('strand') == strand)
        
        if not strand_introns.is_empty():
            # Join introns with exons to find adjacent exons
            intron_exon_annotation = strand_introns.join(
                strand_exons.select([transcript_id_column, 'start', 'end', 'exon_number']),
                on=transcript_id_column,
                how='left',
                suffix='_exon'
            )
            if strand == '+':
                # For positive strand, introns inherit exon_number of preceding exon
                intron_exon_annotation = intron_exon_annotation.filter(pl.col('end_exon') <= pl.col('start'))
                intron_exon_annotation = intron_exon_annotation.group_by([transcript_id_column, 'start', 'end', 'strand', 'type']).agg(
                    pl.col('exon_number').filter(pl.col('end_exon') == pl.col('end_exon').max()).first().alias('exon_number')
                )
            else:
                # For negative strand, introns inherit exon_number of following exon
                intron_exon_annotation = intron_exon_annotation.filter(pl.col('start_exon') >= pl.col('end'))
                intron_exon_annotation = intron_exon_annotation.group_by([transcript_id_column, 'start', 'end', 'strand', 'type']).agg(
                    pl.col('exon_number').filter(pl.col('start_exon') == pl.col('start_exon').min()).first().alias('exon_number')
                )
            intron_exon_annotation_list.append(intron_exon_annotation)
    
    if intron_exon_annotation_list:
        # Concatenate all intron annotations with assigned exon numbers
        intron_exon_annotation = pl.concat(intron_exon_annotation_list)
    else:
        # If no intron entries, assign None to exon_number
        intron_exon_annotation = intron_annotation.with_columns(pl.lit(None).cast(pl.Int64).alias('exon_number'))
    
    # Combine exons, CDS, and introns with their assigned exon numbers
    result_annotation = pl.concat([
        exon_annotation.select([transcript_id_column, 'start', 'end', 'strand', 'type', 'exon_number']),
        cds_exon_annotation.select([transcript_id_column, 'start', 'end', 'strand', 'type', 'exon_number']),
        intron_exon_annotation.select([transcript_id_column, 'start', 'end', 'strand', 'type', 'exon_number'])
    ]).sort([transcript_id_column, 'start'])
    
    return result_annotation
