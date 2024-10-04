import polars as pl
import warnings

def gene_filtering(
    gene_name_to_filter: str,
    annotation: pl.DataFrame,
    counts_matrix: pl.DataFrame = None,
    transcript_feature_column: str = "transcript_id"
):
    """
    Filters genomic annotations and, optionally, a counts matrix based on a specific gene name.

    This function filters the provided annotation DataFrame to include only the specified gene, identified by the `gene_name`.
    If a counts matrix is provided, it will also filter the counts matrix to retain only the transcripts corresponding to the
    filtered gene. The function checks for missing transcripts between the annotation and counts matrix and issues warnings
    if discrepancies are found.

    Parameters
    ----------
    gene_name_to_filter : str
        The gene name to filter in the annotation DataFrame.
    annotation : pl.DataFrame
        A Polars DataFrame containing genomic annotations. Must include 'gene_name' and 'transcript_id' columns.
    counts_matrix : pl.DataFrame, optional
        A Polars DataFrame containing transcript counts data. If provided, it will be filtered to match the filtered annotation.
        Default is None.
    transcript_feature_column : str, optional
        The column name in the counts matrix representing transcript IDs. Default is 'transcript_id'.

    Returns
    -------
    tuple
        If `counts_matrix` is provided, returns a tuple of (filtered_annotation, filtered_counts_matrix).
        If `counts_matrix` is None, returns only the `filtered_annotation`.

    Raises
    ------
    ValueError
        If required columns ('gene_name', 'transcript_id') are missing in the annotation DataFrame, or if the counts matrix
        does not contain the `transcript_feature_column`.
    ValueError
        If the filtered counts matrix is empty after filtering.
    Warning
        If there are discrepancies between transcripts in the annotation and counts matrix.

    Examples
    --------
    Filter an annotation DataFrame by a specific gene:

    >>> import polars as pl
    >>> annotation_df = pl.DataFrame({
    ...     "gene_name": ["BRCA1", "BRCA1", "TP53"],
    ...     "transcript_id": ["tx1", "tx2", "tx3"]
    ... })
    >>> filtered_annotation = gene_filtering("BRCA1", annotation_df)
    >>> print(filtered_annotation)

    Filter an annotation and counts matrix by a specific gene:

    >>> counts_matrix_df = pl.DataFrame({
    ...     "transcript_id": ["tx1", "tx2", "tx4"],
    ...     "sample1": [100, 200, 300],
    ...     "sample2": [150, 250, 350]
    ... })
    >>> filtered_annotation, filtered_counts_matrix = gene_filtering(
    ...     "BRCA1", annotation_df, counts_matrix=counts_matrix_df
    ... )
    >>> print(filtered_annotation)
    >>> print(filtered_counts_matrix)

    Notes
    -----
    - The function assumes the `annotation` DataFrame contains the columns 'gene_name' and 'transcript_id'.
    - If a `counts_matrix` is provided, the function checks for discrepancies between transcripts in the annotation and
      counts matrix and raises warnings if there are differences.
    - If no matching transcripts are found after filtering the counts matrix, the function raises an error.
    """


    # Check if annotation is a polars DataFrame
    if not isinstance(annotation, pl.DataFrame):
        raise ValueError("Annotation must be a polars DataFrame.")

    # Check if annotation has 'gene_name' and 'transcript_id' columns
    required_columns = {"gene_name", "transcript_id"}
    missing_columns = required_columns - set(annotation.columns)
    if missing_columns:
        raise ValueError(
            f"Annotation dataframe must contain the following columns: {', '.join(missing_columns)}."
        )

    # Filter annotation based on 'gene_name'
    filtered_annotation = annotation.filter(pl.col("gene_name") == gene_name_to_filter)

    if counts_matrix is not None:
        # Check if counts_matrix is a polars DataFrame
        if not isinstance(counts_matrix, pl.DataFrame):
            raise ValueError("Counts matrix must be a polars DataFrame.")

        # Check if counts_matrix has the transcript_feature_column
        if transcript_feature_column not in counts_matrix.columns:
            raise ValueError(
                f"Counts matrix must contain the '{transcript_feature_column}' column."
            )

        # Filter counts_matrix based on transcripts in filtered_annotation
        filtered_counts_matrix = counts_matrix.filter(
            pl.col(transcript_feature_column).is_in(filtered_annotation["transcript_id"])
        )

        # If filtered counts matrix is empty, throw informative error
        if filtered_counts_matrix.is_empty():
            raise ValueError(
                "Filtered counts matrix is empty after filtering. No matching transcripts found."
            )

        # Check for discrepancies between transcripts in annotation and counts_matrix
        annotation_transcripts = set(filtered_annotation["transcript_id"].unique())
        counts_transcripts = set(filtered_counts_matrix[transcript_feature_column].unique())

        missing_transcripts = annotation_transcripts - counts_transcripts
        if missing_transcripts:
            warnings.warn(
                f"{len(missing_transcripts)} transcripts are named differently in counts matrix. "
                f"Missing transcripts: {', '.join(missing_transcripts)}"
            )

        return filtered_annotation, filtered_counts_matrix
    else:
        return filtered_annotation
