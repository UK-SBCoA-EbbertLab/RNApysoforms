import polars as pl
import warnings
from typing import Union
from RNApysoforms.utils import check_df

def gene_filtering(
    target_gene: str,
    annotation: pl.DataFrame,
    expression_matrix: pl.DataFrame = None,
    transcript_id_column: str = "transcript_id",
    gene_id_column: str = "gene_name",
    expression_column: str = "counts",
    order_by_expression: bool = True,
    keep_top_expressed_transcripts: Union[str, int] = "all"
) -> Union[pl.DataFrame, tuple]:
    """
    Filters genomic annotations and, optionally, an expression matrix based on a specific gene identifier,
    with options to order and select top expressed transcripts.

    This function filters the provided annotation DataFrame to include only entries corresponding to the specified
    `target_gene`, identified using the column specified by `gene_id_column`. If an expression matrix is provided,
    it will also be filtered to retain only the entries corresponding to the filtered transcripts based on the
    `transcript_id_column`. Additionally, it provides options to order transcripts by their total expression levels
    and to keep only the top expressed transcripts, specified by `keep_top_expressed_transcripts`.

    **Required Columns in `annotation` DataFrame:**
    - `gene_id_column` (default `"gene_name"`): Column containing gene identifiers used for filtering.
    - `transcript_id_column` (default `"transcript_id"`): Column containing transcript identifiers.

    **Required Columns in `expression_matrix` DataFrame (if provided):**
    - `transcript_id_column` (same as in `annotation`): Column containing transcript identifiers matching those in `annotation`.
    - `expression_column` (default `"counts"`): Column containing expression values used for ordering and filtering.

    Parameters
    ----------
    target_gene : str
        The gene identifier to filter in the annotation DataFrame.
    annotation : pl.DataFrame
        A Polars DataFrame containing genomic annotations. Must include the columns specified by `gene_id_column`
        and `transcript_id_column`.
    expression_matrix : pl.DataFrame, optional
        A Polars DataFrame containing expression data. If provided, it will be filtered to match the filtered
        annotation based on `transcript_id_column`. Default is None.
    transcript_id_column : str, optional
        The column name representing transcript identifiers in both the annotation and expression matrix.
        Default is 'transcript_id'.
    gene_id_column : str, optional
        The column name in the annotation DataFrame that contains gene identifiers used for filtering.
        Default is 'gene_name'.
    expression_column : str, optional
        The column name in the expression matrix that contains expression values used for ordering and filtering.
        Default is 'counts'.
    order_by_expression : bool, optional
        If True, transcripts will be ordered by their total expression levels in descending order.
        Default is True.
    keep_top_expressed_transcripts : Union[str, int], optional
        Determines the number of top expressed transcripts to keep after ordering by expression levels.
        Can be 'all' to keep all transcripts or an integer to keep the top N transcripts. Default is 'all'.

    Returns
    -------
    pl.DataFrame or tuple
        - If `expression_matrix` is provided, returns a tuple of (filtered_annotation, filtered_expression_matrix).
        - If `expression_matrix` is None, returns only the `filtered_annotation`.

    Raises
    ------
    TypeError
        If `annotation` or `expression_matrix` are not Polars DataFrames.
    ValueError
        If required columns are missing in the `annotation` or `expression_matrix` DataFrames.
    ValueError
        If the filtered expression matrix is empty after filtering.
    ValueError
        If `keep_top_expressed_transcripts` is not 'all' or a positive integer.
    Warning
        If there are transcripts present in the annotation but missing in the expression matrix.

    Examples
    --------
    Filter an annotation DataFrame by a specific gene:

    >>> import polars as pl
    >>> from RNApysoforms.annotation import gene_filtering
    >>> annotation_df = pl.DataFrame({
    ...     "gene_name": ["BRCA1", "BRCA1", "TP53"],
    ...     "transcript_id": ["tx1", "tx2", "tx3"],
    ...     "counts": [100, 200, 150]
    ... })
    >>> filtered_annotation = gene_filtering("BRCA1", annotation_df)
    >>> print(filtered_annotation)

    Filter an annotation and expression matrix by a specific gene, keeping top 1 expressed transcript:

    >>> expression_matrix_df = pl.DataFrame({
    ...     "transcript_id": ["tx1", "tx2", "tx4"],
    ...     "counts": [100, 200, 300],
    ...     "sample_id": ["sample1", "sample2", "sample3"]
    ... })
    >>> filtered_annotation, filtered_expression_matrix = gene_filtering(
    ...     "BRCA1",
    ...     annotation_df,
    ...     expression_matrix=expression_matrix_df,
    ...     expression_column="counts",
    ...     order_by_expression=True,
    ...     keep_top_expressed_transcripts=1
    ... )
    >>> print(filtered_annotation)
    >>> print(filtered_expression_matrix)

    Notes
    -----
    - The function filters the `annotation` DataFrame to include only entries where `gene_id_column` matches `target_gene`.
    - If an `expression_matrix` is provided, the function filters it to include only transcripts present in the filtered annotation.
    - The function checks for transcripts present in the annotation but missing in the expression matrix and issues a warning for such discrepancies.
    - If `order_by_expression` is True, transcripts are ordered by their total expression levels computed from the `expression_column` in the expression matrix.
    - If `keep_top_expressed_transcripts` is an integer, only the top N expressed transcripts are kept after ordering.
    - If `keep_top_expressed_transcripts` is 'all', all transcripts are kept.
    - If transcripts are present in the expression matrix but not in the annotation, they are silently ignored, and only overlapping transcripts are returned without a warning.

    """
    # Check if 'annotation' is a Polars DataFrame
    if not isinstance(annotation, pl.DataFrame):
        raise TypeError(
            f"Expected 'annotation' to be of type pl.DataFrame, got {type(annotation)}."
            "\nYou can convert a pandas DataFrame to Polars using pl.from_pandas(pandas_df)."
        )

    # Ensure required columns are present in the annotation DataFrame
    check_df(annotation, [gene_id_column, transcript_id_column])

    # Filter annotation based on 'target_gene'
    filtered_annotation = annotation.filter(pl.col(gene_id_column) == target_gene)

    # Check if filtered_annotation is empty and raise an error if true
    if filtered_annotation.is_empty():
        raise ValueError(f"No annotation found for gene: {target_gene} in the '{gene_id_column}' column")

    if expression_matrix is not None:
        # Check if 'expression_matrix' is a Polars DataFrame
        if not isinstance(expression_matrix, pl.DataFrame):
            raise TypeError(
                f"Expected 'expression_matrix' to be of type pl.DataFrame, got {type(expression_matrix)}."
                "\nYou can convert a pandas DataFrame to Polars using pl.from_pandas(pandas_df)."
            )

        # Ensure required columns are present in the expression matrix
        check_df(expression_matrix, [transcript_id_column, expression_column])

        # Filter expression_matrix to include only transcripts present in filtered_annotation
        filtered_expression_matrix = expression_matrix.filter(
            pl.col(transcript_id_column).is_in(filtered_annotation[transcript_id_column])
        )

        # If filtered expression matrix is empty, raise an error
        if filtered_expression_matrix.is_empty():
            raise ValueError(
                f"Expression matrix is empty after filtering. No matching '{transcript_id_column}' entries "
                f"between expression matrix and annotation found for gene '{target_gene}'."
            )

        # Get sets of transcripts in annotation and expression matrix
        annotation_transcripts = set(filtered_annotation[transcript_id_column].unique())
        expression_transcripts = set(filtered_expression_matrix[transcript_id_column].unique())

        # Identify transcripts present in annotation but missing in expression matrix
        missing_in_expression = annotation_transcripts - expression_transcripts

        # Transcripts in expression matrix but not in annotation are silently ignored

        # Warn about transcripts missing in the expression matrix
        if missing_in_expression:
            warnings.warn(
                f"{len(missing_in_expression)} transcript(s) are present in the annotation but missing in the expression matrix. "
                f"Missing transcripts: {', '.join(sorted(missing_in_expression))}. "
                "Only transcripts present in both will be returned."
            )

        # Ensure both filtered_annotation and filtered_expression_matrix contain only common transcripts
        common_transcripts = annotation_transcripts & expression_transcripts
        filtered_annotation = filtered_annotation.filter(
            pl.col(transcript_id_column).is_in(common_transcripts)
        )
        filtered_expression_matrix = filtered_expression_matrix.filter(
            pl.col(transcript_id_column).is_in(common_transcripts)
        )

        if order_by_expression:
            # Aggregate expression data to compute total expression per transcript
            aggregated_df = filtered_expression_matrix.group_by(transcript_id_column).agg(
                pl.col(expression_column).sum().alias("total_expression")
            )

            # Sort transcripts by total expression in descending order
            sorted_transcripts = aggregated_df.sort("total_expression", descending=True)

            # Determine transcripts to keep based on 'keep_top_expressed_transcripts'
            if isinstance(keep_top_expressed_transcripts, int) and keep_top_expressed_transcripts > 0:
                # Keep only the top N expressed transcripts
                if keep_top_expressed_transcripts < len(sorted_transcripts):
                    transcripts_to_keep = sorted_transcripts.head(keep_top_expressed_transcripts)[transcript_id_column]
                else:
                    # If requested number exceeds available transcripts, keep all and issue a warning
                    transcripts_to_keep = sorted_transcripts[transcript_id_column]
                    warnings.warn(
                        "The number specified in 'keep_top_expressed_transcripts' exceeds the total number of transcripts. "
                        "All transcripts will be kept."
                    )
            elif keep_top_expressed_transcripts == "all":
                # Keep all transcripts
                transcripts_to_keep = sorted_transcripts[transcript_id_column]
            else:
                # Raise error if 'keep_top_expressed_transcripts' is invalid
                raise ValueError(
                    f"'keep_top_expressed_transcripts' must be 'all' or a positive integer, got {keep_top_expressed_transcripts}."
                )

            # Filter annotation and expression matrix to include only the selected transcripts
            filtered_annotation = filtered_annotation.filter(
                pl.col(transcript_id_column).is_in(transcripts_to_keep)
            )
            filtered_expression_matrix = filtered_expression_matrix.filter(
                pl.col(transcript_id_column).is_in(transcripts_to_keep)
            )

            # Order annotation and expression matrix by total expression
            filtered_annotation = filtered_annotation.join(
                sorted_transcripts.select([transcript_id_column, "total_expression"]),
                on=transcript_id_column,
                how="inner"
            ).sort("total_expression", descending=False).drop("total_expression")

            filtered_expression_matrix = filtered_expression_matrix.join(
                sorted_transcripts.select([transcript_id_column, "total_expression"]),
                on=transcript_id_column,
                how="inner"
            ).sort("total_expression", descending=False).drop("total_expression")

        return filtered_annotation, filtered_expression_matrix

    else:
        # If no expression_matrix is provided, return only the filtered annotation
        return filtered_annotation
