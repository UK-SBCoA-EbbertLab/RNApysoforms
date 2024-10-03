import polars as pl
import warnings

def gene_filtering(
    gene_name_to_filter: str,
    annotation: pl.DataFrame,
    counts_matrix: pl.DataFrame = None,
    transcript_feature_column: str = "transcript_id"
):
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
