# tests/test_make_traces.py

import pytest
import polars as pl
import plotly.graph_objects as go
from RNApysoforms import make_traces
import warnings

def test_make_traces_basic():
    """
    Test the basic functionality with both annotation and expression_matrix provided.
    """
    # Create sample annotation DataFrame
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
        "start": [100, 200, 150, 250],
        "end": [150, 250, 200, 300],
        "type": ["exon", "CDS", "exon", "CDS"],
        "strand": ["+", "+", "-", "-"],
        "seqnames": ["chr1", "chr1", "chr2", "chr2"]
    })

    # Create sample expression matrix
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
        "sample_id": ["sample1", "sample2", "sample1", "sample2"],
        "counts": [100, 200, 150, 250]
    })

    # Call the function
    traces = make_traces(annotation=annotation_df, expression_matrix=expression_df)

    # Verify that traces are returned
    assert isinstance(traces, list)
    # The last element should be y_dict
    y_dict = traces[-1]
    assert isinstance(y_dict, dict)
    # Verify that traces contain expected number of items
    # Expected: [transcript_traces], [expression_traces], y_dict
    assert len(traces) == 3  # transcript_traces, expression_traces, y_dict

    # Verify that transcript_traces is a list containing dicts
    transcript_traces = traces[0]
    assert isinstance(transcript_traces, list)
    assert all(isinstance(trace, dict) for trace in transcript_traces)

    # Verify that expression_traces is a list containing go.Box objects
    expression_traces = traces[1]
    assert isinstance(expression_traces, list)
    assert all(isinstance(trace, go.Box) for trace in expression_traces)

def test_make_traces_annotation_only():
    """
    Test the function with only the annotation DataFrame provided.
    """
    # Create sample annotation DataFrame
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "start": [100, 200],
        "end": [150, 250],
        "type": ["exon", "CDS"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"]
    })

    # Call the function
    traces = make_traces(annotation=annotation_df)

    # Verify that traces are returned
    assert isinstance(traces, list)
    y_dict = traces[-1]
    assert isinstance(y_dict, dict)
    assert len(traces) == 2  # transcript_traces, y_dict

    transcript_traces = traces[0]
    assert isinstance(transcript_traces, list)
    assert all(isinstance(trace, dict) for trace in transcript_traces)

    # Since expression_matrix is None, there should be no expression traces
    assert len(traces) == 2  # transcript_traces and y_dict only

def test_make_traces_expression_only():
    """
    Test the function with only the expression_matrix DataFrame provided.
    """
    # Create sample expression matrix
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "sample_id": ["sample1", "sample2"],
        "counts": [100, 200]
    })

    # Call the function
    traces = make_traces(expression_matrix=expression_df)

    # Verify that traces are returned
    assert isinstance(traces, list)
    y_dict = traces[-1]
    assert isinstance(y_dict, dict)
    assert len(traces) == 2  # expression_traces, y_dict

    expression_traces = traces[0]
    assert isinstance(expression_traces, list)
    assert all(isinstance(trace, go.Box) for trace in expression_traces)

    # Since annotation is None, there should be no transcript traces
    assert len(traces) == 2  # expression_traces and y_dict only

def test_make_traces_no_input():
    """
    Test that ValueError is raised when neither annotation nor expression_matrix is provided.
    """
    with pytest.raises(ValueError) as excinfo:
        make_traces()
    assert "At least one of 'annotation' or 'expression_matrix' must be provided." in str(excinfo.value)

def test_make_traces_invalid_annotation_type():
    """
    Test that TypeError is raised when annotation is not a Polars DataFrame.
    """
    annotation_df = {"transcript_id": ["tx1"], "start": [100], "end": [150]}
    with pytest.raises(TypeError) as excinfo:
        make_traces(annotation=annotation_df)
    assert "Expected 'annotation' to be of type pl.DataFrame" in str(excinfo.value)

def test_make_traces_invalid_expression_matrix_type():
    """
    Test that TypeError is raised when expression_matrix is not a Polars DataFrame.
    """
    expression_df = {"transcript_id": ["tx1"], "sample_id": ["sample1"], "counts": [100]}
    with pytest.raises(TypeError) as excinfo:
        make_traces(expression_matrix=expression_df)
    assert "Expected 'expression_matrix' to be of type pl.DataFrame" in str(excinfo.value)

def test_make_traces_missing_required_columns_annotation():
    """
    Test that ValueError is raised when required columns are missing in annotation.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [150]
        # Missing 'strand', 'seqnames', 'type'
    })
    with pytest.raises(ValueError) as excinfo:
        make_traces(annotation=annotation_df)
    assert "The DataFrame is missing the following required columns:" in str(excinfo.value)

def test_make_traces_missing_required_columns_expression():
    """
    Test that ValueError is raised when required columns are missing in expression_matrix.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "counts": [100]
        # Missing 'sample_id'
    })
    with pytest.raises(ValueError) as excinfo:
        make_traces(expression_matrix=expression_df)
    assert "The DataFrame is missing the following required columns:" in str(excinfo.value)

def test_make_traces_no_common_transcripts():
    """
    Test that ValueError is raised when there are no matching transcripts between annotation and expression_matrix.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [150],
        "type": ["exon"],
        "strand": ["+"],
        "seqnames": ["chr1"]
    })
    expression_df = pl.DataFrame({
        "transcript_id": ["tx2"],
        "sample_id": ["sample1"],
        "counts": [100]
    })
    with pytest.raises(ValueError) as excinfo:
        make_traces(annotation=annotation_df, expression_matrix=expression_df)
    assert "No matching 'transcript_id' entries between annotation and expression matrix." in str(excinfo.value)

def test_make_traces_warning_missing_in_expression():
    """
    Test that a warning is issued when transcripts are missing in the expression_matrix.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx2"],
        "start": [100, 200],
        "end": [150, 250],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"]
    })
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "sample_id": ["sample1"],
        "counts": [100]
    })
    with pytest.warns(UserWarning) as record:
        traces = make_traces(annotation=annotation_df, expression_matrix=expression_df)
    assert len(record) == 1
    assert "transcript(s) are present in the annotation but missing in the expression matrix" in str(record[0].message)

def test_make_traces_warning_missing_in_annotation():
    """
    Test that a warning is issued when transcripts are missing in the annotation.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [150],
        "type": ["exon"],
        "strand": ["+"],
        "seqnames": ["chr1"]
    })
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx2"],
        "sample_id": ["sample1", "sample2"],
        "counts": [100, 200]
    })
    with pytest.warns(UserWarning) as record:
        traces = make_traces(annotation=annotation_df, expression_matrix=expression_df)
    assert len(record) == 1
    assert "transcript(s) are present in the expression matrix but missing in the annotation" in str(record[0].message)

def test_make_traces_order_transcripts_by_expression_matrix():
    """
    Test that transcripts are ordered by expression_matrix when order_transcripts_by_expression_matrix is True.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx2", "tx1"],
        "start": [200, 100],
        "end": [250, 150],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"]
    })
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx2"],
        "sample_id": ["sample1", "sample2"],
        "counts": [100, 200]
    })
    traces = make_traces(annotation=annotation_df, expression_matrix=expression_df)
    y_dict = traces[-1]
    # Since order_transcripts_by_expression_matrix is True by default, y_dict should order transcripts as in expression_df
    expected_order = {"tx1": 0, "tx2": 1}
    assert y_dict == expected_order

def test_make_traces_order_transcripts_by_annotation():
    """
    Test that transcripts are ordered by annotation when order_transcripts_by_expression_matrix is False.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx2", "tx1"],
        "start": [200, 100],
        "end": [250, 150],
        "type": ["exon", "exon"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"]
    })
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx2"],
        "sample_id": ["sample1", "sample2"],
        "counts": [100, 200]
    })
    traces = make_traces(annotation=annotation_df, expression_matrix=expression_df, order_transcripts_by_expression_matrix=False)
    y_dict = traces[-1]
    # Since order_transcripts_by_expression_matrix is False, y_dict should order transcripts as in annotation_df
    expected_order = {"tx2": 0, "tx1": 1}
    assert y_dict == expected_order

def test_make_traces_expression_plot_style_violin():
    """
    Test that violin plots are created when expression_plot_style is 'violin'.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
        "sample_id": ["sample1", "sample2", "sample1", "sample2"],
        "counts": [100, 200, 150, 250]
    })
    traces = make_traces(expression_matrix=expression_df, expression_plot_style="violin")
    expression_traces = traces[0]
    assert all(isinstance(trace, go.Violin) for trace in expression_traces)

def test_make_traces_invalid_expression_plot_style():
    """
    Test that ValueError is raised when an invalid expression_plot_style is provided.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "sample_id": ["sample1"],
        "counts": [100]
    })
    with pytest.raises(ValueError) as excinfo:
        make_traces(expression_matrix=expression_df, expression_plot_style="invalid_style")
    assert "Invalid expression_plot_style: invalid_style" in str(excinfo.value)

def test_make_traces_annotation_hue():
    """
    Test that annotation_hue correctly colors transcript features.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
        "start": [100, 200, 150, 250],
        "end": [150, 250, 200, 300],
        "type": ["exon", "CDS", "exon", "CDS"],
        "strand": ["+", "+", "-", "-"],
        "seqnames": ["chr1", "chr1", "chr2", "chr2"],
        "feature": ["feature1", "feature2", "feature1", "feature2"]
    })
    traces = make_traces(annotation=annotation_df, annotation_hue="feature")
    transcript_traces = traces[0]
    # Extract colors used in the traces
    colors = [trace['fillcolor'] for trace in transcript_traces]
    # Since there are two features, colors should be assigned accordingly
    assert len(set(colors)) == 2

def test_make_traces_expression_hue():
    """
    Test that expression_hue correctly colors expression plots.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx2", "tx2"],
        "sample_id": ["sample1", "sample2", "sample1", "sample2"],
        "counts": [100, 200, 150, 250],
        "group": ["A", "A", "B", "B"]
    })
    traces = make_traces(expression_matrix=expression_df, expression_hue="group")
    expression_traces = traces[0]
    # Since there are two groups, there should be traces for each group
    assert len(expression_traces) == 2  # Two groups

def test_make_traces_custom_column_names():
    """
    Test that custom column names are handled correctly.
    """
    annotation_df = pl.DataFrame({
        "trans_id": ["tx1"],
        "start_pos": [100],
        "end_pos": [150],
        "type": ["exon"],
        "strand": ["+"],
        "seqnames": ["chr1"]
    })
    expression_df = pl.DataFrame({
        "trans_id": ["tx1"],
        "sample": ["sample1"],
        "counts": [100]
    })
    traces = make_traces(
        annotation=annotation_df,
        expression_matrix=expression_df,
        y="trans_id",
        x_start="start_pos",
        x_end="end_pos",
        sample_id_column="sample",
        hover_start="start_pos", 
        hover_end="end_pos"
    )
    y_dict = traces[-1]
    assert y_dict == {"tx1": 0}

def test_make_traces_expression_columns_list():
    """
    Test that multiple expression_columns are handled correctly.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "sample_id": ["sample1", "sample2"],
        "counts": [100, 200],
        "tpm": [10, 20]
    })
    traces = make_traces(expression_matrix=expression_df, expression_columns=["counts", "tpm"])
    # There should be traces for each expression column
    expression_traces = traces[0:-1]
    assert len(expression_traces) == 2  # Two expression columns

def test_make_traces_expression_columns_string():
    """
    Test that a single string for expression_columns is handled correctly.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "sample_id": ["sample1", "sample2"],
        "counts": [100, 200],
    })
    traces = make_traces(expression_matrix=expression_df, expression_columns="counts")
    expression_traces = traces[0]
    assert len(expression_traces) == 1  # One expression column

def test_make_traces_annotation_color_map():
    """
    Test that annotation_color_map is used correctly.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [150],
        "type": ["exon"],
        "strand": ["+"],
        "seqnames": ["chr1"],
        "feature": ["feature1"]
    })
    color_map = {"feature1": "red"}
    traces = make_traces(annotation=annotation_df, annotation_hue="feature", annotation_color_map=color_map)
    transcript_traces = traces[0]
    # The fillcolor should be 'red' as specified in the color_map
    assert transcript_traces[0]['fillcolor'] == "red"

def test_make_traces_expression_color_map():
    """
    Test that expression_color_map is used correctly.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "sample_id": ["sample1", "sample2"],
        "counts": [100, 200],
        "group": ["A", "A"]
    })
    color_map = {"A": "blue"}
    traces = make_traces(expression_matrix=expression_df, expression_hue="group", expression_color_map=color_map)
    expression_traces = traces[0]
    # The fillcolor should be 'blue' as specified in the color_map
    for trace in expression_traces:
        assert trace.fillcolor == "blue"

def test_make_traces_missing_hue_column_annotation():
    """
    Test that ValueError is raised when annotation_hue column is missing in annotation.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [150],
        "type": ["exon"],
        "strand": ["+"],
        "seqnames": ["chr1"]
        # Missing 'feature' column
    })
    with pytest.raises(ValueError) as excinfo:
        make_traces(annotation=annotation_df, annotation_hue="feature")
    assert "The DataFrame is missing the following required columns:" in str(excinfo.value)

def test_make_traces_missing_hue_column_expression():
    """
    Test that ValueError is raised when expression_hue column is missing in expression_matrix.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "sample_id": ["sample1"],
        "counts": [100]
        # Missing 'group' column
    })
    with pytest.raises(ValueError) as excinfo:
        make_traces(expression_matrix=expression_df, expression_hue="group")
    assert "The DataFrame is missing the following required columns:" in str(excinfo.value)

def test_make_traces_expression_plot_opacity():
    """
    Test that expression_plot_opacity parameter is applied correctly.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "sample_id": ["sample1"],
        "counts": [100]
    })
    traces = make_traces(expression_matrix=expression_df, expression_plot_opacity=0.5)
    expression_traces = traces[0]
    for trace in expression_traces:
        assert trace.opacity == 0.5

def test_make_traces_transcript_plot_opacity():
    """
    Test that transcript_plot_opacity parameter is applied correctly.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [150],
        "type": ["exon"],
        "strand": ["+"],
        "seqnames": ["chr1"]
    })
    traces = make_traces(annotation=annotation_df, transcript_plot_opacity=0.5)
    transcript_traces = traces[0]
    for trace in transcript_traces:
        assert trace['opacity'] == 0.5

def test_make_traces_marker_size():
    """
    Test that marker_size parameter is applied correctly.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "sample_id": ["sample1"],
        "counts": [100]
    })
    traces = make_traces(expression_matrix=expression_df, marker_size=10)
    expression_traces = traces[0]
    for trace in expression_traces:
        assert trace.marker.size == 10

def test_make_traces_custom_hover_info():
    """
    Test that custom hover_start and hover_end are used correctly.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "custom_start": [100],
        "custom_end": [150],
        "start": [100],
        "end": [150],
        "type": ["exon"],
        "strand": ["+"],
        "seqnames": ["chr1"]
    })
    traces = make_traces(annotation=annotation_df, hover_start="custom_start", hover_end="custom_end")
    transcript_traces = traces[0]
    # Verify that hovertemplate contains custom_start and custom_end
    hovertemplate = transcript_traces[0]['hovertemplate']
    assert "Start:</b> 100" in hovertemplate
    assert "End:</b> 150" in hovertemplate

def test_make_traces_arrow_size():
    """
    Test that arrow_size parameter is applied correctly.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [500],  # Long intron to trigger arrow
        "type": ["intron"],
        "strand": ["+"],
        "seqnames": ["chr1"]
    })
    traces = make_traces(annotation=annotation_df, arrow_size=20)
    transcript_traces = traces[0]
    # Find the arrow trace
    arrow_trace = next((trace for trace in transcript_traces if 'marker' in trace and trace['marker']['symbol'] == 'arrow-right'), None)
    assert arrow_trace is not None
    assert arrow_trace['marker']['size'] == 20

def test_make_traces_no_exons():
    """
    Test when there are no exons in the annotation DataFrame.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [200],
        "end": [250],
        "type": ["CDS"],
        "strand": ["+"],
        "seqnames": ["chr1"]
    })
    traces = make_traces(annotation=annotation_df)
    transcript_traces = traces[0]
    # Since there are no exons, only CDS traces should be present
    assert len(transcript_traces) == 1  # Only one CDS trace

def test_make_traces_exon_and_cds_heights():
    """
    Test that exon_height and cds_height parameters are applied correctly.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "start": [100, 200],
        "end": [150, 250],
        "type": ["exon", "CDS"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"]
    })
    traces = make_traces(annotation=annotation_df, exon_height=0.1, cds_height=0.2)
    transcript_traces = traces[0]
    # Verify y coordinates of exon and CDS traces
    exon_trace = next((trace for trace in transcript_traces if trace['hovertemplate'].find('Feature Type:</b> exon') != -1), None)
    cds_trace = next((trace for trace in transcript_traces if trace['hovertemplate'].find('Feature Type:</b> CDS') != -1), None)
    assert exon_trace is not None
    assert cds_trace is not None
    # Exon height should be 0.1
    y_coords_exon = exon_trace['y']
    height_exon = y_coords_exon[2] - y_coords_exon[0]
    assert height_exon == 0.1
    # CDS height should be 0.2
    y_coords_cds = cds_trace['y']
    height_cds = y_coords_cds[2] - y_coords_cds[0]
    assert height_cds == 0.2

def test_make_traces_line_color():
    """
    Test that line_color parameter is applied correctly.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [150],
        "type": ["exon"],
        "strand": ["+"],
        "seqnames": ["chr1"]
    })
    traces = make_traces(annotation=annotation_df, line_color="green")
    transcript_traces = traces[0]
    assert transcript_traces[0]['line']['color'] == "green"

def test_make_traces_intron_line_width():
    """
    Test that intron_line_width parameter is applied correctly.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [150],
        "end": [300],
        "type": ["intron"],
        "strand": ["+"],
        "seqnames": ["chr1"]
    })
    traces = make_traces(annotation=annotation_df, intron_line_width=2)
    transcript_traces = traces[0]
    intron_trace = next((trace for trace in transcript_traces if trace['mode'] == 'lines'), None)
    assert intron_trace is not None
    assert intron_trace['line']['width'] == 2

def test_make_traces_exon_line_width():
    """
    Test that exon_line_width parameter is applied correctly.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [150],
        "type": ["exon"],
        "strand": ["+"],
        "seqnames": ["chr1"]
    })
    traces = make_traces(annotation=annotation_df, exon_line_width=2)
    transcript_traces = traces[0]
    assert transcript_traces[0]['line']['width'] == 2

def test_make_traces_marker_opacity():
    """
    Test that marker_opacity parameter is applied correctly.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "sample_id": ["sample1"],
        "counts": [100]
    })
    traces = make_traces(expression_matrix=expression_df, marker_opacity=0.5)
    expression_traces = traces[0]
    for trace in expression_traces:
        assert trace.marker.opacity == 0.5

def test_make_traces_box_points():
    """
    Test that box_points parameter is applied correctly.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "sample_id": ["sample1"],
        "counts": [100]
    })
    traces = make_traces(expression_matrix=expression_df, box_points=False)
    expression_traces = traces[0]
    for trace in expression_traces:
        assert trace.boxpoints == False

def test_make_traces_show_box_mean():
    """
    Test that show_box_mean parameter is applied correctly.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "sample_id": ["sample1", "sample2"],
        "counts": [100, 200]
    })
    traces = make_traces(expression_matrix=expression_df, show_box_mean=True)
    expression_traces = traces[0]
    for trace in expression_traces:
        assert trace.boxmean == True

def test_make_traces_expression_plot_legend_title():
    """
    Test that expression_plot_legend_title parameter is applied correctly.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "sample_id": ["sample1", "sample2"],
        "counts": [100, 200],
        "group": ["A", "B"]
    })
    traces = make_traces(expression_matrix=expression_df, expression_hue="group", expression_plot_legend_title="Custom Legend Title")
    expression_traces = traces[0]
    # Check if any trace has the custom legend title
    has_custom_title = any(
        getattr(trace, 'legendgrouptitle', None) is not None and 
        getattr(trace.legendgrouptitle, 'text', None) == "Custom Legend Title" 
        for trace in expression_traces
    )
    assert has_custom_title

def test_make_traces_transcript_plot_legend_title():
    """
    Test that transcript_plot_legend_title parameter is applied correctly.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "start": [100, 200],
        "end": [150, 250],
        "type": ["exon", "CDS"],
        "strand": ["+", "+"],
        "seqnames": ["chr1", "chr1"],
        "feature": ["A", "B"]
    })
    traces = make_traces(annotation=annotation_df, annotation_hue="feature", transcript_plot_legend_title="Custom Legend Title")
    transcript_traces = traces[0]
    for trace in transcript_traces:
        if trace.get('legendgrouptitle_text'):
            assert trace['legendgrouptitle_text'] == "Custom Legend Title"
            break

def test_make_traces_spanmode():
    """
    Test that spanmode parameter is applied correctly in violin plots.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "sample_id": ["sample1", "sample2"],
        "counts": [100, 200]
    })
    traces = make_traces(expression_matrix=expression_df, expression_plot_style="violin", spanmode="soft")
    expression_traces = traces[0]
    for trace in expression_traces:
        assert trace.spanmode == "soft"

def test_make_traces_marker_color():
    """
    Test that marker_color parameter is applied correctly.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "sample_id": ["sample1", "sample2"],
        "counts": [100, 200]
    })
    traces = make_traces(expression_matrix=expression_df, marker_color="red")
    expression_traces = traces[0]
    for trace in expression_traces:
        assert trace.marker.color == "red"

def test_make_traces_marker_jitter():
    """
    Test that marker_jitter parameter is applied correctly.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1", "tx1", "tx1"],
        "sample_id": ["sample1", "sample2", "sample3", "sample4"],
        "counts": [100, 150, 200, 250]
    })
    traces = make_traces(expression_matrix=expression_df, marker_jitter=0.5)
    expression_traces = traces[0]
    for trace in expression_traces:
        assert trace.jitter == 0.5

def test_make_traces_expression_fill_color():
    """
    Test that expression_fill_color parameter is applied correctly when no expression_hue is provided.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx1"],
        "sample_id": ["sample1", "sample2"],
        "counts": [100, 200]
    })
    traces = make_traces(expression_matrix=expression_df, expression_fill_color="purple")
    expression_traces = traces[0]
    for trace in expression_traces:
        assert trace.fillcolor == "purple"

def test_make_traces_annotation_fill_color():
    """
    Test that annotation_fill_color parameter is applied correctly when no annotation_hue is provided.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [150],
        "type": ["exon"],
        "strand": ["+"],
        "seqnames": ["chr1"]
    })
    traces = make_traces(annotation=annotation_df, annotation_fill_color="orange")
    transcript_traces = traces[0]
    assert transcript_traces[0]['fillcolor'] == "orange"

def test_make_traces_intron_arrows_negative_strand():
    """
    Test that intron arrows point in the correct direction for negative strand.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [500],
        "type": ["intron"],
        "strand": ["-"],
        "seqnames": ["chr1"]
    })
    traces = make_traces(annotation=annotation_df)
    transcript_traces = traces[0]
    arrow_trace = next((trace for trace in transcript_traces if 'marker' in trace and trace['marker']['symbol'] == 'arrow-left'), None)
    assert arrow_trace is not None

def test_make_traces_intron_arrows_positive_strand():
    """
    Test that intron arrows point in the correct direction for positive strand.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [500],
        "type": ["intron"],
        "strand": ["+"],
        "seqnames": ["chr1"]
    })
    traces = make_traces(annotation=annotation_df)
    transcript_traces = traces[0]
    arrow_trace = next((trace for trace in transcript_traces if 'marker' in trace and trace['marker']['symbol'] == 'arrow-right'), None)
    assert arrow_trace is not None

def test_make_traces_expression_plot_style_default():
    """
    Test that the default expression_plot_style is 'boxplot'.
    """
    expression_df = pl.DataFrame({
        "transcript_id": ["tx1", "tx2"],
        "sample_id": ["sample1", "sample2"],
        "counts": [100, 200]
    })
    traces = make_traces(expression_matrix=expression_df)
    expression_traces = traces[0]
    assert all(isinstance(trace, go.Box) for trace in expression_traces)

def test_make_traces_annotation_missing_strand():
    """
    Test that ValueError is raised when 'strand' column is missing in annotation.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [150],
        "type": ["intron"],
        "seqnames": ["chr1"]
        # Missing 'strand'
    })
    with pytest.raises(ValueError) as excinfo:
        make_traces(annotation=annotation_df)
    assert "The DataFrame is missing the following required columns:" in str(excinfo.value)

def test_make_traces_annotation_missing_seqnames():
    """
    Test that ValueError is raised when 'seqnames' column is missing in annotation.
    """
    annotation_df = pl.DataFrame({
        "transcript_id": ["tx1"],
        "start": [100],
        "end": [150],
        "type": ["intron"],
        "strand": ["+"]
        # Missing 'seqnames'
    })
    with pytest.raises(ValueError) as excinfo:
        make_traces(annotation=annotation_df)
    assert "The DataFrame is missing the following required columns:" in str(excinfo.value)
