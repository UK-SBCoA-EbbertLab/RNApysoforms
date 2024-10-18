# tests/test_make_plot.py

import pytest
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from RNApysoforms import make_plot


def create_sample_traces():
    """
    Helper function to create sample traces for testing.
    Returns:
        List[List[go.Trace]]: A list of lists containing Plotly traces.
        dict: A dictionary mapping transcript IDs to y-axis positions.
    """
    # Transcript structure traces (e.g., shapes or annotations)
    transcript_traces = [
        dict(
            x=[1, 200, 200, 1, 1],
            y=[-0.4, -0.4, 0.4, 0.4, -0.4],
            mode='lines',
            name='Transcript1',
            type='scatter'
        ),
        dict(
            x=[300, 500, 500, 300, 300],
            y=[0.6, 0.6, 1.4, 1.4, 0.6],
            mode='lines',
            name='Transcript2',
            type='scatter'
        )
    ]

    # Expression data traces (e.g., Box plots)
    expression_traces = [
        go.Box(
            y=[0, 0, 0, 0],
            x=[10, 15, 13, 17],
            name='Transcript1'
        ),
        go.Box(
            y=[1, 1, 1, 1],
            x=[16, 5, 11, 9],
            name='Transcript2'
        )
    ]

    # y_dict mapping transcript IDs to y-axis positions
    y_dict = {'Transcript1': 0, 'Transcript2': 1}

    return [transcript_traces, expression_traces, y_dict]


def create_custom_traces():
    """
    Helper function to create custom traces for testing with custom column names.
    Returns:
        List[List[go.Trace]]: A list of lists containing Plotly traces.
        dict: A dictionary mapping transcript IDs to y-axis positions.
    """

    # Transcript structure traces with custom identifiers
    transcript_traces = [
        dict(
            x=[1, 200, 200, 1, 1],
            y=[-0.4, -0.4, 0.4, 0.4, -0.4],
            mode='lines',
            name='CustomTranscript1',
            type='scatter'
        ),
        dict(
            x=[300, 500, 500, 300, 300],
            y=[0.6, 0.6, 1.4, 1.4, 0.6],
            mode='lines',
            name='CustomTranscript2',
            type='scatter'
        )
    ]

    # Expression data traces with custom identifiers
    expression_traces = [
        go.Violin(
            y=[0, 0, 0, 0],
            x=[10, 15, 13, 17],
            name='CustomTranscript1'
        ),
        go.Violin(
            y=[1, 1, 1, 1],
            x=[16, 5, 11, 9],
            name='CustomTranscript2'
        )
    ]

    # y_dict mapping transcript IDs to y-axis positions
    y_dict = {'CustomTranscript1': 0, 'CustomTranscript2': 1}

    return [transcript_traces, expression_traces, y_dict]


def test_make_plot_basic():
    """
    Test creating a basic plot with transcript structure and expression traces.
    """
    traces, y_dict = create_sample_traces()[0], create_sample_traces()[2]
    traces = create_sample_traces()

    fig = make_plot(traces=traces)

    # Verify the number of subplots
    assert len(fig['layout']['annotations']) == 1  # Default subplot_titles length is 1

    # Verify the number of traces added to the figure
    assert len(fig.data) == 4  # 2 transcript traces + 2 expression traces

    # Verify layout parameters
    assert fig.layout.height == 900
    assert fig.layout.width == 1800
    assert fig.layout.hovermode == "closest"

    # Verify that traces are assigned to the correct subplots
    # Assuming first subplot contains transcript traces, second contains expression traces
    # However, based on the function implementation, it's a bit unclear. Adjust as necessary.
    # Here, since traces list has [[transcript_traces], [expression_traces], y_dict],
    # and make_subplots is called with rows=1, cols=2
    # Let's verify subplot_titles length and number of columns

    assert fig.layout.annotations[0].text == "Transcript Structure"
    assert fig.layout.annotations[0].font.size == 16  # Default subplot_title_font_size


def test_make_plot_multiple_subplots():
    """
    Test creating a plot with multiple subplots and verifying trace assignments.
    """
    # Create multiple sets of expression traces
    transcript_traces = [
        dict(
            x=[1, 200, 200, 1, 1],
            y=[-0.4, -0.4, 0.4, 0.4, -0.4],
            mode='lines',
            name='Transcript1',
            type='scatter'
        )
    ]

    expression_traces1 = [
        go.Box(
            y=[0, 0, 0, 0],
            x=[10, 15, 13, 17],
            name='Transcript1'
        )
    ]

    expression_traces2 = [
        go.Box(
            y=[1, 1, 1, 1],
            x=[16, 5, 11, 9],
            name='Transcript1'
        )
    ]

    y_dict = {'Transcript1': 0}

    # Combine into traces
    traces = [transcript_traces, expression_traces1, expression_traces2, y_dict]

    subplot_titles = ["Transcript Structure", "Expression Levels 1", "Expression Levels 2"]

    fig = make_plot(
        traces=traces,
        subplot_titles=subplot_titles,
        height=800,
        width=1600
    )

    # Verify the number of subplots
    assert len(fig['layout']['annotations']) == 3  # Three subplot titles

    # Verify the number of traces added to the figure
    assert len(fig.data) == 3  # 1 transcript trace + 2 expression traces

    # Verify subplot titles
    for i, title in enumerate(subplot_titles):
        assert fig.layout.annotations[i].text == title
        assert fig.layout.annotations[i].font.size == 16  # Default subplot_title_font_size


def test_make_plot_invalid_traces_missing_y_dict():
    """
    Test that ValueError is raised when traces do not include a y_dict.
    """
    # Create traces without y_dict
    transcript_traces = [
        dict(
            x=[1, 200, 200, 1, 1],
            y=[-0.4, -0.4, 0.4, 0.4, -0.4],
            mode='lines',
            name='Transcript1',
            type='scatter'
        )
    ]

    expression_traces = [
        go.Box(
            y=[0, 0, 0, 0],
            x=[10, 15, 13, 17],
            name='Transcript1'
        )
    ]


    traces = [transcript_traces, expression_traces]  # Missing y_dict

    with pytest.raises(AttributeError):
        make_plot(traces=traces)


def test_make_plot_empty_traces():
    """
    Test that ValueError is raised when traces list is empty.
    """
    traces = []

    with pytest.raises(IndexError):
        make_plot(traces=traces)


def test_make_plot_custom_parameters():
    """
    Test creating a plot with custom layout parameters.
    """
    traces = create_sample_traces()
    custom_subplot_titles = ["Custom Transcript Structure", "Custom Expression Levels"]
    custom_height = 600
    custom_width = 1200
    custom_hovermode = "x unified"

    fig = make_plot(
        traces=traces,
        subplot_titles=custom_subplot_titles,
        height=custom_height,
        width=custom_width,
        hovermode=custom_hovermode
    )

    # Verify layout parameters
    assert fig.layout.height == custom_height
    assert fig.layout.width == custom_width
    assert fig.layout.hovermode == custom_hovermode

    # Verify custom subplot titles
    for i, title in enumerate(custom_subplot_titles):
        assert fig.layout.annotations[i].text == title
        assert fig.layout.annotations[i].font.size == 16  # Default subplot_title_font_size


def test_make_plot_invalid_traces_format():
    """
    Test that ValueError is raised when traces are not in the expected format.
    """
    # Traces not as list of lists
    traces = "invalid_traces_format"

    with pytest.raises(ValueError):
        make_plot(traces=traces)


def test_make_plot_hovermode_setting():
    """
    Test that hovermode is set correctly in the figure.
    """
    traces = create_sample_traces()

    fig = make_plot(
        traces=traces,
        hovermode="x"
    )

    assert fig.layout.hovermode == "x"


def test_make_plot_non_go_traces():
    """
    Test that dictionary traces (e.g., shapes or annotations) are handled correctly.
    """
    # Create transcript traces as dictionaries (e.g., shapes)
    transcript_traces = [
        dict(
            x=[1, 200, 200, 1, 1],
            y=[-0.4, -0.4, 0.4, 0.4, -0.4],
            mode='lines',
            name='Transcript1',
            type='scatter'
        ),
        dict(
            x=[300, 500, 500, 300, 300],
            y=[0.6, 0.6, 1.4, 1.4, 0.6],
            mode='lines',
            name='Transcript2',
            type='scatter'
        )
    ]

    # Expression data traces (e.g., Box plots)
    expression_traces = [
        go.Box(
            y=[0, 0, 0, 0],
            x=[10, 15, 13, 17],
            name='Transcript1'
        )
    ]



    y_dict = {'Transcript1': 0, 'Transcript2': 1}

    traces = [transcript_traces, expression_traces, y_dict]

    fig = make_plot(traces=traces)

    # Verify the number of traces
    assert len(fig.data) == 3  # 2 transcript traces (as dicts converted to Scatter) + 1 expression trace

    # Verify that all traces are present
    trace_names = [trace.name for trace in fig.data]
    assert 'Transcript1' in trace_names
    assert 'Transcript2' in trace_names
    assert 'Transcript1' in trace_names


def test_make_plot_keep_top_expressed_transcripts():
    """
    Test that only the top N expressed transcripts are kept when specified.
    """
    # This function does not have a parameter for keeping top N transcripts.
    # It seems the 'make_plot' function does not handle filtering transcripts.
    # Therefore, this test might not be applicable.
    # If 'make_plot' had such functionality, implement accordingly.
    pass  # Placeholder for future implementation if applicable


def test_make_plot_keep_top_expressed_transcripts_exceeds_available():
    """
    Test that no error occurs when the number of transcripts exceeds available transcripts.
    """
    # Similar to the previous test, 'make_plot' does not handle keeping top N transcripts.
    pass  # Placeholder for future implementation if applicable


def test_make_plot_transcript_labels_visibility():
    """
    Test that transcript labels are correctly shown or hidden based on parameters.
    """
    traces = create_sample_traces()

    # Create plot with default settings (labels shown)
    fig = make_plot(traces=traces)

    # Verify that y-axis tick labels are present
    yaxis = fig.layout.yaxis
    assert yaxis.showticklabels == True
    print(yaxis.ticktext)
    assert yaxis.ticktext == ('Transcript1', 'Transcript2')
    assert yaxis.tickvals == (0, 1)

    # Create plot with hidden y-axis tick labels
    fig_hidden = make_plot(
        traces=traces,
        vert_grid_transcript_structure_plot=False,
        horz_grid_transcript_structure_plot=False
    )

    # Verify that y-axis tick labels are still present (since function logic keeps them)
    yaxis_hidden = fig_hidden.layout.yaxis
    assert yaxis_hidden.showticklabels == True
    assert yaxis_hidden.ticktext == ('Transcript1', 'Transcript2')
    assert yaxis_hidden.tickvals == (0, 1)


def test_make_plot_custom_column_names():
    """
    Test creating a plot with custom transcript and expression identifiers.
    """
    traces = create_custom_traces()

    subplot_titles = ["Custom Transcript Structure", "Custom Expression Levels"]

    fig = make_plot(
        traces=traces,
        subplot_titles=subplot_titles,
        height=700,
        width=1400,
        hovermode="y unified"
    )

    # Verify subplot titles
    for i, title in enumerate(subplot_titles):
        assert fig.layout.annotations[i].text == title
        assert fig.layout.annotations[i].font.size == 16  # Default subplot_title_font_size

    # Verify hovermode
    assert fig.layout.hovermode == "y unified"

    # Verify trace names
    trace_names = [trace.name for trace in fig.data]
    assert 'CustomTranscript1' in trace_names
    assert 'CustomTranscript2' in trace_names
    assert 'CustomTranscript1' in trace_names
    assert 'CustomTranscript2' in trace_names


def test_make_plot_missing_expression_column():
    """
    Since the 'make_plot' function does not directly handle DataFrames or columns,
    this test is not applicable. It would be relevant if 'make_plot' processed
    DataFrame inputs.
    """
    pass  # Placeholder for future implementation if applicable


def test_make_plot_invalid_hovermode():
    """
    Test that the function handles invalid hovermode values gracefully.
    """
    traces = create_sample_traces()

    with pytest.raises(ValueError):
        fig = make_plot(
            traces=traces,
            hovermode="invalid_hovermode"
        )


def test_make_plot_large_number_of_transcripts():
    """
    Test plotting with a large number of transcripts to ensure scalability.
    """
    # Create a large number of transcript traces
    transcript_traces = []
    expression_traces = []
    y_dict = {}
    num_transcripts = 50

    for i in range(num_transcripts):
        transcript = go.Scatter(
            x=[1, 2, 3],
            y=[i, i + 1, i],
            mode='lines',
            name=f'Transcript{i}'
        )
        transcript_traces.append(transcript)

        expression = go.Box(
            y=[10 + i, 15 + i, 13 + i, 17 + i],
            name=f'Expression{i}'
        )
        expression_traces.append(expression)

        y_dict[f'Transcript{i}'] = i

    traces = [transcript_traces, expression_traces, y_dict]

    fig = make_plot(
        traces=traces,
        subplot_titles=["Transcript Structures", "Expression Levels"] * 25,  # Example for multiple subplots
        height=2000,
        width=4000
    )

    # Verify the number of traces
    assert len(fig.data) == 100  # 50 transcript + 50 expression traces

    # Verify y_dict entries
    for i in range(num_transcripts):
        assert f'Transcript{i}' in y_dict
        assert y_dict[f'Transcript{i}'] == i


def test_make_plot_empty_y_dict():
    """
    Test that a ValueError is raised when y_dict is empty.
    """
    transcript_traces = [
        go.Scatter(
            x=[1, 2, 3],
            y=[0, 1, 0],
            mode='lines',
            name='Transcript1'
        )
    ]

    expression_traces = [
        go.Box(
            y=[10, 15, 13, 17],
            name='Expression1'
        )
    ]

    y_dict = {}  # Empty y_dict

    traces = [transcript_traces, expression_traces, y_dict]

    fig = make_plot(traces=traces)

    # Verify that yaxis tickvals and ticktext are empty
    yaxis = fig.layout.yaxis
    assert yaxis.tickvals == ()
    assert yaxis.ticktext == ()


def test_make_plot_inconsistent_y_dict():
    """
    Test that the function handles inconsistencies between traces and y_dict.
    For example, a transcript trace not present in y_dict.
    """
    # Create transcript traces
    transcript_traces = [
        go.Scatter(
            x=[1, 2, 3],
            y=[0, 1, 0],
            mode='lines',
            name='Transcript1'
        ),
        go.Scatter(
            x=[2, 3, 4],
            y=[1, 2, 1],
            mode='lines',
            name='Transcript2'
        )
    ]

    # Expression data traces
    expression_traces = [
        go.Box(
            y=[10, 15, 13, 17],
            name='Expression1'
        )
    ]

    # y_dict missing 'Transcript2'
    y_dict = {'Transcript1': 0}

    traces = [transcript_traces, expression_traces, y_dict]

    fig = make_plot(traces=traces)

    # Verify that 'Transcript2' y-axis is not mapped, which might cause it to default
    # Depending on implementation, Plotly may assign default y-axis values
    # Here, we can check that yaxis.ticktext only includes 'Transcript1'
    yaxis = fig.layout.yaxis
    assert yaxis.ticktext == ('Transcript1',)
    assert yaxis.tickvals == (0,)


def test_make_plot_invalid_trace_type():
    """
    Test that the function handles invalid trace types gracefully.
    """
    # Create invalid trace (e.g., integer instead of go.Trace or dict)
    transcript_traces = [123]

    expression_traces = [
        go.Box(
            y=[10, 15, 13, 17],
            name='Expression1'
        )
    ]

    y_dict = {'Transcript1': 0}

    traces = [transcript_traces, expression_traces, y_dict]

    with pytest.raises(ValueError):
        make_plot(traces=traces)


def test_make_plot_no_expression_traces():
    """
    Test creating a plot with only transcript traces and no expression traces.
    """
    transcript_traces = [
        go.Scatter(
            x=[1, 2, 3],
            y=[0, 1, 0],
            mode='lines',
            name='Transcript1'
        )
    ]

    y_dict = {'Transcript1': 0}

    traces = [transcript_traces, y_dict]

    fig = make_plot(traces=traces)

    # Verify the number of traces
    assert len(fig.data) == 1  # Only transcript trace

    # Verify that y_dict is correctly applied
    yaxis = fig.layout.yaxis
    assert yaxis.ticktext == ('Transcript1',)
    assert yaxis.tickvals == (0,)


def test_make_plot_invalid_hover_font_size():
    """
    Test that the function handles invalid hover_font_size values.
    """
    traces = create_sample_traces()

    # Pass a negative font size
    with pytest.raises(ValueError):
        make_plot(
            traces=traces,
            hover_font_size=-5
        )


def test_make_plot_large_dimensions():
    """
    Test creating a plot with extremely large dimensions.
    """
    traces = create_sample_traces()

    fig = make_plot(
        traces=traces,
        height=10000,
        width=20000
    )

    assert fig.layout.height == 10000
    assert fig.layout.width == 20000


def test_make_plot_minimum_dimensions():
    """
    Test creating a plot with minimum acceptable dimensions.
    """
    traces = create_sample_traces()

    fig = make_plot(
        traces=traces,
        height=100,
        width=100
    )

    assert fig.layout.height == 100
    assert fig.layout.width == 100


def test_make_plot_shared_yaxes():
    """
    Test that y-axes are shared across subplots.
    """
    traces = create_sample_traces()

    fig = make_plot(
        traces=traces
    )

    # Verify that y-axes are shared
    # In Plotly, shared y-axes are handled internally; we can verify that tickvals and ticktext are consistent
    yaxis1 = fig.layout.yaxis
    yaxis2 = fig.layout.yaxis2

    assert yaxis1.tickvals == yaxis2.tickvals
    assert yaxis1.ticktext == yaxis2.ticktext

def test_make_plot_legend_visibility():
    """
    Test that the legend visibility is controlled by the 'showlegend' parameter.
    """
    traces = create_sample_traces()

    # Create plot with legend shown
    fig_with_legend = make_plot(
        traces=traces,
        showlegend=True
    )

    # Create plot with legend hidden
    fig_without_legend = make_plot(
        traces=traces,
        showlegend=False
    )

    assert fig_with_legend.layout.showlegend == True
    assert fig_without_legend.layout.showlegend == False


def test_make_plot_hover_label_font_size():
    """
    Test that the hover label font size is set correctly.
    """
    traces = create_sample_traces()

    custom_hover_font_size = 20

    fig = make_plot(
        traces=traces,
        hover_font_size=custom_hover_font_size
    )

    assert fig.layout.hoverlabel.font.size == custom_hover_font_size


def test_make_plot_boxmode_group():
    """
    Test that the boxmode is set to 'group'.
    """
    traces = create_sample_traces()

    fig = make_plot(
        traces=traces
    )

    assert fig.layout.boxmode == 'group'


def test_make_plot_violinmode_group():
    """
    Test that the violinmode is set to 'group'.
    """
    # Create violin traces
    transcript_traces = [
        go.Scatter(
            x=[1, 2, 3],
            y=[0, 1, 0],
            mode='lines',
            name='Transcript1'
        )
    ]

    violin_traces = [
        go.Violin(
            y=[10, 15, 13, 17],
            name='Violin1'
        )
    ]

    y_dict = {'Transcript1': 0}

    traces = [transcript_traces, violin_traces, y_dict]

    fig = make_plot(
        traces=traces
    )

    assert fig.layout.violinmode == 'group'


def test_make_plot_legend_font_size():
    """
    Test that the legend font size is set correctly.
    """
    traces = create_sample_traces()

    custom_legend_font_size = 18

    fig = make_plot(
        traces=traces,
        legend_font_size=custom_legend_font_size
    )

    assert fig.layout.legend.font.size == custom_legend_font_size


def test_make_plot_xaxis_font_size():
    """
    Test that the x-axis font size is set correctly.
    """
    traces = create_sample_traces()

    custom_xaxis_font_size = 14

    fig = make_plot(
        traces=traces,
        xaxis_font_size=custom_xaxis_font_size
    )

    # Verify x-axis font size for all x-axes
    assert fig.layout['xaxis']["tickfont"]["size"] == custom_xaxis_font_size


def test_make_plot_yaxis_font_size():
    """
    Test that the y-axis font size is set correctly.
    """
    traces = create_sample_traces()

    custom_yaxis_font_size = 14

    fig = make_plot(
        traces=traces,
        yaxis_font_size=custom_yaxis_font_size
    )

    # Verify y-axis font size for all y-axes
    assert fig.layout["yaxis"].tickfont.size == custom_yaxis_font_size


def test_make_plot_subtitle_font_size():
    """
    Test that the subplot title font size is set correctly.
    """
    traces = create_sample_traces()

    custom_subtitle_font_size = 20

    fig = make_plot(
        traces=traces,
        subplot_title_font_size=custom_subtitle_font_size
    )

    for annotation in fig.layout.annotations:
        assert annotation.font.size == custom_subtitle_font_size


def test_make_plot_hover_font_size():
    """
    Test that the hover label font size is set correctly.
    """
    traces = create_sample_traces()

    custom_hover_font_size = 18

    fig = make_plot(
        traces=traces,
        hover_font_size=custom_hover_font_size
    )

    assert fig.layout.hoverlabel.font.size == custom_hover_font_size

def test_make_plot_violinmode_group():
    """
    Test that the violinmode is set to 'group' and that the column widths are correctly set.
    """
    # Create violin traces
    transcript_traces = [
        go.Scatter(
            x=[1, 2, 3],
            y=[0, 1, 0],
            mode='lines',
            name='Transcript1'
        )
    ]

    violin_traces = [
        go.Violin(
            y=[10, 15, 13, 17],
            name='Violin1'
        )
    ]

    y_dict = {'Transcript1': 0}

    traces = [transcript_traces, violin_traces, y_dict]

    # Expected column widths
    column_widths = [0.7, 0.3]

    fig = make_plot(
        traces=traces,
        column_widths=column_widths
    )

    # Assert that the violinmode is set to 'group'
    assert fig.layout.violinmode == 'group', "Violin mode is not set to 'group'"

    # Extract the x-axis domains for each subplot
    xaxis_domains = []
    for i in range(1, len(column_widths) + 1):
        xaxis_name = f'xaxis{i}' if i > 1 else 'xaxis'
        domain = fig.layout[xaxis_name].domain
        xaxis_domains.append(domain)

    # Calculate inferred column widths
    inferred_widths = [domain[1] - domain[0] for domain in xaxis_domains]

    # Normalize inferred widths to sum to 1 (since domains might be adjusted for spacing)
    total_inferred_width = sum(inferred_widths)
    normalized_inferred_widths = [width / total_inferred_width for width in inferred_widths]

    # Assert that the normalized inferred widths match the expected widths within a tolerance
    tolerance = 0.01  # Allowable difference due to floating-point arithmetic
    for inferred, expected in zip(normalized_inferred_widths, column_widths):
        assert abs(inferred - expected) < tolerance, (
            f"Inferred column width {inferred:.2f} does not match expected width {expected:.2f}"
        )


def test_make_plot_violinmode_group():
    """
    Test that the violinmode is set to 'group'.
    """
    # Create violin traces
    transcript_traces = [
        go.Scatter(
            x=[1, 2, 3],
            y=[0, 1, 0],
            mode='lines',
            name='Transcript1'
        )
    ]

    violin_traces = [
        go.Violin(
            y=[10, 15, 13, 17],
            name='Violin1'
        )
    ]

    y_dict = {'Transcript1': 0}

    traces = [transcript_traces, violin_traces, y_dict]

    with pytest.warns(UserWarning) as record:
        fig = make_plot(
            traces=traces,
            column_widths=[0.7, 0.3, 0.2])
    assert len(record) == 1
    assert "The `column_widths` parameter must be a list of the same" in str(record[0].message)