import plotly.graph_objects as go

fig = go.Figure()

# Box trace 1
fig.add_trace(go.Box(
    y=[1, 2, 3],
    name='Group A',
    legendgroup='group1',
    legendgrouptitle_text='First Group',
    showlegend=True
))

# Box trace 2
fig.add_trace(go.Box(
    y=[2, 3, 4],
    name='Group B',
    legendgroup='group1',
    showlegend=True
))

# Box trace 3
fig.add_trace(go.Box(
    y=[3, 4, 5],
    name='Group C',
    legendgroup='group2',
    legendgrouptitle_text='Second Group',
    showlegend=True
))

# Box trace 4
fig.add_trace(go.Box(
    y=[4, 5, 6],
    name='Group D',
    legendgroup='group3',
    legendgrouptitle_text='Second Group',
    showlegend=False
))

fig.show()
