import plotly.graph_objects as go
from plotly.subplots import make_subplots

file = open("Tasks/Records/Human_OGs by CDS records.txt", "r")

data_types_description = {"types": "<b>CDS/CDS overlaps by types</b>",
                          "ranks": "<b>CDS/CDS overlaps by maximum<br>expression ranks between pairs</b>",
                          "lengths": "<b>CDS/CDS overlaps by<br>overlapped lengths</b>"}

labels_by_data_types = {"ranks": ['Max rank = 1', 'Max rank = 2', 'Max rank = 3', 'Max rank = 4', 'Max rank > 4'],
                        "types": ['Convergent', 'Divergent', '→[←]*', '→[→]**', 'Tandem'],
                        "lengths": ['<60 bp', '60-120 bp', '120-180 bp', '180-254 bp']}
values_by_data_types = {"types": [0, 0, 0, 0, 0],
                        "ranks": [0, 0, 0, 0, 0],
                        "lengths": [0, 0, 0, 0]}

records_count = 0
for line in file.readlines():
    records_count += 1
    arr = line.split('\t')
    ov_type = arr[12]
    rank1 = arr[5]
    rank2 = arr[11]
    length = int(arr[13])
    if rank1 == '-': rank1 = 10
    if rank2 == '-': rank2 = 10
    rank = max(int(rank1), int(rank2))
    if ov_type == 'Convergent': values_by_data_types["types"][0] += 1
    if ov_type == 'Divergent': values_by_data_types["types"][1] += 1
    if ov_type == 'Nested (anti parallel)': values_by_data_types["types"][2] += 1
    if ov_type == 'Nested (parallel)': values_by_data_types["types"][3] += 1
    if ov_type == 'Tandem': values_by_data_types["types"][4] += 1
    if rank == 1: values_by_data_types["ranks"][0] += 1
    if rank == 2: values_by_data_types["ranks"][1] += 1
    if rank == 3: values_by_data_types["ranks"][2] += 1
    if rank == 4: values_by_data_types["ranks"][3] += 1
    if rank > 4: values_by_data_types["ranks"][4] += 1

    if length < 60: values_by_data_types["lengths"][0] += 1
    if 60 <= length < 120: values_by_data_types["lengths"][1] += 1
    if 120 <= length < 180: values_by_data_types["lengths"][2] += 1
    if length >= 180: values_by_data_types["lengths"][3] += 1

# Create subplots with a 1x3 grid


colors = ['#08519c', '#3182bd', '#6baed6', '#9ecae1', '#c6dbef']

fig = go.Figure()

annotations = []
for row in range(len(labels_by_data_types.keys())):
    data_type = list(labels_by_data_types.keys())[row]
    labels = labels_by_data_types[data_type]
    values = values_by_data_types[data_type]
    # # Sort arrays based on values in reverse order
    sorted_values, sorted_labels = zip(*sorted(zip(values, labels), reverse=True))

    # Print sorted arrays
    labels = sorted_labels
    values = sorted_values

    # Convert to percent values
    percent_values = []
    for index in range(len(values)):
        percent_values.append(values[index] * 100 // records_count)
    values = percent_values

    # Draw horizontal bars
    fig.add_trace(go.Bar(x=values, y=[data_type] * len(values), orientation='h',
                         marker=dict(color=colors, line=dict(color='rgb(248, 248, 249)', width=1))))
    space = 0
    for col in range(len(values)):
        value = values[col]
        label = labels[col]
        # Draw label for each bar
        annotations.append(dict(xref='x', yref='y',
                                x=space + value / 2, y=row + 0.4,
                                text=label,
                                font=dict(family='Arial', size=10,
                                          color='rgb(67, 67, 67)'),
                                showarrow=False))
        # Draw value inside each bar
        annotations.append(dict(xref='x', yref='y',
                                x=space + value / 2, y=data_type,
                                text=str(value) + "%",
                                font=dict(family='Arial', size=14,
                                          color='rgb(248, 248, 255)'),
                                showarrow=False))
        space += value

    # Draw left annotation for each horizontal bars
    annotations.append(dict(xref='paper', yref='y',
                            x=0.14, y=data_type,
                            xanchor='right',
                            text=data_types_description[data_type],
                            font=dict(family='Arial', size=10,
                                      color='rgb(67, 67, 67)'),
                            showarrow=False, align='right'))

fig.update_layout(
    xaxis=dict(
        showgrid=False,
        showline=False,
        showticklabels=False,
        zeroline=False,
        domain=[0.15, 1]
    ),
    yaxis=dict(
        showgrid=False,
        showline=False,
        showticklabels=False,
        zeroline=False,
    ),
    barmode='stack',
    paper_bgcolor='rgb(248, 248, 255)',
    plot_bgcolor='rgb(248, 248, 255)',
    margin=dict(l=120, r=10, t=140, b=80),
    showlegend=False,
)
fig.update_layout(bargap=0.6)
fig.update_layout(
    font_family="Times New Roman",
    font_color="black",
    legend_title_font_color="black",
    legend_tracegroupgap=185,
    height=400, width=1000,
)

# Double column (full width)	190 mm (539 pt)	2244 (300 dpi)
width = 700
#scale = (190 / 25.4) / (width / 300)

fig.update_layout(annotations=annotations)
fig.write_image('./Tasks/Records/figure.jpg',scale=3)

#fig.show()
