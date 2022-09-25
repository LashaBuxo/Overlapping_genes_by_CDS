import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.figure_factory as ff

types_dict = {"Convergent": "Convergent", "Divergent": "Divergent", "Nested (diff strand)": "Nested (diff)",
              "Tandem": "Tandem", "Nested (same strand)": "Nested (same)",
              "MULTI TYPE": "Multi Type<b><sup>1</sup></b>"}

transcripts_count_per_type = {}
exons_count_per_type = {}
for type in list(types_dict.values()):
    transcripts_count_per_type[type] = 0
    exons_count_per_type[type] = 0

bin_size = 25
Overlaps_length = []
Conserved_overlaps_length = []

file = open("generated_data/OGs by CDS records.txt", "r")
lines = file.readlines()

for index in range(1, len(lines)):
    arr = lines[index].split('\t')
    transcriptic_overlap_type, length, exonic_overlap_typ, cons1, cons2 = arr[4], int(arr[5]), arr[7], float(
        arr[9]), float(arr[10])
    transcriptic_overlap_type = types_dict[transcriptic_overlap_type]
    exonic_overlap_type = types_dict[exonic_overlap_typ]

    if not transcripts_count_per_type.__contains__(transcriptic_overlap_type):
        transcripts_count_per_type[transcriptic_overlap_type] = 0
    if not exons_count_per_type.__contains__(exonic_overlap_type):
        exons_count_per_type[exonic_overlap_type] = 0

    transcripts_count_per_type[transcriptic_overlap_type] += 1
    exons_count_per_type[exonic_overlap_type] += 1


    # if min(cons1, cons2) > 0:
    if transcriptic_overlap_type == 'Divergent' or transcriptic_overlap_type == 'Convergent' or transcriptic_overlap_type == 'Nested (diff strand)':

        Overlaps_length.append(length)
    else:
        Conserved_overlaps_length.append(length)

fig = make_subplots(
    rows=2, cols=1,
    subplot_titles=("(a)", "(b)"),
    vertical_spacing=0.1)

# Make traces for graph
trace1 = go.Bar(x=list(types_dict.values()), y=list(transcripts_count_per_type.values()), xaxis='x2', yaxis='y2',
                marker=dict(color='#0099ff'),
                name='Transcripts', legendgroup='1')
trace2 = go.Bar(x=list(types_dict.values()), y=list(exons_count_per_type.values()), xaxis='x2', yaxis='y2',
                marker=dict(color='#404040'),
                name='Exons', legendgroup='1')

fig.add_trace(trace1, row=1, col=1)
fig.add_trace(trace2, row=1, col=1)

fig.update_yaxes(title_text="Quantity", row=1, col=1)
fig.update_yaxes(title_text="Quantity", row=2, col=1)
fig.update_xaxes(title_text="Overlapped CDS Length", row=2, col=1)
fig.update_xaxes(dict(tickmode='linear', tick0=20, dtick=20), row=2, col=1)

fig.update_xaxes(
        tickmode = 'array',
        tickvals = [15,45,75,105,135,165,195,225,255],
        ticktext = ['0-30', '30-60', '60-90', '90-120', '120-150', '150-180', '180-210', '210-240','240-254']
, row=2, col=1
)



fig.add_trace(go.Histogram(x=Overlaps_length, autobinx = False,xbins=dict(start=0,end=260, size=30), legendgroup='2', name='Diff. strands<br>overlaps', ), row=2, col=1)
fig.add_trace(go.Histogram(x=Conserved_overlaps_length, autobinx = False,xbins=dict(start=0,end=260, size=30), legendgroup='2', name='Same strands<br>overlaps', ), row=2,
              col=1)

fig.layout.annotations[0].update(x=0.025)
fig.layout.annotations[1].update(x=0.025)
fig.update_layout(
    font_family="Times New Roman",
    font_color="black",
    legend_title_font_color="black",
    legend_tracegroupgap=185,
    height=400, width=600,
    margin=go.layout.Margin(l=0, r=0, b=40, t=25, )
)

# Double column (full width)	190 mm (539 pt)	2244 (300 dpi)
width = 700
scale = (190 / 25.4) / (width / 300)
fig.write_image('/generated_data/figure.jpg', scale=scale)
