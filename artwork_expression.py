import math

import numpy as np
import scipy
from scipy.stats import spearmanr, pearsonr
import random
import numpy as np
import plotly.graph_objects as go
from scipy.stats import mannwhitneyu
from GenomePackage.Enums import SPECIES

# Brain, Heart, Kidney, Testis, Thyroid, Adipose Tissue, Ovary, Blood Vessel, Lung, Nerve, Muscle, Vagina, Breast, Liver Skin
target_tissues = ["Brain", "Heart", "Kidney", "Testis", "Thyroid",
                  "Ovary", "Lung", "Nerve", "Muscle", "Vagina", "Breast",
                  "Liver", "Skin"]


# region Helper Methods
def get_arrays_correlation(arr1, arr2):
    arr1 = np.array(arr1)
    arr2 = np.array(arr2)

    if len(arr1) == 0 or len(arr2) == 0:
        return 'x', 0, 0
    assert len(arr1) == len(arr2)
    r, p = pearsonr(arr1, arr2)
    r = float("{:.2f}".format(r))
    p = float("{:.5f}".format(p))

    if abs(r) < 0.5 or p > 0.05 or np.isnan(r) or np.isnan(p):
        return 'x', r, p, arr1 / arr2

    if r > 0:
        return '+', r, p, arr1 / arr2
    else:
        return '-', r, p, arr1 / arr2


file = open("./Tasks/Expression/NR1D1_THRA_expression.txt", "r")

data = {}
thra_transcript_expression_samples = {}
nr1d1_transcript_expression_samples = {}

index = 0
lines = file.readlines()
while index < len(lines):
    arr = lines[index].split('\t')
    if len(arr) == 4:
        tissue, val_thra_gene, val_thra_transcript, val_nr1d1_gene = arr[0], float(arr[1]), float(arr[2]), float(arr[3])
        if target_tissues.__contains__(tissue):
            # val_thra_gene=math.log2(val_thra_gene)
            # val_thra_transcript = math.log2(val_thra_transcript)
            # val_nr1d1_gene = math.log2(val_nr1d1_gene)
            data[tissue] = (val_thra_gene, val_thra_transcript, val_nr1d1_gene)
    else:
        arr = lines[index].split('\t')
        tissue = arr[0]
        thra_transcript_expression_samples[tissue] = []
        nr1d1_transcript_expression_samples[tissue] = []
        for i in range(1, len(arr), 2):
            thra_transcript_expression_samples[tissue].append(float(arr[i]))
            nr1d1_transcript_expression_samples[tissue].append(float(arr[i + 1]))

        status, r, p, ratio = get_arrays_correlation(thra_transcript_expression_samples[tissue],
                                                     nr1d1_transcript_expression_samples[tissue])

        # if status != 'x':
        #     print(f"{tissue} r={r} p={p} ratio={ratio}")
    index += 1
file.close()

# Calculate the ratio between expression levels for each tissue
ratio_data = {tissue: gene2 / (transcript_alpa_1 + transcript_alpa_2) for
              tissue, (transcript_alpa_1, transcript_alpa_2, gene2) in data.items()}

# Sort tissues by the ratio between expression levels
sorted_tissues = sorted(ratio_data, key=ratio_data.get, reverse=True)

# Print sorted tissues and corresponding ratios
for tissue in sorted_tissues:
    print(f"{tissue}: {ratio_data[tissue]}")

# Separate gene expression levels for sorted tissues
gene1_expression = []
transcript_alpa1_expression = [data[tissue][0] for tissue in sorted_tissues]
transcript_alpa2_expression = [data[tissue][1] for tissue in sorted_tissues]

ratios=[]
for index in range(len(transcript_alpa1_expression)):
    gene1_expression.append(transcript_alpa1_expression[index] + transcript_alpa2_expression[index])
    ratios.append(transcript_alpa1_expression[index]+transcript_alpa2_expression[index])
gene2_expression = [data[tissue][2] for tissue in sorted_tissues]

print(ratios)
print(gene2_expression)
correlation_coefficient, p_value = pearsonr(ratios, gene2_expression)

print(f"Pearson correlation coefficient: {correlation_coefficient:.3f}")
print(f"P-value: {p_value:.6f}")
# Create horizontal bar plot with side-by-side bars as described
fig = go.Figure()

fig.add_trace(go.Bar(y=sorted_tissues, x=[-val for val in gene1_expression],
                     marker_color='#d6f0fa', name='THRA (TRα1)', orientation='h', offset=-0.45))
fig.add_trace(go.Bar(y=sorted_tissues, x=[-val for val in transcript_alpa2_expression],
                     marker_color='#b5e3f5', name='THRA (TRα2)', orientation='h', offset=-0.45))
fig.add_trace(go.Bar(y=sorted_tissues, x=gene2_expression,
                     marker_color='#e6f1d3', name='NR1D1 (Rev-erbα)', orientation='h', offset=-0.45))

fig.add_trace(go.Scatter(y=sorted_tissues, x=[-val - 10 for val in gene1_expression], textposition="middle left",
                         text=[
                             f'{(gene1_expression[index]) * 100 / (gene1_expression[index] + gene2_expression[index]):.0f}%'
                             for index in range(len(transcript_alpa1_expression))],
                         textfont=dict(size=14, color='black'),
                         mode='text', showlegend=False, ))
fig.add_trace(go.Scatter(y=sorted_tissues, x=[val + 10 for val in gene2_expression], textposition="middle right",
                         text=[
                             f'{gene2_expression[index] * 100 / (gene1_expression[index] + gene2_expression[index]):.0f}%'
                             for index in range(len(transcript_alpa1_expression))],
                         textfont=dict(size=14, color='black'),
                         mode='text', showlegend=False))

fig.update_layout(
    plot_bgcolor='rgba(0,0,0,0)',
    font_family="Times New Roman",
    xaxis_title='TPM (mean)',
    barmode='group',  # Side-by-side bars
    bargap=0.1,  # Gap between bars
    xaxis=dict(range=[-200, 350], tickfont=dict(size=14, color='black' ),titlefont=dict(size=18, color='black'),),
    yaxis=dict( tickfont=dict(size=14, color='black' ),titlefont=dict(size=14, color='black'),),

    height=440, width=960,
    legend=dict(x=0.81, y=1, bgcolor='rgba(0,0,0,0)',font=dict(size=14, color='black'))
)
fig.update_xaxes(showline=False, linewidth=1, linecolor='#cad3d6', showgrid=True, gridcolor='#cad3d6')
# Add lines between texts and bars
for i, tissue in enumerate(sorted_tissues):
    fig.add_shape(type='line', x0=-gene1_expression[i], x1=-gene1_expression[i] - 10, y0=i, y1=i, yref='y', xref='x',
                  line=dict(color='black', width=1))
    fig.add_shape(type='line', x0=gene2_expression[i], x1=gene2_expression[i] + 10, y0=i, y1=i, yref='y', xref='x',
                  line=dict(color='black', width=1))

width = 700
fig.write_image('./Tasks/Expression/figure.jpg', scale=3)
