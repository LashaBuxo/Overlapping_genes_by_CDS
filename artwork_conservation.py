import numpy as np
import random
import numpy as np
import plotly.graph_objects as go
from scipy.stats import mannwhitneyu, ttest_1samp
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px
from GenomePackage.Enums import SPECIES

RANDOM_SAMPLES_PER_RECORD = 100
RANDOM_SAMPLES_PER_NR1D1_THRA = 1000

species_list = [SPECIES.Pan_troglodytes.short_name(), SPECIES.Sus_scrofa.short_name(),
                SPECIES.Canis_lupus_familiaris.short_name(),
                SPECIES.Mus_musculus.short_name(),
                SPECIES.Monodelphis_domestica.short_name(),
                SPECIES.Gallus_gallus.short_name(),
                SPECIES.Danio_rerio.short_name(),
                SPECIES.Drosophila_melanogaster.short_name()]

records_count = 0
OGs_in_species = {}
Control_OGs_in_species = []

NR1D1THRA_values = {}
NR1D1exon7_values = {}
THRAexon9_values = {}

NR1D1THRA_control = {}
NR1D1exon7_control = {}
THRAexon9_control = {}

for i in range(RANDOM_SAMPLES_PER_RECORD):
    Control_OGs_in_species.append({})
for species in species_list:
    OGs_in_species[species] = 0
    for i in range(RANDOM_SAMPLES_PER_RECORD):
        Control_OGs_in_species[i][species] = 0

    NR1D1THRA_values[species] = 0
    NR1D1exon7_values[species] = 0
    THRAexon9_values[species] = 0
    for i in range(RANDOM_SAMPLES_PER_NR1D1_THRA):
        NR1D1THRA_control[species] = []
        NR1D1exon7_control[species] = []
        THRAexon9_control[species] = []

file = open("./Tasks/Conservation/conservation_overlapping.txt", "r")

for line in file.readlines():
    arr = line.split('\t')

    if arr[0].__contains__("CDS_"):
        records_count += 1
    for i in range(1, len(arr), 2):
        species = arr[i]
        identity = float(arr[i + 1]) / 100
        if species_list.__contains__(species):
            if arr[0].__contains__("RANDOM"):
                if arr[0].__contains__("NR1D1THRA"):
                    NR1D1THRA_control[species].append(identity)
                elif arr[0].__contains__("NR1D1exon7"):
                    NR1D1exon7_control[species].append(identity)
                elif arr[0].__contains__("THRAexon9"):
                    THRAexon9_control[species].append(identity)
                else:
                    rand_group_index = int(arr[0].split('_')[0].replace('RANDOM', ''))
                    Control_OGs_in_species[rand_group_index][species] += 1
            else:
                if arr[0].__contains__("NR1D1THRA"):
                    NR1D1THRA_values[species] = identity
                elif arr[0].__contains__("NR1D1exon7"):
                    NR1D1exon7_values[species] = identity
                elif arr[0].__contains__("THRAexon9"):
                    THRAexon9_values[species] = identity
                else:
                    OGs_in_species[species] += 1
file.close()

control_values = {}
for species in species_list:
    arr = []
    for i in range(RANDOM_SAMPLES_PER_RECORD):
        arr.append(Control_OGs_in_species[i][species])
    control_values[species] = np.mean(arr)

for species in species_list:
    control_scores = []
    for i in range(100):
        cons_score = Control_OGs_in_species[i][species]/ records_count
        control_scores.append(cons_score)

    print(f"{species} {OGs_in_species[species] / records_count} {np.mean(control_scores)} {np.std(control_scores)}")

    # Perform a one-sample t-test
    t_statistic, p_value = ttest_1samp(control_scores, OGs_in_species[species] / records_count)

    # Set your significance level (e.g., 0.05 for a 95% confidence level)
    alpha = 0.05

    # Test the hypothesis
    if p_value > alpha:
        print(f"not significantly different p={p_value}")
    else:
        if t_statistic > 0:
            print(f"significantly lower; p={p_value}")
        else:
            print(f"significantly greater; p={p_value}")

# region Create Plot for general conservations
fig = go.Figure()

fig.add_trace(go.Bar(
    x=species_list, y=[OGs_in_species[species] * 100 // records_count for species in species_list],
    name='overlapping CDSs', marker_color='#6baed6',textfont=dict(size=14,color='black'), textposition='auto',
    text=[f"{p:.0f}%" for p in [OGs_in_species[species] * 100 // records_count for species in species_list]],
))

fig.add_trace(go.Bar(
    x=species_list, y=[control_values[species] * 100 // records_count for species in species_list],
    name='control CDSs', marker_color='#c6dbef',textfont=dict(size=14,color='black'), textposition='auto',
    text=[f"{p:.0f}%" for p in [control_values[species] * 100 // records_count for species in species_list]],
))

# Customize the layout
fig.update_layout(
    font_family="Times New Roman", yaxis=dict(showticklabels=False),xaxis=dict(
        tickfont=dict(size=14,color='black')
    ),
    legend_title="Conserved amount (%) of", barmode='group',
    legend=dict(x=0.75, y=1, bgcolor='rgba(0,0,0,0)',
                font=dict(size=14,color='black')),
    height=440, width=960,
)



fig.write_image('./Tasks/Conservation/conservation_overlapping_CDSs.jpg', scale=5)
# endregion

control_genes_df = pd.DataFrame([(k, v) for k, vs in NR1D1THRA_control.items() for v in vs],
                                columns=['species', 'conservation_identity'])
control_genes_df['type'] = 'control'

NR1D1THRA_df = pd.DataFrame([(k, v) for k, v in NR1D1THRA_values.items()], columns=['species', 'conservation_identity'])
NR1D1THRA_df['type'] = 'NR1D1THRA'

NR1D1exon7_df = pd.DataFrame([(k, v) for k, v in NR1D1exon7_values.items()],
                             columns=['species', 'conservation_identity'])
NR1D1exon7_df['type'] = 'NR1D1exon7'

THRAexon9_df = pd.DataFrame([(k, v) for k, v in THRAexon9_values.items()], columns=['species', 'conservation_identity'])
THRAexon9_df['type'] = 'THRAexon9'
data = pd.concat([control_genes_df, NR1D1THRA_df, NR1D1exon7_df, THRAexon9_df], ignore_index=True)

annotations = []


def create_annotations(target_gene_df, off_val):
    for index, row in target_gene_df.iterrows():
        if row['species'] == "Fruit fly": continue
        val = f"{float(row['conservation_identity'] * 100):.1f}%"
        annotations.append(
            go.layout.Annotation(
                x=row['species'], y=row['conservation_identity'],
                xshift=5, yshift=0, text=val, showarrow=True,
                arrowhead=2, arrowsize=1,
                ax=40, ay=off_val, font=dict(size=12),
            )
        )


fig = go.Figure()

fig.add_trace(go.Box(
    x=data.loc[data['type'] == 'control', 'species'],
    y=data.loc[data['type'] == 'control', 'conservation_identity'],
    name='control CDSs',
    marker_color='#6baed6',
))


fig.add_trace(go.Scatter(
    x=list(data.loc[data['type'] == 'NR1D1exon7', 'species'])[0:-1],
    y=list(data.loc[data['type'] == 'NR1D1exon7', 'conservation_identity'])[0:-1],
    mode='markers',
    marker=dict(size=7, color='red', line=dict(width=1, color='DarkSlateGrey')),
    name='NR1D1 7th exon',
))
fig.add_trace(go.Scatter(
    x=list(data.loc[data['type'] == 'THRAexon9', 'species'])[0:-1],
    y=list(data.loc[data['type'] == 'THRAexon9', 'conservation_identity'])[0:-1],

    mode='markers',
    marker=dict(size=7, color='blue', line=dict(width=1, color='DarkSlateGrey')),
    name='THRA 8th exon',
))




fig.add_trace(go.Scatter(
    x=list(data.loc[data['type'] == 'NR1D1THRA', 'species'])[0:-1],
    y=list(data.loc[data['type'] == 'NR1D1THRA', 'conservation_identity'])[0:-1],
    mode='markers',
    marker=dict(size=7, color='#fdd0a2', line=dict(width=1, color='DarkSlateGrey')),
    name='NR1D1-THRA overlapping CDS',
))
#create_annotations(NR1D1THRA_df, 0)
# create_annotations(NR1D1exon7_df, 0)
# create_annotations(THRAexon9_df, 0)

fig.update_layout(annotations=annotations)

fig.update_layout(
    font_family="Times New Roman",
    xaxis=dict( tickfont=dict(size=14,color='black')),
    yaxis=dict(tickformat=".0%", tickfont=dict(size=14,color='black')),
    legend_title="Blast Identities (%) of",
    legend=dict(x=1.01, y=1, bgcolor='rgba(0,0,0,0)', font=dict(size=14,color='black')),
    height=440, width=960,
)

fig.write_image('./Tasks/Conservation/conservation_NR1D1_THRA.jpg', scale=5)
