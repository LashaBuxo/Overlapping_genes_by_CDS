import random
import numpy as np
import plotly.graph_objects as go
from scipy.stats import mannwhitneyu

from GenomePackage.Enums import SPECIES

species_list = [
    SPECIES.Gallus_gallus.short_name(),
    SPECIES.Danio_rerio.short_name(),
]

records_count = 0
records_in_target_species = {}

file = open("./Tasks/Conservation/conservation_overlapping.txt", "r")

for line in file.readlines():
    arr = line.split('\t')

    if arr[0].__contains__("CDS_"):
        records_count += 1
    for i in range(1, len(arr), 2):
        species = arr[i]
        identity = float(arr[i + 1]) / 100
        if species_list.__contains__(species):
            if not arr[0].__contains__("RANDOM"):
                if not arr[0].__contains__("NR1D1THRA") and not arr[0].__contains__("NR1D1exon7") and not arr[
                    0].__contains__("THRAexon9"):
                    records_in_target_species[arr[0].split('_')[1]] = True
                    records_in_target_species[arr[0].split('_')[2]] = True
file.close()

file_disorder_data = open("Tasks/Disorder/disorder_data.txt", "r")

entire_transcripts_CDS_data = []
overlapping_regions_data = []
non_overlapping_regions_data = []

overlapping_regions_data_old = []
overlapping_regions_data_new = []


def Bootstrapping(data, n=10000):
    freq_values = []
    total = 0
    target = 0
    for i in range(len(data)):
        total += data[i][1]
        target += data[i][0]
        freq_values.append(data[i][0] / data[i][1])
    data_mean = target / total

    replicates_means = []
    while n > 0:
        n -= 1
        total_length = 0
        target_length = 0
        for i in range(len(data)):
            index = random.randint(0, len(data) - 1)
            total_length += data[index][1]
            target_length += data[index][0]
        replicates_means.append(target_length / total_length)

    return float(f"{data_mean:.2f}"), np.percentile(replicates_means, [2.5, 97.5]), freq_values


for line in file_disorder_data:
    arr = line.replace('\n', '').split('\t')
    if arr[0].startswith("Gene"): continue
    gene_sym, transcript_id, entire_length, entire_value = arr[0], arr[1], int(arr[2]), int(arr[3])
    non_overlapping_length, non_overlapping_values = int(arr[4]), int(arr[5])
    overlapping_length, overlapping_values = int(arr[6]), int(arr[7])

    entire_transcripts_CDS_data.append((entire_value, entire_length))
    non_overlapping_regions_data.append((non_overlapping_values, non_overlapping_length))
    overlapping_regions_data.append((overlapping_values, overlapping_length))

    if records_in_target_species.__contains__(transcript_id):
        overlapping_regions_data_old.append((overlapping_values, overlapping_length))
    else:
        overlapping_regions_data_new.append((overlapping_values, overlapping_length))

entire_mean, entire_conf, freq_entire = Bootstrapping(entire_transcripts_CDS_data)
non_overlapping_mean, non_overlapping_conf, freq_nonORs = Bootstrapping(non_overlapping_regions_data)
overlapping_mean, overlapping_conf, freq_ORs = Bootstrapping(overlapping_regions_data)

overlapping_mean_old, overlapping_conf_old, freq_ORs_old = Bootstrapping(overlapping_regions_data_old)
overlapping_mean_new, overlapping_conf_new, freq_ORs_new = Bootstrapping(overlapping_regions_data_new)

# Replace these values with your data
values = [entire_mean, non_overlapping_mean, overlapping_mean]

# Replace these confidence intervals with your data
conf_intervals = [entire_conf, non_overlapping_conf, overlapping_conf]
# Separate the lower and upper bounds
lower_bounds = [v - lower for v, (lower, _) in zip(values, conf_intervals)]
upper_bounds = [upper - v for v, (_, upper) in zip(values, conf_intervals)]

fig = go.Figure()

# Create custom box plots using go.Bar
fig.add_trace(go.Bar(
    x=['Entire dataset', 'Non-overlapping regions', 'Overlapping regions'],
    y=values,
    width=[0.3, 0.3, 0.3],  # Adjust the width of the bars
    marker_color='#6baed6',  # Adjust the color of the bars
    error_y=dict(
        type='data',
        symmetric=False,
        array=upper_bounds,
        arrayminus=lower_bounds,
        color="#fdd0a2",
        visible=True,
        thickness=2,  # Adjust the thickness of the error bars
        width=8  # Adjust the width of the error bars
    )
))
# Add text labels for the mean values
fig.add_trace(go.Scatter(
    x=['Entire dataset', 'Non-overlapping regions', 'Overlapping regions'],
    y=[values[0] / 2, values[1] / 2, values[2] / 2],  # Adjust y-coordinate for the text labels
    mode='text',
    text=[f"{int(values[0] * 100)}%", f"{int(values[1] * 100)}%", f"{int(values[2] * 100)}%"],
    textposition='middle center',
    textfont=dict(size=16, color='black')
))

fig.update_layout(
    # xaxis_title='Values',
    yaxis_title='Disorder Content',
    showlegend=False,
    xaxis=dict(
        tickangle=0,  # Set the angle of y-axis tick labels to 0 degrees
        showgrid=False,
        zeroline=False,
        tickfont=dict(size=15, color='black'),
    ),

    yaxis=dict(
        tickangle=0,  # Set the angle of y-axis tick labels to 0 degrees
        titlefont=dict(size=20, color='black'),
        tickfont=dict(size=16, color='black'),
        showgrid=True,
        zeroline=False,
        tickformat=".0%"  # Format the y-axis values as percentages
    ),
    font_family="Times New Roman",
    height=400, width=600,
)

# Double column (full width)	190 mm (539 pt)	2244 (300 dpi)
width = 700
# scale = (190 / 25.4) / (width / 300)
fig.write_image('./Tasks/Disorder/figure.jpg', scale=3)

# Replace these values with your data
values = [overlapping_mean_new, overlapping_mean_old]

# Replace these confidence intervals with your data
conf_intervals = [overlapping_conf_new, overlapping_conf_old]
# Separate the lower and upper bounds
lower_bounds = [v - lower for v, (lower, _) in zip(values, conf_intervals)]
upper_bounds = [upper - v for v, (_, upper) in zip(values, conf_intervals)]

fig = go.Figure()

# Create custom box plots using go.Bar
fig.add_trace(go.Bar(
    x=['Mammalian overlap regions', 'Older overlap regions'],
    y=values,
    width=[0.3, 0.3, 0.3],  # Adjust the width of the bars
    marker_color='#6baed6',  # Adjust the color of the bars
    error_y=dict(
        type='data',
        symmetric=False,
        array=upper_bounds,
        arrayminus=lower_bounds,
        color="#fdd0a2",
        visible=True,
        thickness=2,  # Adjust the thickness of the error bars
        width=8  # Adjust the width of the error bars
    )
))
# Add text labels for the mean values
fig.add_trace(go.Scatter(
    x=['Mammalian overlap regions', 'Older overlap regions'],
    y=[values[0] / 2-0.05, values[1] / 2-0.05],  # Adjust y-coordinate for the text labels
    mode='text',
    text=[f"{int(values[0] * 100)}%", f"{int(values[1] * 100)}%"],
    textposition='middle center',
    textfont=dict(size=16, color='black'),
))

fig.update_layout(
    # xaxis_title='Values',
    yaxis_title='Disorder Content',
    showlegend=False,
    xaxis=dict(
        showgrid=False,
        zeroline=False,
        tickfont=dict(size=16, color='black')
    ),

    yaxis=dict(
        showgrid=True,
        zeroline=False,
        titlefont=dict(size=20, color='black'),
        tickfont=dict(size=16, color='black'),
        tickformat=".0%"  # Format the y-axis values as percentages
    ),
    font_family="Times New Roman",
    height=400, width=600,
)

# Double column (full width)	190 mm (539 pt)	2244 (300 dpi)
width = 700
# scale = (190 / 25.4) / (width / 300)
fig.write_image('./Tasks/Disorder/figure_old_vs_new.jpg', scale=3)

from scipy import stats

print(np.mean(freq_ORs))
print(np.mean(freq_nonORs))
# Perform Mann-Whitney U test
u_value, p_value = stats.mannwhitneyu(freq_ORs, freq_nonORs)

# Check if p-value is less than significance level (0.05)
if p_value < 0.05:
    print(f"Reject the null hypothesis: U-value = {u_value}, one-tailed p-value = {p_value}")
else:
    print(f"Fail to reject the null hypothesis: U-value = {u_value}, one-tailed p-value = {p_value}")

print(np.mean(freq_ORs_new))
print(np.mean(freq_ORs_old))
# Perform Mann-Whitney U test
u_value, p_value = stats.mannwhitneyu(freq_ORs_new, freq_ORs_old)

# Check if p-value is less than significance level (0.05)
if p_value < 0.05:
    print(f"Reject the null hypothesis: U-value = {u_value}, one-tailed p-value = {p_value}")
else:
    print(f"Fail to reject the null hypothesis: U-value = {u_value}, one-tailed p-value = {p_value}")
