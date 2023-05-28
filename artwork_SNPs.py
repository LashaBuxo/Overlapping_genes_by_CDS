import random
import numpy as np
import plotly.graph_objects as go
from scipy.stats import stats, chi2_contingency, fisher_exact, mannwhitneyu

file_SNPs_data = open("Tasks/SNPs/SNPs_data.txt", "r")

dataset_names = ["entire", "nonORs", "ORs"]
dataset_full_names = {"entire": "Entire dataset",
                      "nonORs": "Non-overlapping regions",
                      "ORs": "Overlapping regions"}
parameter_names = ["length", "SNPs", "SSNPs", "SNPs_pos1", "SNPs_pos2", "SNPs_pos3"]
parameter_indexes = [0, 1, 2, 3, 4, 5]
parameter_full_names = {"SNPs": "Total SNPs",
                        "SSNPs": "Synonymous SNPs",
                        "SNPs_pos1": "SNPs at 1st position",
                        "SNPs_pos2": "SNPs at 2nd position",
                        "SNPs_pos3": "SNPs at 3rd position"}

colors = ['#6baed6', '#9ecae1', '#c6dbef', '#c6dbef', '#c6dbef']

data = {}
for dataset_name in dataset_names:
    data[dataset_name] = {}
    for parameter_name in parameter_names:
        data[dataset_name][parameter_name] = []

for line in file_SNPs_data:
    arr = line.replace('\n', '').split('\t')
    if arr[0].startswith("Gene"): continue
    # gene_sym, transcript_id = arr[0], arr[1]

    for i in range(len(parameter_indexes)):
        index = parameter_indexes[i]
        parameter_name = parameter_names[i]
        data["entire"][parameter_name].append(int(arr[2 + index]))

    for i in range(len(parameter_indexes)):
        index = parameter_indexes[i]
        parameter_name = parameter_names[i]
        data["nonORs"][parameter_name].append(int(arr[8 + index]))
    for i in range(len(parameter_indexes)):
        index = parameter_indexes[i]
        parameter_name = parameter_names[i]
        data["ORs"][parameter_name].append(int(arr[14 + index]))


def Bootstrapping(data_values, data_lengths, n=10000):
    freq_values = []
    total = 0
    target = 0
    for i in range(len(data_values)):
        total += data_lengths[i]
        target += data_values[i]
        freq_values.append(data_values[i] / data_lengths[i])
    data_mean = target / total

    replicates_means = []
    while n > 0:
        n -= 1
        total_length = 0
        target_length = 0
        for i in range(len(data_values)):
            index = random.randint(0, len(data_values) - 1)
            total_length += data_lengths[index]
            target_length += data_values[index]
        replicates_means.append(target_length / total_length)

    return float(f"{data_mean:.4f}"), np.percentile(replicates_means, [2.5, 97.5]), freq_values


fig = go.Figure()
annotations = []
bar_width = 0.8 / (len(parameter_names) - 1)

ORs_freqs_by_parameters = {}
NonORs_freqs_by_parameters = {}
i = 0
for parameter_name in parameter_names:
    if parameter_name == "length": continue
    means, upper_bounds, lower_bounds = [], [], []
    for dataset_name in dataset_names:
        mean, confidence, freqs = Bootstrapping(data[dataset_name][parameter_name], data[dataset_name]["length"])
        means.append(mean)
        lower_bounds.append(mean - confidence[0])
        upper_bounds.append(confidence[1] - mean)
        if dataset_name == "nonORs":
            NonORs_freqs_by_parameters[parameter_name] = freqs
        if dataset_name == "ORs":
            ORs_freqs_by_parameters[parameter_name] = freqs

    for j, (param, value) in enumerate(zip(dataset_names, means)):
        fig.add_annotation(
            x=j + (i - (len(parameter_names) - 2) / 2) * bar_width, y=(value - lower_bounds[j]) / 2,
            text=f"{float(f'{value * 100:.2f}')}"[0:], showarrow=False,
            font=dict(size=14, color='black'),
            xanchor='center', yanchor='middle'
        )

    fig.add_trace(go.Bar(
        marker=dict(color=colors[i]),
        legendgroup=parameter_name, name=parameter_full_names[parameter_name],
        x=list(dataset_full_names.values()), y=means,
        error_y=dict(
            type='data', symmetric=False, array=upper_bounds, arrayminus=lower_bounds, visible=True, thickness=2,
            width=8, color="#fdd0a2",
        )))
    i += 1

fig.update_layout(
    showlegend=True,
    legend=dict( font=dict(size=14, color='black')),
    xaxis=dict(showgrid=False, zeroline=False,  tickfont=dict(size=14,color='black')),
    yaxis=dict(showgrid=True,   zeroline=False, tickformat=".1%", tickfont=dict(size=14,color='black')),
    font_family="Times New Roman",
    height=540, width=960
)

width = 700
fig.write_image('./Tasks/SNPs/figure.jpg', scale=3)

for parameter_name in parameter_names:
    if parameter_name == "length": continue
    group1_frequencies = ORs_freqs_by_parameters[parameter_name]
    group2_frequencies = NonORs_freqs_by_parameters[parameter_name]

    # Perform the Mann-Whitney U test
    u_stat, u_p_value = mannwhitneyu(group1_frequencies, group2_frequencies)
    print(parameter_name)
    # Print the results
    # print("Mann-Whitney U test:")
    # print("U statistic:", u_stat)
    # print("P-value:", u_p_value)

    # Determine if the SNP frequency difference is significant
    alpha = 0.05
    if u_p_value < alpha:
        print("The difference in SNP frequencies between the two groups is significant (Mann-Whitney U test).")
    else:
        print("The difference in SNP frequencies between the two groups is not significant (Mann-Whitney U test).")
    print(u_p_value)
