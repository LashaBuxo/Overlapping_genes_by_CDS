import numpy as np
import scipy.stats as stats

import plotly.graph_objects as go

file = open("./Tasks/GC/overlapping_GCs.txt", "r")

overlapping_values = []
random_mean_values = []
semi_random_mean_values = []

for line in file:
    arr = line.replace('\n', '').split('\t')
    if arr[0].startswith("Gene"): continue
    gc = float(arr[2])
    overlapping_values.append(gc)

file = open("./Tasks/GC/control_groups_mean_GCs.txt", "r")
for line in file:
    arr = line.replace('\n', '').split('\t')
    if arr[0].startswith("Random"): continue
    gc_random = float(arr[0])
    gc_semi_random = float(arr[1])
    random_mean_values.append(gc_random)
    semi_random_mean_values.append(gc_semi_random)

mean = np.mean(overlapping_values)
std_dev = np.std(overlapping_values, ddof=1)
sample_size = len(overlapping_values)
sem = std_dev / np.sqrt(sample_size)
confidence_level = 0.95
degrees_of_freedom = sample_size - 1
critical_value = stats.t.ppf((1 + confidence_level) / 2, degrees_of_freedom)
# Calculate the confidence interval error bars
error_bars = critical_value * sem
print(error_bars)

gc_content_target_sequences = overlapping_values  # Replace with your GC content values
gc_content_control_group = random_mean_values  # Replace with control group's GC content values
gc_content_control_group2 = semi_random_mean_values  # Replace with control group 2's GC content values

# Perform one-tailed two-sample t-test
t_stat, p_value_two_tailed = stats.ttest_ind(gc_content_target_sequences, gc_content_control_group2)
p_value_one_tailed = p_value_two_tailed / 2

# Print results
print(f"t-statistic: {t_stat:.3f}")
print(f"one-tailed p-value: {p_value_one_tailed:.6f}")

# Interpret p-value
alpha = 0.05
if p_value_one_tailed < alpha and t_stat > 0:
    print("Reject the null hypothesis: The GC content of target sequences is significantly higher.")
else:
    print("Do not reject the null hypothesis: The GC content of target sequences is not significantly higher.")
# Calculate means
means = [
    np.mean(gc_content_target_sequences),
    np.mean(gc_content_control_group),
    np.mean(gc_content_control_group2)
]

# Calculate sample sizes
sample_sizes = [
    len(gc_content_target_sequences),
    len(gc_content_control_group),
    len(gc_content_control_group2)
]

# Calculate standard errors of the mean (SEM)
sem = [
    stats.sem(gc_content_target_sequences),
    stats.sem(gc_content_control_group),
    stats.sem(gc_content_control_group2)
]

# Calculate t-values for a 95% confidence interval
t_values = [stats.t.ppf(1 - 0.025, n - 1) for n in sample_sizes]

# Calculate 95% confidence intervals
ci_95 = [t * s for t, s in zip(t_values, sem)]
print(ci_95)

fig = go.Figure()

# Create custom box plots using go.Bar
fig.add_trace(go.Bar(
    x=['Overlapping regions', 'Control dataset', "'Dual-Stop-Free' regions"],
    y=means,
    width=[0.3, 0.3, 0.3],  # Adjust the width of the bars
    marker_color='#6baed6',  # Adjust the color of the bars
    error_y=dict(
        type='data',
        color="#fdd0a2",
        symmetric=False,
        array=ci_95,
        arrayminus=ci_95,
        visible=True,
        thickness=2,  # Adjust the thickness of the error bars
        width=8  # Adjust the width of the error bars
    )
))
# Add text labels for the mean values
fig.add_trace(go.Scatter(
    x=['Overlapping regions', 'Control dataset', "'Dual-Stop-Free' regions"],
    y=[0.45+(value - 0.45) / 2 for value in means],
    mode='text',
    text=[f"{int(means[0] * 100)}%", f"{int(means[1] * 100)}%", f"{int(means[2] * 100)}%"],
    textposition='middle center',
    textfont=dict(size=16, color='black' ),
))

fig.update_layout(
    yaxis_title='GC Content',
    showlegend=False,
    xaxis=dict(
        zeroline=False,
        tickfont=dict(size=16, color='black' )
    ),
    yaxis=dict(
        range=[0.45, 0.65],
        zeroline=False,
        titlefont=dict(size=20, color='black'),
        tickfont=dict(size=16, color='black'),
        tickformat=".0%"  # Format the y-axis values as percentages
    ),
    font_family="Times New Roman",
    height=400, width=600,
)

width = 700
fig.write_image('./Tasks/GC/figure.jpg', scale=3)
