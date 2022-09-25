import random
import statistics
import plotly.graph_objects as go
import pandas as pd
import plotly.express as px

from graph import *
from worker_genome import *
from CDS_Overlap_Record import *

genome = GenomeWorker(SPECIES.Homo_sapiens, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS, SEQUENCE_LOAD.LOAD)

records, ATI_overlaps, clusters_graph = get_CDS_records_from_genome(genome, with_clustered_graph=True)

# region output all OGs by CDS records
records.sort(key=lambda x: x.overlapped_length, reverse=True)

file = open("generated_data/OGs by CDS records.txt", "w")

OGs_by_CDS = {}
for record in records:
    OGs_by_CDS[record.gene1_sym] = True
    OGs_by_CDS[record.gene2_sym] = True

clusters = clusters_graph.get_connected_clusters()
OGs_by_cluster_size = {}
for cluster in clusters:
    OGs_by_cluster_size[len(cluster)] = OGs_by_cluster_size.get(len(cluster), 0) + 1
OGs_count_formatted = ''
for size in sorted(list(OGs_by_cluster_size.keys())):
    OGs_count_formatted += f'{" +" if len(OGs_count_formatted) > 0 else ""} {OGs_by_cluster_size[size]}x[{size}]'

file.write(f"{len(OGs_by_CDS)} genes {OGs_count_formatted} overlapped by CDS (without ATI):\n")
for record in records:
    file.write(f"{record.get_tabbed_data_string()}\n")
file.close()

file = open("generated_data/OGs by CDS (ATI genes).txt", "w")
file.write(f"{len(ATI_overlaps)} genes are overlapped by CDS and are ATI:\n")
for gene_sym in ATI_overlaps.keys():
    file.write(f"{gene_sym}\n")
file.close()
# endregion


transcripts_by_cds_length = []
for chr_id in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_id)
    for i in range(0, genes_cnt):
        gene = genome.gene_by_indexes(chr_id, i)
        # if overlapped.__contains__(gene.id): continue
        transcripts = genome.get_transcripts_from_gene(gene.id)
        for transcript in transcripts:
            CDS_length = genome.get_transcript_CDS_length(transcript.id)
            transcripts_by_cds_length.append((transcript.id, CDS_length))
transcripts_by_cds_length.sort(key=lambda x: x[1])


def get_control_transcript(transcript_id):
    observed_length = genome.get_transcript_CDS_length(transcript_id)
    desired_lower_limit, desired_upper_limit = observed_length, observed_length * 1.1
    l = 0
    r = len(transcripts_by_cds_length) - 1
    while l < r:
        mid = (l + r) // 2
        if transcripts_by_cds_length[mid][1] < desired_lower_limit:
            l = mid + 1
        else:
            r = mid
    index_l = l
    l = 0
    r = len(transcripts_by_cds_length) - 1
    while l < r:
        mid = (l + r) // 2
        if transcripts_by_cds_length[mid][1] < desired_upper_limit:
            l = mid + 1
        else:
            r = mid
    index_r = l

    control_transcript_id = transcripts_by_cds_length[random.randint(index_l, index_r)][0]
    return control_transcript_id


# region plot GC (overlapped vs control)

RECORD_MIN_LENGTH = 60
RUNS = 1000

overlapped_gc_contents = {}
for ov_type in OVERLAP_TYPE.get_overlap_types():
    overlapped_gc_contents[ov_type] = []

for record in records:
    if record.overlapped_length < RECORD_MIN_LENGTH: continue
    overlapped_gc_contents[OVERLAP_TYPE.MULTI].append(record.overlapped_GC_content)
    overlapped_gc_contents[OVERLAP_TYPE.MULTI].append(record.overlapped_GC_content)
    for ov_type in OVERLAP_TYPE.get_overlap_types():
        if record.transcripts_overlap_type == ov_type.short_name():
            overlapped_gc_contents[ov_type].append(record.overlapped_GC_content)
            overlapped_gc_contents[ov_type].append(record.overlapped_GC_content)

overlapped_gc_mean = {}
for ov_type in OVERLAP_TYPE.get_overlap_types():
    overlapped_gc_mean[ov_type] = statistics.mean(overlapped_gc_contents[ov_type])

runs_control_gc_means = {}
under_h0 = {}
for ov_type in OVERLAP_TYPE.get_overlap_types():
    runs_control_gc_means[ov_type] = []
    under_h0[ov_type] = 0

for run in range(RUNS):
    control_GC_contents = {}
    for ov_type in OVERLAP_TYPE.get_overlap_types():
        control_GC_contents[ov_type] = []

    for record in records:
        if record.overlapped_length < RECORD_MIN_LENGTH: continue
        transcript_ids_and_cds_start_indexes = [(record.transcript1.id, record.overlap_start_index_on_CDS1),
                                                (record.transcript2.id, record.overlap_start_index_on_CDS2)]
        for transcript_id, CDS_start_index in transcript_ids_and_cds_start_indexes:
            control_transcript_id = get_control_transcript(transcript_id)
            control_cds = genome.get_transcript_CDS(control_transcript_id)
            control_GC = genome.sequence_GC(control_cds[CDS_start_index:])

            control_GC_contents[OVERLAP_TYPE.MULTI].append(control_GC)
            for ov_type in OVERLAP_TYPE.get_overlap_types():
                if record.transcripts_overlap_type == ov_type.short_name():
                    control_GC_contents[ov_type].append(control_GC)

    for ov_type in OVERLAP_TYPE.get_overlap_types():
        runs_control_gc_means[ov_type].append(statistics.mean(control_GC_contents[ov_type]))
        under_h0[ov_type] += 1 if statistics.mean(control_GC_contents[ov_type]) > overlapped_gc_mean[ov_type] else 0

file = open("generated_data/GC control vs observed table.txt", "w")

for ov_type in OVERLAP_TYPE.get_overlap_types():
    file.write(f"{len(overlapped_gc_contents[ov_type])}\t{overlapped_gc_mean[ov_type]}\t"
               f"{statistics.mean(runs_control_gc_means[ov_type])}\t{statistics.stdev(runs_control_gc_means[ov_type])}\t"
               f"{under_h0[ov_type]}\t.<{under_h0[ov_type] + 1}\n")

file.close()


def build_histogram(data, mean_x, legend_title, annot_title, l, r, ax, ay, path):
    fig = go.Figure(data=[
        go.Histogram(x=data, name=legend_title, histnorm='probability', nbinsx=10,
                     marker=dict(color='lightblue'))])
    fig.update_layout(legend=dict(font=dict(size=12), yanchor="top", y=0.99, xanchor="left", x=0.01),
                      showlegend=True, margin=go.layout.Margin(l=0, r=0, b=13, t=5), height=250, width=450)
    fig.update_yaxes(showticklabels=False, )
    fig.update_xaxes(tickfont=dict(size=10), range=[l, r])
    fig.add_vline(x=mean_x, line_color='red', line_width=1)

    fig.add_annotation(
        x=mean_x, y=0.18, xref="x", yref="y", text=annot_title, showarrow=True,
        font=dict(family="Times New Roman", size=14, color="black"),
        align="center", arrowhead=2, arrowsize=1, arrowwidth=2, arrowcolor="black", ax=ax,
        ay=ay, bordercolor="black", borderwidth=2, borderpad=4, opacity=0.8
    )
    width = 700
    scale = (190 / 25.4) / (width / 300)

    fig.write_image(path, scale=scale)


build_histogram(runs_control_gc_means[OVERLAP_TYPE.MULTI], overlapped_gc_mean[OVERLAP_TYPE.MULTI],
                "Control transcripts mean GC contents", "Overlapped genes<br>mean GC content",
                0.48, 0.58, -100, +30, "generated_data/GC (overlapped vs control).jpg")
