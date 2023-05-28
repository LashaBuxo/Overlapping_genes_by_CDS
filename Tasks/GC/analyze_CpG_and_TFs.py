import numpy as np

from GenomePackage.Enums import ANNOTATION_LOAD, OVERLAP_TYPE
from GenomePackage.GenomeData import GenomeData
from GenomePackage.OverlapsFinder import GetCDSOverlapsFromGenome
from Paths import *

MIN_LENGTH_OF_CDS_RECORD = 60
TF_OFFSET = 200

request_species = SPECIES.Homo_sapiens
genome = GenomeData(SPECIES.Homo_sapiens, ANNOTATIONS_PATH[request_species],
                    ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS,
                    sequence_path=SEQUENCES_PATH[request_species],
                    GTEx_multi_tissue_multi_sample_file_path="../../CommonRequiredData/Human Expressions/short_required_expressions.txt"
                    )

records, ATI_overlaps = GetCDSOverlapsFromGenome(genome)
records.sort(key=lambda x: -x.overlapped_length)

# region CpG islands Finder and comparing to Random Groups
RANDOM_RUNS_COUNT = 20
total_reviewed = 0
total_reviewed_convergent = 0
total_reviewed_divergent = 0

CpG_islands_OVERLAPPED = 0
CpG_islands_OVERLAPPED_convergent = 0
CpG_islands_OVERLAPPED_divergent = 0
for record in records:
    l, r = record.get_max_overlapped_segment()
    if r - l + 1 < MIN_LENGTH_OF_CDS_RECORD: continue
    result = genome.find_CGIs_overlapped_on_interval(record.gene1.chrom, l, r)
    CpG_islands_OVERLAPPED += 1 if result != '-' else 0
    total_reviewed += 1

    if record.transcripts_overlap_type == OVERLAP_TYPE.CONVERGENT:
        total_reviewed_convergent += 1
        CpG_islands_OVERLAPPED_convergent += 1 if result != '-' else 0
    if record.transcripts_overlap_type == OVERLAP_TYPE.DIVERGENT:
        total_reviewed_divergent += 1
        CpG_islands_OVERLAPPED_divergent += 1 if result != '-' else 0

CpG_islands_RANDOMS = []
for i in range(RANDOM_RUNS_COUNT):
    print(i)
    random_CGIs = 0
    for record in records:
        l, r = record.get_max_overlapped_segment()
        if r - l + 1 < MIN_LENGTH_OF_CDS_RECORD: continue
        _, chr, cds_l, cds_r = genome.get_random_sub_cds_location(sub_cds_length=r - l + 1)
        result = genome.find_CGIs_overlapped_on_interval(chr, cds_l, cds_r)
        random_CGIs += 1 if result != '-' else 0
    CpG_islands_RANDOMS.append(random_CGIs / total_reviewed)

file = open("./OutputData/CpG_data.txt", "w")
file.write(f"{total_reviewed}\t{total_reviewed_convergent}\t{total_reviewed_divergent}\n")
file.write(f"{CpG_islands_OVERLAPPED}\t{CpG_islands_OVERLAPPED_convergent}\t{CpG_islands_OVERLAPPED_divergent}\n")

file.write(f"Random Freq.: mean - {np.mean(CpG_islands_RANDOMS)}\n")
file.write(f"Random Freq.: std - {np.std(CpG_islands_RANDOMS)}\n")
file.close()
# endregion
# region ChIP-Atlas Query generator
file = open("./OutputData/TFs_enrich_query_convergent_promotors.txt", "w")
for record in records:
    l, r = record.get_max_overlapped_segment()
    if r - l + 1 < MIN_LENGTH_OF_CDS_RECORD: continue
    result = genome.find_CGIs_overlapped_on_interval(record.gene1.chrom, l, r)
    if result == '-': continue
    if record.transcripts_overlap_type == OVERLAP_TYPE.CONVERGENT:
        file.write(f"chr{record.gene1.chrom}\t{l - TF_OFFSET}\t{r + TF_OFFSET}\n")
file.close()

file = open("./OutputData/TFs_enrich_query_random_promotors.txt", "w")
for record in records:
    l, r = record.get_max_overlapped_segment()
    if r - l + 1 < MIN_LENGTH_OF_CDS_RECORD: continue
    result = genome.find_CGIs_overlapped_on_interval(record.gene1.chrom, l, r)
    if result == '-': continue
    if record.transcripts_overlap_type == OVERLAP_TYPE.CONVERGENT:
        _, chr, cds_l, cds_r = genome.get_random_sub_cds_location(sub_cds_length=r - l + 1)
        file.write(f"chr{chr}\t{cds_l - TF_OFFSET}\t{cds_r + TF_OFFSET}\n")
file.close()
# endregion
