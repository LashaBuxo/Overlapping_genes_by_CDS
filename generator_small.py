import random
from CDS_Overlap_Record import *

genome = GenomeWorker(SPECIES.Homo_sapiens, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS, SEQUENCE_LOAD.LOAD)
records,_ = get_CDS_records_from_genome(genome)

total = 0
cnt = 0
OGs={}
for chr_index in range(1, genome.chromosomes_count() + 1):
    genes_cnt = genome.genes_count_on_chr(chr_index)
    total += genes_cnt
    for i in range(0, genes_cnt):
        gene1 = genome.gene_by_indexes(chr_index, i)
        for j in range(i + 1, genes_cnt):
            gene2 = genome.gene_by_indexes(chr_index, j)
            overlap_type = genome.get_features_overlap_type(gene1, gene2)
            # if genes don't overlap, there is no point in searching overlaps by fragments (CDS in this case)
            if overlap_type is OVERLAP_TYPE.NONE: continue
            OGs[gene1.id]=True
            OGs[gene2.id]=True

print(f"{len(OGs) * 100 / total}% ({len(OGs)}) genes overlapped by boundaries")
