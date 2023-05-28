from GenomePackage.Enums import ANNOTATION_LOAD, OVERLAP_TYPE
from GenomePackage.GenomeData import GenomeData
from GenomePackage.Helper import GetFeaturesOverlapType
from GenomePackage.OverlapsFinder import GetCDSOverlapsFromGenome
from Paths import *

request_species = SPECIES.Homo_sapiens
genome = GenomeData(SPECIES.Homo_sapiens, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS)

records, ATI_overlaps = GetCDSOverlapsFromGenome(genome,  exclude_zero_expressed_transrcipts=False,
                                                 exclude_incomplete_utr_transcripts=True)

# region output human CDS/CDS records
records.sort(key=lambda x: -x.overlapped_length)

file = open(f"Tasks/Records/{request_species.short_name()}_OGs by CDS records.txt", "w")
for record in records:
    file.write(f"{record.get_tabbed_data_string()}\n")
file.close()

file = open(f"Tasks/Records/{request_species.short_name()}_OGs by CDS (ATI genes).txt", "w")
file.write(f"{len(ATI_overlaps)} genes are overlapped by CDS and are ATI:\n")
for gene_sym in ATI_overlaps.keys():
    file.write(f"{gene_sym}\n")
file.close()
# endregion

file = open("./Tasks/Records/stats.txt", "w")

OGs_set = {}
OGs_by_types = {}
for chrom in genome.chromosomes:
    genes = genome.genes_on_chr[chrom]
    for i in range(0, len(genes)):
        for j in range(i + 1, len(genes)):
            ov_type = GetFeaturesOverlapType(genes[i], genes[j])
            if ov_type != OVERLAP_TYPE.NONE:
                if not OGs_by_types.__contains__(ov_type):
                    OGs_by_types[ov_type] = 0
                OGs_by_types[ov_type] += 1
                OGs_set[genes[i].id] = True
                OGs_set[genes[j].id] = True

OGs_by_CDS_set = {}
OGs_by_CDS_by_types = {}
for record in records:
    OGs_by_CDS_set[record.gene1.id] = True
    OGs_by_CDS_set[record.gene2.id] = True
    if not OGs_by_CDS_by_types.__contains__(record.transcripts_overlap_type):
        OGs_by_CDS_by_types[record.transcripts_overlap_type] = 0
    OGs_by_CDS_by_types[record.transcripts_overlap_type] += 1
file.write(f"Total OGs: {len(list(OGs_set.keys()))}\n")
file.write(f"OGs by types: {OGs_by_types}\n")
file.write(f"Total OGs by CDS: {len(list(OGs_by_CDS_set.keys()))}\n")
file.write(f"OGs by CDS by types: {OGs_by_CDS_by_types}\n")
# FAHD1 ENSG00000180185
# MEIOB ENSG00000162039

# ZNF575 ENSG00000176472
# ETHE1 ENSG00000105755
file.write(
    f"FAHD1 (overall TPM - {genome.get_gene_expression(genome.get_feature_by_id('gene:ENSG00000180185'))}) transcripts:\n")
for transcript in genome.gene_transcripts["gene:ENSG00000180185"]:
    file.write(
        f"{transcript.id} - TPM: {genome.get_transcript_expression(transcript)} Rank: {genome.get_transcript_expression_rank(transcript)}\n")

file.write(
    f"MEIOB (overall TPM - {genome.get_gene_expression(genome.get_feature_by_id('gene:ENSG00000162039'))}) transcripts:\n")
for transcript in genome.gene_transcripts["gene:ENSG00000162039"]:
    file.write(
        f"{transcript.id} - TPM: {genome.get_transcript_expression(transcript)} Rank: {genome.get_transcript_expression_rank(transcript)}\n")

file.write(
    f"ZNF575 (overall TPM - {genome.get_gene_expression(genome.get_feature_by_id('gene:ENSG00000176472'))}) transcripts:\n")
for transcript in genome.gene_transcripts["gene:ENSG00000176472"]:
    file.write(
        f"{transcript.id} - TPM: {genome.get_transcript_expression(transcript)} Rank: {genome.get_transcript_expression_rank(transcript)}\n")

file.write(
    f"ETHE1 (overall TPM - {genome.get_gene_expression(genome.get_feature_by_id('gene:ENSG00000105755'))}) transcripts:\n")
for transcript in genome.gene_transcripts["gene:ENSG00000105755"]:
    file.write(
        f"{transcript.id} - TPM: {genome.get_transcript_expression(transcript)} Rank: {genome.get_transcript_expression_rank(transcript)}\n")
