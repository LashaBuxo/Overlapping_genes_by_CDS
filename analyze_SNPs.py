# SNPs, codon frequency, GC content and conservation
from gffutils import Feature

from GenomePackage.Enums import ANNOTATION_LOAD
from GenomePackage.GenomeData import GenomeData
from GenomePackage.Helper import GetTranscriptEnsemblId
from GenomePackage.OverlapsFinder import GetCDSOverlapsFromGenome
from Paths import *

MIN_LENGTH_OF_CDS_RECORD = 60

request_species = SPECIES.Homo_sapiens
genome = GenomeData(SPECIES.Homo_sapiens, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS, with_sequence=True,
                    SNPs_path="./Tasks/SNPs/1000GenomesData/short_required_SNPs.txt")

records, ATI_overlaps = GetCDSOverlapsFromGenome(genome,
                                                 exclude_zero_expressed_transrcipts=False,
                                                 exclude_incomplete_utr_transcripts=True,
                                                 )
records.sort(key=lambda x: -x.overlapped_length)

file_SNPs_data = open("Tasks/SNPs/SNPs_data.txt", "w")
file_SNPs_data.write(f"Gene\tTranscript\t"
                     f"Entire Length\tEntire SNPs count\tSynonymous SNPs count\tpos1\tpos2\tpos3\t"
                     f"Non-overlapped regions length\tNon-overlapped regions SNPs count\tSynonymous SNPs count\tpos1\tpos2\tpos3\t"
                     f"Overlapping regions length\tOverlapping regions SNPs count\tSynonymous SNPs count\tpos1\tpos2\tpos3\n")


def print_transcript_data(transcript: Feature, geneSym, l, r):
    overlapped_length = r - l + 1
    overlapped_SNPs_syn, overlapped_SNPs_cnt, overlapped_SNPs_pos = genome.get_SNPs_count_on_transcript_sub_cds(
        transcript, l, r)
    entire_SNPs_syn, entire_SNPs_cnt, entire_length, entire_SNPs_pos = genome.get_transcript_SNPs_ratio(
        transcript)

    file_SNPs_data.write(f"{geneSym}\t{GetTranscriptEnsemblId(transcript)}\t")
    file_SNPs_data.write(
        f"{entire_length}\t"
        f"{entire_SNPs_cnt}\t"
        f"{entire_SNPs_syn}\t"
        f"{entire_SNPs_pos[0]}\t"
        f"{entire_SNPs_pos[1]}\t"
        f"{entire_SNPs_pos[2]}\t"
        f"{entire_length - overlapped_length}\t"
        f"{entire_SNPs_cnt - overlapped_SNPs_cnt}\t"
        f"{entire_SNPs_syn - overlapped_SNPs_syn}\t"
        f"{entire_SNPs_pos[0] - overlapped_SNPs_pos[0]}\t"
        f"{entire_SNPs_pos[1] - overlapped_SNPs_pos[1]}\t"
        f"{entire_SNPs_pos[2] - overlapped_SNPs_pos[2]}\t"
        f"{overlapped_length}\t"
        f"{overlapped_SNPs_cnt}\t"
        f"{overlapped_SNPs_syn}\t"
        f"{overlapped_SNPs_pos[0]}\t"
        f"{overlapped_SNPs_pos[1]}\t"
        f"{overlapped_SNPs_pos[2]}\n"
    )


for record in records:
    l, r = record.get_max_overlapped_segment()
    if r - l + 1 < MIN_LENGTH_OF_CDS_RECORD: continue
    print_transcript_data(record.transcript1, record.sym1, l, r)
    print_transcript_data(record.transcript2, record.sym2, l, r)

file_SNPs_data.close()

#
# val1 = "{:.2f}".format(SNPs_total_on_overlaps / total_on_overlaps * 100)
# val2 = "{:.2f}".format(SNPs_total_on_Non_overlaps / total_on_Non_overlaps * 100)
#
# print(f"SNPs on overlaps {SNPs_total_on_overlaps} = {val1}% (total: {total_on_overlaps})")
# print(f"SNPs on non-overlaps {SNPs_total_on_Non_overlaps} = {val2}% (total: {total_on_Non_overlaps})")
#
# print(
#     f"SNPs on overlaps by index: 1 - {SNPs_index_on_overlaps[0]}, 2 - {SNPs_index_on_overlaps[1]}, 3 - {SNPs_index_on_overlaps[2]}")
# print(
#     f"SNPs on non-overlaps by index: 1 - {SNPs_index_on_Non_overlaps[0]}, 2 - {SNPs_index_on_Non_overlaps[1]}, 3 - {SNPs_index_on_Non_overlaps[2]}")
#
# val1 = "{:.2f}".format(SNPs_syn_on_overlaps / SNPs_total_on_overlaps * 100)
# val2 = "{:.2f}".format(SNPs_syn_on_Non_overlaps / SNPs_total_on_Non_overlaps * 100)
#
# print(f"Synonymous SNPs on overlaps {SNPs_syn_on_overlaps} = {val1}% (total: {SNPs_total_on_overlaps})")
# print(f"Synonymous SNPs on non-overlaps {SNPs_syn_on_Non_overlaps} = {val2}% (total: {SNPs_total_on_Non_overlaps})")
#
# overlapped_CDSs = []
# for record in records:
#     l, r = record.get_max_overlapped_segment()
#     if r - l + 1 < MIN_LENGTH_OF_CDS_RECORD: continue
#     overlapped_CDSs.append(genome.get_transcript_sub_cds(record.transcript1, l, r))
