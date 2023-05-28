import os

from gffutils import Feature
from Bio.Seq import Seq
from GenomePackage.Enums import ANNOTATION_LOAD
from GenomePackage.GenomeData import GenomeData
from GenomePackage.Helper import GetTranscriptEnsemblId
from GenomePackage.OverlapsFinder import GetCDSOverlapsFromGenome
from Paths import *

MIN_LENGTH_OF_CDS_RECORD = 60

request_species = SPECIES.Homo_sapiens
genome = GenomeData(SPECIES.Homo_sapiens, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS, with_sequence=True)
# flDPnn
# (putative function- and linker based Disorder Prediction using deep neural network)

records, ATI_overlaps = GetCDSOverlapsFromGenome(genome, exclude_incomplete_utr_transcripts=True,
                                                 exclude_zero_expressed_transrcipts=False)
records.sort(key=lambda x: -x.overlapped_length)

directory = './Tasks/Disorder/FlDPnnData/Results/'
result_files = []

for filename in os.listdir(directory):
    path = os.path.join(directory, filename)
    if os.path.isfile(path):
        result_files.append(path)

# region query generation
if len(result_files) == 0:
    file = open("Tasks/Disorder/FlDPnnData/Queries/query.txt", "w")

    for record in records:
        l, r = record.get_max_overlapped_segment()
        if r - l + 1 < MIN_LENGTH_OF_CDS_RECORD: continue
        cds1 = genome.get_transcript_cds(record.transcript1)
        cds2 = genome.get_transcript_cds(record.transcript2)

        aa_seq1 = Seq(cds1).translate(to_stop=False)
        aa_seq2 = Seq(cds2).translate(to_stop=False)

        if not aa_seq1.endswith("*"):
            print(f"transcript cds for {record.sym1}: {record.transcript1.id} not ends with stop codon!")
        if not aa_seq2.endswith("*"):
            print(f"transcript cds for {record.sym2}: {record.transcript2.id} not ends with stop codon!")
        aa_seq1 = aa_seq1[:-1]
        aa_seq2 = aa_seq2[:-1]

        if "*" in aa_seq1:
            print(f"transcript cds for {record.sym1}: {record.transcript1.id} contains intra stop codon!")
        if "*" in aa_seq2:
            print(f"transcript cds for {record.sym2}: {record.transcript2.id} contains intra stop codon!")

        file.write(f">{record.sym1}_{GetTranscriptEnsemblId(record.transcript1)}\n")
        file.write(f"{aa_seq1}\n")

        file.write(f">{record.sym2}_{GetTranscriptEnsemblId(record.transcript2)}\n")
        file.write(f"{aa_seq2}\n")
    file.close()


# endregion

def get_disorder_rate_from_transcript(transcript: Feature, sub_aa_seq):
    cds = genome.get_transcript_cds(transcript)
    transcript_aa_seq = Seq(cds).translate(to_stop=False)[:-1]
    assert not transcript_aa_seq.__contains__("*")
    data = disorder_data_by_transcript_id[transcript.id]

    full_seq = ""
    for aa, state in data:
        full_seq += aa
    assert full_seq == transcript_aa_seq

    assert transcript_aa_seq.index(sub_aa_seq) != -1
    if transcript_aa_seq.count(sub_aa_seq) > 1:
        print(transcript_aa_seq)
        print(sub_aa_seq)

    start_index = transcript_aa_seq.index(sub_aa_seq)

    disorder_count = 0
    for i in range(start_index, start_index + len(sub_aa_seq)):
        disorder_count += 1 if float(data[i][1]) > 0.5 else 0
    return disorder_count, len(sub_aa_seq)


# region read results
if len(result_files) != 0:

    disorder_data_by_transcript_id = {}
    for result_file in result_files:
        file = open(result_file, "r")
        lines = file.readlines()
        index = 0
        while index < len(lines):
            if lines[index].startswith(">"):
                transcript_id = lines[index].replace('\n', '').split('_')[1]
                transcript_id = f"transcript:{transcript_id}"
                disorder_data_by_transcript_id[transcript_id] = []
                arr_aa_seq = lines[index + 1].replace('\n', '').split(',')
                arr_states = lines[index + 2].replace('\n', '').split(',')
                assert len(arr_aa_seq) == len(arr_states)
                for i in range(len(arr_aa_seq)):
                    disorder_data_by_transcript_id[transcript_id].append((arr_aa_seq[i], arr_states[i]))
                index += 10
            else:
                index += 1
        file.close()

    file_supplementary=open("Tasks/Disorder/supplementary_disorder.txt","w")
    file_supplementary.write(f"Gene1\tGene2\tOverlapped amino acid seq1\tOverlapped amino acid seq2\n")

    file_disorder_data = open("Tasks/Disorder/disorder_data.txt", "w")
    file_disorder_data.write(f"Gene\tTranscript\tEntire Length\tEntire disorder value\t"
                             f"Non-overlapped regions length\tNon-overlapped regions disorder value\t"
                             f"Overlapping regions length\tOverlapping regions disorder value\n")


    def print_transcript_data(transcript: Feature, geneSym, l, r):
        cds = genome.get_transcript_cds(transcript)
        aa_seq = Seq(cds).translate(to_stop=False)[:-1]
        assert not aa_seq.__contains__("*")
        sub_cds = genome.get_transcript_sub_cds(transcript, l, r, phased=True)
        sub_aa_seq = Seq(sub_cds).translate(to_stop=False)
        if sub_aa_seq.count("*") > 0:
            sub_aa_seq = sub_aa_seq[:-1]
            assert sub_aa_seq.count("*") == 0

        transcript_rate = get_disorder_rate_from_transcript(transcript, aa_seq)
        parts_rate = get_disorder_rate_from_transcript(transcript, sub_aa_seq)

        file_disorder_data.write(
            f"{geneSym}\t{GetTranscriptEnsemblId(transcript)}\t"
            f"{transcript_rate[1]}\t{transcript_rate[0]}\t"
            f"{transcript_rate[1] - parts_rate[1]}\t{transcript_rate[0] - parts_rate[0]}\t"
            f"{parts_rate[1]}\t{parts_rate[0]}\n")

        return sub_aa_seq
    for record in records:
        l, r = record.get_max_overlapped_segment()
        if r - l + 1 < MIN_LENGTH_OF_CDS_RECORD: continue
        seq1=print_transcript_data(record.transcript1, record.sym1, l, r)
        seq2=print_transcript_data(record.transcript2, record.sym2, l, r)
        file_supplementary.write(f"{record.sym1}\t{record.sym2}\t{seq1}\t{seq2}\n")
    file_supplementary.close()
    file_disorder_data.close()
# endregion
