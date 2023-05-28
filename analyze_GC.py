import numpy as np
from Bio.Seq import Seq
from GenomePackage.Enums import ANNOTATION_LOAD, OVERLAP_TYPE
from GenomePackage.GenomeData import GenomeData
from GenomePackage.OverlapsFinder import GetCDSOverlapsFromGenome
from Paths import *

MIN_LENGTH_OF_CDS_RECORD = 60
RANDOM_RUNS_COUNT = 100

request_species = SPECIES.Homo_sapiens
genome = GenomeData(SPECIES.Homo_sapiens,
                    ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS, with_sequence=True)

records, ATI_overlaps = GetCDSOverlapsFromGenome(genome, exclude_zero_expressed_transrcipts=False,
                                                 exclude_incomplete_utr_transcripts=True)
records.sort(key=lambda x: -x.overlapped_length)

target_group_CDSs = []
random_group_CDSs = [[] for _ in range(RANDOM_RUNS_COUNT)]
semi_random_group = [[] for _ in range(RANDOM_RUNS_COUNT)]


def GC_for_seq(cds_seq: str):
    gc = 0
    for ch in cds_seq.upper():
        if ch == 'C' or ch == 'G': gc += 1
    return gc / len(cds_seq)


def find_random_CDSs_for_target(target_cds):
    for i in range(RANDOM_RUNS_COUNT):
        random_cds = genome.get_random_sub_cds(len(target_cds))
        random_group_CDSs[i].append(GC_for_seq(random_cds))


def does_cds_contains_intra_stop_codon(cds_seq):
    aa_seq = Seq(cds_seq).translate(to_stop=False)
    return aa_seq.count("*") > 1 or (aa_seq.count("*") == 1 and not aa_seq.endswith("*"))


def get_cds_other_ORFs(cds_seq):
    ORFs = []
    for in_strand_frame in range(1, 3):
        other_orf_cds = cds_seq[in_strand_frame:]
        other_orf_cds = other_orf_cds[:len(other_orf_cds) - (len(other_orf_cds) % 3)]
        assert len(other_orf_cds) == len(cds_seq) - 3
        ORFs.append(other_orf_cds)

    random_cds = str(Seq(cds_seq).reverse_complement())
    for in_strand_frame in range(0, 3):
        other_orf_cds = random_cds[in_strand_frame:]
        other_orf_cds = other_orf_cds[:len(other_orf_cds) - (len(other_orf_cds) % 3)]
        assert len(other_orf_cds) == len(cds_seq) or len(other_orf_cds) == len(cds_seq) - 3
        ORFs.append(other_orf_cds)
    return ORFs


def find_semi_random_CDSs_for_target(target_cds, ov_type: OVERLAP_TYPE):
    for i in range(RANDOM_RUNS_COUNT):
        count = 0
        while True:
            count += 1
            random_cds = genome.get_random_sub_cds(len(target_cds), phased=True)
            if does_cds_contains_intra_stop_codon(random_cds): continue

            ORFs = get_cds_other_ORFs(random_cds)
            flag = False
            for orf in ORFs:
                if not does_cds_contains_intra_stop_codon(orf):
                    flag = True
            if flag: break
        print(count)
        semi_random_group[i].append(GC_for_seq(random_cds))


file = open("./Tasks/GC/overlapping_GCs.txt", "w")
file.write(f"Gene1\tGene2\tGC content\n")

for record in records:
    l, r = record.get_max_overlapped_segment()
    if r - l + 1 < MIN_LENGTH_OF_CDS_RECORD: continue

    cds = genome.get_transcript_sub_cds(record.transcript1, l, r)
    file.write(f"{record.sym1}\t{record.sym2}\t{GC_for_seq(cds)}\n")

    find_random_CDSs_for_target(cds)
    find_semi_random_CDSs_for_target(cds, record.transcripts_overlap_type)
file.close()

file_control = open("./Tasks/GC/control_groups_mean_GCs.txt", "w")
file_control.write(f"Random Groups\tSemi-random Groups\n")
for i in range(RANDOM_RUNS_COUNT):
    for j in range(len(random_group_CDSs[i])):
        file_control.write(f"{random_group_CDSs[i][j]}\t{semi_random_group[i][j]}\n")
file_control.close()
