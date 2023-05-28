import subprocess
import numpy as np
from GenomePackage.Enums import ANNOTATION_LOAD, OVERLAP_TYPE
from GenomePackage.GenomeData import GenomeData
from GenomePackage.Helper import GetTranscriptEnsemblId
from GenomePackage.OverlapsFinder import GetCDSOverlapsFromGenome, CDSRecord
from Paths import *

MIN_LENGTH_OF_CDS_RECORD = 60
RANDOM_SAMPLES_PER_RECORD = 100
RANDOM_SAMPLES_PER_NR1D1_THRA = 1000
request_species = SPECIES.Homo_sapiens
genome = GenomeData(SPECIES.Homo_sapiens,
                    ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS, with_sequence=True
                    )

records, ATI_overlaps = GetCDSOverlapsFromGenome(genome, exclude_incomplete_utr_transcripts=True,
                                                 exclude_zero_expressed_transrcipts=False)
records.sort(key=lambda x: -x.overlapped_length)

# region Blast query file generator
CDS_entries_file = open(f"./Tasks/Conservation/BlastData/Queries/CDS_entries_{request_species.tax_id()}.fa", "w")

conserved_species_by_rec_ids = {}
query_rec_ids = []


def add_query_entry(seq_id, sequence):
    query_rec_ids.append(seq_id)
    CDS_entries_file.write(f">{seq_id}\n")
    CDS_entries_file.write(f"{sequence}\n")


def deep_query_generator_for_NR1D1_THRA(transcript1, transcript2, cds):
    thra_last_exon = ''
    for frag in genome.transcript_fragments[transcript1.id]:
        if frag.featuretype == 'exon' and frag.attributes['rank'][0] == '8':
            thra_last_exon = genome.get_transcript_sub_cds(transcript1, frag.start, frag.end, phased=False)
    thra_last_exon = thra_last_exon[-len(cds):]

    nr1d1_last_exon = ''
    for frag in genome.transcript_fragments[record.transcript2.id]:
        if frag.featuretype == 'exon' and frag.attributes['rank'][0] == '7':
            nr1d1_last_exon = genome.get_transcript_sub_cds(record.transcript2, frag.start, frag.end, phased=False)
    nr1d1_last_exon = nr1d1_last_exon[-len(cds):]

    add_query_entry(f"NR1D1THRA", cds)
    add_query_entry(f"THRAexon9", thra_last_exon)
    add_query_entry(f"NR1D1exon7", nr1d1_last_exon)

    for i in range(RANDOM_SAMPLES_PER_NR1D1_THRA):
        add_query_entry(f"NR1D1THRA_RANDOM{i}", genome.get_random_sub_cds(len(cds)))
        add_query_entry(f"NR1D1exon7_RANDOM{i}", genome.get_random_sub_cds(len(nr1d1_last_exon)))
        add_query_entry(f"THRAexon9_RANDOM{i}", genome.get_random_sub_cds(len(thra_last_exon)))


for record in records:
    l, r = record.get_max_overlapped_segment()
    if r - l + 1 < MIN_LENGTH_OF_CDS_RECORD: continue
    cds = genome.get_transcript_sub_cds(record.transcript1, l, r, phased=False)
    add_query_entry(
        f"CDS_{record.transcript1.id.replace('transcript:', '')}_{record.transcript2.id.replace('transcript:', '')}",
        cds)

    if record.sym1 == "THRA":
        deep_query_generator_for_NR1D1_THRA(record.transcript1, record.transcript2, cds)

    for i in range(RANDOM_SAMPLES_PER_RECORD):
        random_cds_plus = genome.get_random_sub_cds(len(cds))
        add_query_entry(
            f"RANDOM{i}_{record.transcript1.id.replace('transcript:', '')}_{record.transcript2.id.replace('transcript:', '')}",
            random_cds_plus)

CDS_entries_file.close()
print(f"Blast query generated with {len(query_rec_ids)} entries for {request_species.short_name()}")

# endregion

# region Blast run and update records count
for rec_id in query_rec_ids:
    conserved_species_by_rec_ids[rec_id] = []

for target_species in ANNOTATIONS_PATH.keys():
    if target_species == request_species: continue
    query_path = f"./Tasks/Conservation/BlastData/Queries/CDS_entries_{request_species.tax_id()}.fa"
    db_path = f"./Tasks/Conservation/BlastData/db/{target_species.tax_id()}_db"
    query_result_path = f"./Tasks/Conservation/BlastData/Results/{request_species.name}_records_hits_in_{target_species.name}.txt"
    cmd_query = f"blastn -db {db_path} -query {query_path} -out {query_result_path} " \
                f"-num_threads 8 -task dc-megablast -outfmt 7 -evalue 1e-10 -dust no -max_target_seqs 5"

    print(f"CDS records from {request_species.short_name()} are being blasted on {target_species.short_name()} dna...")
    subprocess.run(cmd_query, shell=True)

    query_result_file = open(query_result_path, "r")
    lines = query_result_file.readlines()

    index = 0
    cur_rec_id = None

    while index < len(lines):
        if lines[index].startswith("# Query:"):
            cur_rec_id = lines[index].replace('\n', '').split(': ')[1]
            index += 1
        elif lines[index].__contains__("hits found"):
            index += 1
            flag = False
            max_identity = 0
            while index < len(lines) and len(lines[index]) > 1 and not lines[index].startswith("#"):
                blast_data = lines[index].split('\t')
                max_identity = max(max_identity, float(blast_data[2]))
                flag = True
                index += 1
            if flag:
                if not conserved_species_by_rec_ids.__contains__(cur_rec_id):
                    conserved_species_by_rec_ids[cur_rec_id] = []
                conserved_species_by_rec_ids[cur_rec_id].append((target_species, max_identity))
        else:
            index += 1

# endregion

file = open("./Tasks/Conservation/conservation_overlapping.txt", "w")

for rec_id in query_rec_ids:
    file.write(f"{rec_id}")
    for s, identity in conserved_species_by_rec_ids[rec_id]:
        file.write(f"\t{s.short_name()}\t{identity}")
    file.write(f"\n")
file.close()
file = open("./Tasks/Conservation/supplementary_conservation.txt", "w")

for rec_id in query_rec_ids:
    if not rec_id.__contains__("CDS"): continue
    file.write(f"{rec_id}")
    for s, identity in conserved_species_by_rec_ids[rec_id]:
        file.write(f"\t{s.short_name()}\t{identity}")
    file.write(f"\n")
file.close()
