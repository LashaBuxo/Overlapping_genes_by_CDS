from CDS_Overlap_Record import *

human_genome = GenomeWorker(SPECIES.Homo_sapiens, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS,
                            SEQUENCE_LOAD.LOAD)
human_records,_ = get_CDS_records_from_genome(human_genome,with_clustered_graph=False)

mouse_genome = GenomeWorker(SPECIES.Mus_musculus, ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_FRAGMENTS,
                            SEQUENCE_LOAD.LOAD)
mouse_records,_ = get_CDS_records_from_genome(mouse_genome,with_clustered_graph=False)

orth_data = {}

file = open("used_data/orthologous_human_mouse.txt", "r")
lines = file.readlines()
for index in range(1, len(lines)):
    arr = lines[index].replace('\n', '').split(',')
    human_gene_id = arr[0]
    mouse_gene_id = arr[2]

    if human_genome.feature_by_id(f"gene:{human_gene_id}") is None or mouse_genome.feature_by_id(
            f"gene:{mouse_gene_id}") is None:
        continue

    if not orth_data.__contains__(human_gene_id):
        orth_data[human_gene_id] = []
    if not orth_data.__contains__(mouse_gene_id):
        orth_data[mouse_gene_id] = []
    orth_data[human_gene_id].append(mouse_gene_id)
    orth_data[mouse_gene_id].append(human_gene_id)


def get_match_score(s1: str, s2: str):
    s1 = s1.replace('.', '').replace('|', '')
    s2 = s2.replace('.', '').replace('|', '')
    dp = []
    for i in range(0, len(s1) + 1):
        arr = [0] * (len(s2) + 1)
        dp.append(arr)
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            dp[i][j] = max(max(dp[i - 1][j], dp[i][j - 1]),
                           dp[i - 1][j - 1] + 1 if s1[i - 1] == s2[j - 1] else dp[i - 1][j - 1])

    return dp[len(s1)][len(s2)]


data = []
for human_record in human_records:
    for mouse_record in mouse_records:
        human_gene_id1 = human_genome.get_transcript_parent(human_record.transcript1.id).id.replace("gene:", "")
        human_gene_id2 = human_genome.get_transcript_parent(human_record.transcript2.id).id.replace("gene:", "")

        mouse_gene_id1 = mouse_genome.get_transcript_parent(mouse_record.transcript1.id).id.replace("gene:", "")
        mouse_gene_id2 = mouse_genome.get_transcript_parent(mouse_record.transcript2.id).id.replace("gene:", "")

        if not orth_data.__contains__(human_gene_id1) or not orth_data.__contains__(human_gene_id2):
            continue

        if (orth_data[human_gene_id1].__contains__(mouse_gene_id1) and orth_data[human_gene_id2].__contains__(
                mouse_gene_id2)) or (
                orth_data[human_gene_id1].__contains__(mouse_gene_id2) and orth_data[human_gene_id2].__contains__(
            mouse_gene_id1)):
            match_score = get_match_score(human_record.overlapped_CDS, mouse_record.overlapped_CDS)
            match_string = f"{match_score} ({match_score * 100 // human_record.overlapped_length}%)"
            data.append((human_record.gene1_sym, human_record.gene2_sym, human_record.transcripts_overlap_type,
                         human_record.exonic_overlap_type, human_record.overlapped_length,
                         mouse_record.overlapped_length,
                         match_string))
data.sort(key=lambda x: int(x[6].split(' ')[0]), reverse=True)
file = open("generated_data/human_mouse_conservation.txt", "w")
for sym1, sym2, type1, type2, ov_length1, ov_length2, match_string in data:
    file.write(f"{sym1}\t{sym2}\t{type1}\t{type2}\t{ov_length1}\t{ov_length2}\t{match_string}\n")
file.close()
