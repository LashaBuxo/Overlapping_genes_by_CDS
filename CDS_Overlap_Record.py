from worker_genome import *
from graph import *


def get_CDS_records_from_genome(genome: GenomeWorker, with_clustered_graph=True):
    graph_clusters = AnalyzerGraph()
    ATI_overlaps = {}
    records = []
    for chr_id in range(1, genome.chromosomes_count() + 1):
        genes_cnt = genome.genes_count_on_chr(chr_id)
        for i in range(0, genes_cnt):
            for j in range(i + 1, genes_cnt):
                gene_A, gene_B = genome.gene_by_indexes(chr_id, i), genome.gene_by_indexes(chr_id, j)

                overlap_type = genome.get_features_overlap_type(gene_A, gene_B)

                # if genes don't overlap, there is no point in searching overlaps by fragments (CDS in this case)
                if overlap_type is OVERLAP_TYPE.NONE: continue

                transcripts_A = genome.get_transcripts_from_gene(gene_A.id)
                transcripts_B = genome.get_transcripts_from_gene(gene_B.id)

                best_record = None
                for transcript1 in transcripts_A:
                    for transcript2 in transcripts_B:
                        transcript_a, transcript_b = (transcript1, transcript2) if transcript1.strand == '+' else (
                            transcript2, transcript1)

                        non_ati_segments, ov_l, ATI = genome.get_overlaps_between_transcripts(transcript_a.id,
                                                                                              transcript_b.id)

                        if ov_l > 0 and not ATI:
                            record = CDSRecord(transcript_a, transcript_b, genome, non_ati_segments)
                            best_record = record.get_better_record(best_record)
                        if ATI:
                            ATI_overlaps[genome.get_gene_symbol(gene_A)] = True
                            ATI_overlaps[genome.get_gene_symbol(gene_B)] = True
                if best_record is not None:
                    edge = AnalyzerGraph.GraphEdge(gene_A.id, gene_B.id)
                    graph_clusters.add_edge(edge)
                    records.append(best_record)
    if with_clustered_graph == False:
        return records, ATI_overlaps
    else:
        return records, ATI_overlaps, graph_clusters


class CDSRecord:
    def __init__(self, transcript1: Feature, transcript2: Feature, genome: GenomeWorker, overlap_segments):
        self.transcript1 = transcript1
        self.transcript2 = transcript2
        self.genome = genome
        self.overlap_segments = overlap_segments

        self.gene1_sym = genome.get_gene_symbol(genome.get_transcript_parent(transcript1.id))
        self.gene2_sym = genome.get_gene_symbol(genome.get_transcript_parent(transcript2.id))

        self.transcripts_overlap_type = genome.get_features_overlap_type(transcript1, transcript2).short_name()
        self.exonic_overlap_type = self.get_exonic_overlap_type().short_name()

        self.overlapped_length = self.get_overlapped_length()
        self.segments = len(overlap_segments)

        self.overlapped_CDS = self.get_overlapped_CDS()

        self.overlapped_GC_content = self.genome.sequence_GC(self.overlapped_CDS)

        self.conservation1 = self.genome.get_transcript_conservation_score(self.transcript1.id)
        self.conservation2 = self.genome.get_transcript_conservation_score(self.transcript2.id)

        self.overlap_start_index_on_CDS1 = self.get_transcript_CDS_overlap_start_index(self.transcript1.id)
        self.overlap_start_index_on_CDS2 = self.get_transcript_CDS_overlap_start_index(self.transcript2.id)

    def get_tabbed_data_string(self):
        return f"{self.gene1_sym}\t{self.gene2_sym}\t" \
               f"{self.transcript1.id.replace('transcript:', '')}\t" \
               f"{self.transcript2.id.replace('transcript:', '')}\t" \
               f"{self.transcripts_overlap_type}\t{self.overlapped_length}\t{self.segments}\t" \
               f"{self.exonic_overlap_type}\t{self.overlapped_GC_content}\t" \
               f"{self.conservation1}\t{self.conservation2}"

    def get_transcript_CDS_overlap_start_index(self, transcript_id):
        transcript = self.genome.feature_by_id(transcript_id)
        frags = self.genome.get_fragments_from_transcript(transcript_id)
        cds_intervals = []
        for frag in frags:
            if frag.featuretype == 'CDS':
                cds_intervals.append((frag.start, frag.end))
        cds_intervals.sort(key=lambda x: x[0])
        index = 0
        for i in range(len(cds_intervals)):
            l, r = cds_intervals[i] if transcript.strand == '+' else cds_intervals[len(cds_intervals) - 1 - i]
            flag = False
            for ov_l, ov_r in self.overlap_segments:
                if self.genome.get_segments_overlap_type((l, r, '+'), (ov_l, ov_r, '+')) != OVERLAP_TYPE.NONE:
                    flag = True
                    if transcript.strand == '+':
                        index += 1 if ov_l < l else ov_l - l + 1
                    else:
                        index += 1 if ov_r > r else r - ov_r + 1
            if not flag:
                index += r - l + 1
            else:
                break
        return index

    def get_better_record(self, alt_record):
        if alt_record is None: return self
        self_score = min(self.conservation1, self.conservation2)
        alt_score = min(alt_record.conservation1, alt_record.conservation2)
        return self if self_score > alt_score else alt_record

    def get_overlapped_CDS(self):
        strand = self.transcript1.strand if self.gene1_sym < self.gene2_sym else self.transcript2.strand
        chr_index = self.genome.chr_name2index(self.transcript1.chrom)
        cds = ""
        for l, r in self.overlap_segments:
            cds += self.genome.retrieve_segment_sequence(chr_index, l, r, strand)
        return cds

    def get_overlapped_length(self):
        length = 0
        for l, r in self.overlap_segments:
            length += r - l + 1
        return length

    def get_exonic_overlap_type(self):
        fragments1 = self.genome.get_fragments_from_transcript(self.transcript1.id)
        fragments2 = self.genome.get_fragments_from_transcript(self.transcript2.id)

        exonic_overlap_type = None
        for frag1 in fragments1:
            for frag2 in fragments2:
                if frag1.featuretype != 'exon' or frag2.featuretype != 'exon': continue
                ov_type = self.genome.get_features_overlap_type(frag1, frag2)
                if ov_type != OVERLAP_TYPE.NONE:
                    if exonic_overlap_type is None:
                        exonic_overlap_type = ov_type
                    elif exonic_overlap_type != ov_type:
                        return OVERLAP_TYPE.MULTI
        return exonic_overlap_type
