import scipy

from GenomePackage.GenomeData import *


def GetCDSOverlapsFromGenome(genome: GenomeData, exclude_zero_expressed_transrcipts=False,
                             exclude_incomplete_utr_transcripts=False):
    ATI_overlaps = {}
    records = []

    for chrom in genome.chromosomes:
        genes = genome.genes_on_chr[chrom]
        for i in range(0, len(genes)):
            for j in range(i + 1, len(genes)):
                overlap_type = GetFeaturesOverlapType(genes[i], genes[j])

                # if genes don't overlap, there is no point in searching overlaps by fragments (CDS in this case)
                if overlap_type is OVERLAP_TYPE.NONE: continue

                transcripts_i = genome.gene_transcripts[genes[i].id]
                transcripts_j = genome.gene_transcripts[genes[j].id]

                best_record = None
                for transcript_i in transcripts_i:
                    if transcript_i is None: continue;
                    if exclude_incomplete_utr_transcripts:
                        if not genome.is_transcript_cds_annotated(transcript_i): continue
                    if exclude_zero_expressed_transrcipts:
                        val = genome.get_transcript_expression(transcript_i)
                        if val == '-' or val == 0: continue
                    for transcript_j in transcripts_j:
                        if transcript_j is None: continue;
                        if exclude_incomplete_utr_transcripts:
                            if not genome.is_transcript_cds_annotated(transcript_j): continue
                        if exclude_zero_expressed_transrcipts:
                            val = genome.get_transcript_expression(transcript_j)
                            if val == '-' or val == 0: continue

                        if transcript_i.strand == '+':
                            transcript_a, transcript_b = (transcript_i, transcript_j)
                            gene_a, gene_b = genes[i], genes[j]
                        else:
                            transcript_a, transcript_b = (transcript_j, transcript_i)
                            gene_a, gene_b = genes[j], genes[i]

                        non_ati_segments, ov_l, ATI = GetOverlaps(genome.transcript_fragments[transcript_a.id],
                                                                  genome.transcript_fragments[transcript_b.id])

                        if ov_l > 0 and not ATI:
                            record = CDSRecord(transcript_a, gene_a, transcript_b, gene_b, genome,
                                               non_ati_segments)

                            best_record = record.get_better_record(best_record)

                        if ATI:
                            ATI_overlaps[(GetGeneEnsemblId(gene_a), GetGeneSymbol(gene_a))] = True
                            ATI_overlaps[(GetGeneEnsemblId(gene_b), GetGeneSymbol(gene_b))] = True
                if best_record is not None:
                    records.append(best_record)
    return records, ATI_overlaps


class CDSRecord:
    def __init__(self, transcript1: Feature, gene1: Feature, transcript2: Feature, gene2: Feature, genome: GenomeData,
                 overlap_segments):
        self.genome = genome

        self.gene1 = gene1
        self.gene2 = gene2

        self.transcript1 = transcript1
        self.transcript2 = transcript2

        self.overlap_segments = overlap_segments

        self.id1 = GetGeneEnsemblId(gene1)
        self.id2 = GetGeneEnsemblId(gene2)

        self.sym1 = GetGeneSymbol(gene1)
        self.sym2 = GetGeneSymbol(gene2)

        self.t_expression1 = genome.get_transcript_expression(transcript1)
        self.t_expression2 = genome.get_transcript_expression(transcript2)
        self.expression_rank1 = genome.get_transcript_expression_rank(transcript1)
        self.expression_rank2 = genome.get_transcript_expression_rank(transcript2)
        self.g_expression1 = genome.get_gene_expression(gene1)
        self.g_expression2 = genome.get_gene_expression(gene2)

        self.transcripts_overlap_type = GetFeaturesOverlapType(transcript1, transcript2)

        self.overlapped_length = self.get_overlapped_length()
        self.segments = f"{len(overlap_segments)} interval on chr{transcript1.chrom} - "
        for l, r in overlap_segments:
            self.segments += f"({l}:{r}) "

    def get_tabbed_data_string(self):
        return f"{self.id1}\t{self.sym1}\t{self.transcript1.id.replace('transcript:', '')}\t" \
               f"{self.g_expression1}\t{self.t_expression1}\t{self.expression_rank1}\t" \
               f"{self.id2}\t{self.sym2}\t{self.transcript2.id.replace('transcript:', '')}\t" \
               f"{self.g_expression2}\t{self.t_expression2}\t{self.expression_rank2}\t" \
               f"{self.transcripts_overlap_type.short_name()}\t" \
               f"{self.overlapped_length}\t{self.segments}\t"

    def get_score(self):
        score = 0
        score += self.expression_rank1 if self.expression_rank1 != '-' else 100
        score += self.expression_rank2 if self.expression_rank2 != '-' else 100
        return score

    def get_better_record(self, alt_record):
        if alt_record is None: return self

        self_score = self.get_score()
        alt_score = alt_record.get_score()
        if self_score != alt_score:
            return self if self_score < alt_score else alt_record

        return self if self.overlapped_length > alt_record.overlapped_length else alt_record

    def get_overlapped_length(self):
        length = 0
        for l, r in self.overlap_segments:
            length += r - l + 1
        return length

    def get_max_overlapped_segment(self):
        max_l = 0
        max_r = 0
        for l, r in self.overlap_segments:
            if r - l > max_r - max_l:
                max_l = l
                max_r = r
        return max_l, max_r
