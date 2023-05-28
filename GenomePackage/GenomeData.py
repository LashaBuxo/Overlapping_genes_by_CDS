import gc
import math
import os
import random
from os.path import exists
import gffutils
import numpy as np
import scipy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from gffutils import *
from numpy import median
from random import randrange
from os import walk

from scipy.stats import wilcoxon, ttest_ind

# Specific tools for working with Ensembl annotation and with sequence data
from GenomePackage.Enums import *
from GenomePackage.Filter import *
from Paths import *


class GenomeData:
    def __init__(self, species: SPECIES, annotation_load_type: ANNOTATION_LOAD, with_sequence=False, SNPs_path=None):

        print(f"loading data for {species}:")

        self.__transcript_expression_values = {}
        self.__transcript_expression_rank = {}

        if species == SPECIES.Homo_sapiens:
            self.__load_GTEx_data_multi_sample_tissues(f"{GTEx_EXPRESSION_DATA_PATH}")

        self.SNPs = {}
        if SNPs_path is not None:
            self.__load_SNPs_data(SNPs_path)

        gc.collect()

        self.species = species
        self.annotation_path = ANNOTATIONS_PATH[species]
        self.annotation_load_type = annotation_load_type
        self.sequence_path = SEQUENCES_PATH[species]

        self.chromosomes = []
        self.genes_on_chr = {}  # {chrom, [list of gene features who is located on chrom]}
        self.gene_transcripts = {}  # {gene_id, [list of mRNA features whose parent is gene_id]}
        self.transcript_fragments = {}  # {transcript_id, [list of Exons/UTRs/CDSs features whose parent is gene_id]}

        self.__sequences = {}  # {chrom, DNA sequence of chromosome with id chr_id}

        # storage for gene or transcript features to access by id
        self.__loaded_feature_by_id = {}

        self.__find_chromosomes()
        self.__load_requested_features()
        if with_sequence:
            self.__load_requested_sequences()

        self.all_transcript_list = list(self.transcript_fragments.keys())

    def is_transcript_cds_annotated(self, transcript):
        frags = self.transcript_fragments[transcript.id]
        flag1 = False
        flag2 = False
        for frag in frags:
            if frag.featuretype == 'five_prime_UTR': flag1 = True
            if frag.featuretype == 'three_prime_UTR': flag2 = True
        return flag1 and flag2

    def __load_SNPs_data(self, path):
        file = open(path, "r")
        lines = file.readlines()
        for line in lines:
            arr = line.split('\t')
            chrom = arr[0]
            ind = int(arr[1])
            from_nt = arr[2]
            to_nt = arr[3]
            if not self.SNPs.__contains__(chrom):
                self.SNPs[chrom] = []
            if len(from_nt) == 1 and len(to_nt) == 1:
                self.SNPs[chrom].append((ind, from_nt, to_nt))
            elif len(from_nt) == 1 and len(to_nt) > 1 and to_nt.__contains__(','):
                to_nt = to_nt.replace(' ', '')
                arr = to_nt.split(',')
                for to_nt in arr:
                    self.SNPs[chrom].append((ind, from_nt, to_nt))

        file.close()
        gc.collect()
        print(f"SNPs data loaded!")

    # region random

    def get_random_sub_cds_location(self, sub_cds_length):
        while True:
            rand_index = random.randint(0, len(self.all_transcript_list) - 1)
            random_transcript_id = self.all_transcript_list[rand_index]

            frags = self.transcript_fragments[random_transcript_id]
            if len(frags) == 0: continue

            valid_frags = []
            for frag in frags:
                frag_len = frag.end - frag.start + 1
                if frag.featuretype != 'CDS': continue
                if sub_cds_length is None or frag_len > sub_cds_length:
                    valid_frags.append(frag)
            if len(valid_frags) == 0: continue
            random_frag = random.choice(valid_frags)
            if sub_cds_length is None: sub_cds_length = random.randint(0, len(random_frag) - 1)

            random_l = random.randint(random_frag.start, random_frag.end - sub_cds_length + 1)
            l, r = random_l, random_l + sub_cds_length - 1

            return self.get_segment_sequence(random_frag.chrom, random_frag.strand, l, r), random_frag.chrom, l, r


    def get_random_sub_cds(self, sub_cds_length=None, phased=False):
        while True:
            rand_index = random.randint(0, len(self.all_transcript_list) - 1)
            random_transcript_id = self.all_transcript_list[rand_index]

            frags = self.transcript_fragments[random_transcript_id]
            if len(frags) == 0: continue

            valid_frags = []
            for frag in frags:
                frag_len = frag.end - frag.start + 1
                if frag.featuretype != 'CDS': continue
                if sub_cds_length is None or frag_len > sub_cds_length:
                    valid_frags.append(frag)
            if len(valid_frags) == 0: continue
            random_frag = random.choice(valid_frags)
            if sub_cds_length is None: sub_cds_length = random.randint(0, len(random_frag) - 1)

            random_l = random.randint(random_frag.start, random_frag.end - sub_cds_length + 1)
            l, r = random_l, random_l + sub_cds_length - 1

            if phased:
                parent_id = GetFeatureParentID(random_frag)
                parent_transcript = self.get_feature_by_id(parent_id)
                return self.get_transcript_sub_cds(parent_transcript, l, r, phased=True)

            return self.get_segment_sequence(random_frag.chrom, random_frag.strand, l, r)

    # endregion
    def get_transcript_cds(self, transcript: Feature):
        frags = self.transcript_fragments[transcript.id]
        CDSs = []
        for frag in frags:
            if frag.featuretype == 'CDS':
                seq = self.get_segment_sequence(transcript.chrom, transcript.strand, frag.start, frag.end)
                CDSs.append((frag.start, frag.end, seq, frag))

        CDSs.sort(key=lambda x: x[0], reverse=(True if transcript.strand == '-' else False))
        res_seq = ""
        for cds in CDSs:
            res_seq += cds[2]
        if int(CDSs[0][3].frame) != 0:
            print(f"\ttranscript:{transcript.id} first CDS not 0-phased!")
        if len(res_seq) % 3 != 0:
            print(f"\ttranscript:{transcript.id} CDS not 3 dividable!")
        return res_seq.upper()

    def get_transcript_sub_cds(self, transcript: Feature, l, r, phased=False):
        frags = self.transcript_fragments[transcript.id]
        target_CDS = None
        for frag in frags:
            if frag.featuretype == 'CDS':
                if frag.start <= l and frag.end >= r:
                    target_CDS = frag
        assert target_CDS is not None
        phase = int(target_CDS.frame)
        cds_seq = self.get_segment_sequence(transcript.chrom, transcript.strand, target_CDS.start, target_CDS.end)
        l -= target_CDS.start
        r -= target_CDS.start

        # reverse target interval, as we retrieved reverse-complement sequence
        if transcript.strand == '-':
            r_temp = len(cds_seq) - 1 - l
            l_temp = len(cds_seq) - 1 - r
            l, r = l_temp, r_temp
        if phased:
            if phase != 0:
                cds_seq = cds_seq[phase:]
                l -= phase
                r -= phase
                if l < 0: l = 0
            if l % 3 != 0:
                l += 3 - (l % 3)
            r -= (r - l + 1) % 3
            if l > r: return ''

        return cds_seq[l:r + 1].upper()

    # region SNPs
    def get_transcript_SNPs_ratio(self, transcript: Feature):
        frags = self.transcript_fragments[transcript.id]
        total = 0
        snps_total = 0
        snps_synonymous = 0
        total_freqs = {0: 0, 1: 0, 2: 0, 3: 0}
        for frag in frags:
            if frag.featuretype == 'CDS':
                total += frag.end - frag.start + 1
                x, y, freqs = self.get_SNPs_count_on_transcript_sub_cds(transcript, frag.start,
                                                                        frag.end)
                snps_synonymous += x
                snps_total += y
                for ind in freqs.keys():
                    total_freqs[ind] += freqs[ind]

        return snps_synonymous, snps_total, total, total_freqs

    def get_SNPs_count_on_transcript_sub_cds(self, transcript: Feature, x, y):
        if not self.SNPs.__contains__(transcript.chrom): return 0, 0, {0: 0, 1: 0, 2: 0, 3: 0}
        arr = self.SNPs[transcript.chrom]
        l_ind = 0
        r_ind = len(arr) - 1
        while l_ind < r_ind:
            mid = (l_ind + r_ind) // 2
            if arr[mid][0] < x:
                l_ind = mid + 1
            else:
                r_ind = mid
        if arr[l_ind][0] < x: return 0, 0, {0: 0, 1: 0, 2: 0, 3: 0}
        start_index = l_ind
        l_ind = 0
        r_ind = len(arr) - 1
        while l_ind < r_ind:
            if l_ind == r_ind - 1:
                if arr[r_ind][0] <= y:
                    l_ind = r_ind
                break
            mid = (l_ind + r_ind) // 2
            if arr[mid][0] > y:
                r_ind = mid - 1
            else:
                l_ind = mid
        end_index = l_ind
        if start_index > end_index: return 0, 0, {0: 0, 1: 0, 2: 0, 3: 0}

        frags = self.transcript_fragments[transcript.id]
        CDSs = []
        for frag in frags:
            if frag.featuretype == 'CDS':
                CDSs.append((frag.start, frag.end))
        CDSs.sort(key=lambda x: x[0], reverse=(True if transcript.strand == '-' else False))

        synonymous = 0
        cds_seq = self.get_transcript_cds(transcript).upper()
        freq = {0: 0, 1: 0, 2: 0, 3: 0}
        for i in range(start_index, end_index + 1):
            status, index = self.__is_SNP_synonymous(arr[i][0], arr[i][1], arr[i][2], transcript.strand, CDSs, cds_seq)
            if status: synonymous += 1
            freq[index] += 1

        return synonymous, end_index - start_index + 1, freq

    def __revcomp(self, ch):
        if ch == 'C':
            return 'G'
        elif ch == 'G':
            return 'C'
        elif ch == 'T':
            return 'A'
        elif ch == 'A':
            return 'T'
        print(ch)
        return 'C'

    def __is_SNP_synonymous(self, index, ref_allele, alt_allele, strand, cds_intervals, cds_seq):
        if ref_allele is None or alt_allele is None: return False, 3
        if len(ref_allele) != 1 or len(alt_allele) != 1: return False, 3

        cur_ind = 0
        for cds_l, cds_r in cds_intervals:
            if cds_l <= index <= cds_r:
                cur_ind += index - cds_l if strand == '+' else cds_r - index
                break
            else:
                cur_ind += cds_r - cds_l + 1

        nashti = cur_ind % 3
        if strand == '-':
            ref_allele = self.__revcomp(ref_allele)
            alt_allele = self.__revcomp(alt_allele)

        if cds_seq[cur_ind] != ref_allele:
            print(f"{index} {ref_allele} {alt_allele}")
            print(cds_seq[cur_ind])
        assert cds_seq[cur_ind] == ref_allele
        ref_codon = cds_seq[cur_ind - nashti] + cds_seq[cur_ind - nashti + 1] + cds_seq[cur_ind - nashti + 2]
        if nashti == 0:
            alt_codon = alt_allele + cds_seq[cur_ind + 1] + cds_seq[cur_ind + 2]
        if nashti == 1:
            alt_codon = cds_seq[cur_ind - 1] + alt_allele + cds_seq[cur_ind + 1]
        if nashti == 2:
            alt_codon = cds_seq[cur_ind - 2] + cds_seq[cur_ind - 1] + alt_allele

        aa_seq_ref = Seq(ref_codon).translate(to_stop=False)
        aa_seq_alt = Seq(alt_codon).translate(to_stop=False)
        return aa_seq_ref == aa_seq_alt, nashti

    # endregion

    def get_feature_by_id(self, id):
        if self.__loaded_feature_by_id.__contains__(id):
            return self.__loaded_feature_by_id[id]
        return None

    def get_segment_protein_range_of_transcript(self, transcript: Feature, l, r):
        frags = self.transcript_fragments[transcript.id]
        CDSs = []
        for frag in frags:
            if frag.featuretype == 'CDS':
                CDSs.append((frag.start, frag.end))
        CDSs.sort(key=lambda x: x[0], reverse=(True if transcript.strand == '-' else False))

        start_cds_ind = -1
        cur_cds_ind = 0

        for index in range(len(CDSs)):
            cds_segment = CDSs[index if transcript.strand == '+' else (len(CDSs) - index - 1)]
            if l >= cds_segment[0] and r <= cds_segment[1]:
                start_cds_ind = cur_cds_ind
                start_cds_ind += (l - cds_segment[0] + 1) if transcript.strand == '+' else cds_segment[1] - r + 1
            cur_cds_ind += cds_segment[1] - cds_segment[0] + 1

        assert start_cds_ind != -1
        start_cds_ind -= 1
        end_cds_ind = start_cds_ind + r - l
        return start_cds_ind // 3, end_cds_ind // 3

    def get_transcript_protein_sequence(self, transcript: Feature):
        frags = self.transcript_fragments[transcript.id]
        CDSs = []
        for frag in frags:
            if frag.featuretype == 'CDS':
                seq = self.get_segment_sequence(transcript.chrom, transcript.strand, frag.start, frag.end)
                CDSs.append((frag.start, frag.end, seq))

        CDSs.sort(key=lambda x: x[0], reverse=(True if transcript.strand == '-' else False))
        res_seq = ""
        for cds in CDSs:
            res_seq += cds[2]

        cds_seq = Seq(res_seq)
        aa_seq = cds_seq.translate(to_stop=False)
        # if not aa_seq.endswith("*"):
        #     print(aa_seq)

        return aa_seq[:-1]

    def get_promotor_BED_line(self, transcript: Feature):
        promotor_range_upstream = 500
        promotor_range_downstream = 100
        frags = self.transcript_fragments[transcript.id]

        l, r = 0, 0
        for frag in frags:
            if frag.featuretype == 'exon':
                l = frag.start if l == 0 else min(l, frag.start)
                r = frag.end if r == 0 else max(r, frag.end)

        if transcript.strand == '+':
            l, r = l - promotor_range_upstream, l + promotor_range_downstream
        else:
            l, r = r - promotor_range_downstream, r + promotor_range_upstream

        return f"chr{transcript.chrom}\t{l}\t{r}"

    def find_CGIs_overlapped_on_interval(self, chrom, l, r):
        cgi_min_length = 200
        max_oe = 0
        for start_index in range(l - cgi_min_length + 1, r):
            seq = self.get_segment_sequence(chrom, '+', start_index, start_index + cgi_min_length - 1).lower()
            g, c, a, t, observed = seq.count('g'), seq.count('c'), seq.count('a'), seq.count(
                't'), seq.count('cg')

            total = g + c + a + t
            if total < len(seq):
                print("lasha")

            if (g + c) > (total / 2):
                observed = observed / total
                expected = (c + g) * (c + g) / 4 / total / total
                max_oe = max(max_oe, observed / expected)
        return '-' if max_oe < 0.6 else "{:.2f}".format(max_oe)

    def get_random_pair_genes(self, count, min_distance, max_distance):
        valid_pairs = []
        for chrom in self.chromosomes:
            genes = self.genes_on_chr[chrom]
            for i in range(0, len(genes)):
                for j in range(i + 1, len(genes)):
                    gene1, gene2 = genes[i], genes[j]
                    l1, r1 = gene1.start, gene1.end
                    l2, r2 = gene2.start, gene2.end
                    dits = 0
                    if GetFeaturesOverlapType(gene1, gene2) != OVERLAP_TYPE.NONE:
                        dist = 0
                    else:
                        if l1 > r2:
                            dist = l1 - r2 + 1
                        else:
                            dist = l2 - r1 + 1
                    if min_distance <= dist <= max_distance:
                        valid_pairs.append((gene1, gene2))
        random.shuffle(valid_pairs)
        out = []
        for i in range(0, count):
            out.append(valid_pairs[i])
        return out

    def get_expression_data_tissues(self):
        return list(self.__tissues_and_samples_count.keys())

    # raw subject values
    def get_transcript_expressions_in_tissue_by_samples(self, transcript: Feature, tissue_name):
        if not self.__transcript_expression_values_by_tissues.__contains__(transcript.id):
            return []
        if not self.__transcript_expression_values_by_tissues[transcript.id].__contains__(tissue_name):
            return []
        return self.__transcript_expression_values_by_tissues[transcript.id][tissue_name]

    # median value from raw subject values in tissue
    def get_transcript_expression_in_tissue(self, transcript: Feature, tissue_name):
        values = self.get_transcript_expressions_in_tissue_by_samples(transcript, tissue_name)
        if len(values) == 0: return '-'
        return median(values)

    # median of medians of expression values in tissue
    def get_transcript_expressions_by_tissues(self, transcript: Feature):
        if self.__transcript_expression_values.__contains__(transcript.id):
            return self.__transcript_expression_values[transcript.id]

        if not self.__transcript_expression_values_by_tissues.__contains__(transcript.id):
            return {}

        data = {}
        exp_by_tissues = self.__transcript_expression_values_by_tissues[transcript.id]
        for tissue in exp_by_tissues.keys():
            data[tissue] = self.get_transcript_expression_in_tissue(transcript, tissue)

        self.__transcript_expression_values[transcript.id] = data
        return data

    def get_transcript_expression(self, transcript: Feature):
        tissue_values = self.get_transcript_expressions_by_tissues(transcript)
        values = []
        for tissue in tissue_values.keys():
            values.append(tissue_values[tissue])
        if len(values) == 0: return '-'
        return median(values)

    def get_transcript_expression_rank(self, query_transcript: Feature):
        if self.__transcript_expression_rank.__contains__(query_transcript.id):
            return self.__transcript_expression_rank[query_transcript.id]

        if self.get_transcript_expression(query_transcript) == '-':
            return '-'

        parent_gene_id = query_transcript.attributes['Parent'][0]

        arr = []
        for transcript in self.gene_transcripts[parent_gene_id]:
            exp_value = self.get_transcript_expression(transcript)
            if exp_value == '-':
                exp_value = -1
            arr.append((transcript.id, exp_value))
        arr.sort(key=lambda x: x[1], reverse=True)

        for i in range(0, len(arr)):
            if arr[i][1] == -1:
                self.__transcript_expression_rank[arr[i][0]] = '-'
            else:
                self.__transcript_expression_rank[arr[i][0]] = i + 1

        return self.get_transcript_expression_rank(query_transcript)

    def get_gene_expression(self, gene: Feature):
        transcripts = self.gene_transcripts[gene.id]
        sum = 0
        for transcript in transcripts:
            val = self.get_transcript_expression(transcript)
            if val != '-': sum += val
        return sum

    def get_gene_expression_in_tissue(self, gene: Feature, tissue):
        transcripts = self.gene_transcripts[gene.id]
        sum = 0
        for transcript in transcripts:
            val = self.get_transcript_expression_in_tissue(transcript, tissue)
            if val != '-': sum += val
        return sum

    def get_gene_expressions_in_tissue_by_samples(self, gene: Feature, tissue):
        transcripts = self.gene_transcripts[gene.id]
        sum_values = []
        for transcript in transcripts:
            values = self.get_transcript_expressions_in_tissue_by_samples(transcript, tissue)
            if len(sum_values) == 0:
                if len(values) > 0:
                    sum_values = values
            else:
                if len(values) > 0:
                    assert len(values) == len(sum_values)
                    for index in range(len(values)):
                        sum_values[index] += values[index]
        return sum_values

    # in subjects
    def feature_expression_values_in_tissue(self, feature: Feature, tissue):
        if not feature.id.__contains__("gene"):
            return self.get_transcript_expressions_in_tissue_by_samples(feature, tissue)
        return self.get_gene_expressions_in_tissue_by_samples(feature, tissue)

    #
    # def debug(self, tf: Feature, target_gene: Feature, tissue):
    #     if tf == None or target_gene == None:
    #         print("genes not found")
    #         return
    #     file = open(f"Debug-{GetGeneSymbol(tf)}-{GetGeneSymbol(target_gene)}.txt", 'w')
    #
    #     tf_positively_regulates_gene = []
    #     tf_expression = self.get_gene_expressions_in_tissue_by_samples(tf, tissue)
    #     gene_expression = self.get_gene_expressions_in_tissue_by_samples(target_gene, tissue)
    #     if len(tf_expression) != 0 and len(gene_expression) != 0:
    #         assert len(tf_expression) == len(gene_expression)
    #         n = len(tf_expression)
    #         arr = []
    #         for i in range(n):
    #             arr.append((gene_expression[i], tf_expression[i]))
    #         arr.sort(key=lambda x: x[1])
    #         file.write(f"{GetGeneSymbol(tf)}\t{tissue}\t")
    #         for i in range(n):
    #             file.write(f"{arr[i][1]}\t")
    #         file.write(f"\n{GetGeneSymbol(target_gene)}\t{tissue}\t")
    #         for i in range(n):
    #             file.write(f"{arr[i][0]}\t")
    #         file.write(f"\n")
    #
    #         # Define the threshold expression level for the transcription factor (TF)
    #         tf_threshold = np.percentile(tf_expression, 75)  # use the upper quartile as the threshold
    #
    #         # Divide the samples into two groups based on the expression level of the TF
    #         tf_high_index = np.where(np.array(tf_expression) >= tf_threshold)[0]
    #         tf_low_index = np.where(np.array(tf_expression) < tf_threshold)[0]
    #
    #         # Extract the gene expression levels for the two groups
    #         gene_high = np.array(gene_expression)[tf_high_index]
    #         gene_low = np.array(gene_expression)[tf_low_index]
    #
    #         # Perform a one-tailed t-test to test for upregulation
    #         t_statistic_up, p_value_up = ttest_ind(gene_high, gene_low, alternative='greater')
    #
    #         # Perform a one-tailed t-test to test for downregulation
    #         t_statistic_down, p_value_down = ttest_ind(gene_high, gene_low, alternative='less')
    #
    #         # Print the results
    #         print("TF threshold:", tf_threshold)
    #         print("Number of samples in high TF group:", len(tf_high_index))
    #         print("Number of samples in low TF group:", len(tf_low_index))
    #         print("Mean gene expression (high TF group):", np.mean(gene_high))
    #         print("Mean gene expression (low TF group):", np.mean(gene_low))
    #         print("t-statistic (upregulation):", t_statistic_up)
    #         print("p-value (upregulation):", p_value_up / 2)  # divide by 2 for one-tailed test
    #         print("t-statistic (downregulation):", t_statistic_down)
    #         print("p-value (downregulation):", p_value_down / 2)  # divide by 2 for one-tailed test
    #     file.close()
    #
    # def Tissues_where_TF_regulates_feature(self, tf: Feature, feature: Feature):
    #     tissues = self.get_expression_data_tissues()
    #
    #     result = ""
    #     for tissue in tissues:
    #         tf_expression = self.get_gene_expressions_in_tissue_by_samples(tf, tissue)
    #         feature_expression = self.__feature_Expression_values_in_tissue(feature, tissue)
    #         if len(tf_expression) != 0 and len(feature_expression) != 0:
    #             assert len(tf_expression) == len(feature_expression)
    #
    #             # Define the threshold expression level for the transcription factor (TF)
    #             tf_threshold = np.percentile(tf_expression, 75)  # use the upper quartile as the threshold
    #
    #             # Divide the samples into two groups based on the expression level of the TF
    #             tf_high_index = np.where(np.array(tf_expression) >= tf_threshold)[0]
    #             tf_low_index = np.where(np.array(tf_expression) < tf_threshold)[0]
    #
    #             # Extract the gene expression levels for the two groups
    #             feature_high = np.array(feature_expression)[tf_high_index]
    #             feature_low = np.array(feature_expression)[tf_low_index]
    #
    #             # Perform a one-tailed t-test to test for upregulation
    #             t_statistic_up, p_value_up = ttest_ind(feature_high, feature_low, alternative='greater')
    #
    #             # Perform a one-tailed t-test to test for downregulation
    #             t_statistic_down, p_value_down = ttest_ind(feature_high, feature_low, alternative='less')
    #
    #             if p_value_up < 0.00001:
    #                 result += f"{tissue}:+ "
    #             if p_value_down < 0.00001:
    #                 result += f"{tissue}:- "
    #     return result

    def get_segment_sequence(self, chrom, strand, start, end) -> str:
        if not self.__sequences.__contains__(chrom):
            print("Sequence must be loaded during processing! Error...")
            return ""

        seq_record = self.__sequences[chrom]
        interval_record = seq_record[max(0, start - 1):min(end, len(seq_record))]
        if strand == '-':
            interval_record = interval_record.reverse_complement()

        return str(interval_record.seq)

    def __find_chromosomes(self):
        annotation = open(self.annotation_path, 'r')
        while True:
            line = annotation.readline()
            if line.__contains__("#!"):
                break
            if line.__contains__("##sequence-region"):
                self.chromosomes.append(line.split(' ')[-3])

    # def __load_TFs_and_motifs_data(self, motifs_directory_path, TFs_motifs_file_path):
    #     print(f"Loading motifs and calculating best scores...")
    #     if not os.path.exists('./temp'): os.makedirs('./temp')
    #     temp_file_path = "./temp/motifs_best_scores.txt"
    #     if exists(temp_file_path):
    #         temp_file = open(temp_file_path, "r")
    #         for line in temp_file.readlines():
    #             arr = line.replace('\n', '').split('\t')
    #             motif_id, best_score = arr[0], float(arr[1])
    #             file_path = f"{motifs_directory_path}{motif_id}.txt"
    #             pwm = GetMotifPwmFromFile(file_path)
    #             self.__motifs_PWMs_and_max_score[motif_id] = (pwm, best_score)
    #         temp_file.close()
    #     else:
    #         temp_file = open(temp_file_path, "w")
    #         for (_, _, filenames) in walk(motifs_directory_path):
    #             cnt = 0
    #             for file_name in filenames:
    #                 cnt += 1
    #                 file_path = f"{motifs_directory_path}{file_name}"
    #                 motif_id = file_name.replace('.txt', '')
    #                 pwm = GetMotifPwmFromFile(file_path)
    #                 # best score from 10^5 variant
    #                 best_score = GetMotifBestMatchScore(pwm, 1000000)
    #                 self.__motifs_PWMs_and_max_score[motif_id] = (pwm, best_score)
    #                 temp_file.write(f"{motif_id}\t{best_score}\n")
    #                 print(cnt)
    #         temp_file.close()
    #     print(f"Calculated BestMatchScore for {len(self.__motifs_PWMs_and_max_score.keys())} Motifs (10^5 random seq)!")
    #
    #     TFs_motifs_file = open(TFs_motifs_file_path, "r")
    #     lines = TFs_motifs_file.readlines()
    #     for index in range(1, len(lines)):
    #         arr = lines[index].replace('\n', '').split('\t')
    #         gene_id, gene_sym, motif_id = arr[0], arr[1], arr[6]
    #
    #         if not self.__motifs_PWMs_and_max_score.__contains__(motif_id):
    #             # print("Maybe License required?")
    #             continue
    #         if not self.__motifs_with_associated_TFs.__contains__(motif_id):
    #             self.__motifs_with_associated_TFs[motif_id] = []
    #         self.__motifs_with_associated_TFs[motif_id].append(f"gene:{gene_id}")
    #     TFs_motifs_file.close()
    #     print(f"{len(self.__motifs_with_associated_TFs.keys())} motifs associated with TFs!")

    def __load_GTEx_data_multi_sample_tissues(self, GTEx_multi_tissue_multi_sample_file_path):
        self.__transcript_expression_values_by_tissues = {}
        self.__transcript_expression_values_by_tissues_and_age = {}
        self.__tissues_and_samples_count = {}

        file = open(GTEx_multi_tissue_multi_sample_file_path, 'r')
        lines = file.readlines()

        header_tissues = lines[0].split('\t')
        for j in range(1, len(header_tissues)):
            if header_tissues[j] == '\n' or header_tissues[j] == '': continue
            tissue = header_tissues[j]
            if not self.__tissues_and_samples_count.__contains__(tissue):
                self.__tissues_and_samples_count[tissue] = 0
            self.__tissues_and_samples_count[tissue] += 1

        print("Loading Tissue Expression Data:")
        for tissue_name in self.__tissues_and_samples_count.keys():
            print(f"{self.__tissues_and_samples_count[tissue_name]} samples per {tissue_name}...")

        for i in range(1, len(lines)):
            data = lines[i].split('\t')
            transcript_id = f"transcript:{data[0]}"
            self.__transcript_expression_values_by_tissues[transcript_id] = {}
            for j in range(1, len(data)):
                if data[j] == '\n' or data[j] == '': continue
                tissue_name = header_tissues[j]
                exp_value = float(data[j])
                if not self.__transcript_expression_values_by_tissues[transcript_id].__contains__(tissue_name):
                    self.__transcript_expression_values_by_tissues[transcript_id][tissue_name] = []
                self.__transcript_expression_values_by_tissues[transcript_id][tissue_name].append(exp_value)

                # if not self.__transcript_expression_values_by_tissues_and_age.__contains__(tissue_name):
                #     self.__transcript_expression_values_by_tissues_and_age[tissue_name] = {}
                #
                # if not self.__transcript_expression_values_by_tissues_and_age[tissue_name].__contains__(age):
                #     self.__transcript_expression_values_by_tissues_and_age[tissue_name][age] = []
                #
                # self.__transcript_expression_values_by_tissues_and_age[tissue_name][age].append(exp_value)
        print("Loaded!")
        file.close()
        gc.collect()

    # if database does not exists for specific chromosome, then
    # builds database from annotation file and fills/loads all
    # necessary features from the database.
    #
    # if database exists, then fills/loads all necessary data structures

    def __load_requested_features(self):
        db_file = f"{'./temp/'}{self.annotation_path.split('/')[-1]}{'.db'}"

        if not os.path.exists('../temp'):
            os.makedirs('../temp')

        if not exists(db_file):
            print("Creating database for " + self.species.name + " only for first RUN it takes that long!")
            gffutils.create_db(self.annotation_path,
                               dbfn=db_file,
                               verbose=True, force=False, keep_order=False, merge_strategy='create_unique',
                               sort_attribute_values=False)

        features_db = gffutils.FeatureDB(db_file, keep_order=False)

        self.__load_genes(features_db)
        if self.annotation_load_type == ANNOTATION_LOAD.GENES:
            return
        self.__load_transcripts_and_utrs(features_db)
        if self.annotation_load_type == ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS:
            return
        self.__load_CDS(features_db)
        if self.annotation_load_type == ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS:
            return
        self.__load_exons(features_db)

    def __load_genes(self, features_db: FeatureDB):
        features_generator = features_db.features_of_type('gene')
        feature_genes = list(features_generator)

        self.genes_on_chr = {}
        for chrom in self.chromosomes:
            self.genes_on_chr[chrom] = []

        # choose only genes who has protein_coding attribute
        for gene in feature_genes:
            if (gene.featuretype == 'gene' and gene.attributes.__contains__('biotype')
                    and gene.attributes['biotype'][0] == 'protein_coding'):
                assert self.genes_on_chr.__contains__(gene.chrom)
                self.genes_on_chr[gene.chrom].append(gene)  # store gene feature on chromosome list

        # filter some unnecessary genes
        self.genes_on_chr, filtered_genes_dict = filter_by_ensembl_attributes(self.genes_on_chr)

        loaded_chromosomes = 0
        loaded_genes = 0
        for chrom in self.chromosomes:
            if self.genes_on_chr[chrom] is not None and len(self.genes_on_chr[chrom]) > 0:
                loaded_chromosomes += 1
                loaded_genes += len(self.genes_on_chr[chrom])
            else:
                print(f'WARNING: Chromosome {chrom} not contains any genes!')
            for gene in self.genes_on_chr[chrom]:
                self.__loaded_feature_by_id[gene.id] = gene
                # self.__gene_symbols_set[GenomeWorker.get_gene_symbol(gene)] = gene.id
                # self.__gene_accessions_set[GenomeWorker.get_gene_accession(gene)] = gene.id

        print(
            f"\t{loaded_genes} protein-coding genes loaded from chromosomes {loaded_chromosomes}/{len(self.chromosomes)}")
        print(f"\t\tFiltered out Genes: {str(filtered_genes_dict)}")

    def __load_transcripts_and_utrs(self, features_db: FeatureDB):
        features_generator = features_db.features_of_type('mRNA')
        loaded_transcripts = 0
        for mRNA in list(features_generator):
            if self.__load_feature_by_type(self.gene_transcripts, mRNA, 'mRNA'): loaded_transcripts += 1
        print(f"\t{loaded_transcripts} mRNAs loaded successfully!")

        features_generator = features_db.features_of_type('three_prime_UTR')
        loaded_utr3s = 0
        for utr3 in list(features_generator):
            if self.__load_feature_by_type(self.transcript_fragments, utr3, 'three_prime_UTR'): loaded_utr3s += 1
        print(f"\t{loaded_utr3s} UTR3s loaded successfully!")

        features_generator = features_db.features_of_type('five_prime_UTR')
        loaded_utr5s = 0
        for utr5 in list(features_generator):
            if self.__load_feature_by_type(self.transcript_fragments, utr5, 'five_prime_UTR'): loaded_utr5s += 1
        print(f"\t{loaded_utr5s} UTR5s loaded successfully!")

    def __load_CDS(self, features_db: FeatureDB):
        # load preferred fragments and link them to transcripts in dictionaries
        features_generator = features_db.features_of_type('CDS')
        loaded_CDSs = 0
        for CDS in list(features_generator):
            if self.__load_feature_by_type(self.transcript_fragments, CDS, 'CDS'): loaded_CDSs += 1
        print(f"\t{loaded_CDSs} CDSs loaded successfully!")

    def __load_exons(self, features_db: FeatureDB):
        features_generator = features_db.features_of_type('exon')
        loaded_exons = 0
        for exon in list(features_generator):
            if self.__load_feature_by_type(self.transcript_fragments, exon, 'exon'): loaded_exons += 1
        print(f"\t{loaded_exons} EXONs loaded successfully!")

    def __load_feature_by_type(self, dict_to_fill, feature: Feature, feature_type):
        if feature.featuretype != feature_type: return False

        assert feature.attributes.__contains__('Parent') and len(feature.attributes['Parent']) == 1
        assert len(feature.attributes['Parent']) == 1

        # parent_gene_id for mRNA and parent_transcript_id for CDS/UTR5'/UTR3'/Exon
        parent_id = feature.attributes['Parent'][0]

        # for mRNA parent is gene, and for CDS,Exon,UTR's it is transcript.
        parent_is_loaded = self.__loaded_feature_by_id.__contains__(parent_id)

        if not parent_is_loaded:
            return False  # it seems feature belongs to transcript or gene, which was filtered out

        if not dict_to_fill.__contains__(parent_id):
            dict_to_fill[parent_id] = []

        dict_to_fill[parent_id].append(feature)

        assert not self.__loaded_feature_by_id.__contains__(feature.id)

        self.__loaded_feature_by_id[feature.id] = feature

        return True

    def __load_requested_sequences(self):
        for record in SeqIO.parse(self.sequence_path, 'fasta'):
            self.__sequences[record.id] = record

        loaded_chromosomes = 0
        for chrom in self.chromosomes:
            loaded_chromosomes += 1 if self.__sequences[chrom] is not None and len(self.__sequences[chrom]) > 0 else 0

        print(f"\tDNA chromosomes loaded {loaded_chromosomes}/{len(self.chromosomes)}!")
