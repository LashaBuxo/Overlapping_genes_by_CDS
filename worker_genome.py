from os.path import exists
import gffutils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from gffutils import Feature

# Specific tools for working with Ensembl annotation and with sequence data
from worker_genome_values import *
from worker_genome_enums import *


class GenomeWorker:
    def __init__(self, species: SPECIES, annotation_load_type: ANNOTATION_LOAD, sequence_load_type: SEQUENCE_LOAD):

        print(f"loading data for {species}:")

        self.species = species
        self.annotation_load_type = annotation_load_type
        self.sequence_load_type = sequence_load_type

        self.imported_protein_coding_genes = 0
        self.ignored_protein_coding_genes = 0
        self.ignored_genes_by_types = {}

        # storage for gene or transcript features to access by id
        self.__loaded_feature_by_id = {}

        # gene features clustered by chromosome, they are on.
        self.__genes_on_chr = [None]
        self.__gene_symbols_set = {}
        self.__gene_accessions_set = {}

        self.__gene_transcripts = {}  # {gene_id, [list of mRNA features/isoforms whose parent is gene_id]}

        self.__transcript_APPRIS_data = {}

        self.__transcript_fragments = {}
        self.__transcript_fragments_is_sorted = {}  # optimizes speed, if we are not sorting it if we don't need

        self.__gene_transcript_by_criteria = {}  # {gene_id, representative mRNA features whose parent is gene_id}

        # {chr_id,DNA sequence of chromosome with id chr_id}
        self.__sequences = [None]

        self.__load_requested_features()
        self.__load_APPRIS_data()
        self.__load_requested_sequences()

        self.__gene_mitoCarta_data = {}

    # region load methods for Annotation, Assembly(sequence), APPRIS and mitoCarta

    def __load_APPRIS_data(self):
        scores_file = APPRIS_DATA[self.species.value]
        assert exists(scores_file)
        x = {}
        file = open(scores_file, 'r')
        lines = file.readlines()
        for index in range(1, len(lines)):
            line_data = lines[index].split('\t')
            transcript_id = line_data[2]

            if self.feature_by_id('transcript:' + transcript_id) is None: continue

            # if self.get_transcript_APPRIS_status(transcript_id) is None: continue
            residues, structure, conservation, domains, helices, signals, trifid_score, mapped_peptides, appris_score, appris_annotation = \
                line_data[10], line_data[11], line_data[12], line_data[13], line_data[14], line_data[15], line_data[16], \
                line_data[17], line_data[18], line_data[19]

            not_found_tag = line_data[6]

            trifid_score = 0 if trifid_score == '' else float(trifid_score)
            if not self.__transcript_APPRIS_data.__contains__(transcript_id):
                self.__transcript_APPRIS_data[transcript_id] = {}
            self.__transcript_APPRIS_data[transcript_id]['functional_residues'] = residues
            self.__transcript_APPRIS_data[transcript_id]['structure_score'] = structure
            self.__transcript_APPRIS_data[transcript_id]['domain_residues'] = domains
            self.__transcript_APPRIS_data[transcript_id]['conservation_score'] = conservation
            self.__transcript_APPRIS_data[transcript_id]['Trifid_Score'] = trifid_score
            self.__transcript_APPRIS_data[transcript_id]['APPRIS_score'] = appris_score
            self.__transcript_APPRIS_data[transcript_id]['APPRIS_annotation'] = appris_annotation
            self.__transcript_APPRIS_data[transcript_id]['transcript_incomplete'] = not_found_tag

            x[not_found_tag] = 1
        print(f"\t{len(self.__transcript_APPRIS_data)}/{len(lines) - 1} transcripts scores loaded from APPRIS!")

        if self.species != SPECIES.Homo_sapiens: return
        assert exists(APPRIS_HOMOLOGS_DATA_PATH)
        file = open(APPRIS_HOMOLOGS_DATA_PATH, 'r')
        lines = file.readlines()
        cnt = 0
        transcript_id = ''
        for index in range(0, len(lines)):
            line = lines[index].replace('"', '')
            if line.startswith('>'):
                transcript_id = line.replace('>', '').split('\t')[0]
                cnt += 1
                if not self.__transcript_APPRIS_data.__contains__(transcript_id):
                    self.__transcript_APPRIS_data[transcript_id] = {'homologue': []}
                elif not self.__transcript_APPRIS_data[transcript_id].__contains__('homologue'):
                    self.__transcript_APPRIS_data[transcript_id]['homologue'] = []
                continue
            if len(line) < 3: continue
            species = line.split('\t')[0]
            self.__transcript_APPRIS_data[transcript_id]['homologue'].append(species)

        print(f"\t{cnt} transcripts homologs (CORSAIR) loaded from APPRIS!")

    def get_gene_principal_transcript_id(self, gene_id):
        transcripts = self.get_transcripts_from_gene(gene_id)
        principal_transcript_ids = []
        for transcript in transcripts:
            transcript_id = transcript.id.replace("transcript:", "")
            if self.__transcript_APPRIS_data.__contains__(transcript_id):
                if self.__transcript_APPRIS_data[transcript_id]['APPRIS_annotation'].__contains__("PRINCIPAL"):
                    principal_transcript_ids.append(transcript.id)
        if len(principal_transcript_ids) == 0 or len(principal_transcript_ids) > 1:
            return None
        return principal_transcript_ids[0]

    def get_transcript_conservation_score(self, transcript_id):
        transcript_id = transcript_id.replace('transcript:', '')
        if not self.__transcript_APPRIS_data.__contains__(transcript_id): return 0
        score = float(self.__transcript_APPRIS_data[transcript_id]['conservation_score'])
        return score

    def get_transcript_homologue_species(self, transcript_id):
        transcript_id = transcript_id.replace('transcript:', '')
        if not self.__transcript_APPRIS_data.__contains__(transcript_id):
            print("APPRIS not contains any info for that transcript!")
            return []
        if not self.__transcript_APPRIS_data[transcript_id].__contains__('homologue'):
            print("APPRIS not contains homolog info for that transcript!")
            return []
        return self.__transcript_APPRIS_data[transcript_id]['homologue']

    def chromosomes_count(self):
        return NUMBER_OF_CHROMOSOMES[self.species.value]

    def genes_count_on_chr(self, chr_id):
        assert chr_id <= self.chromosomes_count()
        return 0 if self.__genes_on_chr[chr_id] is None else len(self.__genes_on_chr[chr_id])

    def gene_by_indexes(self, chr_id, index) -> Feature:
        assert chr_id <= self.chromosomes_count()
        return self.__genes_on_chr[chr_id][index]

    def gene_by_symbol(self, symbol) -> Feature:
        if self.__gene_symbols_set.__contains__(symbol):
            return self.feature_by_id(self.__gene_symbols_set[symbol])
        return None

    def get_transcript_parent(self, transcript_id):
        feature = self.feature_by_id(transcript_id)
        return self.feature_by_id(feature.attributes['Parent'][0])

    def feature_by_id(self, feature_id) -> Feature:
        return self.__loaded_feature_by_id[feature_id] if self.__loaded_feature_by_id.__contains__(feature_id) else None

    def __get_annotation_file_path(self):
        return ENSEMBL_ANNOTATIONS[self.species.value]

    def __get_annotation_db_path(self):
        added_name = self.__get_annotation_file_path().replace('/', '_')
        annotation_name = "ENSEMBL"
        generated_path = GENOME_DATABASES_DIRECTORY + self.species.name + "_" + annotation_name + "_" + added_name
        return generated_path + '.db'

    def __get_sequence_file_path(self):
        return ENSEMBL_SEQUENCES[self.species.value]

    def get_transcripts_from_gene(self, gene_id) -> list[Feature]:
        return self.__gene_transcripts[gene_id] if self.__gene_transcripts.__contains__(gene_id) else []

    def get_transcript_CDS(self, transcript_id):
        transcript = self.feature_by_id(transcript_id)
        frags = self.__transcript_fragments[transcript_id]

        chr_index = self.chr_name2index(transcript.chrom)
        cds = ""
        for index in range(len(frags)):
            frag = frags[index] if transcript.strand == '+' else frags[len(frags) - 1 - index]
            if frag.featuretype == 'CDS':
                cds += self.retrieve_segment_sequence(chr_index, frag.start, frag.end, frag.strand)
        return cds

    def get_transcript_CDS_length(self, transcript_id):
        frags = self.__transcript_fragments[transcript_id]
        CDS_length = 0
        for frag in frags:
            if frag.featuretype == 'CDS':
                CDS_length += frag.end - frag.start + 1
        return CDS_length

    def get_fragments_from_transcript(self, transcript_id) -> list[Feature]:
        if transcript_id not in self.__transcript_fragments: return []

        fragments = self.__transcript_fragments[transcript_id]
        if self.__transcript_fragments_is_sorted.__contains__(transcript_id):
            return fragments

        fragments.sort(key=lambda x: x.start)

        self.__transcript_fragments_is_sorted[transcript_id] = True
        self.__transcript_fragments[transcript_id] = fragments

        return fragments

    # endregion

    # region load methods for Annotation, Assembly(sequence)

    def __load_requested_sequences(self):
        if self.sequence_load_type != SEQUENCE_LOAD.LOAD:
            return

        sequence_file_path = self.__get_sequence_file_path()
        for record in SeqIO.parse(sequence_file_path, 'fasta'):
            chr_index = self.chr_name2index(record.id)
            if chr_index == -1: continue
            self.__sequences[chr_index] = record

        loaded_chromosomes = 0
        for chr_index in range(1, NUMBER_OF_CHROMOSOMES[self.species.value] + 1):
            loaded_chromosomes += 1 if self.__sequences[chr_index] is not None and len(
                self.__sequences[chr_index]) > 0 else 0

        print(f"\tDNA chromosomes loaded {loaded_chromosomes}/{self.chromosomes_count()}!")

    def __check_gene_for_Ensembl_filters(self, gene: Feature, feature_ids, feature_syms, feature_accs):
        # https: // www.biostars.org / p / 5304 /  # 9521037
        # ignore genes with miscellaneous chromosome/scaffold names
        chr_index = self.chr_name2index(gene.chrom)
        if chr_index == -1: return False, "located_on_miscellaneous_scaffold"

        gene_symbol = GenomeWorker.get_gene_symbol(gene)
        gene_accession = GenomeWorker.get_gene_accession(gene)
        description = GenomeWorker.get_gene_description(gene)

        # ignore genes without description or gene_symbol
        if description == "no_desc": return False, "no_desc"
        if gene_symbol == "no_sym": return False, "no_symbol"
        if gene_accession == "no_acc": return False, "no_accession"

        # ignore genes if they are pseudogene, novel or predicted, readthrough
        if description.__contains__('readthrough'):
            return False, "is_readthrough"
        if description.__contains__('pseudogene') or description.__contains__('pseudogene'):
            return False, "is_pseudogene"
        if description.__contains__('novel') or description.__contains__('Novel'):
            return False, "is_novel"
        if description.__contains__('predicted') or description.__contains__('Predicted'):
            return False, "is_predicted"

        # ignore genes with duplicate Names or accessions
        if gene_symbol != "no_sym" and feature_syms.__contains__(gene_symbol):
            return False, "name_duplicated"
        # ignore genes if gene with same NCBI accession already imported
        if gene_accession != 'no_acc' and feature_accs.__contains__(gene_accession):
            return False, "acc_duplicated"

        # for any annotation, there must not be 2 gene record with same id
        assert not feature_ids.__contains__(gene.id)

        return True, "passed_filter"

    def filter_by_ensembl_attributes(self):
        feature_ids = {}
        feature_syms = {}
        feature_accs = {}
        for chr_index in range(1, self.chromosomes_count() + 1):
            genes_on_chr = self.__genes_on_chr[chr_index]
            new_genes_on_chr = []
            cnt = len(genes_on_chr)
            for i in range(0, cnt):
                gene = genes_on_chr[i]
                is_valid, status = self.__check_gene_for_Ensembl_filters(gene, feature_ids, feature_syms, feature_accs)
                if is_valid:
                    new_genes_on_chr.append(gene)
                    feature_ids[gene.id] = True
                    feature_syms[GenomeWorker.get_gene_symbol(gene)] = True
                    feature_accs[GenomeWorker.get_gene_accession(gene)] = True
                else:
                    self.filter_gene_by_reason(status)

            self.__genes_on_chr[chr_index] = new_genes_on_chr

    def filter_gene_by_reason(self, status):
        if not self.ignored_genes_by_types.__contains__(status):
            self.ignored_genes_by_types[status] = 0
        self.ignored_genes_by_types[status] += 1
        self.ignored_protein_coding_genes += 1
        self.imported_protein_coding_genes -= 1

    # if database does not exists for specific chromosome, then
    # builds database from annotation file and fills/loads all
    # necessary features from the database.
    #
    # if database exists, then fills/loads all necessary data structures

    def __load_requested_features(self):
        if not exists(self.__get_annotation_db_path()):
            print("Creating database for " + self.species.name + " only for first RUN it takes that long!")
            gffutils.create_db(self.__get_annotation_file_path(),
                               dbfn=self.__get_annotation_db_path(),
                               verbose=True, force=False, keep_order=False, merge_strategy='create_unique',
                               sort_attribute_values=False)

        features_db = gffutils.FeatureDB(self.__get_annotation_db_path(), keep_order=False)
        features_generator = features_db.features_of_type('gene')
        feature_genes = list(features_generator)

        self.__genes_on_chr = [None] * (self.chromosomes_count() + 1)
        self.__sequences = [None] * (self.chromosomes_count() + 1)

        # choose genes who has protein_coding attribute and additional filter values

        for gene in feature_genes:
            if (gene.featuretype == 'gene' and gene.attributes.__contains__('biotype')
                    and gene.attributes['biotype'][0] == 'protein_coding'):
                chr_id = self.chr_name2index(gene.chrom)
                if self.__genes_on_chr[chr_id] is None:
                    self.__genes_on_chr[chr_id] = []

                self.imported_protein_coding_genes += 1
                self.__genes_on_chr[chr_id].append(gene)  # store gene feature on chromosome list

        self.filter_by_ensembl_attributes()

        # store gene feature->id, id->symbol, id->accession
        for chr_index in range(1, self.chromosomes_count() + 1):
            genes_on_chr = self.__genes_on_chr[chr_index]
            for i in range(0, len(genes_on_chr)):
                gene = genes_on_chr[i]
                self.__loaded_feature_by_id[gene.id] = gene
                self.__gene_symbols_set[GenomeWorker.get_gene_symbol(gene)] = gene.id
                self.__gene_accessions_set[GenomeWorker.get_gene_accession(gene)] = gene.id

        loaded_chromosomes = 0
        for chr_index in range(1, NUMBER_OF_CHROMOSOMES[self.species.value] + 1):
            if self.__genes_on_chr[chr_index] is not None and len(self.__genes_on_chr[chr_index]) > 0:
                loaded_chromosomes += 1
            else:
                print(f'WARNING: Chromosome {chr_index} not contains any genes!')

        print(f"\t{self.imported_protein_coding_genes} genes loaded "
              f"({self.ignored_protein_coding_genes} filtered out) "
              f"from chromosomes {loaded_chromosomes}/{self.chromosomes_count()}")
        print(f"\t\tFiltered out Genes: {str(self.ignored_genes_by_types)}")

        if self.annotation_load_type == ANNOTATION_LOAD.GENES:
            return

        features_generator = features_db.features_of_type('mRNA')
        features_mRNA = list(features_generator)
        loaded_transcripts = 0
        for mRNA in features_mRNA:
            if self.__load_feature_by_type(self.__gene_transcripts, mRNA, 'mRNA'): loaded_transcripts += 1
        print(f"\t{loaded_transcripts} mRNAs loaded successfully!")

        if self.annotation_load_type == ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS:
            return

        # load preferred fragments and link them to transcripts in dictionaries
        features_generator = features_db.features_of_type('CDS')
        features_CDSs = list(features_generator)
        loaded_CDSs = 0
        for CDS in features_CDSs:
            if self.__load_feature_by_type(self.__transcript_fragments, CDS, 'CDS'): loaded_CDSs += 1
        print(f"\t{loaded_CDSs} CDSs loaded successfully!")

        if self.annotation_load_type == ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS:
            return

        features_generator = features_db.features_of_type('exon')
        features_exons = list(features_generator)
        loaded_exons = 0
        for exon in features_exons:
            if self.__load_feature_by_type(self.__transcript_fragments, exon, 'exon'): loaded_exons += 1
        print(f"\t{loaded_exons} EXONs loaded successfully!")

        features_generator = features_db.features_of_type('three_prime_UTR')
        features_UTR3s = list(features_generator)
        loaded_utr3s = 0
        for utr3 in features_UTR3s:
            if self.__load_feature_by_type(self.__transcript_fragments, utr3, 'three_prime_UTR'): loaded_utr3s += 1
        print(f"\t{loaded_utr3s} UTR3s loaded successfully!")

        features_generator = features_db.features_of_type('five_prime_UTR')
        features_UTR5s = list(features_generator)
        loaded_utr5s = 0
        for utr5 in features_UTR5s:
            if self.__load_feature_by_type(self.__transcript_fragments, utr5, 'five_prime_UTR'): loaded_utr5s += 1
        print(f"\t{loaded_utr5s} UTR5s loaded successfully!")

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

    def chr_name2index(self, chr_name):
        if chr_name == 'MT': return NUMBER_OF_CHROMOSOMES[self.species.value]
        if chr_name == 'Y': return NUMBER_OF_CHROMOSOMES[self.species.value] - 1
        if chr_name == 'X': return NUMBER_OF_CHROMOSOMES[self.species.value] - 2

        try:
            index = int(chr_name)
        except ValueError:
            index = -1

        return index

    def chr_index2name(self, index):
        if index == NUMBER_OF_CHROMOSOMES[self.species.value]: return 'MT'
        if index == NUMBER_OF_CHROMOSOMES[self.species.value] - 1: return 'Y'
        if index == NUMBER_OF_CHROMOSOMES[self.species.value] - 2: return 'X'
        assert 1 <= index <= NUMBER_OF_CHROMOSOMES[self.species.value]
        return f"{index}"

    # endregion

    # region overlap findings between transcript coding sequences (CDS)

    def __calculate_frame_of_fragment_interval(self, fragment: Feature, interval):
        if fragment.strand == '+':
            return (3 + int(fragment.frame) - (interval[0] - fragment.start) % 3) % 3
        return (3 + int(fragment.frame) - (fragment.end - interval[1]) % 3) % 3

    # between each pairs of fragments from each transcript (mRNA) finds fragments overlapped by CDS

    def _are_features_same_framed(self, fragment_a: Feature, fragment_b: Feature) -> bool:
        l = fragment_a.start + int(fragment_a.frame) if fragment_a.strand == '+' else fragment_a.end - int(
            fragment_a.frame)
        r = fragment_b.start + int(fragment_b.frame) if fragment_b.strand == '+' else fragment_b.end - int(
            fragment_b.frame)
        return l % 3 == r % 3

    def get_overlaps_between_transcripts(self, transcript1_id, transcript2_id):
        non_ati_overlaps = []

        fragments_a = self.get_fragments_from_transcript(transcript1_id)
        fragments_b = self.get_fragments_from_transcript(transcript2_id)

        contains_in_phase_overlaps = False

        ov_length = 0
        for fragment_a in fragments_a:
            if fragment_a.featuretype != 'CDS': continue
            for fragment_b in fragments_b:
                if fragment_b.featuretype != 'CDS': continue
                if fragment_b.end < fragment_a.start or fragment_b.start > fragment_a.end: continue
                if fragment_b.end <= fragment_a.end:
                    if fragment_b.start >= fragment_a.start:
                        overlap = (fragment_b.start, fragment_b.end)
                    else:
                        overlap = (fragment_a.start, fragment_b.end,)
                else:
                    if fragment_b.start <= fragment_a.start:
                        overlap = (fragment_a.start, fragment_a.end)
                    else:
                        overlap = (fragment_b.start, fragment_a.end)

                if fragment_a.strand == fragment_b.strand and self._are_features_same_framed(fragment_a, fragment_b):
                    contains_in_phase_overlaps = True
                    continue

                ov_length += overlap[1] - overlap[0] + 1
                non_ati_overlaps.append(overlap)

        return non_ati_overlaps, ov_length, contains_in_phase_overlaps

    # endregion

    # region sequence retriever methods

    # if sequence for specific chromosome is not loaded then loads it.
    # if its already loaded retrieves it.
    def retrieve_sequence_record(self, chr_index) -> SeqRecord:
        if self.__sequences[chr_index] is None:
            print("Sequence must be loaded during processing! Error...")
        return self.__sequences[chr_index]

    def retrieve_feature_sequence(self, chr_id, feature: Feature) -> str:
        seq_record = self.retrieve_sequence_record(chr_id)
        feature_record = seq_record[feature.start - 1:feature.end]
        if feature.strand == '-':
            feature_record = feature_record.reverse_complement()
        return str(feature_record.seq)

    def retrieve_segment_sequence(self, chr_id, start, end, strand) -> str:
        seq_record = self.retrieve_sequence_record(chr_id)
        interval_record = seq_record[max(0, start - 1):min(end, len(seq_record))]
        if strand == '-': interval_record = interval_record.reverse_complement()
        return str(interval_record.seq)

    # endregion

    # region static Methods

    @staticmethod
    def get_gene_accession(gene: Feature):
        if not gene.attributes.__contains__('description'): return "no_acc"
        desc_parts = gene.attributes['description'][0].split("Acc:")
        if len(desc_parts) != 2: return "no_acc"  # it seems record does not have NCBI ID
        ncbi_id = desc_parts[1][0:(len(desc_parts[1]) - 1)]
        return ncbi_id

    @staticmethod
    def get_gene_symbol(gene: Feature):
        if not gene.attributes.__contains__('Name'): return "no_sym"  # loaded genes must be filtered
        assert len(gene.attributes['Name']) == 1
        return gene.attributes['Name'][0].upper()

    @staticmethod
    def get_gene_description(gene: Feature):
        if not gene.attributes.__contains__('description'): return "no_desc"  # loaded genes must be filtered
        return gene.attributes['description'][0]

    ########################################Distance################################################

    def get_features_overlap_type(self, segment1: Feature, segment2: Feature):
        if segment1.chrom != segment2.chrom: return OVERLAP_TYPE.NONE
        return self.get_segments_overlap_type((segment1.start, segment1.end, segment1.strand),
                                              (segment2.start, segment2.end, segment2.strand))

    # must be assumption that segments located on same chromosome
    def get_segments_overlap_type(self, segment1, segment2):
        assert segment1 is not None and segment2 is not None
        l_a, r_a, strand_a = segment1
        l_b, r_b, strand_b = segment2

        if (l_a <= l_b <= r_a and l_a <= r_b <= r_a) or (l_b <= l_a <= r_b and l_b <= r_a <= r_b):
            return OVERLAP_TYPE.SAME_NESTED if strand_a == strand_b else OVERLAP_TYPE.DIFF_NESTED

        if strand_a == strand_b:
            return OVERLAP_TYPE.TANDEM if l_a <= l_b <= r_a or l_a <= r_b <= r_a else OVERLAP_TYPE.NONE

        # if it comes there, it is DIFF strand
        if strand_a == '-':
            l_a, r_a, strand_a = segment2
            l_b, r_b, strand_b = segment1

        if l_a <= l_b <= r_a:
            return OVERLAP_TYPE.CONVERGENT

        if l_a <= r_b <= r_a:
            return OVERLAP_TYPE.DIVERGENT

        return OVERLAP_TYPE.NONE

    # retrieves nucleotide composition of sequence
    # output: (C_count,G_count,A_count,T_count)
    @staticmethod
    def sequence_GC(sequence):
        sequence = sequence.upper()
        gc = 0
        for char in sequence:
            assert 'CGAT'.__contains__(char)
            gc += 1 if char == 'C' or char == 'G' else 0

        return gc / len(sequence)

    # endregion
