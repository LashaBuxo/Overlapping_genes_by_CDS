import enum
from enum import IntEnum


class OVERLAP_TYPE(enum.Enum):
    NONE = 0
    # Specific classification of overlaps
    CONVERGENT = 1  # anti parallel + convergent
    NESTED_ANTI_PARALLEL = 2  # anti parallel + nested
    DIVERGENT = 3  # anti parallel + divergent
    NESTED_PARALLEL = 4  # one of the gene entirely located into others boundaries on same strand
    TANDEM = 5  # genes located on same strand and neither is nested
    MULTI = 6

    def short_name(self):
        if self == OVERLAP_TYPE.CONVERGENT:
            return "Convergent"
        elif self == OVERLAP_TYPE.NESTED_ANTI_PARALLEL:
            return "Nested (anti parallel)"
        elif self == OVERLAP_TYPE.DIVERGENT:
            return "Divergent"
        elif self == OVERLAP_TYPE.NESTED_PARALLEL:
            return "Nested (parallel)"
        elif self == OVERLAP_TYPE.TANDEM:
            return "Tandem"
        elif self == OVERLAP_TYPE.MULTI:
            return "MULTI TYPE"
        else:
            assert False

    @staticmethod
    def get_overlap_types():
        return [OVERLAP_TYPE.CONVERGENT, OVERLAP_TYPE.NESTED_ANTI_PARALLEL, OVERLAP_TYPE.DIVERGENT,
                OVERLAP_TYPE.NESTED_PARALLEL, OVERLAP_TYPE.TANDEM, OVERLAP_TYPE.MULTI]


class ANNOTATION_LOAD(enum.Enum):
    GENES = 1
    GENES_AND_TRANSCRIPTS = 2
    GENES_AND_TRANSCRIPTS_AND_CDS = 3
    GENES_AND_TRANSCRIPTS_AND_FRAGMENTS = 4


class SPECIES(enum.Enum):
    Homo_sapiens = 0
    Pan_troglodytes = 1,
    Mus_musculus = 2,
    Sus_scrofa = 3,
    Canis_lupus_familiaris = 4,
    Monodelphis_domestica = 5,
    Gallus_gallus = 6,
    Danio_rerio = 7,
    Drosophila_melanogaster = 8,

    def tax_id(self):
        if self == SPECIES.Homo_sapiens:
            return "9606"
        elif self == SPECIES.Pan_troglodytes:
            return "9598"
        elif self == SPECIES.Mus_musculus:
            return "10090"
        elif self == SPECIES.Sus_scrofa:
            return "9823"
        elif self == SPECIES.Canis_lupus_familiaris:
            return "9615"
        elif self == SPECIES.Monodelphis_domestica:
            return "13616"
        elif self == SPECIES.Gallus_gallus:
            return "9031"
        elif self == SPECIES.Danio_rerio:
            return "7955"
        elif self == SPECIES.Drosophila_melanogaster:
            return "7227"
        else:
            assert False

    def short_name(self):
        if self == SPECIES.Homo_sapiens:
            return "Human"
        elif self == SPECIES.Pan_troglodytes:
            return "Chimpanzee"
        elif self == SPECIES.Mus_musculus:
            return "Mouse"
        elif self == SPECIES.Sus_scrofa:
            return "Pig"
        elif self == SPECIES.Canis_lupus_familiaris:
            return "Dog"
        elif self == SPECIES.Monodelphis_domestica:
            return "Opossum"
        elif self == SPECIES.Gallus_gallus:
            return "Chicken"
        elif self == SPECIES.Danio_rerio:
            return "Zebrafish"
        elif self == SPECIES.Drosophila_melanogaster:
            return "Fruit fly"
        else:
            assert False

    def __str__(self):
        return str(self.name).replace('_', ' ')


class TRANSCRIPT_CHOOSE_CRITERIA(enum.Enum):
    NONE = 0
    LONGEST = 1
    LONGEST_CDS = 2
    LONGEST_CDS_AND_UTRs = 3  # similar to exons
    RANDOM = 4
