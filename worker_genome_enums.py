# region ENUM classes
import enum
from enum import IntEnum


class OVERLAP_TYPE(IntEnum):
    NONE = 0
    # Specific classification of overlaps
    CONVERGENT = 1  # different strand + convergent
    DIFF_NESTED = 2  # different strand + nested
    DIVERGENT = 3  # different strand + divergent
    SAME_NESTED = 4  # one of the gene entirely located into others boundaries on same strand
    TANDEM = 5  # genes located on same strand and neither is nested
    MULTI = 6

    def short_name(self):
        if self == OVERLAP_TYPE.CONVERGENT:
            return "Convergent"
        elif self == OVERLAP_TYPE.DIFF_NESTED:
            return "Nested (diff strand)"
        elif self == OVERLAP_TYPE.DIVERGENT:
            return "Divergent"
        elif self == OVERLAP_TYPE.SAME_NESTED:
            return "Nested (same strand)"
        elif self == OVERLAP_TYPE.TANDEM:
            return "Tandem"
        elif self == OVERLAP_TYPE.MULTI:
            return "MULTI TYPE"
        else:
            assert False

    @staticmethod
    def get_overlap_types():
        return [OVERLAP_TYPE.CONVERGENT, OVERLAP_TYPE.DIFF_NESTED, OVERLAP_TYPE.DIVERGENT,
                OVERLAP_TYPE.SAME_NESTED, OVERLAP_TYPE.TANDEM, OVERLAP_TYPE.MULTI]


class ANNOTATION_LOAD(enum.Enum):
    GENES = 1
    GENES_AND_TRANSCRIPTS = 2
    GENES_AND_TRANSCRIPTS_AND_CDS = 3
    GENES_AND_TRANSCRIPTS_AND_FRAGMENTS = 4


class SEQUENCE_LOAD(enum.Enum):
    LOAD = 1
    NOT_LOAD = 2


class SPECIES(enum.Enum):
    Homo_sapiens = 0
    Mus_musculus = 1

    @staticmethod
    def from_string(str):
        if str == "Homo sapiens":
            return SPECIES.Homo_sapiens
        elif str == "Mus musculus":
            return SPECIES.Mus_musculus
        else:
            if not str.__contains__('_'): assert False
            return SPECIES.from_string(str.replace('_', ' '))

    def shortest_name(self):
        if self == SPECIES.Homo_sapiens:
            return "Hs"
        elif self == SPECIES.Mus_musculus:
            return "Ms"
        else:
            assert False

    def short_name(self):
        if self == SPECIES.Homo_sapiens:
            return "Human"
        elif self == SPECIES.Mus_musculus:
            return "Mouse"
        else:
            assert False

    def __str__(self):
        return str(self.name).replace('_', ' ')


class TRANSCRIPT_CRITERIA(enum.Enum):
    NONE = 0
    LONGEST = 1
    LONGEST_CDS = 2
    LONGEST_CDS_AND_UTRs = 3  # similar to exons
    RANDOM = 4
