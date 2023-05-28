import numpy as np
from gffutils import Feature

from GenomePackage.Enums import *


def GetGeneAccession(gene: Feature):
    if not gene.attributes.__contains__('description'): return "no_acc"
    desc_parts = gene.attributes['description'][0].split("Acc:")
    if len(desc_parts) != 2: return "no_acc"  # it seems record does not have NCBI ID
    ncbi_id = desc_parts[1][0:(len(desc_parts[1]) - 1)]
    return ncbi_id


def GetGeneSymbol(gene: Feature):
    if not gene.attributes.__contains__('Name'): return "no_sym"  # loaded genes must be filtered
    assert len(gene.attributes['Name']) == 1
    return gene.attributes['Name'][0].upper()



def GetFeatureParentID(feature: Feature):
    if not feature.attributes.__contains__('Parent'): return None
    assert len(feature.attributes['Parent']) == 1
    return feature.attributes['Parent'][0]


def GetGeneEnsemblId(gene: Feature):
    return gene.id.replace('gene:', '')

def GetTranscriptEnsemblId(transcript: Feature):
    return transcript.id.replace('transcript:', '')

def GetGeneDescription(gene: Feature):
    if not gene.attributes.__contains__('description'): return "no_desc"  # loaded genes must be filtered
    return gene.attributes['description'][0]


def GetGeneFullName(gene: Feature):
    id = gene.id.replace('gene:', '')
    return f"{id}({GetGeneSymbol(gene)})" if GetGeneSymbol(gene) != "no_sym" else f"{id}"


def GetFeaturesOverlapType(feature1: Feature, feature2: Feature):
    if feature1.chrom != feature2.chrom: return OVERLAP_TYPE.NONE
    return GetSegmentsOverlapType((feature1.start, feature1.end, feature1.strand),
                                  (feature2.start, feature2.end, feature2.strand))


# must be assumption that segments located on same chromosome
def GetSegmentsOverlapType(segment1, segment2):
    assert segment1 is not None and segment2 is not None
    l_a, r_a, strand_a = segment1
    l_b, r_b, strand_b = segment2

    if (l_a <= l_b <= r_a and l_a <= r_b <= r_a) or (l_b <= l_a <= r_b and l_b <= r_a <= r_b):
        return OVERLAP_TYPE.NESTED_PARALLEL if strand_a == strand_b else OVERLAP_TYPE.NESTED_ANTI_PARALLEL

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


def AreFeaturesSameFramed(fragment_a: Feature, fragment_b: Feature) -> bool:
    l = fragment_a.start + int(fragment_a.frame) if fragment_a.strand == '+' else fragment_a.end - int(
        fragment_a.frame)
    r = fragment_b.start + int(fragment_b.frame) if fragment_b.strand == '+' else fragment_b.end - int(
        fragment_b.frame)
    return l % 3 == r % 3


def GetOverlaps(fragments_a, fragments_b):
    non_ati_overlaps = []

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

            if fragment_a.strand == fragment_b.strand and AreFeaturesSameFramed(fragment_a, fragment_b):
                contains_in_phase_overlaps = True
                continue

            ov_length += overlap[1] - overlap[0] + 1
            non_ati_overlaps.append(overlap)
    non_ati_overlaps.sort(key=lambda x: x[0])
    return non_ati_overlaps, ov_length, contains_in_phase_overlaps

    # retrieves nucleotide composition of sequence
    # output: (C_count,G_count,A_count,T_count)


def GetMotifPwmFromFile(file_path):
    file = open(file_path, "r")
    lines = file.readlines()
    pwm = []
    for index in range(1, len(lines)):
        arr = lines[index].replace('\n', '').split('\t')
        pwm.append((float(arr[1]), float(arr[2]), float(arr[3]), float(arr[4])))
    file.close()
    return pwm


def GetMotifBestMatchScore(pwm, n):
    best_match_score = 0
    for i in range(n):
        score = 0
        arr = np.random.randint(low=0, high=4, size=(len(pwm),))
        for index in range(len(pwm)):
            score += pwm[index][arr[index]]
        best_match_score = max(best_match_score, score)
    return best_match_score
