import random

import numpy as np
import scipy
from gffutils import Feature

from GenomePackage.Enums import ANNOTATION_LOAD
from GenomePackage.GenomeData import GenomeData
from GenomePackage.OverlapsFinder import GetCDSOverlapsFromGenome
from Paths import *

request_species = SPECIES.Homo_sapiens
genome = GenomeData(SPECIES.Homo_sapiens,
                    ANNOTATION_LOAD.GENES_AND_TRANSCRIPTS_AND_CDS,
                    )

records, ATI_overlaps = GetCDSOverlapsFromGenome(genome, exclude_incomplete_utr_transcripts=True,
                                                 exclude_zero_expressed_transrcipts=True)
records.sort(key=lambda x: -x.overlapped_length)

file = open("Tasks/Expression/NR1D1_THRA_expression.txt", "w")
for record in records:
    if record.sym1 != "THRA" and record.sym2 != "NR1D1": continue
    for tissue in genome.get_expression_data_tissues():

        val1 = genome.get_transcript_expression_in_tissue(genome.get_feature_by_id("transcript:ENST00000450525"),tissue)
        val1_0 = genome.get_transcript_expression_in_tissue(record.transcript1, tissue)
        val2 = genome.get_gene_expression_in_tissue(record.gene2, tissue)
        file.write(f"{tissue}\t{val1}\t{val1_0}\t{val2}\n")

    for tissue in genome.get_expression_data_tissues():
        samples1 = genome.get_transcript_expressions_in_tissue_by_samples(record.transcript1, tissue)
        samples2 = genome.get_transcript_expressions_in_tissue_by_samples(record.transcript2, tissue)
        assert len(samples1) == len(samples2)
        file.write(f"{tissue}")

        for index in range(len(samples1)):
            file.write(f"\t{samples1[index]}\t{samples2[index]}")
        file.write(f"\n")
file.close()
