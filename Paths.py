from GenomePackage.Enums import SPECIES

GTEx_EXPRESSION_DATA_PATH = "./CommonRequiredData/Human Expressions/short_required_expressions.txt"

ANNOTATIONS_PATH = {
    SPECIES.Homo_sapiens: "./CommonRequiredData/Annotations/Homo_sapiens.GRCh38.109.chr.gff3",
    SPECIES.Pan_troglodytes: "./CommonRequiredData/Annotations/Pan_troglodytes.Pan_tro_3.0.109.chr.gff3",
    SPECIES.Mus_musculus: "./CommonRequiredData/Annotations/Mus_musculus.GRCm39.109.chr.gff3",
    SPECIES.Sus_scrofa: "./CommonRequiredData/Annotations/Sus_scrofa.Sscrofa11.1.109.chr.gff3",
    SPECIES.Canis_lupus_familiaris: "./CommonRequiredData/Annotations/Canis_lupus_familiaris.ROS_Cfam_1.0.109.chr.gff3",
    SPECIES.Monodelphis_domestica: "./CommonRequiredData/Annotations/Monodelphis_domestica.ASM229v1.109.chr.gff3",
    SPECIES.Gallus_gallus: "./CommonRequiredData/Annotations/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.109.chr.gff3",
    SPECIES.Danio_rerio: "./CommonRequiredData/Annotations/Danio_rerio.GRCz11.109.chr.gff3",
    SPECIES.Drosophila_melanogaster: "./CommonRequiredData/Annotations/Drosophila_melanogaster.BDGP6.32.109.chr.gff3",
}

# DNA sequences (downloaded from ensembl)
SEQUENCES_PATH = {
    SPECIES.Homo_sapiens: "./CommonRequiredData/Sequences/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa",
    SPECIES.Pan_troglodytes: "./CommonRequiredData/Sequences/Pan_troglodytes.Pan_tro_3.0.dna_sm.toplevel.fa",
    SPECIES.Mus_musculus: "./CommonRequiredData/Sequences/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa",
    SPECIES.Sus_scrofa: "./CommonRequiredData/Sequences/Sus_scrofa.Sscrofa11.1.dna_sm.toplevel.fa",
    SPECIES.Canis_lupus_familiaris: "./CommonRequiredData/Sequences/Canis_lupus_familiaris.ROS_Cfam_1.0.dna_sm.toplevel.fa",
    SPECIES.Monodelphis_domestica: "./CommonRequiredData/Sequences/Monodelphis_domestica.ASM229v1.dna_sm.toplevel.fa",
    SPECIES.Gallus_gallus: "./CommonRequiredData/Sequences/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna_sm.toplevel.fa",
    SPECIES.Danio_rerio: "./CommonRequiredData/Sequences/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa",
    SPECIES.Drosophila_melanogaster: "./CommonRequiredData/Sequences/Drosophila_melanogaster.BDGP6.32.dna_sm.toplevel.fa",
}
