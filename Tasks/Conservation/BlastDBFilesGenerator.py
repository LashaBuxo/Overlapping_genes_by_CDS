# DNA sequences (downloaded from ensembl)
import os

from GenomePackage.Enums import SPECIES

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

for species in SEQUENCES_PATH.keys():
    print(f"Building Blast database for {species.short_name()} DNA!")
    fasta_file = SEQUENCES_PATH[species]
    db_path = f".\Tasks\Conservation\BlastData\db\{species.tax_id()}_db"
    db_title = f"{species.tax_id()}"
    cmd_query = f"makeblastdb -in {fasta_file}  -parse_seqids -blastdb_version 5 -title {db_title} -dbtype nucl -out {db_path}"
    os.system(cmd_query)
print("Success!!!")
