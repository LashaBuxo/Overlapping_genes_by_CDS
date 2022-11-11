
Project was executed on Python 3.
In order to generate all output files 3 scripts must be executed:
```console
python generator_main_data.py
python generator_human_mouse.py
python generator_small.py
```
But generation of output files in `./generated_data` requires installing some libraries and downloading genome data.
### External required libraries:
    pip install plotly
    pip install pandas
    pip install gffutils
    pip install biopython
    pip install kaleido
### Required Files:

Some files must be downloaded to run all parts of the project and must be named and placed in directory according to this table or parameters declared in `worker_genome_values.py` file.

| **File Name**                               | **File Description**                  | **Download URL**                                                   | **Directory to save** | **Already located?**                 |
|---------------------------------------------|---------------------------------------|--------------------------------------------------------------------|-----------------------|--------------------------------------|
| Homo_sapiens.GRCh38.107.chr.gff3            | Annotation of Human Genome            | [link](https://ftp.ensembl.org/pub/release-108/gff3/homo_sapiens/) | used_data/genomes/    | <span style="color:red">No</span>    |
| Mus_musculus.GRCm39.107.chr.gff3.gz         | Annotation of Mouse Genome            | [link](https://ftp.ensembl.org/pub/release-108/gff3/mus_musculus/)   | used_data/genomes/    | <span style="color:red">No</span>    |
| Homo_sapiens.GRCh38.dna.primary_assembly.fa | DNA sequence of Human Genome          | [link](https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/) | used_data/genomes/    | <span style="color:red">No</span>    |
| Mus_musculus.GRCm39.dna.primary_assembly.fa | DNA sequence of Mouse Genome          | [link](https://ftp.ensembl.org/pub/release-108/fasta/mus_musculus/dna/) | used_data/genomes/    | <span style="color:red">No</span>    |
| human_appris_data.appris.txt                | Isoforms scores for genes | [link](https://appris.bioinfo.cnio.es/#/downloads)                    | used_data/APPRIS/     | <span style="color:green">Yes</span> |
| mouse_appris_data.appris.txt                | Isoforms scores for genes                                      | [link](https://appris.bioinfo.cnio.es/#/downloads)                    | used_data/APPRIS/     | <span style="color:green">Yes</span>    |
| orthologous_human_mouse.txt                 | BioMart othologoue genes              | [link](http://asia.ensembl.org/biomart/martview/)                     | used_data/APPRIS/     | <span style="color:green">Yes</span>    |