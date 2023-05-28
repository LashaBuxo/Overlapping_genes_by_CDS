# install https://bitsilla.com/wiki/doku.php?id=howto:pymol_install_on_windows
from os.path import exists

import pymol
from pymol import cmd
import requests


def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i + 2], 16) / 255.0 for i in (0, 2, 4))


def get_residues(obj, chain):
    stored = {'residues': set()}
    cmd.iterate(f"{obj} and chain {chain} and name CA", "stored['residues'].add(resi)", space=locals())
    return sorted(list(map(int, stored['residues'])))


def highlight_protein_part(protein_name, sequence, highlight_id):
    chains = cmd.get_chains(protein_name)
    if (protein_name == "RANBP1"):
        print(f"{chains}")
    for chain in chains:
        fasta = cmd.get_fastastr(f"{protein_name} and chain {chain}")
        seq = "".join(line.strip() for line in fasta.split("\n") if not line.startswith(">"))

        if seq.__contains__(sequence):
            start_index = seq.index(sequence)
            residues = get_residues(protein_name, chain)

            highlight_range = f"resi {residues[start_index]}-{residues[start_index] + len(ov_seq1)}"
            cmd.select(highlight_id, f"{protein_name} and chain {chain} and {highlight_range}")

            if (protein_name == "RANBP1"):
                print(start_index)
                print(residues)
                print(f"resi {residues[start_index]}-{residues[start_index] + len(ov_seq1)}")


pymol.finish_launching()

cmd.set("ray_trace_gain", 0.001)
cmd.set("ray_trace_disco_factor", 0)
cmd.set("antialias", 2)
cmd.set("ray_trace_mode", 0)
cmd.set("specular", 0.2)
cmd.set("ambient", 0.5)

cmd.bg_color("white")
cmd.set_color("highlight_color1", (181, 227, 245))
cmd.set_color("highlight_color2", (230, 241, 211))
cmd.set_color("protein_color", hex_to_rgb("#999999ff"))

ind = 0
file = open("./Tasks/Disorder/supplementary_detailed.txt", 'r')
for line in file.readlines():
    arr = line.replace('\n', '').split('\t')
    gene1, gene2 = arr[0], arr[1]
    link1, link2 = arr[9], arr[14]

    ov_seq1, ov_seq2 = arr[2], arr[3]

    if link1 != "" and link2 != "":
        if not exists(f"./Tasks/Disorder/Structures/{gene1}.cif"):
            response = requests.get(link1)
            pdb_data = response.text

            file1 = open(f"./Tasks/Disorder/Structures/{gene1}.cif", "w")
            file1.write(pdb_data)
            file1.close()

            response = requests.get(link2)
            pdb_data = response.text

            file2 = open(f"./Tasks/Disorder/Structures/{gene2}.cif", "w")
            file2.write(pdb_data)
            file2.close()

        cmd.load(f"./Tasks/Disorder/Structures/{gene1}.cif", gene1)
        cmd.load(f"./Tasks/Disorder/Structures/{gene2}.cif", gene2)
        cmd.color("protein_color", gene1)
        cmd.color("protein_color", gene2)
        highlight_protein_part(gene1, ov_seq1, f"highlight_{gene1}")
        highlight_protein_part(gene2, ov_seq2, f"highlight_{gene2}")
        cmd.color("highlight_color1", f"highlight_{gene1}")
        cmd.color("highlight_color2", f"highlight_{gene2}")

# Show as cartoon representation
cmd.show("cartoon")
cmd.zoom("all")
