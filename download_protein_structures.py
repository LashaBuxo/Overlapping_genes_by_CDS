import requests


def search_pdb_by_sequence(sequence, evalue_cutoff=0.1, identity_cutoff=0, force_predicted=False):
    api_url = "https://search.rcsb.org/rcsbsearch/v2/query?json="
    query_data = {
        "query": {
            "type": "terminal",
            "service": "sequence",
            "parameters": {
                "target": "pdb_protein_sequence",
                "value": sequence
            }
        },
        "request_options": {
            "return_all_hits": False,
            "results_content_type": ["computational", "experimental"] if not force_predicted else ["computational"]
        },
        "return_type": "entry"
    }

    response = requests.post(api_url, json=query_data)
    response.raise_for_status()

    if len(response.text) == 0:
        return []
    pdb_ids = [[item['identifier'], item['score']] for item in response.json()['result_set']]
    return pdb_ids


def search_pdb_and_download_structure(protein_sequence):
    force_predicted = False
    # For MRPL43 gene crystal data is absent for this segment
    if protein_sequence == "DTGLRLSAVAPQILLPGWPDP":
        force_predicted = True
    # For DUOXA1 gene crystal data is absent for this segment
    if protein_sequence == "PWRFGPEGCEERWAEHTGDSPRPLRGRGTGR":
        force_predicted = True
    # For RANBP1 gene crystal data is absent for this segment
    if protein_sequence == "MGSALGRARRTLSGRPFQRAPCKTRRALSLSAALRNVTKAQGGCPKSLVLWGCRPKRPRKRRTSLKLAWR":
        force_predicted = True
    entry_ids_and_scores = search_pdb_by_sequence(protein_sequence, force_predicted=force_predicted)
    if len(entry_ids_and_scores) == 0:
        return "", "", "", "", ""

    # Choose best entry
    entry_id, score = entry_ids_and_scores[0]
    base_url = "https://data.rcsb.org/rest/v1/core"

    entry_url = f"{base_url}/entry/{entry_id}"
    entry_response = requests.get(entry_url)
    entry_data = entry_response.json()

    method = entry_data.get("rcsb_entry_info", {}).get("structure_determination_methodology")
    entry_title = entry_data.get("struct", {}).get("title")
    citation_title = entry_data.get("citation", {})[0].get('title')
    citation_doi = entry_data.get("citation", {})[0].get('pdbx_database_id_doi')
    # for data in entry_data.get("citation", {}):
    #     print(f"\tCITATION: '{data.get('title')}' ({data.get('pdbx_database_id_doi')})")

    structure_file_url = ""
    if method == "experimental":
        structure_file_url = f"https://files.rcsb.org/download/{entry_id}.cif"
    elif method == "computational":
        structure_file_url = entry_data.get("rcsb_comp_model_provenance", {}).get('source_url')

    return method, entry_title, f"{citation_title} ({citation_doi})", score, structure_file_url


out_file = open("./Tasks/Disorder/supplementary_detailed.txt", "w")

in_file = open("./Tasks/Disorder/supplementary_disorder.txt", "r")
lines = in_file.readlines()
for i in range(1, len(lines)):
    arr = lines[i].replace('\n', '').split('\t')
    seq1 = arr[2]
    seq2 = arr[3]
    out_file.write(f"{arr[0]}\t{arr[1]}\t{seq1}\t{seq2}\t")
    data = [search_pdb_and_download_structure(seq1), search_pdb_and_download_structure(seq2)]

    if data[0][0] != "" and data[1][0] != "":
        out_file.write("Yes\t")
    else:
        out_file.write("No\t")

    for method, entry_title, entry_citation, entry_similarity, entry_url in data:
        out_file.write(f"{method}\t")
        out_file.write(f"{entry_title}\t")
        out_file.write(f"{entry_citation}\t")
        out_file.write(f"{entry_similarity}\t")
        out_file.write(f"{entry_url}\t")
    out_file.write(f"\n")
    print(f"Detailed structures search finished for overlap: {arr[0]}-{arr[1]}")

out_file.close()
in_file.close()
