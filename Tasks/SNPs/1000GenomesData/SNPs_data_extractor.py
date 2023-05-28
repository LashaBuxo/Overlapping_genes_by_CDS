file_path = f"../../InputData/1000GenomesData/00-common_all.vcf"
out_file_path = f"../../InputData/1000GenomesData/short_required_SNPs.txt"

out_file = open(out_file_path, "w")
last_ind = -1
with open(file_path, 'r') as file:
    index = 0
    while True:
        line = file.readline()
        if not line:
            break
        if line.startswith("#"): continue
        # if line.replace('\n', '').endswith("SNP"):
        arr = line.split("\t")
        chrom = arr[0]
        ind = arr[1]
        from_nt = arr[3]
        to_nt = arr[4]
        if ind == last_ind: continue
        last_ind = ind
        out_file.write(f"{chrom}\t{ind}\t{from_nt}\t{to_nt}\t\n")
        index += 1
        if index % 100000 == 0:
            print(f"{chrom}\t{ind}\t{from_nt}\t{to_nt}\t\n")
        # process the line here
out_file.close()
