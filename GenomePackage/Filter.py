from GenomePackage.Helper import *


def __check_gene_for_Ensembl_filters(gene: Feature, chroms):
    # https: // www.biostars.org / p / 5304 /  # 9521037
    # ignore genes with miscellaneous chromosome/scaffold names
    if not chroms.__contains__(gene.chrom): return False, "located_on_miscellaneous_scaffold"

    # gene_symbol = GetGeneSymbol(gene)
    # gene_accession = GetGeneAccession(gene)
    description = GetGeneDescription(gene)

    # ignore genes without description or gene_symbol
    # if description == "no_desc": return False, "no_desc"
    # if gene_symbol == "no_sym": return False, "no_symbol"
    # if gene_accession == "no_acc": return False, "no_accession"

    # ignore genes if they are pseudogene, novel or predicted, readthrough
    if description.__contains__('readthrough'):
        return False, "is_readthrough"
    if description.__contains__('pseudogene') or description.__contains__('pseudogene'):
        return False, "is_pseudogene"
    if description.__contains__('novel') or description.__contains__('Novel'):
        return False, "is_novel"
    if description.__contains__('predicted') or description.__contains__('Predicted'):
        return False, "is_predicted"

    # ignore genes with duplicate Names or accessions
    # if feature_syms[gene_symbol] > 1: return False, "name_duplicated"
    # ignore genes if gene with same NCBI accession already imported
    # if feature_accs[gene_accession] > 1: return False, "acc_duplicated"
    # ignore genes if gene with same NCBI accession already imported
    # if feature_ids[gene.id] > 1: return False, "id_duplicated"

    return True, "passed_filter"


def filter_by_ensembl_attributes(genes_on_chr):
    filtered_genes_on_chr = {}
    ignored_genes_by_types = {}

    chroms = genes_on_chr.keys()

    # feature_ids = {}
    # feature_syms = {}
    # feature_accs = {}
    # for chrom in chroms:
    #     for gene in genes_on_chr[chrom]:
    #         feature_ids[gene.id] = feature_ids.get(gene.id, 0) + 1
            # feature_syms[GetGeneSymbol(gene)] = feature_syms.get(GetGeneSymbol(gene), 0) + 1
            # feature_accs[GetGeneAccession(gene)] = feature_accs.get(GetGeneAccession(gene), 0) + 1

    for chrom in chroms:
        new_genes_on_chr = []
        for gene in genes_on_chr[chrom]:
            is_valid, status = __check_gene_for_Ensembl_filters(gene, chroms)
            if is_valid:
                new_genes_on_chr.append(gene)
            else:
                if not ignored_genes_by_types.__contains__(status):
                    ignored_genes_by_types[status] = 0
                ignored_genes_by_types[status] += 1
        filtered_genes_on_chr[chrom] = new_genes_on_chr

    return filtered_genes_on_chr, ignored_genes_by_types
