# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 17:18:37 2022

@author: dinah
"""

## This script takes the output of DefenseFinder and considers only systems
## found in finished genomes. These are then mapped to E. coli K-12 MG1655 to
## find their core flanking genes, i.e. where exactly they have been
## integrated in the genome.

import sys
from ast import literal_eval
from datetime import datetime
import pandas as pd

date = datetime.now().strftime('%Y-%m-%d')

Ecoli_allgenes_path = sys.argv[1]
all_islands_path = sys.argv[2]
DefenseFinder_path = sys.argv[3]
finished_genomes_path = sys.argv[4]
genome_name_path = sys.argv[5]
K12_symbols_path = sys.argv[6]


defence_name_dict = {
    "2TM-1TM-TIR" : "Rst_2TM_1TM_TIR",
    "3HP" : "Rst_3HP",
    "AbiD" : "Abi2",
    "AbiE" : "AbiEii",
    "ApeA" : "Gao_Ape",
    "AVAST_II" : "AVAST_II",
    "AVAST_III" : "AVAST_III",
    "AVAST_IV" : "AVAST_IV",
    "AVAST_type_II" : "AVAST_II",
    "AVAST_type_III" : "AVAST_III",
    "AVAST_type_IV" : "AVAST_IV",
    "AVAST_V" : "AVAST_V",
    "BREX_I" : "BREX_I",
    "brex_type_I" : "BREX_II",
    "CAS_Class1-Subtype-I-E" : "CRISPR_IE",
    "CAS_Class1-Subtype-I-F" : "CRISPR_IF",
    "CAS_Class1-Subtype-IV-A" : "CRISPR_IVA",
    "CAS_Class1-Type-I" : "CRISPR_I",
    "CAS_Cluster" : "CRISPR_cluster",
    "cas_type_I-E" : "CRISPR_IE",
    "cas_type_I-F1" : "CRISPR_IF",
    "cas_type_other" : "CRISPR_other",
    "CBASS_I" : "CBASS_I",
    "CBASS_II" : "CBASS_II",
    "CBASS_III" : "CBASS_III",
    "cbass_type_I" : "CBASS_I",
    "cbass_type_II" : "CBASS_II",
    "cbass_type_III" : "CBASS_III",
    "darTG" : "DarTG",
    "Dnd_ABCDE" : "Dnd",
    "Dnd_ABCDEFGH" : "Dnd",
    "DprA-PRTase" : "Rst_DprA-PRTase",
    "DRT_1" : "DRT_I",
    "DRT_2" : "DRT_II",
    "DRT_3" : "DRT_III",
    "DRT_4" : "DRT_IV",
    "DRT_5" : "DRT_V",
    "DRT_class_I" : "DRT_I",
    "DRT_class_II" : "DRT_II",
    "DRT_class_III" : "DRT_III",
    "DRT_other" : "DRT_other",
    "DRT_type_I" : "DRT_I",
    "DRT_type_III" : "DRT_III",
    "DRT_type_IV" : "DRT_IV",
    "DRT_type_V" : "DRT_V",
    "Druantia_I" : "Druantia_I",
    "Druantia_III" : "Druantia_III",
    "druantia_type_I" : "Druantia_I",
    "druantia_type_III" : "Druantia_III",
    "Dsr_I" : "Dsr_I",
    "dsr1" : "Dsr_I",
    "DUF4238" : "Rst_DUF4238",
    "dXTPase" : "dGTPase",
    "gabija" : "Gabija",
    "GAO_19" : "Gao_Her_SIR",
    "GAO_20" : "Gao_Her_DUF",
    "GAO_29" : "Gao_RL",
    "gop_beta_cll" : "Rst_gop_beta_cll",
    "hachiman_type_I" : "Hachiman_I",
    "Helicase-DUF2290" : "Rst_HelicaseDUF2290",
    "hhe" : "Gao_Hhe",
    "Hydrolase-TM" : "Rst_Hydrolase-Tm",
    "ietAS" : "Gao_Iet",
    "kiwa" : "Kiwa",
    "Lamassu_like" : "Lamassu_like",
    "lamassu_type_I" : "Lamassu_I",
    "lamassu_type_II" : "Lamassu_II",
    "Lamassu-Fam" : "Lamassu_like",
    "Mokosh_TypeI" : "Mokosh_I",
    "Mokosh_TypeII" : "Mokosh_II",
    "mza" : "Gao_Mza",
    "mza_other" : "Gao_Mza",
    "NMD_typeIV" : "NMD",
    "Paris" : "PARIS",
    "Paris_fused" : "PARIS_fused",
    "PARIS_I" : "PARIS_I",
    "PARIS_II" : "PARIS_II",
    "PARIS_II_merge" : "PARIS_merge",
    "PifA" : "Pif",
    "ppl" : "Gao_Ppl",
    "PT_DndABCDE" : "Dnd",
    "PT_DndFGH" : "Dnd",
    "PT_SspABCD" : "Ssp",
    "PT_SspE" : "Ssp",
    "PT_SspFGH" : "Ssp",
    "pycsar_effector" : "CBASS",
    "qatABCD" : "Gao_Qat",
    "radar_I" : "RADAR_I",
    "Retron_I_A" : "Retron_IA",
    "Retron_I_B" : "Retron_IB",
    "Retron_I_C" : "Retron_IC",
    "retron_I-A" : "Retron_IA",
    "retron_I-B" : "Retron_IB",
    "retron_I-C" : "Retron_IC",
    "Retron_II" : "Retron_II",
    "retron_II-A" : "Retron_IIA",
    "Retron_III" : "Retron_III",
    "retron_III-A" : "Retron_IIIA",
    "retron_IV" : "Retron_IV",
    "Retron_IV" : "Retron_IV",
    "retron_other" : "Retron_other",
    "Retron_VI" : "Retron_VI",
    "retron_VI" : "Retron_VI",
    "retron_XII" : "Retron_XII",
    "Retron_XII" : "Retron_XII",
    "RM_Type_I" : "RM_type_I",
    "RM_type_IIG" : "RM_type_IIG",
    "RM_Type_II" : "RM_type_II",
    "RM_Type_IIG" : "RM_type_IIG",
    "RM_Type_III" : "RM_type_III",
    "RM_Type_IV" : "RM_type_IV",
    "septu_type_I" : "Septu_I",
    "septu_type_II" : "Septu_II",
    "shedu" : "Shedu",
    "SspBCDE" : "Ssp",
    "TerY-P" : "Gao_TerY",
    "Thoeris_I" : "Thoeris_I",
    "Thoeris_II" : "Thoeris_II",
    "thoeris_type_I" : "Thoeris_I",
    "TIR-NLR" : "Rst_TIR-NLR",
    "Wadjet_I" : "Wadjet_I",
    "Wadjet_III" : "Wadjet_III",
    "zorya_type_I" : "Zorya_I",
    "Zorya_TypeI" : "Zorya_I",
    "zorya_type_II" : "Zorya_II",
    "Zorya_TypeII" : "Zorya_II"
}


# This dictionary has the format { genome : genome_name }
genome_name_dict = {}
with open(genome_name_path) as f:
    next(f)
    for line in f:
        fields = line.strip().split("\t")
        genome = fields[0]
        name = fields[9]
        genome_name_dict[genome] = name


# This dictionary has the format { geneID : gene_symbol}
gene_symbol_dict = {}
with open(K12_symbols_path) as f:
    next(f)
    for line in f:
        fields = line.strip().split("\t")

        if len(fields) == 7:
            geneID, locus_tag, gene_product, genomeID, genome_name, batch, gene_symbol = fields
            gene_symbol_dict[int(geneID)] = gene_symbol


# List of finished genomes
finished_genomes_list = []
with open(finished_genomes_path) as f:
    next(f)
    for line in f:
        genome = line.strip()
        finished_genomes_list.append(int(genome))


# Dictionary of original clusters and super-clusters
origcl_SC1_dict = {}

# Dictionary of genes and their original clusters
gene2cluster_dict = {}

# Dictionary of K-12 clusters and their genes
K12_cluster_genes_dict = {}

# Dictionary of genes and their contigs
gene_contig_dict = {}

# Dictionary of genes and their strands
gene_strand_dict = {}

# Dictionary of genes and their coordinates
gene_coord_dict = {}

# Dictionary of contigIDs and their contig names
contigID_contig_dict = {}

with open(Ecoli_allgenes_path) as f:
    for line in f:
        SC2, SC1, orig_cl, geneID, contig, genome, domain, deflated, from_coord, to_coord, strand, contig_ID, contig_pos, num_genes_contig = line.split()
        origcl_SC1_dict[int(orig_cl)] = int(SC1)
        gene2cluster_dict[int(geneID)] = int(orig_cl)

        gene_coord_dict[int(geneID)] = [from_coord, to_coord]
        contigID_contig_dict[contig_ID] = contig
        
        # If the genes are in K-12
        if contig_ID == "2687453745":
            if int(orig_cl) in K12_cluster_genes_dict:
                K12_cluster_genes_dict[int(orig_cl)].append(int(geneID))
            else:
                K12_cluster_genes_dict[int(orig_cl)] = [int(geneID)]

        gene_contig_dict[int(geneID)] = contig_ID
        gene_strand_dict[int(geneID)] = strand


# Dictionary of all genes found in islands at hotspots, as previously detected
island_gene_ID_dict = {}
with open(all_islands_path) as f:
    for line in f:
        fields = line.strip().split("\t")
        hotspot, K12_up_gene, K12_down_gene, K12_link, hotspot_found, occupied, genome, contig, contigID, contig_genes, island_up_gene, island_down_gene, start_coord, end_coord, size, num_genes, island_link, island_genes, island_clusters, island_SC1s, strand = fields
        island_genes = literal_eval(island_genes)

        ID = hotspot + "_" + genome

        for gene in island_genes:
            island_gene_ID_dict[gene] = ID


# DefenseFinder output
DefenseFinder_df = pd.read_csv(DefenseFinder_path, sep="\t")

# Dataframe for defence systems not previously detected
DF_not_islands_df = pd.DataFrame()

# Dictionary of each system and the genes in the system
system_genes_dict = {}

# Dictionary of all information found about each gene
all_info_dict = {}

# Iterating through DefenseFinder output, we check if the genes are found at
# the hotspots previously identified.
for index, row in DefenseFinder_df.iterrows():

    gene = int(row["gene_id"])
    sys_id = row["sys_id"]
    genome = row["replicon"]
    contigID = gene_contig_dict[gene]
    strand = gene_strand_dict[gene]

    if int(contigID) in finished_genomes_list:

        DF_not_islands_df = DF_not_islands_df.append(row)

        up_gene_found = False
        down_gene_found = False

        up_K12_gene = False
        down_K12_gene = False

        up_K12_gene_found = False
        down_K12_gene_found = False

        up_coord = False
        down_coord = False
        link = False

        up_gene = gene
        down_gene = gene

        # Tracking upstream of the defence system to try to find a core
        # gene. We do this by considering genes with a gene symbol as a
        # proxy for genes found in K-12. Later this will be checked
        # manually.
        while up_K12_gene_found == False and abs(gene - up_gene) < 300:
            up_gene -= 1

            try:
                up_orig_cluster = gene2cluster_dict[up_gene]
                up_K12_genes_in_cluster = K12_cluster_genes_dict[up_orig_cluster]
                
                # Only continue if there is only one gene in this cluster
                if len(up_K12_genes_in_cluster) == 1:

                    up_gene_strand = gene_strand_dict[up_gene]
                    up_K12_gene_strand = gene_strand_dict[up_K12_genes_in_cluster[0]]

                    up_5_gene_clusters = []
                    up_5_genes = []

                    K12_up_5_genes = []
                    K12_up_5_clusters = []

                    for i in range(10):
                        if up_gene_strand == up_K12_gene_strand:

                            try:
                                add_up_gene_cluster = str(gene2cluster_dict[up_gene - i])
                                up_5_genes.append(up_gene - i)
                                up_5_gene_clusters.append(add_up_gene_cluster)
                            except KeyError:
                                pass

                            try:
                                add_up_K12_gene_cluster = str(K12_cluster_genes_dict[up_orig_cluster][0] - i)
                                K12_up_5_genes.append(str(K12_cluster_genes_dict[up_orig_cluster][0] - i))
                                K12_up_5_clusters.append(str(gene2cluster_dict[K12_cluster_genes_dict[up_orig_cluster][0] - i]))
                            except KeyError:
                                pass

                        else:
                            try:
                                add_up_gene_cluster = str(gene2cluster_dict[up_gene - i])
                                up_5_genes.append(up_gene - i)
                                up_5_gene_clusters.append(add_up_gene_cluster)
                            except KeyError:
                                pass

                            try:
                                add_up_K12_gene_cluster = str(K12_cluster_genes_dict[up_orig_cluster][0] + i)
                                K12_up_5_genes.append(str(K12_cluster_genes_dict[up_orig_cluster][0] + i))
                                K12_up_5_clusters.append(str(gene2cluster_dict[K12_cluster_genes_dict[up_orig_cluster][0] + i]))
                            except KeyError:
                                pass
                    
                    # Checking that the genes upstream are in the same cluster
                    # as those in K-12
                    up_5_gene_clusters = up_5_gene_clusters[:5]
                    K12_up_5_clusters = K12_up_5_clusters[:5]

                    if (up_5_gene_clusters == K12_up_5_clusters or up_5_gene_clusters == K12_up_5_clusters[::-1]) and gene_contig_dict[up_gene] == contigID:

                        up_K12_gene_found = True
                        up_K12_gene = up_K12_genes_in_cluster[0]
                        up_gene_found = up_gene

            except KeyError:
                pass
            
            
        # Now doing the same to find the downstream gene
        while down_K12_gene_found == False and abs(gene - down_gene) < 300:
            down_gene += 1

            try:
                down_orig_cluster = gene2cluster_dict[down_gene]
                down_K12_genes_in_cluster = K12_cluster_genes_dict[down_orig_cluster]

                if len(down_K12_genes_in_cluster) == 1:

                    down_gene_strand = gene_strand_dict[down_gene]
                    down_K12_gene_strand = gene_strand_dict[down_K12_genes_in_cluster[0]]

                    down_5_gene_clusters = []
                    down_5_genes = []

                    K12_down_5_genes = []
                    K12_down_5_clusters = []

                    for i in range(10):

                        if down_gene_strand == down_K12_gene_strand:
                            try:
                                add_down_gene_cluster = str(gene2cluster_dict[down_gene + i])
                                down_5_genes.append(down_gene + i)
                                down_5_gene_clusters.append(add_down_gene_cluster)
                            except KeyError:
                                pass

                            try:
                                add_down_K12_gene_cluster = str(K12_cluster_genes_dict[down_orig_cluster][0] + i)
                                K12_down_5_genes.append(str(K12_cluster_genes_dict[down_orig_cluster][0] + i))
                                K12_down_5_clusters.append(str(gene2cluster_dict[K12_cluster_genes_dict[down_orig_cluster][0] + i]))
                            except KeyError:
                                pass

                        else:

                            try:
                                add_down_gene_cluster = str(gene2cluster_dict[down_gene + i])
                                down_5_genes.append(down_gene + i)
                                down_5_gene_clusters.append(add_down_gene_cluster)
                            except KeyError:
                                pass

                            try:
                                add_down_K12_gene_cluster = str(K12_cluster_genes_dict[down_orig_cluster][0] - i)
                                K12_down_5_genes.append(str(K12_cluster_genes_dict[down_orig_cluster][0] - i))
                                K12_down_5_clusters.append(str(gene2cluster_dict[K12_cluster_genes_dict[down_orig_cluster][0] - i]))
                            except KeyError:
                                pass
                    
                    # Checking that the genes downstream are in the same
                    # cluster as those in K-12
                    down_5_gene_clusters = down_5_gene_clusters[:5]
                    K12_down_5_clusters = K12_down_5_clusters[:5]

                    if (down_5_gene_clusters == K12_down_5_clusters or down_5_gene_clusters == K12_down_5_clusters[::-1]) and gene_contig_dict[down_gene] == contigID:

                        down_K12_gene_found = True
                        down_K12_gene = down_K12_genes_in_cluster[0]
                        down_gene_found = down_gene

            except KeyError:
                pass

        try:
            up_gene_symbol = gene_symbol_dict[up_K12_gene]
        except KeyError:
            up_gene_symbol = False

        try:
            down_gene_symbol = gene_symbol_dict[down_K12_gene]
        except KeyError:
            down_gene_symbol = False

        # Have we found both core flanking genes?
        if up_gene_found and down_gene_found:
            try:
                up_coord = gene_coord_dict[up_gene_found][1]
                down_coord = gene_coord_dict[down_gene_found][0]
                link = "https://img.jgi.doe.gov/cgi-bin/mer/main.cgi?section=ScaffoldGraph&page=scaffoldGraph&scaffold_oid={}&start_coord={}&end_coord={}".format(contigID, up_coord, down_coord)
            except KeyError:
                pass

            if sys_id in system_genes_dict:
                system_genes_dict[sys_id].append(gene)
            else:
                system_genes_dict[sys_id] = [gene]

        else:
            print("COULD NOT BE MAPPED", gene)
        
        
        # Saving all information
        all_info_dict[gene] = [gene, up_gene_found, up_K12_gene, up_gene_symbol, down_gene_found, down_K12_gene, down_gene_symbol, contigID, genome, up_coord, down_coord, link]
        print(gene, up_gene_found, up_K12_gene, up_gene_symbol, down_gene_found, down_K12_gene, down_gene_symbol, contigID, genome, up_coord, down_coord, link)

output_path = r"\\data.wexac.weizmann.ac.il\sorek\dinah\2022_Ecoli-defence-islands\Work\2022-08Aug_revisions\Oct\{}_DF-outside-islands_formatted2.output".format(date)

output_file = open(output_path, "a+")
output_file.write("\t".join(["Representative island genome name", "Representative island genome IMG ID", "Representative island contig IMG accession", "Representative island contig IMG ID", "Representative island upstream flanking gene (IMG ID)", "Representative island downstream flanking gene (IMG ID)", "Representative island start coordinate", "Representative island end coordinate", "Island size (bp)", "Number of genes in representative island", "IMG Browser link to representative island", "IMG gene IDs in representative island", "Mobile genetic element?", "DF_defence_genes", "DF_defence_systems", "DF_defence_simple_systems", "Found in hotspot?"])+"\n")

# Converting this into the desired format with all additional information
# about the island
for key, value in system_genes_dict.items():

    first_gene = value[0]
    all_info = all_info_dict[first_gene]
    gene, up_gene_found, up_K12_gene, up_gene_symbol, down_gene_found, down_K12_gene, down_gene_symbol, contigID, genome, up_coord, down_coord, link = all_info
    genome_name = genome_name_dict[str(genome)]
    genome = genome
    contig = contigID_contig_dict[contigID]
    contigID = contigID
    up_gene = up_gene_found
    down_gene = down_gene_found
    start_coord = up_coord
    end_coord = down_coord
    size = int(end_coord) - int(start_coord) + 1
    num_genes = down_gene - up_gene - 1
    link = link
    gene_list = list(range(up_gene, down_gene + 1))
    MGE = 0
    defence_genes = system_genes_dict[key]
    defence_system = key

    try:
        simple_defence_system = defence_name_dict["_".join(key.split("_")[1:-1])]
    except KeyError:
        simple_defence_system = "_".join(key.split("_")[1:-1])

    # Is the island found at a hotspot?
    try:
        hotspot = island_gene_ID_dict[first_gene]

        up_gene = up_gene_found
        down_gene = down_gene_found
        start_coord = up_coord
        end_coord = down_coord
        size = int(end_coord) - int(start_coord) + 1
        num_genes = down_gene - up_gene + 1
        link = link
        gene_list = list(range(up_gene, down_gene + 1))
        MGE = 0

    except KeyError:
        hotspot = 0

    output_file.write("\t".join([genome_name, str(genome), contig, contigID, str(up_gene), str(down_gene), start_coord, end_coord, str(size), str(num_genes), link, str(gene_list), str(MGE), str(defence_genes), defence_system, simple_defence_system, str(hotspot)]) + "\n")

output_file.close()

DF_not_islands_df["gene"] = pd.to_numeric(DF_not_islands_df["gene"])
DF_not_islands_df["all_info"] = DF_not_islands_df["gene"].map(all_info_dict)
DF_not_islands_df[["gene", "up_gene_found", "up_K12_gene", "up_gene_symbol", " down_gene_found", "down_K12_gene", "down_gene_symbol", "contigID", "genome", "up_coord", "down_coord", "link"]] = pd.DataFrame(DF_not_islands_df.all_info.tolist(), index=DF_not_islands_df.index)
DF_not_islands_df = DF_not_islands_df.reset_index(drop=True)
DF_not_islands_df.to_excel(r"\\data.wexac.weizmann.ac.il\sorek\dinah\2022_Ecoli-defence-islands\Work\2022-08Aug_revisions\{}_DefenseFinder-in-finished-genomes.xlsx".format(date))
