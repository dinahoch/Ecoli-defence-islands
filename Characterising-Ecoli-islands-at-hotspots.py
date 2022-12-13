# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 15:45:49 2022

@author: dinah
"""

import numpy as np
import os
import pandas as pd
import re
from datetime import datetime
from ast import literal_eval
from Bio import SeqIO, Entrez

date = datetime.now().strftime('%Y-%m-%d')


## Detecting genes of mobile genetic element found in each sequence

## 1. Prophages: Nucleotide sequences were submitted to Phaster
## 2. Phage satellites: Protein sequences were submitted to SatelliteFinder
## 3. Integrated conjugative elements and integrative mobilisable elements:
## Nucleotide sequences were submitted to ConjScan

## and genes of defence systems

## 1. DefenseFinder: Amino acid sequences were submitted to DefenseFinder
## 2. PADLOC: Amino acid sequences and gff3 files were submitted to PADLOC


S2_path = sys.argv[1]
Phaster_path = sys.argv[2]
DefenseFinder_path = sys.argv[3]
PADLOC_path = sys.argv[4]
Ecoli_allgenes_path = sys.argv[5]
output_path = sys.argv[6]


# Creating dictionaries of E. coli genes and their coordinates
gene_coords_dict = {}
with open(Ecoli_allgenes_path) as f:
    for line in f:
        fields = line.strip().split("\t")
        gene = fields[3]
        start_coord = fields[8]
        end_coord = fields[9]
        contig = fields[11]
        gene_coords_dict[gene] = start_coord + "_" + end_coord


# Dataframe of information from Supplementary Table S2
S2_df = pd.read_csv(S2_path, sep="\t")
S2_df["ID"] = S2_df["Hotspot"].astype(str) + "_" + S2_df["Representative island genome IMG ID"].astype(str)


# This is a dictionary with the format { ID : (start_coord, end_coord) }
island_coords_dict = dict(zip(S2_df["ID"], (S2_df["Representative island start coordinate"].astype(int), S2_df["Representative island end coordinate"].astype(int))))
island_coords_dict = {k: (df["Representative island start coordinate"].item(), df["Representative island end coordinate"].item()) for k, df in S2_df.groupby(["ID"])}


# Making defence system names consistent between DefenseFinder and PADLOC
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


# Dictionary of the simple versions of each system name
simple_defence_dict = {
    "AVAST_II" : "AVAST",
    "AVAST_III" : "AVAST",
    "AVAST_IV" : "AVAST",
    "AVAST_V" : "AVAST",
    "BREX_I" : "BREX",
    "brex_type_I" : "BREX",
    "CRISPR_IE" : "CRISPR",
    "CRISPR_IF" : "CRISPR",
    "CRISPR_IVA" : "CRISPR",
    "CRISPR_other" : "CRISPR",
    "CBASS_I" : "CBASS",
    "CBASS_II" : "CBASS",
    "CBASS_III" : "CBASS",
    "DRT_I" : "DRT",
    "DRT_II" : "DRT",
    "DRT_III" : "DRT",
    "DRT_IV" : "DRT",
    "DRT_V" : "DRT",
    "DRT_other" : "DRT",
    "Druantia_I" : "Druantia",
    "Druantia_III" : "Druantia",
    "Dsr_I" : "Dsr",
    "Hachiman_I" : "Hachiman",
    "Lamassu_like" : "Lamassu",
    "Lamassu_I" : "Lamassu",
    "Lamassu_II" : "Lamassu",
    "Mokosh_I" : "Mokosh",
    "Mokosh_II" : "Mokosh",
    "Paris_fused" : "PARIS",
    "PARIS_I" : "PARIS",
    "PARIS_II" : "PARIS",
    "PARIS_merge" : "PARIS",
    "RADAR_I" : "RADAR",
    "Retron_IA" : "Retron",
    "Retron_IB" : "Retron",
    "Retron_IC" : "Retron",
    "Retron_II" : "Retron",
    "Retron_IIA" : "Retron",
    "Retron_III" : "Retron",
    "Retron_IIIA" : "Retron",
    "Retron_IV" : "Retron",
    "Retron_other" : "Retron",
    "Retron_VI" : "Retron",
    "Retron_XII" : "Retron",
    "RM_type_IIG" : "RM_type_II",
    "Septu_I" : "Septu",
    "Septu_II" : "Septu",
    "Thoeris_I" : "Thoeris",
    "Thoeris_II" : "Thoeris",
    "Wadjet_I" : "Wadjet",
    "Wadjet_III" : "Wadjet",
    "Zorya_I" : "Zorya",
    "Zorya_II" : "Zorya",
    }


# Adding phage satellite predictions from SatelliteFinder. In my dataset I
# only see P4-like and PICI satellites.

# Also adding information from CONJScan on ICEs and IMEs.

P4_dict = {}
P4_full_type_dict = {}
P4_simple_type_dict = {}
P4_full_set_dict = {}

PICI_dict = {}
PICI_full_type_dict = {}
PICI_simple_type_dict = {}
PICI_full_set_dict = {}

P4_coord_dict = {}

ConjScan_dict = {}

for index, row in S2_df.iterrows():

    hotspot = row["Hotspot"]
    ID = row["ID"]

    P4_path = r"\\data.wexac.weizmann.ac.il\sorek\dinah\2022_Ecoli-defence-islands\Work\2022-08Aug_revisions\MGEs\Prophages\SatelliteFinder\P4\{}\{}.faa.csv".format(hotspot, ID)

    if os.path.exists(P4_path):
        with open(P4_path) as f:
            next(f)
            for line in f:
                fields = line.split(",")
                full_type = fields[4]
                simple_type = fields[5]

                fields_split_quote = line.split('"')

                full_set = literal_eval(fields_split_quote[1])

                gene_list = [full_set[i][0] for i in range(len(full_set))]
                gene_list = sorted(gene_list)
                first_gene = gene_list[0]
                last_gene = gene_list[-1]

                start_coord = gene_coords_dict[str(first_gene)].split("_")[0]
                end_coord = gene_coords_dict[str(last_gene)].split("_")[1]

                P4_coord_dict[ID] = range(int(start_coord), int(end_coord))

                P4_dict[ID] = [full_type, simple_type, full_set]
                P4_full_type_dict[ID] = full_type
                P4_simple_type_dict[ID] = simple_type
                P4_full_set_dict[ID] = full_set


    PICI_path = r"\\data.wexac.weizmann.ac.il\sorek\dinah\2022_Ecoli-defence-islands\Work\2022-08Aug_revisions\MGEs\Prophages\SatelliteFinder\PICI\{}\{}.faa.csv".format(hotspot, ID)

    if os.path.exists(PICI_path):
        with open(PICI_path) as f:
            next(f)
            for line in f:
                fields = line.split(",")
                full_type = fields[4]
                simple_type = fields[5]

                fields_split_quote = line.split('"')

                full_set = literal_eval(fields_split_quote[1])
                gene_list = [full_set[i][0] for i in range(len(full_set))]
                gene_list = sorted(gene_list)
                first_gene = gene_list[0]
                last_gene = gene_list[-1]

                PICI_dict[first_gene] = [full_type, simple_type, full_set]
                PICI_full_type_dict[ID] = full_type
                PICI_simple_type_dict[ID] = simple_type
                PICI_full_set_dict[ID] = full_set


    ConjScan_path = r"\\data.wexac.weizmann.ac.il\sorek\dinah\2022_Ecoli-defence-islands\Work\2022-08Aug_revisions\MGEs\ConjScan\macsyfinder\output\{}.output\all_best_solutions.tsv".format(ID)

    if os.path.exists(ConjScan_path):
        with open(ConjScan_path) as f:
            for line in f:
                if line != "\n" and line[0] != "#" and line[0] != "s":
                    fields = line.strip().split("\t")
                    gene = fields[2]
                    hit = fields[3]
                    system = fields[5].split("/")[-1]

                    ConjScan_dict[ID] = [gene, hit, system]


# Adding prophage annotation from Phaster for intact, questionable and
# incomplete prophages. If the Phaster prediction overlaps a phage satellite
# prediction, take the satellite prediction, which is more accurate.

Phaster_df = pd.read_csv(Phaster_path, header=None, sep=r"\s+")
Phaster_df.columns = ["region", "region_length", "completeness(score)", "specific_keyword", "region_position", "trna_num", "total_protein_num", "phage_hit_protein_num", "hypothetical_protein_num", "phage+hypo_protein_percentage", "bacterial_protein_num", "att_site_showup", "phage_species_num", "most_common_phage_name(hit_genes_count)", "first_most_common_phage_num", "first_most_common_phage_percentage", "gc_percentage"]

intact_island_phage_dict = {}
quest_island_phage_dict = {}
incom_island_phage_dict = {}

for index, row in Phaster_df.iterrows():

    island = row["region_position"].split(":")[0]
    position = row["region_position"].split(":")[1]
    phage_start = position.split("-")[0]
    phage_end = position.split("-")[1]
    phage_coords = (int(phage_start), int(phage_end))
    try:
        shifted_phage_coords = tuple(np.add(island_coords_dict[island], phage_coords))
    except KeyError:
        pass

    if island in P4_coord_dict:

        overlap = range(max(P4_coord_dict[island][0], shifted_phage_coords[0]), min(P4_coord_dict[island][-1], shifted_phage_coords[-1]))

        if P4_coord_dict[island][0] >= shifted_phage_coords[0] - 3000 and P4_coord_dict[island][1] <= shifted_phage_coords[1] + 3000:
            print("Island ID:", island, "P4 predicted in Phaster phage")

        else:
            print("Island ID:", island, ". P4 coordinates:", P4_coord_dict[island], ". Phaster coordinates:", shifted_phage_coords, ". Overlap:", overlap)

    else:
        score = row["completeness(score)"].split("(")[0]

        # For one island, we split the phages into a list
        phages = row["most_common_phage_name(hit_genes_count)"]
        phages = phages.split(",")

        intact_tax_num_dict = {}
        quest_tax_num_dict = {}
        incom_tax_num_dict = {}

        tax_genus_pos = False
        tax_order_pos = False
        tax_species_pos = False

        # There are multiple hits for each island. We only consider hits with more
        # than one gene mapped.
        for phage in phages:

            print(phage)

            NC = "_".join(phage.split("_")[-2:]).split("(")[0]
            num = int(re.findall('\((.*?)\)', phage)[-1])

            if num > 1:
                protein_id = NC
                protein_record_raw = Entrez.efetch(id=protein_id, db="nuccore", rettype="gb", retmode="text")
                protein_record = SeqIO.read(protein_record_raw, "genbank")
                taxonomy = protein_record.annotations["taxonomy"]
                print(taxonomy)

                taxonomy = [i for i in taxonomy if " " not in i]

                try:
                    tax_genus_pos = next(i for i in reversed(range(len(taxonomy))) if "virus" in taxonomy[i])
                except StopIteration:
                    pass
                try:
                    tax_order_pos = next(i for i in reversed(range(len(taxonomy))) if "virales" in taxonomy[i])
                except StopIteration:
                    pass
                try:
                    tax_class_pos = next(i for i in reversed(range(len(taxonomy))) if "viricetes" in taxonomy[i])
                except StopIteration:
                    pass
                try:
                    tax_species_pos = next(i for i in reversed(range(len(taxonomy))) if " " in taxonomy[i])
                except StopIteration:
                    pass

                if tax_species_pos:
                    if tax_genus_pos and tax_order_pos:
                        tax = "_".join(taxonomy[tax_order_pos:tax_species_pos])
                    elif tax_genus_pos and tax_class_pos:
                        tax = "_".join(taxonomy[tax_class_pos:tax_species_pos])
                    elif tax_class_pos:
                        tax = "_".join(taxonomy[tax_class_pos:tax_species_pos])
                    else:
                        print("NO", taxonomy)
                else:
                    if tax_genus_pos and tax_order_pos:
                        tax = "_".join(taxonomy[tax_order_pos:tax_genus_pos+1])
                    elif tax_genus_pos and tax_class_pos:
                        tax = "_".join(taxonomy[tax_class_pos:tax_genus_pos+1])
                    elif tax_class_pos:
                        tax = "_".join(taxonomy[tax_class_pos:])
                    else:
                        print("NO", taxonomy)


                if score == "intact":
                    if tax not in intact_tax_num_dict:
                        intact_tax_num_dict[tax] = num
                    else:
                        intact_tax_num_dict[tax] += num

                if score == "questionable":
                    if tax not in quest_tax_num_dict:
                        quest_tax_num_dict[tax] = num
                    else:
                        quest_tax_num_dict[tax] += num

                if score == "incomplete":
                    if tax not in incom_tax_num_dict:
                        incom_tax_num_dict[tax] = num
                    else:
                        incom_tax_num_dict[tax] += num


        intact_tax_num_dict = dict(sorted(intact_tax_num_dict.items(), reverse=True, key=lambda item: item[1]))
        quest_tax_num_dict = dict(sorted(quest_tax_num_dict.items(), reverse=True, key=lambda item: item[1]))
        incom_tax_num_dict = dict(sorted(incom_tax_num_dict.items(), reverse=True, key=lambda item: item[1]))

        intact_entry = False
        quest_entry = False
        incom_entry = False

        if len(list(intact_tax_num_dict.keys())) > 1:
            intact_entry = list(intact_tax_num_dict.keys())[0]

        if len(list(quest_tax_num_dict.keys())) > 1:
            quest_entry = list(quest_tax_num_dict.keys())[0]

        if len(list(incom_tax_num_dict.keys())) > 1:
            incom_entry = list(incom_tax_num_dict.keys())[0]

        if len(list(intact_tax_num_dict.keys())) == 1:
            intact_entry = list(intact_tax_num_dict.keys())[0]

        if len(list(quest_tax_num_dict.keys())) == 1:
            quest_entry = list(quest_tax_num_dict.keys())[0]

        if len(list(incom_tax_num_dict.keys())) == 1:
            incom_entry = list(incom_tax_num_dict.keys())[0]

        if intact_entry:
            if island in intact_island_phage_dict:
                intact_island_phage_dict[island].append(intact_entry)
            else:
                intact_island_phage_dict[island] = [intact_entry]

        if quest_entry:
            if island in quest_island_phage_dict:
                quest_island_phage_dict[island].append(quest_entry)
            else:
                quest_island_phage_dict[island] = [quest_entry]

        if incom_entry:
            if island in incom_island_phage_dict:
                incom_island_phage_dict[island].append(incom_entry)
            else:
                incom_island_phage_dict[island] = [incom_entry]



# Making a dictionary of the DefenseFinder output in the format
# { gene : gene_defence }
DF_defence_gene_dict = {}
DF_defence_system_dict = {}
DF_defence_system_genes_dict = {}
DF_defence_gene_IDs_dict = {}
DF_defence_simple_system_dict = {}
DF_gene_simple_system_dict = {}
DF_defence_system_only_genes_dict = {}

with open(DefenseFinder_path) as f:
    next(f)
    for line in f:
        fields = line.strip().split("\t")
        if len(fields) == 23:
            replicon, hit_id, gene, gene_name, hit_pos, model_fqn, sys_id, sys_loci, locus_num, sys_wholeness, sys_score, sys_occ, hit_gene_ref, hit_status, hit_seq_len, hit_i_eval, hit_score, hit_profile_cov, hit_seq_cov, hit_begin_match, hit_end_match, counterpart, used_in = fields
        if len(fields) == 21:
            replicon, hit_id, gene, gene_name, hit_pos, model_fqn, sys_id, sys_loci, locus_num, sys_wholeness, sys_score, sys_occ, hit_gene_ref, hit_status, hit_seq_len, hit_i_eval, hit_score, hit_profile_cov, hit_seq_cov, hit_begin_match, hit_end_match = fields

        genome = replicon

        defence = sys_id
        system = "_".join(sys_id.split("_")[1:])
        simple_system = "_".join(system.split("_")[:-1])
        system_id = system.split("_")[-1]

        if simple_system in defence_name_dict:
            simple_system = defence_name_dict[simple_system]
            defence_system = "_".join([simple_system, system_id])
        else:
            defence_system = system

        if simple_system in simple_defence_dict:
            very_simple_system = simple_defence_dict[simple_system]
        else:
            very_simple_system = simple_system

        if gene in DF_defence_gene_dict:
            print("DF WARNING: Gene {} already in dictionary".format(str(gene)))
            print(DF_defence_gene_dict[gene], simple_system)

        DF_defence_gene_dict[gene] = "_".join([str(gene), gene_name])
        DF_gene_simple_system_dict[gene] = simple_system
        DF_defence_gene_IDs_dict[gene] = gene
        DF_defence_simple_system_dict[gene] = "_".join(defence_system.split("_")[:-1])
        DF_defence_system_dict[gene] = defence_system

        if defence_system in DF_defence_system_genes_dict:
            DF_defence_system_genes_dict[defence_system].append("_".join([str(gene), gene_name]))
            DF_defence_system_only_genes_dict[defence_system].append("_".join([str(gene), very_simple_system]))
        else:
            DF_defence_system_genes_dict[defence_system] = ["_".join([str(gene), gene_name])]
            DF_defence_system_only_genes_dict[defence_system] = ["_".join([str(gene), very_simple_system])]



# Doing the same for the PADLOC output
other_list = ['zorya_other', 'DMS_other', 'cbass_other', 'druantia_other', 'wadjet_other']

PADLOC_defence_gene_dict = {}
PADLOC_defence_system_genes_dict = {}
PADLOC_defence_system_dict = {}
PADLOC_defence_gene_IDs_dict = {}
PADLOC_defence_simple_system_dict = {}
PADLOC_gene_simple_system_dict = {}
PADLOC_defence_system_only_genes_dict = {}

with open(PADLOC_path) as f:
    next(f)
    for line in f:
        fields = line.strip("\n").split(",")

        fields = fields[:17]
        DI, system_number, seqid, system, target_name, hmm_accession, hmm_name, protein_name, full_seq_E_value, domain_iE_value, target_coverage, hmm_coverage, start, end, strand, target_description, relative_position = fields

        if system in other_list:
            continue

        if system in defence_name_dict:
            simple_system = defence_name_dict[system]
            defence_system = "_".join([simple_system, system_number])
        else:
            simple_system = system
            defence_system = "_".join([system, system_number])

        if simple_system in simple_defence_dict:
            very_simple_system = simple_defence_dict[simple_system]
        else:
            very_simple_system = simple_system

        gene = target_name

        if gene in PADLOC_defence_gene_dict:
            print("PADLOC WARNING: Gene {} already in dictionary".format(str(gene)))
            print(PADLOC_defence_gene_dict[gene], simple_system)

        PADLOC_defence_gene_dict[gene] = "_".join([str(gene), protein_name])
        PADLOC_gene_simple_system_dict[gene] = simple_system
        PADLOC_defence_gene_IDs_dict[gene] = gene
        PADLOC_defence_simple_system_dict[gene] = "_".join(defence_system.split("_")[:-1])
        PADLOC_defence_system_dict[gene] = defence_system

        if defence_system in PADLOC_defence_system_genes_dict:
            PADLOC_defence_system_genes_dict[defence_system].append("_".join([str(gene), protein_name]))
            PADLOC_defence_system_only_genes_dict[defence_system].append("_".join([str(gene), very_simple_system]))
        else:
            PADLOC_defence_system_genes_dict[defence_system] = ["_".join([str(gene), protein_name])]
            PADLOC_defence_system_only_genes_dict[defence_system] = ["_".join([str(gene), very_simple_system])]


# Adding information to dataframd and exporting
S2_df["Intact phage hit"] = S2_df["ID"].map(intact_island_phage_dict).fillna("")
S2_df["Questionable phage hit"] = S2_df["ID"].map(quest_island_phage_dict)
S2_df["Incomplete phage hit"] = S2_df["ID"].map(incom_island_phage_dict)

S2_df["P4_full_type"] = S2_df["ID"].map(P4_full_type_dict).fillna("")
S2_df["P4_simple_type"] = S2_df["ID"].map(P4_simple_type_dict).fillna("")
S2_df["P4_full_set"] = S2_df["ID"].map(P4_full_set_dict).fillna("")

S2_df["PICI_full_type"] = S2_df["ID"].map(PICI_full_type_dict).fillna("")
S2_df["PICI_simple_type"] = S2_df["ID"].map(PICI_simple_type_dict).fillna("")
S2_df["PICI_full_set"] = S2_df["ID"].map(PICI_full_set_dict).fillna("")

S2_df["ConjScan_info"] = S2_df["ID"].map(ConjScan_dict).fillna("")
S2_df["ConjScan_gene"] = [k[0] if isinstance(k, list) else k for k in S2_df["ConjScan_info"]]
S2_df["ConjScan_full_system"] = [k[1] if isinstance(k, list) else k for k in S2_df["ConjScan_info"]]
S2_df["ConjScan_simple_system"] = [k[2] if isinstance(k, list) else k for k in S2_df["ConjScan_info"]]


# Treating the gene lists as lists rather than string representations of lists
S2_df["IMG gene IDs in representative island"] = S2_df["IMG gene IDs in representative island"].apply(lambda x: x.split(", "))

S2_df["DF_defence_genes"] = S2_df["IMG gene IDs in representative island"].apply(lambda x: [DF_defence_gene_dict[y] for y in x if y in DF_defence_gene_dict])
S2_df["DF_defence_gene_IDs"] = S2_df["IMG gene IDs in representative island"].apply(lambda x: [DF_defence_gene_IDs_dict[y] for y in x if y in DF_defence_gene_IDs_dict])
S2_df["DF_defence_systems"] = S2_df["IMG gene IDs in representative island"].apply(lambda x: [DF_defence_system_dict[y] for y in x if y in DF_defence_system_dict]).apply(lambda x: list(set(x)))
S2_df["DF_defence_simple_systems"] = S2_df["DF_defence_systems"].apply(lambda x: ["_".join(y.split("_")[:-1]) for y in x])

S2_df.to_csv(output_path, header=True, sep="\t")
