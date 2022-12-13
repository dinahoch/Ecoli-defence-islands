# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 11:25:37 2020
Modified on Wed Nov 17 09:02:31 2021

@author: dinah
"""

## This script will islands in 1351 E.coli genomes based on the core flanking
## genes identified in E. coli K-12, defined as genes present in >=80% of E.
## coli genomes

import sys
from pathlib import Path
from itertools import product
from Bio import SeqIO
from Bio import Align
import pandas as pd
aligner = Align.PairwiseAligner()

# File paths
Ecoli_allgenes_path = sys.argv[1]
lib_path = sys.argv[2]
SC1_path = sys.argv[3]
output_path = sys.argv[4]

hotspot = sys.argv[5]


# This is a dictionary with the format { orig_cluster : super_cluster }
origcl_SC1_dict = {}
with open(SC1_path) as f:
    for line in f:
        SC1, orig_cl = line.strip().split()
        origcl_SC1_dict[int(orig_cl)] = int(SC1)


# Dataframe of all Ecoli genes
Ecoli_allgenes_df = pd.read_csv(Ecoli_allgenes_path, dtype=str, sep="\t", header=None, engine="python")
Ecoli_allgenes_df.columns = ["SC2", "SC1", "orig_cl", "gene", "contig", "genome", "domain", "deflated", "from", "to", "strand", "contig_ID"]


# Treating all dataframe values as numbers if possible
Ecoli_allgenes_df = Ecoli_allgenes_df.apply(pd.to_numeric, downcast="integer", errors="coerce").fillna(Ecoli_allgenes_df)
Ecoli_allgenes_df["contig_ID"] = Ecoli_allgenes_df["contig_ID"].astype('Int64')

# Assigning each gene in each contig a position
Ecoli_allgenes_df = Ecoli_allgenes_df.assign(pos=Ecoli_allgenes_df.groupby(Ecoli_allgenes_df.contig.ne(Ecoli_allgenes_df.contig.shift()).cumsum()).cumcount())

# Counting the number of genes per contig
contig_genes_dict = Ecoli_allgenes_df.groupby("contig")["gene"].nunique().to_dict()
Ecoli_allgenes_df["num_contig_genes"] = Ecoli_allgenes_df["contig"].map(contig_genes_dict)


# This is a dictionary of dictionaries with the format
# { contig : { pos : gene } }
contig_pos_gene_dict = Ecoli_allgenes_df.groupby("contig").apply(lambda x: x.set_index("pos")["gene"].to_dict()).to_dict()


# This is a dictionary of dictionaries with the format
# { contig : { gene : pos } }
contig_gene_pos_dict = Ecoli_allgenes_df.groupby("contig").apply(lambda x: x.set_index("gene")["pos"].to_dict()).to_dict()


# Filtering the contig_pos_gene dictionary to get only contigs with more than
# 20 genes
contig_pos_gene_dict_over20 = {key:value for (key, value) in contig_pos_gene_dict.items() if len(value) > 20}
contig_pos_gene_over20_list = list(contig_pos_gene_dict_over20.keys())


# Subsetting the original dataframe to include only contigs with more than 20
# genes
Ecoli_allgenes_over20_df = Ecoli_allgenes_df[Ecoli_allgenes_df.contig.isin(contig_pos_gene_over20_list)]


# This function returns the original cluster for a given gene
def gene_to_cluster(gene):
    cluster = Ecoli_allgenes_over20_df.loc[Ecoli_allgenes_over20_df["gene"] == gene, "orig_cl"].iloc[0]
    return cluster


# Flanking core genes in E. coli K-12
hotspot_dict = {
    "1" : "2688075859_2688075861",
    "2" : "2688075873_2688075874",
    "3" : "2688076039_2688076040",
    "4" : "2688076467_2688076469",
    "5" : "2688076569_2688076571",
    "6" : "2688076586_2688076587",
    "7" : "2688076778_2688076788",
    "8" : "2688076809_2688076811",
    "9" : "2688076891_2688076921",
    "10" : "2688077039_2688077040",
    "11" : "2688077186_2688077203",
    "12" : "2688077475_2688077478",
    "13" : "2688077576_2688077578",
    "14" : "2688077579_2688077581",
    "15" : "2688077581_2688077583",
    "16" : "2688077587_2688077590",
    "17" : "2688077590_2688077593",
    "18" : "2688077652_2688077656",
    "19" : "2688077856_2688077857",
    "20" : "2688078071_2688078072",
    "21" : "2688078208_2688078238",
    "22" : "2688078395_2688078397",
    "23" : "2688078416_2688078417",
    "24" : "2688078439_2688078464",
    "25" : "2688078567_2688078572",
    "26" : "2688078596_2688078597",
    "27" : "2688078700_2688078701",
    "28" : "2688078719_2688078721",
    "29" : "2688078757_2688078759",
    "30" : "2688078785_2688078786",
    "31" : "2688078832_2688078833",
    "32" : "2688079080_2688079112",
    "33" : "2688079363_2688079408",
    "34" : "2688079435_2688079437",
    "35" : "2688079682_2688079683",
    "36" : "2688079707_2688079723",
    "37" : "2688079754_2688079798",
    "38" : "2688079931_2688079933",
    "39" : "2688080010_2688080011",
    "40" : "2688080164_2688080165",
    "41" : "2688080373_2688080374"
    }

K12_up_gene, K12_down_gene = hotspot_dict[hotspot].split("_")

K12_up_gene = int(K12_up_gene)
K12_down_gene = int(K12_down_gene)

K12_contig = "Ga0123766_11"
K12_genome = 2687453259
K12_contig_ID = 2687453745


# Starting the analysis for a given hotspot
gff_path = Path(lib_path + "/" + str(K12_genome) + "/" + str(K12_genome) + ".gff")
contig_fasta_path = Path(lib_path + "/" + str(K12_genome) + "/" + str(K12_genome) + ".fna")
K12_gene_fasta_path = Path(lib_path + "/" + str(K12_genome) + "/" + str(K12_genome) + ".genes.faa")

gff_df = pd.read_csv(gff_path, dtype=str, sep="\t", skiprows=1, header=None)
gff_df.columns = ["seqid", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "ID"]

K12_gene_fasta_file = open(K12_gene_fasta_path)
K12_gene_fasta_dict = SeqIO.to_dict(SeqIO.parse(K12_gene_fasta_file, "fasta"))

K12_up_seq = str(K12_gene_fasta_dict[str(K12_up_gene)].seq)
K12_down_seq = str(K12_gene_fasta_dict[str(K12_down_gene)].seq)

K12_start_coord = gff_df[gff_df["ID"].str.contains(str(K12_up_gene), na=False)]["start_coord"].iloc[0]
K12_end_coord = gff_df[gff_df["ID"].str.contains(str(K12_down_gene), na=False)]["end_coord"].iloc[0]

K12_link = "https://img.jgi.doe.gov/cgi-bin/mer/main.cgi?section=ScaffoldGraph&page=scaffoldGraph&scaffold_oid={}&start_coord={}&end_coord={}".format(str(K12_contig_ID), str(K12_start_coord), str(K12_end_coord))

# This is a dataframe with all K12 genes and their information
K12_contig = "Ga0123766_11"
K12_allgenes_df = Ecoli_allgenes_over20_df.loc[Ecoli_allgenes_over20_df["contig"] == K12_contig]

# Obtaining the information for these genes
K12_up_cluster = gene_to_cluster(K12_up_gene)
K12_up_strand = Ecoli_allgenes_over20_df.loc[Ecoli_allgenes_over20_df["gene"] == K12_up_gene, "strand"].iloc[0]
K12_up_pos = Ecoli_allgenes_over20_df.loc[Ecoli_allgenes_over20_df["gene"] == K12_up_gene, "pos"].iloc[0]

K12_down_cluster = gene_to_cluster(K12_down_gene)
K12_down_strand = Ecoli_allgenes_over20_df.loc[Ecoli_allgenes_over20_df["gene"] == K12_down_gene, "strand"].iloc[0]
K12_down_pos = Ecoli_allgenes_over20_df.loc[Ecoli_allgenes_over20_df["gene"] == K12_down_gene, "pos"].iloc[0]


# Storing ten genes upstream of the K12_up_gene and their information
K12_up10_df = K12_allgenes_df[K12_up_pos-9:K12_up_pos+1]
K12_up_gene_list = K12_up10_df["gene"].tolist()
K12_up_cluster_list = K12_up10_df["orig_cl"].tolist()
K12_up_cluster_set = set(K12_up_cluster_list)

# and downstream
K12_down10_df = K12_allgenes_df[K12_down_pos:K12_down_pos+10]
K12_down_gene_list = K12_down10_df["gene"].tolist()
K12_down_cluster_list = K12_down10_df["orig_cl"].tolist()
K12_down_cluster_set = set(K12_down_cluster_list)

# minimum intersections
min_up_int = 0.5*len(K12_up_cluster_set)
min_down_int = 0.5*len(K12_down_cluster_set)


# Looping through each contig, the logic will be:
# 1. Check if both K12 flanking clusters are found in the contig within 200
# genes. If yes, these are the flanking genes and we have our island hotspot.
# If no, continue to (2).
# 2. Go through the contig in chunks of ten genes. Compare each set of ten
# genes to the upstream and downstream K12 sets, requiring an intersection of
# 5/10 between the two. First check for upstream. If found, repeat for
# downstream within 200 genes. If no upstream found, this contig does not
# contain the hotspot.

mode = "run"

output_file = open(output_path, "a+")

# Looping through each contig
for contig in contig_pos_gene_dict_over20.keys():

    # print(contig)

    contig_ID, genome = Ecoli_allgenes_df.loc[Ecoli_allgenes_df["contig"] == contig, ["contig_ID", "genome"]].iloc[0]

    gff_path = Path(lib_path + "/" + str(genome) + "/" + str(genome) + ".gff")
    contig_fasta_path = Path(lib_path + "/" + str(genome) + "/" + str(genome) + ".fna")
    gene_fasta_path = Path(lib_path + "/" + str(genome) + "/" + str(genome) + ".genes.faa")

    gff_df = pd.read_csv(gff_path, dtype=str, sep="\t", skiprows=1, header=None)
    gff_df.columns = ["seqid", "source", "type", "start_coord", "end_coord", "score", "strand", "phase", "ID"]

    gene_fasta_file = open(gene_fasta_path, "r")
    gene_fasta_dict = SeqIO.to_dict(SeqIO.parse(gene_fasta_file, "fasta"))

    # Initial parameters
    up_gene_found = 0
    up_uniq = 0
    up_pos = 0
    up_cluster_found = 0
    real_up_gene = 0
    up_int_dict = {}
    real_up_found = 0
    up_gene_int = 0

    down_gene_found = 0
    down_uniq = 0
    down_pos = 0
    down_cluster_found = 0
    real_down_gene = 0
    down_int_dict = {}
    real_down_found = 0
    down_gene_int = 0

    # Storing information for this contig in a dataframe
    contig_df = Ecoli_allgenes_over20_df.loc[Ecoli_allgenes_over20_df["contig"] == contig]

    # Number of genes on contig
    num_contig_genes = Ecoli_allgenes_over20_df.loc[Ecoli_allgenes_over20_df["contig"] == contig, "num_contig_genes"].iloc[0]

    # Creating a dictionary in the format
    # { gene : cluster }
    gene_cluster_dict = dict(zip(contig_df.gene, contig_df.orig_cl))

    # Creating a dictionary in the format
    # { cluster : gene(s) } where gene(s) are in the form of a list
    cluster_gene_dict = {}

    for key, value in gene_cluster_dict.items():
        if value in cluster_gene_dict:
            cluster_gene_dict[value].append(key)
        else:
            cluster_gene_dict[value] = [key]

    # Before we go through contig_df in chunks, let's see if the K12_up_gene
    # and K12_down gene clusters are found in the genome.

    # First looking for the upstream gene
    try:
        up_gene_found = cluster_gene_dict[K12_up_cluster]
        up_cluster_found = 1
    except KeyError:
        pass

    # Now downstream
    try:
        down_gene_found = cluster_gene_dict[K12_down_cluster]
        down_cluster_found = 1
    except KeyError:
        pass

    if mode == "check":
        print("up_cluster_found, down_cluster_found", up_gene_found, down_gene_found)


    # We have three scenarios:
        # 1. Both up- and downstream genes are found
        # 2. Only one is found
        # 3. Neither are found

    # Let's deal with (1) and (2) first. If both are found, consider their
    # strand(s) compared to the K12 strand(s) to determine the direction of the
    # contig, strand_dir. Based on this, take either ten genes up- or
    # downstream and intersect their set with that of the K12 genes. If the
    # intersection is greater than a minimum parameter, min_int, then add to
    # the dictionary up_int_dict or down_int_dict. Take the dictionary entry
    # that has the largest intersection with the K12 set as the real gene.

    if up_cluster_found:

        for up in up_gene_found:
            up_pos = contig_gene_pos_dict[contig][up]
            up_strand = gff_df[gff_df["ID"].str.contains(str(up), na=False)]["strand"].iloc[0]

            if up_strand == K12_up_strand:
                strand_dir = True
            else:
                strand_dir = False

            n = 10 # chunk size

            if strand_dir:

                if up_pos-n+1 < 0:
                    start_pos = 0
                else:
                    start_pos = up_pos-n+1

                if up_pos+1 > num_contig_genes:
                    end_pos = num_contig_genes
                else:
                    end_pos = up_pos+1

            else:
                start_pos = up_pos

                if up_pos+n > num_contig_genes:
                    end_pos = num_contig_genes
                else:
                    end_pos = up_pos+n

            up_chunk_df = contig_df[start_pos:end_pos]

            up_gene_list = up_chunk_df["gene"].tolist()
            up_cluster_list = up_chunk_df["orig_cl"].tolist()
            up_cluster_set = set(up_cluster_list)

            up_int = up_cluster_set.intersection(K12_up_cluster_set)

            if len(up_int) >= min_up_int:
                real_up_gene = up

    # Now we repeat for the downstream gene
    if down_cluster_found:

        for down in down_gene_found:
            down_pos = contig_gene_pos_dict[contig][down]
            down_strand = gff_df[gff_df["ID"].str.contains(str(down), na=False)]["strand"].iloc[0]

            if down_strand == K12_down_strand:
                strand_dir = True
            else:
                strand_dir = False

            n = 10 # chunk size

            if strand_dir:

                if down_pos-n+1 < 0:
                    start_pos = 0
                else:
                    start_pos = down_pos-n+1

                if down_pos+1 > num_contig_genes:
                    end_pos = num_contig_genes
                else:
                    end_pos = down_pos+1

            else:

                start_pos = down_pos

                if down_pos+n > num_contig_genes:
                    end_pos = num_contig_genes
                else:
                    end_pos = down_pos+n

            down_chunk_df = contig_df[start_pos:end_pos]

            down_gene_list = down_chunk_df["gene"].tolist()
            down_cluster_list = down_chunk_df["orig_cl"].tolist()
            down_cluster_set = set(down_cluster_list)

            down_int = down_cluster_set.intersection(K12_down_cluster_set)

            if len(down_int) >= min_down_int:
                real_down_gene = down


        if mode == "check":
            print("1:", real_up_gene, real_down_gene)

    # In case (1), both genes are found. In this case, we are done and can
    # record these genes as the real flanking genes.

    if (real_up_gene and real_down_gene):
        pass

    # Now we deal with case (2): only one of the genes is found. Its
    # intersection was calculated above. Now we look for the other gene that
    # was not found previously. We do this by scanning the contig in the
    # correct direction until we find a set of ten genes with a sufficient
    # overlap with the K12 set.

    elif (real_up_gene or real_down_gene):

        # If the up gene was found, we look for the down gene
        if (real_up_gene and not real_down_gene):

            down_int_dict = {}

            for pos in range(0, num_contig_genes):

                # Subsetting contig_df in ten-gene chunks
                n = 10  # chunk row size

                down_chunk_df = contig_df[pos:pos+n]

                # Getting the information for this ten-gene chunk
                down_gene_list = down_chunk_df["gene"].tolist()
                down_cluster_list = down_chunk_df["orig_cl"].tolist()
                down_cluster_set = set(down_cluster_list)

                down_int = down_cluster_set.intersection(K12_down_cluster_set)
                down_int_index = [down_cluster_list.index(i) for i in down_int]
                down_int_genes = [down_gene_list[i] for i in down_int_index]


                if len(down_int) >= min_down_int:
                    len_down_int = len(down_int)
                    if len_down_int in down_int_dict:
                            down_int_dict[len_down_int] += down_int_genes
                    else:
                        down_int_dict[len_down_int] = down_int_genes

            if down_int_dict:
                real_down_gene = down_int_dict[max(down_int_dict, key=float)]


        elif (real_down_gene and not real_up_gene):

            for pos in range(10, num_contig_genes):

                # Subsetting contig_df in ten-gene chunks
                n = 10  # chunk row size

                up_chunk_df = contig_df[pos-10:pos]

                # Getting the information for this ten-gene chunk
                up_gene_list = up_chunk_df["gene"].tolist()
                up_cluster_list = up_chunk_df["orig_cl"].tolist()
                up_cluster_set = set(up_cluster_list)

                up_int = up_cluster_set.intersection(K12_up_cluster_set)
                up_int_index = [up_cluster_list.index(i) for i in up_int]
                up_int_genes = [up_gene_list[i] for i in up_int_index]

                if len(up_int) >= min_up_int:
                    len_up_int = len(up_int)
                    if len_up_int in up_int_dict:
                        up_int_dict[len_up_int] += up_int_genes
                    else:
                        up_int_dict[len_up_int] = up_int_genes

            if up_int_dict:
                real_up_gene = up_int_dict[max(up_int_dict, key=float)]


    if mode == "check":
        print("2:", real_up_gene, real_down_gene)

    # At this point, we may have found the real_up_gene and real_down_gene
    # through either finding their clusters or using one as an anchor and
    # finding the second.
    if (real_up_gene and real_down_gene):
        pass


    # If this is not the case and only one was found, we take the cluster of
    # the other (if found) as the second flanking gene. If there are multiple
    # genes corresponding to this cluster, we take the one that is most
    # similar to the K12 gene. If there are multiple similar genes, we take
    # the gene closest to the up or down gene previously found.
    elif (real_up_gene or real_down_gene):

        if (real_up_gene and not real_down_gene):

            if down_cluster_found:

                if len(down_gene_found) == 1:
                    real_down_gene = down_gene_found[0]

                else:

                    down_score_dict = {}

                    for down in down_gene_found:
                        down_seq = str(gene_fasta_dict[str(down)].seq)
                        score = aligner.score(down_seq, K12_down_seq)
                        down_score_dict[down] = score

                    real_down_gene = max(down_score_dict, key=down_score_dict.get)


        elif (real_down_gene and not real_up_gene):

            if up_cluster_found:

                if len(up_gene_found) == 1:
                    real_up_gene = up_gene_found[0]

                else:

                    up_score_dict = {}

                    for up in up_gene_found:
                        up_seq = str(gene_fasta_dict[str(up)].seq)
                        score = aligner.score(up_seq, K12_up_seq)
                        up_score_dict[up] = score

                    real_up_gene = max(up_score_dict, key=up_score_dict.get)

        if mode == "check":
            print("3:", real_up_gene, real_down_gene)


    # This is case (3), in which neither gene was found. In this case, we scan
    # the entire contig to find the up and down genes using their
    # intersections. If this does not result in the real up and down genes, if
    # the clusters were found, we take these as the real genes.
    else:

        up_int_dict = {}
        down_int_dict = {}

        # Scanning the whole contig
        for pos in range(0, num_contig_genes):

            # Subsetting contig_df in ten-gene chunks
            n = 10  # chunk row size

            up_chunk_df = contig_df[pos:pos+n]

            # Getting the information for this ten-gene chunk
            up_gene_list = up_chunk_df["gene"].tolist()
            up_cluster_list = up_chunk_df["orig_cl"].tolist()
            up_cluster_set = set(up_cluster_list)

            # Intersecting the clusters of these ten genes with those in the
            # K12 upstream set
            up_int = up_cluster_set.intersection(K12_up_cluster_set)
            up_int_index = [up_cluster_list.index(i) for i in up_int]
            up_int_genes = [up_gene_list[i] for i in up_int_index]

            # If the intersection is more than a specified value, add to
            # dictionary, with the position included in the dictionary key
            if len(up_int) >= min_up_int:
                len_up_int = len(up_int)
                if len_up_int in up_int_dict:
                    up_int_dict[len_up_int] += up_int_genes
                else:
                    up_int_dict[len_up_int] = up_int_genes

            down_chunk_df = contig_df[pos:pos+n]

            # Getting the information for this ten-gene chunk
            down_gene_list = down_chunk_df["gene"].tolist()
            down_cluster_list = down_chunk_df["orig_cl"].tolist()
            down_cluster_set = set(down_cluster_list)

            # Intersecting the clusters of these ten genes with those in the
            # K12 downstream set
            down_int = down_cluster_set.intersection(K12_down_cluster_set)
            down_int_index = [down_cluster_list.index(i) for i in down_int]
            down_int_genes = [down_gene_list[i] for i in down_int_index]

            # If the intersection is more than a specified value, add to
            # dictionary, with the position included in the dictionary key
            if len(down_int) >= min_down_int:
                len_down_int = len(down_int)
                if len_down_int in down_int_dict:
                    down_int_dict[len_down_int] += down_int_genes
                else:
                    down_int_dict[len_down_int] = down_int_genes

        if up_int_dict and down_int_dict:
            real_gene_list = [min(up_int_dict[max(up_int_dict, key=float)]), max(up_int_dict[max(up_int_dict, key=float)]), min(down_int_dict[max(down_int_dict, key=float)]), max(down_int_dict[max(down_int_dict, key=float)])]
            real_gene_list.sort()
            real_up_gene = real_gene_list[1]
            real_down_gene = real_gene_list[2]

        else:
            if up_int_dict:
                real_up_gene = up_int_dict[max(up_int_dict, key=float)] # not sure this is right
            if down_int_dict:
                real_down_gene = down_int_dict[max(down_int_dict, key=float)] # or this

    if mode == "check":
        print("4:", real_up_gene, real_down_gene)

    # Have we now found both genes?
    if (real_up_gene and real_down_gene):
        pass

    # If not, check if we've found one of them
    elif (real_up_gene or real_down_gene):

        if (real_up_gene and not real_down_gene):

            if down_cluster_found:

                if len(down_gene_found) == 1:
                    real_down_gene = down_gene_found[0]

                else:

                    down_score_dict = {}

                    for down in down_gene_found:
                        down_seq = str(gene_fasta_dict[str(down)].seq)
                        score = aligner.score(down_seq, K12_down_seq)
                        down_score_dict[down] = score

                    real_down_gene = max(down_score_dict, key=down_score_dict.get)


        elif (real_down_gene and not real_up_gene):
            if up_cluster_found:

                if len(up_gene_found) == 1:
                    real_up_gene = up_gene_found[0]

                else:

                    up_score_dict = {}

                    for up in up_gene_found:
                        up_seq = str(gene_fasta_dict[str(up)].seq)
                        score = aligner.score(up_seq, K12_up_seq)
                        up_score_dict[up] = score

                    real_up_gene = max(up_score_dict, key=up_score_dict.get)


        if mode == "check":
            print("5:", real_up_gene, real_down_gene)

    # Neither genes were found. If the clusters were found, take these.
    # Otherwise, the contig does not contain the hotspot.
    else:

        if up_cluster_found:

            if len(up_gene_found) == 1:
                real_up_gene = up_gene_found[0]

            else:

                up_score_dict = {}

                for up in up_gene_found:
                    up_seq = str(gene_fasta_dict[str(up)].seq)
                    score = aligner.score(up_seq, K12_up_seq)
                    up_score_dict[up] = score

                real_up_gene = max(up_score_dict, key=up_score_dict.get)


        if down_cluster_found:

            if len(down_gene_found) == 1:
                real_down_gene = down_gene_found[0]

            else:

                down_score_dict = {}

                for down in down_gene_found:
                    down_seq = str(gene_fasta_dict[str(down)].seq)
                    score = aligner.score(down_seq, K12_down_seq)
                    down_score_dict[down] = score

                real_down_gene = max(down_score_dict, key=down_score_dict.get)

    if mode == "check":
        print("6:", real_up_gene, real_down_gene)

    # Now we summarise and format the output information
    if (real_up_gene and real_down_gene):

        if (type(real_up_gene) == list or type(real_down_gene) == list):
            if (type(real_up_gene) == list and type(real_down_gene) == int):
                if up_gene_found:
                    real_up_found = set(real_up_gene).intersection(set(up_gene_found))
                if real_up_found:
                    real_up_gene = list(real_up_found)[0]
                else:
                    real_up_gene = min(real_up_gene, key=lambda x: abs(x-real_down_gene) if x < real_down_gene else 10000000000000)

            elif (type(real_up_gene) == int and type(real_down_gene) == list):
                if down_gene_found:
                    real_down_found = set(real_down_gene).intersection(set(down_gene_found))
                if real_down_found:
                    real_down_gene = list(real_down_found)[0]
                else:
                    real_down_gene = min(real_down_gene, key=lambda x: abs(x-real_up_gene) if x > real_up_gene else 10000000000000)

            elif (type(real_up_gene) == list and type(real_down_gene) == list):
                if up_gene_found:
                    real_up_found = set(real_up_gene).intersection(set(up_gene_found))
                if down_gene_found:
                    real_down_found = set(real_down_gene).intersection(set(down_gene_found))
                if real_up_found:
                    real_up_gene = list(real_up_found)[0]
                if real_down_found:
                    real_down_gene = list(real_down_found)[0]

                if (real_up_found and real_down_found):
                    pass

                else:
                    if (real_up_found and not real_down_found):
                        real_down_gene = min(real_down_gene, key=lambda x: abs(x-real_up_gene) if x > real_up_gene else 10000000000000)

                    elif (real_down_found and not real_up_found):
                        real_up_gene = min(real_up_gene, key=lambda x: abs(x-real_down_gene) if x < real_down_gene else 10000000000000)

                    else:
                        up_gene, down_gene = min(product(real_up_gene, real_down_gene), key=lambda x: abs(x[0]-x[1]))

                        if up_gene == down_gene:
                            real_up_gene = up_gene
                            real_down_gene = min(real_down_gene, key=lambda x: abs(x-up_gene) if x > up_gene else 10000000000000)
                        else:
                            real_up_gene = up_gene
                            real_down_gene = down_gene

        if mode == "check":
            print("7:", real_up_gene, real_down_gene)



        if real_up_gene < real_down_gene:
            hotspot_genes = list(range(real_up_gene, real_down_gene+1))
        elif real_down_gene < real_up_gene:
            hotspot_genes = list(range(real_down_gene, real_up_gene+1))
        elif real_up_gene == real_down_gene:
            print("Error: Up and down genes are the same")


        if up_gene_found and down_gene_found:
            up_gene_int = set(hotspot_genes).intersection(set(up_gene_found))
            down_gene_int = set(hotspot_genes).intersection(set(down_gene_found))

        if up_gene_int and down_gene_int:
            if real_up_gene < real_down_gene:
                real_up_gene = list(up_gene_int)[0]
                real_down_gene = list(down_gene_int)[0]
            elif real_down_gene < real_up_gene:
                real_up_gene = list(down_gene_int)[0]
                real_down_gene = list(up_gene_int)[0]
            else:
                print("Error: Up and down genes are the same")

        if mode == "check":
            print("8:", real_up_gene, real_down_gene)

        up_start_coord = int(gff_df[gff_df["ID"].str.contains(str(real_up_gene), na=False)]["start_coord"].iloc[0])
        up_end_coord = int(gff_df[gff_df["ID"].str.contains(str(real_up_gene), na=False)]["end_coord"].iloc[0])
        down_start_coord = int(gff_df[gff_df["ID"].str.contains(str(real_down_gene), na=False)]["start_coord"].iloc[0])
        down_end_coord = int(gff_df[gff_df["ID"].str.contains(str(real_down_gene), na=False)]["end_coord"].iloc[0])

        if real_up_gene < real_down_gene:
            strand = "+"
            real_start_coord = up_end_coord
            real_end_coord = down_start_coord
            up_pos = contig_gene_pos_dict[contig][real_up_gene]
            down_pos = contig_gene_pos_dict[contig][real_down_gene]

        elif real_down_gene < real_up_gene:
            strand = "-"
            real_start_coord = down_end_coord
            real_end_coord = up_start_coord
            up_pos = contig_gene_pos_dict[contig][real_down_gene]
            down_pos = contig_gene_pos_dict[contig][real_up_gene]

        else:
            print("Error: Up and down genes are the same")

        size = abs(real_end_coord - real_start_coord)
        num_genes = abs(down_pos - up_pos) - 1

        updown_chunk_df = contig_df[up_pos:down_pos-1]
        updown_gene_list = updown_chunk_df["gene"].tolist()
        updown_cluster_list = updown_chunk_df["orig_cl"].tolist()


        if updown_cluster_list:
            updown_SC1_list = list(map(origcl_SC1_dict.get, updown_cluster_list))
            empty = "island"
        else:
            updown_SC1_list = []
            empty = "empty"


        with open(contig_fasta_path, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                if record.id == contig:
                    if strand == "+":
                        sequence = str(record.seq[int(real_start_coord):int(real_end_coord)])
                    if strand == "-":
                        sequence = str(record.seq[int(real_start_coord):int(real_end_coord)].reverse_complement())

        hotspot_link = "https://img.jgi.doe.gov/cgi-bin/mer/main.cgi?section=ScaffoldGraph&page=scaffoldGraph&scaffold_oid={}&start_coord={}&end_coord={}".format(str(contig_ID), str(real_start_coord), str(real_end_coord))

        output_file.write("\t".join([str(hotspot), str(K12_up_gene), str(K12_down_gene), str(K12_link), "hotspot found", str(empty), str(genome), str(contig), str(contig_ID), str(num_contig_genes), str(real_up_gene), str(real_down_gene), str(real_start_coord), str(real_end_coord), str(size), str(num_genes), str(hotspot_link), str(updown_gene_list), str(updown_cluster_list), str(updown_SC1_list), strand, str(sequence)])+"\n")

    elif (real_up_gene or real_down_gene):

        if real_up_gene:
            output_file.write("\t".join([str(hotspot), str(K12_up_gene), str(K12_down_gene), str(K12_link), "up_found", "0", str(genome), str(contig), str(contig_ID), str(num_contig_genes), str(real_up_gene), "0"])+"\n")

            end_contig = False

            if type(real_up_gene) == list:
                up_pos = [contig_gene_pos_dict[contig][i] for i in real_up_gene]
                for pos in up_pos:
                    if pos <= 10 or num_contig_genes - pos <= 10:
                        end_contig = True

            elif type(real_up_gene) == int:
                up_pos = contig_gene_pos_dict[contig][real_up_gene]
                if up_pos <= 10 or num_contig_genes - up_pos <= 10:
                    end_contig = True


        if real_down_gene:
            output_file.write("\t".join([str(hotspot), str(K12_up_gene), str(K12_down_gene), str(K12_link), "down_found", "0", str(genome), str(contig), str(contig_ID), str(num_contig_genes), "0", str(real_down_gene)])+"\n")

            end_contig = False

            if type(real_down_gene) == list:
                down_pos = [contig_gene_pos_dict[contig][i] for i in real_down_gene]
                for pos in down_pos:
                    if pos <= 10 or num_contig_genes - pos <= 10:
                        end_contig = True

            elif type(real_down_gene) == int:
                down_pos = contig_gene_pos_dict[contig][real_down_gene]
                if down_pos <= 10 or num_contig_genes - down_pos <= 10:
                    end_contig = True

    else:
        output_file.write("\t".join([str(hotspot), str(K12_up_gene), str(K12_down_gene), str(K12_link), "no hotspot", "0", str(genome), str(contig), str(contig_ID), str(num_contig_genes)])+"\n")


output_file.close()
