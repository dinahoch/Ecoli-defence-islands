## Ecoli-defence-islands

This is the source code for the article <a href="https://www.biorxiv.org/content/10.1101/2022.06.09.495481v1"><b>The defence island repertoire of the *Escherichia coli* pan-genome</b></a>.

### 1. Finding the hotspots identified in *E. coli* K-12 MG1655 in all genomes of *E. coli* using their core flanking genes
`Finding-all-hotspots-in-Ecoli.py`
###### Each hotspot was searched for in all analyzed E. coli genomes using the two genes immediately flanking the hotspot in the K-12 genome (Table S1). To exclude fragmented contigs, only contigs with more than 20 genes were considered. For cases where the immediately flanking genes were not found in the target genome, or in cases where multiple instances of immediately flanking genes were found, a window of ten genes on either side of the flanking genes in K-12 was searched for in the target genome, requiring at least five of the genes matching in gene order and cluster identity. If multiple matches were found, the closest sets of upstream and downstream flanking genes were selected to define the hotspot. Both flanking regions were required in the same contig to declare a hotspot. The resulting islands at these hotspots were filtered for those containing 200 or fewer genes and an individual hotspot was defined as “empty” if it comprised three or fewer genes.

### 2. Characterising the islands found at these hotspots
`Characterising-Ecoli-islands-at-hotspots.py`
#### This includes detecting MGEs using <a href="https://pubmed.ncbi.nlm.nih.gov/27141966/">PHASTER</a>, <a href="https://pubmed.ncbi.nlm.nih.gov/31584169/">CONJScan</a> and <a href="https://www.biorxiv.org/content/10.1101/2022.09.14.508007v1">SatelliteFinder</a>
###### Prophages: Prophages were identified using PHASTER analysis of the nucleotide sequences of each island. A phage hit was only considered if the island had more than one gene that matched the phage. When PHASTER identified intact prophages, the taxonomy of the phage hit was recorded using NCBI classification.
###### Phage satellites: Phage satellites were detected using SatelliteFinder (Galaxy Version 0.9) analysis of the amino acid sequences of genes in each island. P4-like satellites were only considered if they were predicted to be of Types A, B or C, and PICI satellites if predicted to be of Types A or B. Manual inspection of islands annotated to contain PICI satellites revealed several of these to be intact Uetakevirus prophages; the annotation was changed accordingly. When PHASTER and SatelliteFinder gave overlapping predictions, the SatelliteFinder prediction was used and the predicted was checked by manual inspection.
###### Conjugative plasmids and integrative conjugative elements: The amino acid sequences of genes in each island were submitted to CONJscan with default parameters (Galaxy Version 2.0+galaxy1). This identified both conjugative type IV secretion systems (T4SS) and also integrative mobilizable elements that can hijack T4SS. In addition, islands were searched for the genes VirB5 and VirB6, which are known to be components of non-canonical integrative mobilizable elements that lack a relaxase gene.
#### and defence systems using <a href="https://academic.oup.com/nar/article/49/19/10868/6381132">PADLOC</a> and <a href="https://www.nature.com/articles/s41467-022-30269-9">DefenseFinder</a>
###### DefenseFinder and PADLOC were utilized to identify known defence systems in each island. Amino acid sequences of genes in all 1,351 E. coli genomes were submitted to DefenseFinder release version 1.2.0. Genes predicted to be part of multiple different defence systems were inspected manually for proper annotation. Amino acid sequences and gff3 files of genes in each island were submitted to the PADLOC web server v1.1.0 with defence systems included in padlocdb v1.4.0. Systems annotated by PADLOC as “[system]_other” were excluded, since these represent partial or separated defence systems. When two overlapping systems of the same type were predicted by the two tools, all constituent genes were considered part of this system.

### 3. Mapping defence systems found in finished genomes by DefenseFinder to *E. coli* K-12
`Mapping-DefenseFinder-systems-to-K12.py`
###### Amino acid sequences of genes in the main chromosomes of all finished E. coli genomes were similarly analyzed using DefenseFinder release version 1.2.0. Defense systems were mapped to the E. coli K-12 reference genome as described above. These were then manually examined to identify exactly where in the genome they were integrated. If the integration position constituted a hotspot but this had not been detected due to deletion of flanking core genes on one or both sides, this was manually recorded in Table S3.
