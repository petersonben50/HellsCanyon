### Phylogeny for Prolixibacter-like bin

I wanted to look more closely at the fine-scale phylogeny of the Prolixibacter genus to determine if our bin is in that genus.
If it is, then we'll look more closely at the metabolism and the *hgcA* genes from those bins.

**GTDB tree**

First we'll look at the GTDB tree to see if there are more genomes to go after.
Download `~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/taxonomy/classify/gtdbtk.bac120.classify.tree` to my local computer: `/Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/bins/binAnalysis/phylogeny/bacteroidetes/prolixibacter`.
I loaded in the GTDB tree, then subsetted the tree 2 levels back from the bin of interest (anvio_hgcA_0130).
This shows that the two genomes that we have in the Bacteroidetes tree are the two closest we have in GTDB as well.


**NCBI taxonomy browser**

I also looked here for more references:
https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi.
Looks like there are more we can include.
We already have GCF_000621705.1 (*Prolixibacter bellariivorans* ATCC BAA-1284) and GCF_003014495.1 (*Prolixibacter denitrificans*).

Others to include:
- Prolixibacter denitrificans MIC1-1: GCF_009617855.1
- Prolixibacter bellariivorans JCM 13498: GCF_009617915.1
- Prolixibacter sp. NT017: GCF_009617875.1
- Prolixibacter sp. SD074: GCF_009617895.1

This is a small enough number that we should just go back to the Bacteroidetes tree and incorporate them into there. 
