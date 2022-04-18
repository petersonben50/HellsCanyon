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



**Tree generation**

Notes on rp16 alignment. Had to redo this, and on the first time, there were several bins with duplicate hits that were pains to fix, so I removed them from analysis (they weren't relevant seqs anyways, based on the first tree):
- GCA_021740345.1
- GCA_020719265.1
- GCA_016788725.1
- GCA_002749385.1
- GCA_016937755.1

These ones didn't have enough hits to make it into the final alignment, so I removed them ahead of time on round two as well:
GCA_016936715.1
GCA_003486765.1
GCA_002749305.1
GCA_007116765.1
GCA_013139955.1
GCA_903933215.1
GCA_913063265.1
GCA_003537645.1
GCA_021648495.1
GCA_003250655.1
GCA_903903255.1
GCA_913061395.1
GCA_016938055.1
GCA_903854615.1

*Notes from second run*

fall2017cluster6_bin_0136 had a duplicate at some point, so after concatenating the proteins I had to remove the duplicate.
List of them: `prolix_bins_to_remove.txt`.
GCA_913063105.1 also had one duplicate that I missed.
fall2017coassembly_bin_0181 only had one of the rp16 genes, despite being 61% complete according to checkM.
I removed this one.


*Notes from the first run on this*

There were duplicate sequences in all of the rp16 gene alignments. I manually fixed this in Geneious:
- rpL2: GCA_020719265.1 was duplicated, deleted one. GCA_021740345.1 was also duplicated, deleted one.
- rpL3: GCA_002749385.1 was duplicated, deleted one. GCA_016788725.1 was also duplicated, deleted one.
- rpL4: GCA_016788725.1 was duplicated, deleted one.
- rpL5: GCA_021740345.1 duplicated, deleted second one.
- rpL6: GCA_021740345.1 duplicated, deleted second one.
- rpl14: GCA_021740345.1 duplicated, deleted second one.
- rpL15: GCA_021740345.1 duplicated, deleted second truncated one. GCA_913063105.1 also duplicated, deleted second one.
- rpL16: GCA_016937755.1 is split. Combined the two sequences. GCA_021740345.1 is duplicated, deleted second truncated one.
- rpL18: GCA_021740345.1 duplicated, deleted second one.
- rpL22: GCA_020719265.1 duplicated, deleted second one. GCA_021740345.1 also duplicated, deleted second truncated one.
- rpL24: GCA_021740345.1 duplicated, deleted second one.
- rpS3: GCA_021740345.1 duplicated, deleted second one.
- rpS8: GCA_021740345.1 duplicated, deleted second one.
- rpS10: GCA_002749385.1 duplicated, deleted second one. GCA_016788725.1 also duplicated, deleted second one.
- rpS17: GCA_021740345.1 duplicated, deleted second one.
- rpS19: GCA_020719265.1 duplicated, deleted second one. GCA_021740345.1 also duplicated, deleted second one.

Alignment had 148 sequences in it, so a lot of the reference bin didn't have any of the rp16 genes.
Removed alignments that had less than 1000 residues in the alignment:

Left with 134 sequences. Mask at 50%.
