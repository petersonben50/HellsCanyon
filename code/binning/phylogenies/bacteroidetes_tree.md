### Scripts to generate Bacteroidetes tree for HCC

I need to generate a tree that focuses on the Bacteroidetes bins from HCC, since the discussion of this bin will be a major talking point for the paper.
I'll base this on the rp16 proteins.
For references, I included the hgcA+ Bacteroidetes bins from Lake Mendota, some reference isolate genomes that I had picked out for the 5M study, and a few other hand-selected ones from Elizabeth's paper.


**Collect references**

The first thing I needed to do was grab the needed references.
I decided to grab the ORFs when I could, and settle for the nucleic acid sequences when I had to.


*Collect 5M Bacteroidetes hgcA+ bins*

I first grabbed the Bacteroidetes bins from the 5M folder.
These already had the ORFs predicted, so I was able to just snag those.

*References selected for 5M project*

I had a bunch of them listed from the 5M paper, and the taxonomy ID's can be found here: `~/5M/dataEdited/binAnalysis/phylogeny/PVC/reference_taxonomy.tsv`.
The ORFs are here: `~/references/genomes/ORFs`.
Let's copy over the Bacteroidetes ORFs, but just the ones from RefSeq.
These are mostly Bacteroidales bins, but include two Flavobacterales as an outgroup.
There was a third outgroup member that we don't have a genome for, so we'll remove that one (sp001027725).

*hgcA+ Bacteroidetes from McDaniel et al, 2020*

I also wanted to include a range of hgcA+ Bacteroidetes genomes that Elizabeth identified in her study.
I manually selected a number of genomes from Supplementary Table 1 from her paper (subset of that data here with a column added to denote genomes I selected: `dataEdited/binning/phylogeny/bacteroidetes/mcdanielEtAl2020_bacteroidetes_bins.xlsx`).
List of accession numbers of selected bins is here: `dataEdited/binning/phylogeny/bacteroidetes/hgcA_bacteroidetes_genome_list.txt`.
I then downloaded the needed genomes using Batch Entrez: https://www.ncbi.nlm.nih.gov/sites/batchentrez.
I retrieved them as assemblies.
I couldn't download them all as protein faa files, so I just got the genomic fna files.
Two had RefSeq entries, which I downloaded, otherwise grabbed the GenBank entries.
I added them to project reference folder here: `references/genomes/bacteroidetes`.
I then unzipped them and cleaned the files (see scripts sheet).
I then uploaded them to the GLBRC.

Once in GLBRC, I predicted the ORFs.

*Add bins from this project*

Only anvio_hgcA_0130 to be included here.

*Add reference bins identified from NCBI*

After this tree was generated the first time, I went back looking for closer references to solidify where the bin was within Prolixibacteraceae (`code/binning/phylogenies/prolixibacter_tree.md`).
I found a few additional reference genomes to use (GCF_009617855.1, GCF_009617915.1, GCF_009617875.1, GCF_009617895.1).
There were two other ones (GCF_000621705.1 and GCF_003014495.1) that were already in the tree, but not accounted for in this document.
So, we'll include those too, for posterity.
I downloaded these six genomes using Batch Entrez: https://www.ncbi.nlm.nih.gov/sites/batchentrez.
File name list here: `references/genomes/bacteroidetes/NCBI/GTDB_Prolix_gene_list.txt`.
Hmm, Batch Entrez not working for this, will need to manually download them.
Downloaded them here: `references/genomes/bacteroidetes/NCBI/`.



*Set up dataframe for renaming tips*

I manually made a tsv file with a name for the bin, the accession number, the phylum, and the name of the tip label: `dataEdited/bins/binAnalysis/phylogeny/bacteroidetes/tip_naming.tsv`.



**Search for rp16 genes**

*Concatenate all bins*

First I concatenated all the bins into a single ORF faa file.
Then I generated a gene to bin file that linked each gene ID to its respective bin.

*Search for all genes*

I searched the concatenated ORF file for each of the rp16 genes using HMMs with hmmsearch (v3.3.1).
I then pulled out the amino acid sequences for the rp16 genes and aligned them.

*Rename fasta headers with bin name*

I then renamed the fasta file headers with the name of the bin so that they could be concatenated later in Geneious.

*Generate total alignment*

I then downloaded all the clean alignments to my local computer and imported them into Geneious.
The rpL5 gene has duplicate hits in GCA_004526055.1.
They seem to be two halves one sequence, so I just concatenated them.
I then concatenated all the alignments and masked them at 50% gaps.
I exported this as `rp16_alignment_masked.afa` and uploaded it to the GLBRC servers: `~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/Bacteroidetes/tree_building`.


**Identify hgcA sequences in bins**

I wanted to know which of my reference sequences had hgcA, so I searched the ORFs using the custom HMM.
I pulled out the ORFs, aligned them, and downloaded them to my computer for manual inspection.
Checked them in Geneious, they all looked good, so I made a list of the bins with hgcA.


**Generate tree**

I first generated a tree using FastTree, just to inspect it.
I checked it out in R: `code/phylogenies/bacteroidetes_tree.R`.
There were a few bins that I had thought might be missing too many of the rp16 proteins, but they all fit within the tree fairly well, so we'll stick with this alignment.

I then generated a maximum likelihood tree using RAxML (v8.2.11).
The R scripts to prepare this tree are also here: `code/phylogenies/bacteroidetes_tree.R`.

I then renamed all the branches using this file: `dataEdited/binning/phylogeny/bacteroidetes/tip_naming.tsv`.
I added the hgcA information to this vector (hgcA+ bins get two asterisks).
I also generated a color vector to color code the bins by source.
Doing this, I realized the Rhodothermia shouldn't be in there, they're far from my bins and kinda messing with the tree.
I removed these two: GCA_003557245.1 and GCA_003566645.1.
