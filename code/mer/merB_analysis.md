### Notes on analyzing merB genes in the HCC assemblies

**Identify merB sequences**

We used the PFAM HMM to identify the putative MerB sequences from each of the assemblies.
We then pulled out all the amino acid sequences and aligned them to the HMM.
I downloaded this alignment and looked at it in Geneious.
It seems that we have different types of sequences here, but that all align fairly well.

First we have sequences that are around 400 bp long.
The second half of the sequence is what seems to line up with the HMM, although we should probably align the sequences against a well-characterized sequence so we can check for the cysteines.

*Concatenate and align all merB seqs with references for tree generation*

Actually, let's make a tree to compare how these sequences fall out.
I concatenated the sequences along with the merB reference sequences, then aligned it to the HMM.
I masked the alignments in Geneious at 50% gaps, then ran FastTree on the masked alignment.

I inspected the tree in R: `code/mer/merB_tree.R`.
Looks like the vast majority of our sequences are on a couple of deep branches.
I pulled out the tip labels from all of the groups that had sequences of interest in them, and saved them in a separate list.

*Inspect individual clusters*

I then looked at the alignments of each of these groups individually.
Downloaded them to my local computer and inspected them in Geneious.
Let's look for those cysteines. The first one (Cys96) is at position 849 in our alignment. Cys117 is at 881. Cys159 is at residue 975.

Let's make a tree with these clusters and add in the alignments.
First, trim the aligment. Removed first 526 residues.
Also removed last 52 residues.
Then mask it at 50% gaps again.
Run FastTree to get a tree.
Load this into R, then use ggtree to plot it together with the alignment: `/Users/benjaminpeterson/Documents/research/HellsCanyon/results/mer/merB_tree_subset_with_alignment.pdf`.


- Divergent cluster 1: This cluster, after plotting just the subset, does cluster in the middle of the other clusters. However, it is missing the Cys96 residue. They've been replaced with G, S, or A. We will not count these as merB.
- merB cluster 1: Cys96 is there. Cys117 has been converted to D. Cys975 is there as well. These include all the long sequences, seem to be some sort of fused merB gene. Interesting, the fused part seems to encode a thioredoxin of some kind. Manual inspection reveals a number of CXXC motifs in the added part of each of the fusions sections.
- merB cluster 2: Many of these are the truncated ones, which are missing the cysteine at Cys159. First I removed the sequences with no cysteine at Cys96. Then I removed the truncated ones with no cysteine at Cys 159. There was a small subset of proteins here with no cysteine at the usual Cys117, but we'll keep them for now, but in a separate group from merB_2
- merB cluster 3: This is the one I have the most confidence in, as it clusters most closely with the references. The two seqs are a bit truncated relative to the references. But, it does have all the requisite cysteines. I think these two seqs will be clustered together, one is from coassembly and one from cluster 3.

Alright, so now we have four groups of merB.
1. merB_fused: These are the long fused sequences
2. merB2: These are the sequences from group 2 that have all cysteines intact.
3. merB2_C117I-L: These are the sequences from group 2 that have had the Cys117 residue replaced with either I or L.
4. merB3: These are the two well-conserved merB sequences.

Upload these lists to GLBRC


**Pull out depth of merB+ scaffolds**

Next we'll extract the depths of the scaffolds with the *merB* genes.
We'll do this for all the genes that were identified using the HMM, rather than just the curated set.
Do this using my usual method, with looking at the coverage of the scaffolds excluding 150bp from either end.



**Genomic context for merB**

Let's pull out the gene neighborhoods of these genes.
Eventually, we should write a script to search all the proteins within 10kbp or so for other mer genes.
When I first did this analysis, I noticed that some scaffolds had multiple hits for *merB*.
I needed to deal with this, so I searched for these scaffolds.
Each of the pairs of hits were next to each other on the scaffolds.
Interesting.
None of these duplicated ones are in the groups we kept.
Of the pairs I checked out, one of the pairs was in the tree, and they were always the heavily truncated ones. I'm willing to bet that the second ones will be the very truncated sequences that align with the C-terminus end.

Interesting, these might require their own analysis.
Let's remove them from the current neighborhood analysis.

Then, I pulled out the GFF files and the scaffolds that had *merB* hits on them.

*Isolate gene neighborhoods*

Then I used the scripts that I modified from Tyler Barnum's scripts to pull out just the gene neighborhoods, 8000bp to either side of the gene.
I then concatenated all the results.

*Pull out amino acid sequences*

Next, I needed to pull out the amino acid sequences of the nearby sequences so that I could identify them.
I saved out these files here: `dataEdited/merB_analysis/scaffolds/friendly_neighborhood_genes.faa`

*Run Kofamscan on ORFs*

I then wanted to annotate all these ORFs using KOFAMscan (v1.3.0).
I did this using the same workflow that I used for HeCaMT.

*Generate gene neighborhood images*

Finally, I generated gene neighborhood images for each of the clusters of sequences, using this script: `extract_viz_table_from_gff.py`
I wanted to add in colors to correspond to the additional *mer* genes on the scaffold.
Got a list of the KO families here: https://www.kegg.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=kegg&dbkey=orthology&keywords=mercuric&page=1
We'll make *merA* genes red,
