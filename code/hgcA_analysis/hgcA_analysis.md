## Walkthrough of assembly-based hgcA analysis

This document outlines the assembly-based analyses I conducted to analyze the *hgcA* gene in the Hells Canyon metagenomes.
The initial data processing was done in the `code/hgcA_analysis/hgcA_analysis.sh` script.
Further processing and data visualization was done on Geneious and in R.


**Identify hgcA sequences**


*Identify putative hgcA genes with HMM*

I used the HMM I built for the 5M project to search through the open reading frames from the assemblies for HgcA-like sequences.
These results were stored in folders by year.
I then pulled out the amino acid sequences for each of the putative HgcA sequences.


*Concatenate and align all hgcA seqs for curation*

I aligned all the identified hgcA sequences to the HMM using hmmalign, then converted the file to the fasta format.
I downloaded this file (`dataEdited/hgcA_analysis/hgcA_raw.afa`) and loaded it into Geneious.
Here are my notes on the curation of those sequences:

- fall2017cluster5_000000264767_1 is truncated at the C-terminus. Let's hang on to it for now. All seqs truncated like this will be added to this list: `dataEdited/hgcA_analysis/identification/truncated_hgcA_seq.txt`
- fall2017cluster6_000000245170_1 also is truncated. Also, let's keep this one around for now.
- fall2017cluster6_000000208658_3 doesn't have the cap helix domain. Cut it!
- fall2017coassembly_000001368876_2 truncated at C-terminus.
- fall2017coassembly_000000596207_3 truncated at C-terminus.
- fall2017coassembly_000000643563_1 truncated at C-terminus.
- fall2017coassembly_000000746684_1 missing cap helix domain
- fall2017coassembly_000000746684_1 truncated at N-terminus
- fall2017coassembly_000001272590_2 truncated at C-terminus, very close to cap helix domain.
- HC18HY300_000000198900_1 truncated at C-terminus, very close to cap helix domain.
- HC18HY300_000000087507_1 truncated at C-terminus, very close to cap helix domain.
- KMBP004E_000000161924_1 truncated at C-terminus, very close to cap helix domain.
- KMBP004F_000000165420_1 truncated at C-terminus
- KMBP004F_000000644814_3 truncated at N-terminus
- KMBP004F_000000253505_6 truncated at N-terminus

Exported `hgcA_good.afa` and uploaded it to the GLBRC.


**Pull out depth of hgcA+ scaffolds**

I then extracted the depths of all the hgcA+ scaffolds from the metagenomes.
I used samtools for this.
When aggregating the depth over a scaffold, I took the mean of the coverage, not including the 150 bp at either end of the scaffold.
Wonder if I should calculate the median instead?


**Genomic context for hgcA**

Next I pulled out the hgcA+ scaffolds and the corresponding GFF files.

*Search downstream genes for hgcB*

I pulled out the gene immediately downstream from ecah *hgcA* using `retrieve_downstream_gene_name.py`.
I ran those ORFs against the *hgcB* HMM I built for the 5M project.
I aligned the hits against the HMM and checked them in Geneious.
These all look good for the most part. fall2017cluster4_000000001920_1 and HC18HY300_000000119163_1 are a little truncated, but the motif we need to see is there (or at least mostly).
We'll stick with the output we have there.

I then manually checked the seven non-hgcB sequences using [MOTIF](https://www.genome.jp/tools/motif/):
  - fall2017cluster6_000000000428_49 - PF02627 is the only hit. Carboxymuconolactone decarboxylase. Not a great hit.
  - fall2017coassembly_000000448578_3: No motif found. Checked it on blast: Closest was hypothetical protein Bacteroidales bacterium.
  - fall2017coassembly_000001334838_1: Voltage-dependent anion channel, super low support for this. Probably hypothetical.
  - HC18HY300_000000013755_6 - this matches to a sodium bile symporter... weird. Pretty good match too, e-value of 5.9e-51. Very weird
  - HC18ME02_000000018262_3 - This one had the identical sequence to the one above. Must be the same sequence, just assembled two different times.
  - KMBP004F_000000216801_1: Nothing with high support. BLAST: hypothetical. Closest hit is 37% identity.
  - KMBP009B_000000084934_2: transposase nuclease. Interesting. Wonder if this is a transposon?

*Isolate gene neighborhoods*

I extracted the region 5000 bp on either side of the hgcA sequences.
I loaded these scaffolds into Geneious and took a closer look at the scaffolds that didn't have a downstream hgcB.

*Manually inspect scaffold for hgcB*

- fall2017cluster6_000000000428_49 - There a ~350 bp gap between *hgcA* and the next predicted ORF. I copied the downstream DNA info and blasted it on NCBI using blastx (starting with ATGAAGATGC). This returned hits to a verrucomicrobial 4Fe-4S binding protein, which is likely *hgcB*, so we'll call this one hgcB+.

- fall2017coassembly_000000448578_3: This one was directly downstream from *hgcA* and in the right spot. No downstream hgcB here.

- fall2017coassembly_000001334838_1: The identified ORF is in the opposite direction as *hgcA*. Ran blastx on the downstream DNA...  No similarity found. No hgcB.

- HC18HY300_000000013755_6 - There's a gap in the ORFs downstream of this *hgcA*, about 300 bp long. Started a few bp down, at ATGGAATGCG. Blast comes back as 4Fe-4S binding protein. Definitely *hgcB*

- HC18ME02_000000018262_3 - Same *hgcA* sequences as HC18HY300_000000013755. Scaffold looks similar too, with the downstream gap. Blast also identified a 4Fe-4S binding protein downstream.

- KMBP004F_000000216801_1: ORF downstream is reverse compared to *hgcA*. Not much of gap. Blasted the downstream nucleotides in all three same-strand frames. Nothing. Likely to be true non-*hgcB* containing

- KMBP009B_000000084934_2: Nothing here either, tried blasting the downstream DNA, which was included in another ORF.

So, to summarize, we have 4 sequences with truly no downstream *hgcB* gene: fall2017coassembly_000000448578_3, fall2017coassembly_000001334838_1, KMBP004F_000000216801_1, KMBP009B_000000084934_2.
The other three have hgcB, they're just not predicted as such.

In the `hgcA_data_aggregation.R` file, I read out a csv file (`hgcB_notes.csv`) with the *hgcA* and *hgcB* sequences linked by their scaffold IDs.
I saved this out as an xlsx file and added a notes column with notes from my manual inspection of the gene neighborhoods.


**Classify hgcA seqs with pplacer workflow**

I used a workflow developed by Caitlin Gionfriddo to classify the *hgcA* genes, so all credit to her for this really nice workflow.
A very helpful tutorial can be found here: https://caitlingio.com/tutorial-for-hgcab-amplicon-sequencing-data/.
It relies on the Hg-MATE database (I used version 1 for this analysis).
First is the generation of an alignment of the HgcA sequences from this study using MUSCLE.
Then, consensus alignment is made between this alignment and the alignment of the Hg-MATE database.
The sequences are placed on a pre-generated tree using pplacer.
The program rppr is used to make a sqlite database of the taxonomy of the reference sequences.
Guppy was then used to classify the sequences and visualize the placements on a tree.
Finally, I used the R script `clean_hgcA_classification.R` (also adapted from Caitlin's well-annotated workflow) to clean up the classification and save out a csv file with that information.


**hgcA dereplication**

Okay, did an initial dereplication with CD-HIT with all the sequences, to see if there was any cross-over across different years.
There is... now how do I account for this?

I am including all the *hgcA* sequences in the final hgcA table, with a column for overall representative and another key with the year that that sequence will be counted in.
Essentially, I picked a representative from each cluster overall and one for each year.
I aggregated all the data here: `hgcA_data_aggregation.R`.
Save out a csv here: `dataEdited/hgcA_analysis/hgcA_dereplication_data.csv`.
Save this as an Excel file `hgcA_dereplication.xlsx`.
Add a column (representative) and assign a TRUE or FALSE depending if this gene is to be a representative of its cluster.
Also add a column (usedForAbundance) to filter out the sequences that will be used for abundance calculations.
These two columns will be similar, except that the usedForAbundance column will include multiples of some clusters, for clusters that span multiple years.
I am not including sequences that are truncated in the abundance calculations though (i.e., cluster 13 has an *hgcA* from 2017 and 2018 samples, but the 2017 sample is truncated).
Also not including the truncated samples in the tree.


**Phylogenetic analysis of hgcA**

Lastly, I wanted to generated a phylogenetic tree with representatives from each of the HgcA clusters.
I first pulled out the sequences that I determined I would use in the tree (1 from each cluster).
I then aligned these sequences to each other, then generated a consensus alignment with the Hg-MATE reference database.
I then used FastTree to generate a rough HgcA tree.
I looked at this tree in R (`clean_hgcA_tree.R`).
I first looked at the unrooted tree to find the branch leading to the paralogs, then rooted the tree there.
I read out a pdf of the tree and identified nodes that led to tips that were not needed.
I removed all the sequences under these nodes and read out a list of the sequence names that I wanted to include.
I then used a custom script to pull out the sequences that I didn't want to include.
I masked the alignment at 50% gaps in Geneious.

*Generate tree in RAxML*

I then used this alignment to generate a maximum-likelihood tree using RAxML.
There were a few duplicated sequences, so I used the `.reduced` file that RAxML generated.


*Generate good tree with subset of references using RAxML*

The HgcA seqs are concentrated enough into a few clades, so I'll just manually select sequences to use for a final RAxML tree.
I'll also include all the hgcA sequences from the 2017 Mendota study, so won't take those ones from Hg-MATE.
Hg-MATE names of references to be used can be found here: `dataEdited/hgcA_analysis/phylogeny/reference_names_to_use.txt`.
I used this list to pull out the sequences I needed, locally.
Needed to do a couple of iterations, since there were a few sequences that were in the refpackage but not in the database.
I think Caitlin had updated the refpackage separately or something.
Either way, good to go now.
I uploaded the references to GLBRC.

I also pulled the references from my study and from Jones et al, 2019, and the paralogs Caitlin uses.
Just snagged them from the HCC folder, and removed a few that I knew would not be relevant.
I concatenated all these sequences and aligned them using MUSCLE.
I downloaded the alignment and inspected it in Geneious to ensure it looks good.
Which it does.
I then masked the alignment at 50% gaps, which seems to work well.
I exported it (`hgcA_for_tree_final_masked.afa`) and uploaded it to the GLBRC.
I then ran RAxML (v8.2.11) on it to generate a ML tree.
I used rapid bootstrapping with automatic detection of limits and autodetection of the mutation model.
