
**Identification of metabolic proteins**

*Run HMMs*

I first ran a set of HMMs from various sources against the ORFs from each assembly.
I'll detail the information about each of the HMMs that I end up using below.
Baseline thing to know is that we're using the internal trusted cut-off for each of the HMMs (`cut_tc`).
After identifying all the genes, I pulled them out and aligned them to the HMM.

*Extract depths of all scaffolds*

I then pulled out the depth of each scaffold that contained one of the identified genes using a submission script (`aggregate_depth_proteins.sub`)

*Dereplicate sequences*

I then used CD-HIT to dereplicate the genes within a single year (since that's how I split up the mapping).
I grouped them at 97% identity with the default alignment fraction.


**Verify and classify narG gene**

I used the reference dataset from Lüke et al, 2016 to generate a phylogenetic tree of the identified *narG* sequences.
First step was to build a database with the information from Luke et al.
I collected the fasta file for the blast reference dataset for the narG/nxrA gene (from the Supplementary Info of Lüke et al, 2016) and saved it as `narG_luke_database.faa`.
I extracted the header info from `narG_luke_database.faa` into `narG_luke_database.tsv`.
I also collected the "blast reference dataset to distinguish the narG gene and the nxrA gene of the Nitrospira/Nitrospina/ anammox group" and the "the blast reference dataset to distinguish the narG gene and the nxrA gene of the Nitrobacter/Nitrococcus/Nitrolancea group" into `nxrA_1.faa` and `nxrA_2.faa`, respectively.
I then pulled out the header information from these two, and combined it here: `nxrA.tsv`
I then manually generated an Excel file (`narG_luke_database.xlsx`) where I combined this information.

Second step was to build a phylogeny of the reference dataset and see if I could place the divide between *narG* and *nxrA*.
I cleaned up the fasta file, aligned the sequences using MUSCLE, then generated a ML tree with FastTree.
I looked at the tree here: `code/metabolic_analyses/narG_tree.R`, using the metadata sheet I generated previously to label the branches.
Hmm, this actually didn't come out very cleanly, the *nxrA* seqs from *Nitrobacter* are close to the nitrate reductases.
Need to further check this out.
I looked at the alignment, which looks fine.
I pulled out the hits against the NarG HMM and added them to the tree.
Also looked at this one in `narG_tree.R`.
All of the sequences aren't clustering with the *nxrA* reference sequences.
Pretty impressive how well the HMM works!
I'll go forward with using all the *narG* genes we identified with the HMM.
For completeness sake, I checked the alignment in Geneious.
There were a number of truncated sequences, but for the most part they looked good.
I trimmed the alignment to mask residues with 50% gaps and ran it through RAxML for a higher quality tree.


**Classify dsrA genes**



Alignment editing:
Upload dsrA alignment to the phylemon2 server.
Use a user-defined method, just mask residues that have 50% gaps.
Window size of 1, no minimum similarity threshold.
Download output file `outfile.out` and rename: `dsrA_phylogeny_masked.afa`.
Generate a tree in FastTree.

Look at unrooted tree in R.
Read out a pdf of it.
Look for root:


**Identify potential PCCs**


*Search adjacent genes for BBOMP*




After aligning the sequences to the HMM, download the sequences and manually inspect them in Geneious.
There are 40 sequences in here, most of the top hits are from the 2019 data.
Nothing really jumps out as not belonging, so might be best to build a tree here.

First I pulled out some reference sequences from the NCBI refseq database.
I then clustered these sequences at 97% identity, then removed ones with duplicated fasta headers.
Then, I dereplicated these refseq references against the reference sequences I had previously accumulated.
I then generated a list of the final set of refseq proteins and downloaded the taxonomic information for these genes.

*Generate a phylogeny*

I then concatenated the references and the study sequences into a single file and aligned them using MUSCLE.
I trimmed the alignment at 50% gaps using `trimal`, then generated a quick ML tree using FastTree.
I then looked at this tree in R: `code/metabolic_analyses/bbomp_tree.R`.

I midpoint rooted the tree, added in the reference names and saved it out into a PDF.
There were three clusters of sequences of interest, all the others were too divergent to argue that they were truly BB-OMP.
I subsetted the tree to divvy up the sequences into three different groups:

1. Cluster 1: ExtE-like
2. Cluster 2: Omb-like
3. Cluster 3: Brocadia-like

One tricky part is that these sequences have already been clustered, and clustered across the years.
To look at the depth, I need to include a sequence from each year that it was detected, since we only mapped within years.
To do this, I'll read in the clustering information and identify the BB-OMP cluster that are to be represented.
There's few enough of these that I'm just going to do it manually.
In the R script, I'll combine the dereplication files and the BB-OMP key to generate a table that has all the genes from each cluster.
Then I'll manually filter out the genes that we do not want to keep around.
