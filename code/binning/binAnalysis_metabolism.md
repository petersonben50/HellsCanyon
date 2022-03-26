**Run METABOLIC**

I wanted to run the METABOLIC pipeline on the hgcA+ bins I collected.
First I transferred the files for the bins that we wanted to use and renamed them to use fasta as the file extension.
Ran METABOLIC using the docker image that Patricia generated with a submission file.

























**Metabolic analyses**

However, the focus will be on the abundances in the assembly, so the binning will just be used for metabolic pathway analysis and taxonomy/phylogeny.
Meaning, I want to just use the bins to inform the assembly-based analysis.
Thus, I ran the metabolic analyses on all the bins to assign a general metabolic strategy to the high matching sets.
I did the same thing with taxonomy.

*Custom set of metabolic HMMs*

I have a custom set of HMMs that I used to search through my bins to examine the metabolic potential of these organisms.

*Search for MHCs*

I was also specifically interested in searching for multiheme cytochrome c proteins.
These seem to be good markers for respiratory vs. obligatory bins.
I used Shaomei's python script to identify the CXXCH motif in the ORF amino acid sequences.

*Search for BBOMPs*

I was also interested in seeing if any of these were part of a porin cytochrome c complex (PCC) that might be used for Mn reduction in the system.
These PCCs have a beta-barrel outer membrane protein (BBOMP) next to a MHC.
I pulled out all the sequences that were adjacent to the MHCs, then searched them with a custom BOMP HMM.
Looks like I only had one substantial hit, from that Geobacterales in KMBP009B.

*Generate a phylogeny*

I wanted to see what type of PCC this might be.
First I BLASTed the refseq database for 5 hits to this sequence.
I then pulled out these sequences and dereplicated them against the PCC reference set that we have.
I then concatenated all the references and our one PCC sequence and aligned them.
I trimmed the alignment at 50% gaps using trimal, then generated a tree using FastTree.
I downloaded this tree and the refseq protein list to my local computer.
I used Entrez to pull out the taxonomic information associated with the reference sequences.
Finally, I checked out the tree in R: `code/binning/metabolism/bbomp_tree.R`.
Looks like this one gene is similar to ExtE.
The corresponding MHC is KMBP009B_000000035879_3, which also has another MHC adjacent to it (KMBP009B_000000035879_2).
I ran both of these through the CELLO localization software here (http://cello.life.nctu.edu.tw/cello2go/alignment.php).
KMBP009B_000000035879_2 is predicted to be extracellular and KMBP009B_000000035879_3 is predicted to be periplasmic, which also support the fact that this is a true PCC.
This is from bin anvio_hgcA_0210.


*Check PCC_porin hits from batch HMMs*

Hmm, all except one have very low scores (29-31) range, e-values around 10^-6, 10^-7.
The one that's higher was predicted to be a PCC porin through our custom workflow.
Let's see if there are MHCs on either side of it, but I doubt there will be.
Yeah, none of them except for the ExtE from Geobacterales has adjacent MHCs, so these are probably not actually PCCs.


*Confirm dsrA phylogeny*

I wanted to confirm that the dsrA sequences we identified were reductive dsrA genes and were not the reverse dsrA sequences.
I only had two sequences, so shouldn't be too hard.
I aligned these two to each other, then aligned that to karthik's *dsrA* database, using MUSCLE for both.
I then used trimal to trim up the alignment to remove residues with over 50% gaps.
I then generated a ML tree using FastTree.
I inspected this tree in R: `code/binning/metabolism/dsr_tree.R`.
Both of the dsrA genes we identified are reductive *dsrA* sequences.
We'll just leave this as is in the spreadsheet then.


*Run through METABOLIC*

I set up a METABOLIC run on the GLBRC servers using a submission script that Charles had set up for me.


*Run kofamscan on ORFs*

It seems that using kofamscan and identifying the present pathways with KEGG-Decoder might be a more complete way to look at the metabolic pathways in our organisms.
So, I went ahead and ran KOFAMscan on all of my bins, individually.


*Run genomes through Charles's workflow*

I sent my genomes to Charles and he ran them through his EET workflow.
I downloaded the output from that.


*Identify CAZymes in each bin*

I realized a lot of the bins we identified were oligate fermenters and METABOLIC predicted that they had a lot of genes used in complex-C degradation.
To investigate this further, I put the genes through the [dbCAN meta server](http://bcb.unl.edu/dbCAN2/blast.php) (accessed 2021-04-26).
I downloaded the ORFs for all the bins here: `dataEdited/binning/bin_orfs/ORFs.faa`.
This is a little too large, so I'll split it up into two files.
I uploaded both these files to dbCAN, and ran with HMMER (E-Value < 1e-15, coverage > 0.35).
I combined the output text into the file here: `dataEdited/binning/metabolism/GHs/cazyme_output.tsv`.
I checked out the output in R (`code/binning/metabolism/GHs.R`) and saved out a csv (`dataEdited/binning/metabolism/GHs/clean_cazyme_data.csv`) with counts of each class of CAZyme within each bin.
I then added this to a separate Excel sheet in the `metabolic_summary.xlsx` file.

*Characterize MoORs in bins*

Next I wanted to characterize the MoORs I found in the bins.
I used an HMM I developed for the 5M project to identify putative MoORs in the bin ORFs.
I then pulled out all the hits and aligned them to the HMM.
I downloaded this alignment (`dataEdited/binning/metabolism/MoORs/putative_MoORs.afa`) and manually inspected it in Geneious.
I cut out about 19 of them that didn't match up well.
The rest we'll keep for now, save out new version: `dataEdited/binning/metabolism/MoORs/putative_MoORs_cut.afa`.
I aligned this to the set of reference MoORs I gathered for the 5M project.
I trimmed at 50% gaps using trimal, then generated a ML tree using FastTree.
I inspected the tree here: `code/binning/metabolism/MoORs.R`.
This doesn't look too bad actually, let's go ahead with a RAxML tree from this same alignment.


*Notes on metabolic assignment*

Finally, I went through my metabolic gene data and manually assigned potential metabolic functions for each of the bins.
R script to aggregate the information is here: `code/binning/metabolism/metabolism_summary.R`.
Saved it out here: `dataEdited/binning/metabolism/metabolic_summary.csv`.
I then saved this as a xlsx file and manually edited that one.
First I made a manual taxonomy column to make this easier to follow.
Then I went through what my workflow pulled out an assigned a metabolic function to each of these.
Notes are included in the sheet.


**hgcA analysis in bins**

First, I'm mostly interested in linking the bins to the HgcA phylogeny.
I searched through the bin ORFs for HgcA sequences.
I identified 19 of them (which is good, that's what I expected) and saved out a list of them.
I then used the G2B file to link the hgcA sequences to the bins I generated.


**Depth analysis**

I extracted the depth of each of the bins using a submission script.
I cleaned this data here: `/Users/benjaminpeterson/Documents/research/HellsCanyon/code/binning/bin_depth_aggregate.R`.
I used the SCG abundance to normalize the bin coverage across metagenomes.
I plotted the coverages of the bins here: `code/binning/depth_plots.R`.