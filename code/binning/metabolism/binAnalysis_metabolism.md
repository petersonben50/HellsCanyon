#### Bin analysis - metabolism

This is the protocol file for the various methods I'll use to analyze the hgcA+ bins I previously generated.


**Run METABOLIC**

I wanted to run the METABOLIC pipeline on the hgcA+ bins I collected.
First I transferred the files for the bins that we wanted to use and renamed them to use fasta as the file extension.
Ran METABOLIC using the docker image that Patricia generated with a submission file.
Github page here: https://github.com/AnantharamanLab/METABOLIC


**Custom set of metabolic HMMs**

I have a custom set of HMMs that I used to search through my bins to examine the metabolic potential of these organisms.
This was done with a custom script I generated to run multiple HMMs at once.

*Confirm dsrA phylogeny*

Most likely, all the *dsrA* genes that were identified were reductive *dsrA* genes, but I did want to confirm this.
I used the data from Karthik, from his sulfate-reduction paper, to generate a ML tree using FastTree.
I downloaded the `dsrA_phylogeny_trimmed.tree` file to my local computer (`dataEdited/bins/binAnalysis/metabolism/batch_HMMS/dsrA`), and analyzed that in R (`code/binning/metabolism/dsrA_tree.R`).
I used the Nitrospirae_bacterium_RBG_19FT_COMBO_42_15 and RBG_16_scaffold_151951 tip labels as markers for the branch leading to reverse *dsrA*.
None of the *dsrA* that we identified were reverse *dsrA* sequences, so we'll just stick with the lists of genes that we have.

**Search for MHCs**

I was also specifically interested in searching for multiheme cytochrome c proteins.
These seem to be good markers for respiratory vs. obligatory bins.
I used Shaomei's python script to identify the CXXCH motif in the ORF amino acid sequences.


**Search for BB-OMPs**

I was also interested in seeing if any of these were part of a porin cytochrome c complex (PCC) that might be used for Mn reduction in the system.
These PCCs have a beta-barrel outer membrane protein (BBOMP) next to a MHC.

*Pull out genes adjacent to MHCs*

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


**Check PCC_porin hits from batch HMMs**

Hmm, all except one have very low scores (29-31) range, e-values around 10^-6, 10^-7.
The one that's higher was predicted to be a PCC porin through our custom workflow.
Let's see if there are MHCs on either side of it, but I doubt there will be.
Yeah, none of them except for the ExtE from Geobacterales has adjacent MHCs, so these are probably not actually PCCs, except for the one identified above.



**Identify CAZymes in each bin**

I realized a lot of the bins we identified were oligate fermenters and METABOLIC predicted that they had a lot of genes used in complex-C degradation.
To investigate this further, I put the genes through the [dbCAN meta server](http://bcb.unl.edu/dbCAN2/blast.php) (accessed 2021-04-26).
I downloaded the ORFs for all the bins here: `dataEdited/bins/binAnalysis/bin_ORFs/ORFs.faa`.
This is a little too large, so I'll split it up into two files.
I uploaded both these files to dbCAN, and ran with HMMER (E-Value < 1e-15, coverage > 0.35).
I combined the output text into the file here: `dataEdited/binning/metabolism/GHs/cazyme_output.tsv`.
I checked out the output in R (`code/binning/metabolism/GHs.R`) and saved out a csv (`dataEdited/bins/binAnalysis/metabolism/GHs/clean_cazyme_data.csv`) with counts of each class of CAZyme within each bin.



**Characterize MoORs in bins**

Next I wanted to characterize the MoORs I found in the bins.
I used an HMM I developed for the 5M project to identify putative MoORs in the bin ORFs.
I then pulled out all the hits and aligned them using MUSCLE.
I downloaded this alignment (`dataEdited/binning/metabolism/MoORs/putative_MoORs.afa`) and manually inspected it in UGENE.
Looked decent, so we'll stick with them for now.
Might need to cut some out later though.
I aligned this to the set of reference MoORs I gathered for the 5M project.
I trimmed at 50% gaps using trimal.
I inspected this in UGENE too, and it also looked good for now.
I then generated a ML tree using FastTree.
I inspected the tree here: `code/binning/metabolism/MoORs.R`.
There are a number of sequences from the metagenome in this tree that are not closely related to the reference sequences.
I think I want to remove them, but I'd like to see how they did on the HMM first (if they did really well, it might be worth keeping them).
This is leading me down a rabbit hole of writing packages, since the `rhmmer` library seems to be breaking my libraries, and I should really be writing some functions to read in specific datasets.
Okay, I generated a stand-in function for this.
There's one branch that leads to sequences with low scores and doesn't include any reference sequences.
This one will be removed, since we cannot garner any additional info from these seqs without more extensive analysis.
I generated a list of the IDs for the sequences to save and uploaded that to GLBRC.

*MoOR tree take 2*

Run through the same process as before: align the seqs from this study, align that to the references gathered previously, and clean/trim the final alignment (trim at 50% gaps).
Let's check this in UGENE before generating a tree.
The alignment looks solid, let's go forward with the tree.


*Run kofamscan on ORFs*

It seems that using kofamscan and identifying the present pathways with KEGG-Decoder might be a more complete way to look at the metabolic pathways in our organisms.
So, I went ahead and ran KOFAMscan on all of my bins, individually.


**Aggregate bin metabolism data**

R script to aggregate the information is here: `code/binning/metabolism/metabolism_summary.R`.

Finally, I went through my metabolic gene data and manually assigned potential metabolic functions for each of the bins.
Saved it out here: `dataEdited/binning/metabolism/metabolism_summary.csv`.
I then saved this as a xlsx file and manually edited that one.
First I made a manual taxonomy column to make this easier to follow.
Then I went through what my workflow pulled out an assigned a metabolic function to each of these.
Notes are included in the sheet.
