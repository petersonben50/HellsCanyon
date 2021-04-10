### Notes from the binning workflow

I will be doing two separate binning workflows for two separate purposes.
First, I will manually bin hgcA+ genomes using anvi'o.
These bins will be used to investigate the metabolic potential of the methylating organisms.
Second, I will generate an uncurated set of bins made through automatic binning algorithms.
I won't look at the detailed metabolic pathways of these organisms, but will use them to look at the broad taxonomic composition of the community and observe general trends in metabolic pathways.









**Prepare anvi'o databases for manual binning**



**Manually bin hgcA+ bins**

After doing the binning, I made a tsv to link the bin names to the assembly they were associated with: `/Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/binning/manualBinning/notes/bin_to_assembly.tsv`.
I then uploaded this to `~/HellsCanyon/dataEdited/binning/manualBinning/binsGood`.


**Summarize/export curated bins**

I summarized our curated bins from anvi'o first.
Then I pulled out the DNA sequences of the scaffolds and made a list of the bins.


**Check quality of bins**

*Completeness/redundancy estimates from anvio*

First I pulled out the completeness/redundancy values that anvi'o generates.
Then I concatenated these all into a single file.
Finally, I identified hgcA bins with completeness over 50% and redundancy less than 10%.
Saved this information locally here: `dataEdited/binning/quality/bins_summary_hgcA.txt`

*Completeness/redundancy estimates from CheckM*

I then ran CheckM on the bins from anvi'o using the lineage_wf and the qa workflows.
I then pulled out the necessary information.
Saved this information locally here: `dataEdited/binning/quality/checkM_stats.csv`
I also then pulled out all the bin IDs that had a completeness of over 50% and a redundancy of less than 10%.

*Final bin list*

I used the list of the bins that met the 50% C, 10% R criteria for either CheckM or anvi'o.
I'll sort out later which of these we actually want to use, but for now I'll include them in all the analyses.
I pulled out the scaffolds for these bins and cleaned them.


**Check out taxonomy of bins with GTDB**

I then ran the classify_wf workflow in GTDB-TK to generate a taxonomic classification for each of the high qualified bins we identified.
I downloaded this information here: `dataEdited/binning/taxonomy/gtdbtk.bac120.summary.tsv`.


**Pull out coverage information from anvio**

I then pulled out the coverage information from anvi'o for use in the differential coverage analysis, making sure that the bins in our HMSs have a similar distribution.
I pulled out that coverage information then concatenated it all together.
These are split up by year, since the mapping was split out by year.
I downloaded the `coverage_goodBins.txt` file to my local computer.


**Get ORFs for bins**

I used Prodigal to predict the open reading frames for each of the genomes.
We saved out the nucleic acid and amino acid sequences of all the genes.
I also saved out a GFF file for each bin.
I then cleaned up each of these files and the scaffold files.
I concatenated them and generated a scaffold to bin and gene to bin file, and concatenated the sequence files.


**Run ANI comparisons on good bins**

I needed to dereplicate across the years and the different assemblies.
I used Sarah's DAGman ANI workflow for this.
I downloaded the `group.sub` file to my computer to edit it: `/Users/benjaminpeterson/Documents/research/HellsCanyon/code/binning/ANI`
I then ran the script to compare all these bins to each other.
I then downloaded `goodBins.all.ani.out.cleaned` here: `dataEdited/binning/ANI/goodBins.all.ani.out.cleaned`


**Aggregate bin data**

The code for this can be found here: `code/binning/aggregate_bin_data.R`
First I took at look at the ANI values.
Looks like we can use 97% ANI and 50% coverage and capture the bins that are close together, so let's go with that.
I got 4 HMSs, two with three bins and two with two bins.
So all told, it looks like I ended up with 13 sets of bins.
Eight were from 2017, 3 from 2018, and 2 from 2019.
Not very many, really, when you consider how much sequencing went into all this.
No HMSs including bins from multiple years.
I was going to look at the differential coverage and ordinate the bins by that, but didn't think we had enough HMSs to warrant it.



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


*Notes on metabolic assignment*

Finally, I went through my metabolic gene data and manually assigned potential metabolic functions for each of the bins.
R script to aggregate the information is here: `code/binning/metabolism/metabolism_summary.R`.
Saved it out here: `dataEdited/binning/metabolism/metabolic_summary.csv`.
I then saved this as a xlsx file and manually edited that one.
First I made a manual taxonomy column to make this easier to follow.
Then I went through what my workflow pulled out an assigned a metabolic function to each of these.
Notes are included in the sheet.


**Depth analysis**

I wanted to check the depth profiles of our bins and was too lazy to write up a script so I used the depth calculations from anvi'o.
I looked at this here: `code/binning/depth_plots.R`.
This is a pretty rough analysis, since the normalization scheme of anvi'o is pretty confusing, so I'll want to extract the depths myself at some point to make sure I've got this right.
Nice way to see the relative abundance of these different groups though.


**hgcA analysis in bins**

First, I'm mostly interested in linking the bins to the HgcA phylogeny.
I searched through the bin ORFs for HgcA sequences.
I identified 19 of them (which is good, that's what I expected) and saved out a list of them.
I then used the G2B file to link the hgcA sequences to the bins I generated.
