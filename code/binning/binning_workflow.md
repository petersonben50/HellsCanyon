## Notes from the binning workflow

I will be doing two separate binning workflows for two separate purposes.
First, I will manually bin hgcA+ genomes using anvi'o.
These bins will be used to investigate the metabolic potential of the methylating organisms.
Second, I will generate an uncurated set of bins made through automatic binning algorithms.
I won't look at the detailed metabolic pathways of these organisms, but will use them to look at the broad taxonomic composition of the community and observe general trends in metabolic pathways.

Obsidian notes: "HCC binning strategy"


**Prepare scaffolds and mapping files**

I only did the binning on scaffolds at least 2000 bp in length.
To do this, I generated a new assembly file to only include these scaffolds, then I remapped the reads to these new assemblies.
Mapping was done with bowtie2 with a submission file.


**Run automatic binning algorithms**

First, I automatically generated bins using Metabat2 and Maxbin2, which were then consolidated using Das Tool.
This was all done in a submission file.
Little messy, since I am changing a previous workflow, but the final DNA sequences of the auto-generated bins should be here: `~/HellsCanyon/dataEdited/binning/autoBinning/completeBinSet/DNA/`.


**Initial bin characterization of auto-generated bins**

I wanted to do some quick bin characterization before moving on.

*Get quality bins*

Das Tool generally removes super low quality bins, but I wanted to check them anyways.
I ran the `checkm lineage_wf` workflow.
All bins that were above 50% complete and less than 10% redundant were included in the next step of the analysis (moved to `hqBinSet`).
This bin set will be used for the remainder of the initial bin characterization.

*Run GTDB*

I used the GTDB-TK `classify_wf` to predict the taxonomy of the bins.


*Run ORF prediction on all bins*

I used Prodigal to predict the ORFs for all the bins, done on single mode.
The `cleanFASTA.py` was then used to clean them up.
Finally, I concatenated all the ORF files and generated a scaffold to bin file.

*Check for hgcA gene in bins*

I used my custom HMM to identify bins with putative hgcA sequences.
This identified 17 bins with *hgcA*, all sequences of high quality (lowest score was 324).
All bins were from 2017 and 2018.

*Run ANI comparisons on good bins*

Then I compared the final set of bins using Sarah Stevens's ANI calculation workflow (written as a DAG script), which is done using JGI's ANI Calculator.
This resulted in the `hqBins.all.ani.out.cleaned`.
*Will not be used in final analysis for this paper*

**Dereplicate bins**:

*Not used in final analysis. Keeping this here for legacy purposes but I switched to dereplicating across both methods*.

I then dereplicated the auto-generated bins.
First I downloaded the ANI, checkM, GTDB, and *hgcA* presence data to `dataEdited/binning/autoBinning`.
These were combined in an R script: `code/binning/sandbox/aggregate_bin_data_autoBinning.R`.


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

*Check out taxonomy of bins with GTDB*

I then ran the classify_wf workflow in GTDB-TK to generate a taxonomic classification for each of the high qualified bins we identified.
I downloaded this information here: `dataEdited/binning/taxonomy/gtdbtk.bac120.summary.tsv`.

*Pull out coverage information from anvio*

I then pulled out the coverage information from anvi'o for use in the differential coverage analysis, making sure that the bins in our HMSs have a similar distribution.
I pulled out that coverage information then concatenated it all together.
These are split up by year, since the mapping was split out by year.
I downloaded the `coverage_goodBins.txt` file to my local computer.

*Get ORFs for bins*

I used Prodigal to predict the open reading frames for each of the genomes.
We saved out the nucleic acid and amino acid sequences of all the genes.
I also saved out a GFF file for each bin.
I then cleaned up each of these files and the scaffold files.
I concatenated them and generated a scaffold to bin and gene to bin file, and concatenated the sequence files.


*Run ANI comparisons on good bins*

I needed to dereplicate across the years and the different assemblies.
I used Sarah's DAGman ANI workflow for this.
I downloaded the `group.sub` file to my computer to edit it: `/Users/benjaminpeterson/Documents/research/HellsCanyon/code/binning/ANI`
I then ran the script to compare all these bins to each other.
I then downloaded `goodBins.all.ani.out.cleaned` here: `dataEdited/binning/ANI/goodBins.all.ani.out.cleaned`
*Will not be used in the final analysis of this paper*

**Aggregate bin data**

*This is not going to be used in the final manuscript, but leaving it here for... well, in case it's needed for some reason.*
The code for this can be found here: `code/binning/sandbox/aggregate_bin_data_manualBinning.R`
First I took at look at the ANI values.
Looks like we can use 97% ANI and 50% coverage and capture the bins that are close together, so let's go with that.
I got 4 HMSs, two with three bins and two with two bins.
So all told, it looks like I ended up with 13 sets of bins.
Eight were from 2017, 3 from 2018, and 2 from 2019.
Not very many, really, when you consider how much sequencing went into all this.
No HMSs including bins from multiple years.
I was going to look at the differential coverage and ordinate the bins by that, but didn't think we had enough HMSs to warrant it.


**Aggregate all hgcA+ bins**

Next I combined the hgcA+ bins from the manual and automatic methods.
Moved them all into `~/HellsCanyon/dataEdited/binning/bins_hgcA`.
I took the DNA, ORFs, taxonomy, and checkM files.

*Run ANI comparisons on hgcA+ bins*

I used Sarah Steven's workflow to check the ANI between all the hgcA+ bins that were generated through this study.
The `binsHgcA.all.ani.out.cleaned` was downloaded to my computer: (`dataEdited/bins/binning/bins_hgcA`).
Also downloaded the checkM and GTDB files.
Aggregated everything here: `code/binning/aggregate_hgcA_bin_data.R`
Save out a data file here: `dataEdited/bins/binning/bins_hgcA/bin_dereplication_data.csv`.
Convert this to a xlsx file for easier viewing.
We are going to want to avoid having duplicate genomes (like, both manual and autobinned from the same assembly).
Made a list of the ones to keep: `bins_hgcA_keepers_list.txt`.
Uploaded that to GLBRC

*Concatenate all ORFs*

We'll run subsequent analyses on all of these bins, so let's just concatenate all the ORFs.
Generate a G2B file as well.


**hgcA analysis in bins**

Lastly, I wanted to do a deeper dive on the *hgcA* genes from these bins.
I searched through the concatenated ORFs for HgcA sequences.
Wow, they all look good, scores over 320.
I pulled out all of the HgcA amino acid sequences, then aligned them.
Checked the alignment out in UGENE, and they all have the conserved cap-helix domain.
Think we're all set to go with that.
Put together tsv files that link the binID to the *hgcA* ORF ID, and the *hgcA* ORF ID to the representative *hgcA* ORF ID for the assembly-based anlayses.
