## Workflow for further examination of hgcA+ PVC bins


#### Collect PVC references

Scripts for this here: `code/binning/PVC_details/collect_PVC_references.sh`.
For all of these, I gathered the genome fna files (here: `~/HellsCanyon/dataEdited/binAnalysis/PVC_details/genome_dereplication/genomes`).
After gathering them all I predicted the ORFs.

**Dereplicate bins from this project**

I wanted to include both hgcA+ and hgcA- bins from HCC.
I needed to dereplicate these, sine the hgcA+ bins were manually curated, but the rest were autogenerated (which includes some hgcA+ ones).
I used Sarah's ANI workflow for this.
R scripts to clean up the data are here: `PVC_hgcA_minus_derep.R`.


**Collect 5M PVC hgcA+ bins**

All hgcA+ bins from the 5M project on Lake Mendota were included.

**Bins from Jones et al, 2019**

PVC bins from Dan Jones's 2020 paper.
I made a preliminary tree, and these were not closely related.
I kept them in, but they aren't very close.

**References selected for 5M project**

Known isolate genomes for both these phyla.
I picked out some to use with the 5M project, which were stored on the GLBRC servers here: `~/references/genomes`.
I copied over the Kiritimatiellaeota and Lentisphaerae ORFs.
Also include two Planctomycetes genomes, for an outgroup (GCF_000255705.1 and GCF_000025185.1).
In our initial tree, there weren't many good RefSeq or even GenBank genomes that were close to our bins.
I looked at the GTDB tree (`PVC_tree.R`)to confirm that there wasn't more genomes to be pulled here, but there weren't many, so I left the reference set as is.

**Additional genomes**

There were a few that I wanted to include that we didn't have previously.
My notes on collecting them are here.
All of them will be added to the reference folder: `~/references/genomes/bins/`

1. *Victivallis vadensis* genome in here, which we know is GCA_003096415.1. Let's download this from NCBI as well, here: https://www.ncbi.nlm.nih.gov/assembly/GCF_003096415.1. Saved it by its accession ID: `GCF_003096415.1_ASM309641v1_genomic.fna`.

2. I also wanted to include the two genomes from van Vliet et al, 2019, in which they isolated two Kiritimatiellaeota organisms. These are listed at ENA under these accession numbers: SAMEA5207384 and SAMEA5207385. I downloaded the genome fasta files from NCBI here:
  - Pontiella desulfatans: https://www.ncbi.nlm.nih.gov/assembly/GCF_900890425.1.
  - Pontiella sulfatireligans: https://www.ncbi.nlm.nih.gov/assembly/GCF_900890705.1

Moved all of these over, cleaned them up.


**Remove scaffolds from HC18HY300_bin_0036**

The tree analysis downstream identified two scaffolds that were present in two genomes (this is why I need to clean up my bin-generating pipelines).
The two bins are anvio_hgcA_0110 and HC18HY300_bin_0036.
We'll remove the two scaffolds (HC18HY300_000000021381 and HC18HY300_000000059556) from HC18HY300_bin_0036.
Now will continue as before.


**Removal of one genome**

Added 2022-06-17.
The genome GCA_001804865.1 is isolated in the tree.
I'm worried that it is interfering with tree generation, so I am going to remove it.
I won't re-run all the metabolic analyses, just re-generate the tree.


**Predict ORFs**

Then I ran the IMMA_ORF_STAN workflow to get the needed ORFs.
After it was done, I concatenated the ORFs and generated the gene-to-bin file.



#### Run genomes through GTDB

I ran GTDB on these genomes, just so we're not going through and manually assigning names.
Scripts here: `code/binning/PVC_details/metabolic_genes_PVC.sh`.
The results were downloaded to my local computer: `dataEdited/bins/binAnalysis/PVC_details/taxonomy/`.


#### hgcA analysis for PVC

**Identify hgcA sequences**

First I used the `protein_identification_and_alignment.py` script in HomeBio to identify potential *hgcA* sequences and saved out an alignment file.
I downloaded this alignment file to my local computer here: `dataEdited/bins/binAnalysis/PVC_details/hgcA_analysis/hgcA_hits/PVC_hgcA.afa`.
I then used the alignment viewer within DECIPHER to look at the alignment with `BrowseSeqs`.
R script is here: `code/binning/PVC_details/hgcA_alignment.R`
All these sequences look like true *hgcA* sequences.
Next I made a list and downloaded it to my local computer (`PVC_hgcA_geneID_list.txt`).
I also pulled out the relevant entries from the G2B file for these *hgcA* sequences (`PVC_hgcA_G2B.tsv`).


**Generate hgcA tree**

I then took the HgcA alignment I generated earlier and masked it at 50% gaps using Trimal.
Then, I used RAxML to generate a maximum likelihood phylogeny.
I downloaded the bipartitions (`RAxML_bipartitions.PVC_hgcA`) and the info (`RAxML_info.PVC_hgcA`) files to my local computer (`dataEdited/bins/binAnalysis/PVC_details/hgcA_analysis/tree_generation`).
This data was aggregated into a tree here: `code/binning/PVC_details/hgcA_tree.R`.
Note, the final version of this includes some bin information that was gathered in the "Tree generation for PVC" step, described below.


**Gene neighborhood analysis**

*Set up gene neighborhood analysis*:
Pull out the gene neighborhood data from around the *hgcA* genes using the `doctor_petersons_neighborhood.py` script.
I grabbed 20,000 residues to either side of the *hgcA* seqs.

*Run KOFAMSCAN on the neighborhood ORFs*:
I then ran kofamscan on the ORFs within the neighborhood.

*Cluster ORFs with CD-HIT*:
I then clustered all these ORFs with CD-HIT, using 0.40 as a cut-off.
Not sure that this is very helpful.

*Search for hgcB*:
One of the main targets was *hgcB*, considering its importance as the partner to *hgcA*.
I wanted to extract the downstream genes and check them with the *hgcB* HMM to see if they looked like *hgcB*.
When I first did this, I was getting some weird results, where the downstream gene was the right size, and clustered with *hgcB* genes, but wasn't annotated as such with the HMM.
Additionally, there were references that didn't hit the *hgcB* HMM despite me knowing they should have it.
But, I realized that the code we use to extract the downstream gene is reliant on the direction that *hgcA* is on the strand.
Since I had been using the gene neighborhood GFF file, which involves manipulating the strand sign, but the ORF numbers aren't changed, this meant I was extracting the wrong ORF.
To fix this, I need to use the original ORF GFFs.
To do this, I added a step to the ORF prediction standard output, in that I also concatenated the GFF files.
Kinda funny, the `retrieve_downstream_genes.py` script was broken anyways, but should be fixed now.

Doing it this way, it looks like there are three genes that are not marked as *hgcB* that are about the right size.
They are the downstream genes for GCA_002408885.1, LEN_0037, and KIR_0001.
Let's manually inspect these ones.
The gene ids for them are KIR_0001_19_36, LEN_0037_131_5, DHTZ01000076.1_16.
Well, right off the bat, DHTZ01000076.1_16 is assigned as a ferredoxin.
Let's align these three to the other hgcB sequences.
I downloaded this to my local computer and checked it in `hgcB_alignment.R`.
This did get pulled out by HMM, just didn't make it into the GN figure for some reason.
The other two are definitely not *hgcB*.
Checking them with MOTIF and BLAST, neither of them returned any hits.
Two other sequences are very long, from GCA_001804865.1 and GCA_001803235.1.
- MHBL01000254.1_3 is the hgcB hit from GCA_001803235.1. No hits to MOTIF in the rest protein (just the ferredoxin similarities)
- MHBK01000043.1_101 is the hgcB hit from GCA_001804865.1. No MOTIF hits to this one either other than the ferredoxin hits.

For both of these, we'll just keep them as is for now.
New figure from this.


#### Tree generation for PVC

Scripts for this are here: `code/binning/PVC_details/tree_generation_for_PVC.sh`


**Run tree generation workflow**

I then ran my tree generation workflow from HomeBio.
I used the rp16 bacteria HMMs requiring hits agains at least 8 of the HMMs.
I masked the alignment at 50% gaps, then ran RAxML on the alignment.
I downloaded these files here: `dataEdited/bins/binAnalysis/PVC_details/tree_generation`.


**Set up dataframe for renaming tips**

I manually made a tsv file with a name for the bin, the accession number, the phylum, and the name of the tip label: `tip_naming.xlsx`.
This is also in the `tree_generation` folder.
I was mostly interested in getting names for the reference sequences here.
The RefSeq genomes (GCF prefix) will be named with their genus and species (when available).
The GenBank genomes (GCA prefix) will be named with their order names, as determined by GTDB.
These two groups will get an accession ID after their name in the tree.
All bins will be prefaced with the family name (as determined by GTDB).

**Visualize tree**

Scripts to generate the tree are here: `code/phylogenies/PVC_tree.R`.

First I read out an unrooted tree, just to check it.
I did a few iterations of it this way.
Settled on using some Planctomycetes as outgroup to root this.

I then renamed all the branches using this file: `dataEdited/binning/phylogeny/PVC/tip_naming.tsv`.
I also generated a color vector to color code the bins by source and hgcA+ content.
I also labeled the Lentisphaerae and Kiritimatiellaeota branches.
PDF is saved out here: `results/bins/binAnalysis/PVC_details/PVC_tree.pdf`.


#### Run all PVC bins through some metabolic analyses

Code here: `code/binning/PVC_details/metabolic_genes_PVC.sh`.

This analysis looked for differences MHC content and the metabolic genes that we included in the batch HMM analysis between the PVC families, particularly the two Kiritimatiellaeota families.

I looked at this data here: `code/binning/PVC_details/metabolic_genes_PVC.R`.
It seems there's a cutoff at about 5 MHCs, where it's more likely that the bin will also have a cydAB, so I went with 5 MHCs as a point to split the bins into high and low MHC classifications.
Save them out here: `dataEdited/bins/binAnalysis/PVC_details/metabolism/high_MHC_bins.rds` and `dataEdited/bins/binAnalysis/PVC_details/metabolism/low_MHC_bins.rds`.
Let's use these to look at the depth profiles of all the bins, see if there are any trends there. See this below.


#### Investigate abundance of all PVC bins

Scripts here: `code/binning/PVC_details/abundance_PVC_bins.R`.

There is no major difference in the depth patterns of all of these bins.
In fact, they nearly all have identical trends in abundance across the water column.
The one exception is one of the Lentisphaerae bins, which is most abundant in the top part of the hypolimnion.
This bin has *hgcA*.
Other than that though, nothing much to write about.
