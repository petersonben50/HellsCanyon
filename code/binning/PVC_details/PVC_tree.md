### Scripts to generate PVC tree for HCC

I need to generate a tree that includes all the bins from the PVC superphylum, including Kiritimatiellaeota and Lentisphaerae.
I'll base this on the rp16 proteins.
For references, I'd like to include the hgcA+ bins from Lake Mendota, known isolate genomes for both these phyla, Dan Jones's references, and a few others.
Actually, the few others is covered, I picked out some to use with the 5M project.
Those are stored here: `~/references/genomes`
First let's collect those references.
We'll get the ORF sequences.


**Collect references**

The first thing I needed to do was grab the needed references.
I decided to grab the ORFs when I could, and settle for the nucleic acid sequences when I had to.


*Collect 5M PVC hgcA+ bins*

I first grabbed the Kiritimatiellaeota and Lentisphaerae bins from the 5M folder.
These already had the ORFs predicted, so I was able to just snag those.

*Bins from Jones et al, 2019*

I also wanted to include the PVC bins from the Jones et al 2019 paper.
The genomes from this study are listed here `~/references/jonesGenomes`.
We'll copy over the ORFs for the three PVC genes too.

*References selected for 5M project*

I had a bunch of them listed from the 5M paper, and the taxonomy ID's can be found here: `~/5M/dataEdited/binAnalysis/phylogeny/PVC/reference_taxonomy.tsv`.
The ORFs are here: `~/references/genomes/ORFs`.
Let's copy over the Kiritimatiellaeota and Lentisphaerae ORFs.
Also include two Planctomycetes genomes, for an outgroup.
Can't seem to find the *Victivallis vadensis* genome in here, which we know is GCA_003096415.1.
Let's download this from NCBI as well, here: https://www.ncbi.nlm.nih.gov/genome/1349.
Saved it by its accession ID: `GCF_003096415.1.faa`, then upload to the ORFs file: `~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/PVC/ORFs`.


*Cultured isolates*

I also wanted to include the two genomes from van Vliet et al, 2019, in which they isolated two Kiritimatiellaeota organisms.
These are listed at ENA under these accession numbers: SAMEA5207384 and SAMEA5207385.
I downloaded the protein fasta (the ORF annotations) from NCBI here:
- Pontiella desulfatans: https://www.ncbi.nlm.nih.gov/assembly/GCF_900890425.1
- Pontiella sulfatireligans: https://www.ncbi.nlm.nih.gov/assembly/GCF_900890705.1

Unzipped them to here: `/Users/benjaminpeterson/Documents/research/HellsCanyon/references/genomes/PVC/`.
Then I renamed them to be accessionNumber.faa
Then I uploaded these to `~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/PVC/ORFs`.


*Add bins from this project*
We'll use anvio_hgcA_0220, anvio_hgcA_0261, anvio_hgcA_0040, and anvio_hgcA_0110.

*Set up dataframe for renaming tips*

I manually made a tsv file with a name for the bin, the accession number, the phylum, and the name of the tip label: `dataEdited/binning/phylogeny/PVC/tip_naming.tsv`.



**Search for rp16 genes**


*Concatenate all bins*

First I concatenated all the bins into a single ORF faa file.
Then I generated a gene to bin file that linked each gene ID to its respective bin.
I searched the concatenated ORF file for each of the rp16 genes using HMMs with hmmsearch (v3.3.1).
I then pulled out the amino acid sequences for the rp16 genes and aligned them.

*Rename fasta headers with bin name*

I then renamed the fasta file headers with the name of the bin so that they could be concatenated later in Geneious.

*Generate alignment*

I then downloaded all the clean alignments to my local computer and imported them into Geneious.
The rpS3 gene has duplicate hits in LEN_0037, so I deleted LEN_0037 2.
I then concatenated all the alignments.
KIR_0040 had less than 1000 ungapped residues and was clearly missing a bunch of the proteins, so I initially removed it from the analysis, but in later iterations kept it, due to its similarity to a few of the genomes from HCC.
I exported this as `rp16_alignment_masked.afa` and uploaded it to the GLBRC servers: `~/HellsCanyon/dataEdited/binning/manualBinning/binsGood/phylogeny/PVC/tree_building`.

**Generate tree**

I then generated a maximum likelihood tree using RAxML (v8.2.11).
The R scripts to prepare this tree are here: `code/phylogenies/PVC_tree.R`.

First I read out an unrooted tree, just to check it.
I did a few iterations of it this way.
Settled on using some Planctomycetes as outgroup to root this.

I then renamed all the branches using this file: `dataEdited/binning/phylogeny/PVC/tip_naming.tsv`.
I also generated a color vector to color code the bins by source.

I also labeled the Lentisphaerae and Kiritimatiellaeota branches.
I always have a hard time making these actually pretty, so I saved it out as "unedited", and manually cleaned it up in Illustrator.
