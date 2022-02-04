### Protocol for assembly mapping to 2017 metagenomes

This is the protocol I'm using for mapping reads to my genomes from the 2017 dataset.


**Pick assemblies and metagenomes for analysis**

First, I set up two tsv files that contained the metagenomes and assemblies that I was going to use.
I added a second column to each one to include the sampling year that they were from, since I'm planning to separate out those MGs and assemblies for abundance analysis.
I saved them out here: `~/HellsCanyon/metadata/lists`.
They're saved as as `metagenome_key.txt` and `assembly_key.txt`

*2017 assemblies*: I decided to use the six coassemblies of the similar redox state sequences, as well as the coassembly of all the 2017 metagenomes.

*2018 assemblies*: I used the seven coassemblies of metagenomes from similar redox states and the coassembly of all the 2018 metagenomes.

*2019 assemblies*: I just used the individual assemblies here.

*Metagenomes*: I used all the metagenomes from the water column for these analyses, separated out by year. That's 34 total metagenomes.

In R, I then joined the `assembly_key.csv` and `metagenome_key.csv` into a `mapping_key.csv` file that I could use to run the mapping scripts.
This R script is here: `code/metagenome_metadata.R`.
It saved out a tsv to the `metadata/lists` folder called `mapping_key.tsv`.
I then uploaded this file to GLBRC.

**Calculate total coverage for each metagenome**

I wanted to calculate the total number of nucleotides from each metagenome, which I could use to normalize the coverage of all our sequences.
For this, I used the [readfq program](https://github.com/billzt/readfq), counting the number of nucleotides in the forward, reverse, single, and merged read files, all of which will be used for mapping.
I downloaded all these files to my local computer (`dataEdited/metagenomes/2017_analysis_assembly/mapping/metagenome_coverage.tsv`).
In the R script `code/2017_analysis_assembly/mapping_normalization.R`, I generated an R object containing a vector that has normalization values for each metagenome.
For normalization, we assumed 100 million 150bp PE sequencing reads.


**Map reads and process output**

First I had to build an index of each assembly.
Then, I checked to see which of the mapping steps I had already done (since this analysis was done in stages to a certain extent).
I saved out this list of the mapping pairs that needed to be completed yet.

Then I mapped the paired-end reads, the single reads, and the merged reads to the indices using bowtie2.
Once the mapping was done, I used samtools to convert the files to BAM files, then sorted and indexed the files.
