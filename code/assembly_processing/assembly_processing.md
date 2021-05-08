### Assembly Processing 2017, 2018, and 2019 metagenomes

This document outlines the protocol I followed for processing the metagenomes from the 2017, 2018, and 2019 intensives.
It will also include my notes about the process, and directions to other files as needed.

**Data retrieval and inspection**

First thing we'll need to do is retrieve the data.
QB3 emailed us directions, I stored them here: `protocols/sequencingQB3/2019_intensive/round_1/191204_NovaSeq- Peterson data is available on the GSL secure FTP server.eml`.
They followed up with an email that the multiplexing got screwed up for some of the samples.
Let's see what we have so far.

First I retrieved the data as outlined in the email, using lftp.
The multiplexing really did get messed up. A, B, C, and D were all ~30million or less, with KMBP003C at 4.8 million.
KMBP004F is 478 million though!

I made a file (`misrun_samples`) to store the small metagenomes (KMBP004A-KMBP004D)

On Dec. 26th, the misrun samples were re-run (KMBP004A-KMBP004D).
I downloaded these and stored them on the GLBRC.
All the good metagenomes are now in a single folder.

*KMBP007, KMBP008, KMBP009*

These sets of metagenomes were submitted to MSU sequencing center in the summer of 2020.
They are mostly 2018 samples, with KMBP009 being additional 2019 samples.
I received an email from MSU on August 4th, 2020 that the sequencing was completed.
I used wget to grab these, as outlined here: https://rtsf.natsci.msu.edu/genomics/data-retrieval/.
I then added these metagenome IDs to `dataEdited/metagenomes/MG_list_all.txt`.


**Trimming metagenomes**

Next I quality trimmed the metagenome.
To do this, I used the fastp program (version 0.20.0), which trims the metagenome to our liking and gives us the stats about the metagenome both before and after trimming.
I download the program using conda, into our bioinformatics virtual environment.

I used this program to save out both the paired and unpaired reads after trimming.
I also stored the failed reads, just so we can check to see what's driving the loss of reads if need be.

I cut the tail end of the reads based on the read quality.
The front end quality of the reads are usually pretty good so I didn't bother with trimming the front or using the `cut_right option`.
I used a 10bp window with a quality score cutoff of 20.
Since I was trimming these down, I wanted to cut out any short reads that get overtrimmed.
Length filtering is set to 15 by default, but I changed it to 100.
I left the default quality filtering in place as well, which is a quality limit of 15, and a maximum of 40% of bases in a read below that threshold.
I also merged the reads.
For some reason, once I use the merge function, the single reads are filled in at the end, so don't get worried if you're checking the size of the files as the program runs and the single reads file is empty.
I then went through and took a look at the output of fastp, the files of which are stored here: `/Users/benjaminpeterson/Documents/research/HellsCanyon/dataRaw/metagenomes/2019_intensive/fastp`.

- KMBP004A: The duplication rate here was relatively high (17.2%), so we may want to deal with that. Most of the reads that have been duplicated have only been duplicated a single time (2 overall) which makes me think it's more likely due to coverage. Might be good to ultimately check to see where these duplicate reads are coming from, but I don't think I can justify cutting them out since it could just be due to depth. 98% of reads passed, lost 1.4% to reads with overall low quality. Most of other ones lost due to overtrimming. ~11% of them are overlapping. Seems to be a good deal of C and G replication.

- KMBP004B: Duplication rate around 11%, mostly of reads found twice. Increase in GC ratio from 20% to ~30% in these duplicated reads. 5% overlapping reads. 97% of reads passed filters, again mostly due to overall low quality of read quality. Seems to be a good deal of C and G replication.

- KMBP004C: Whoa, high duplication here, 28.8%. Similar story, increasing GC content in duplicated reads. 98% passing.

- KMBP004D: 13% duplication. 97.2% of reads are kept.

- KMBP004E: 98.4% of reads passed filters, ones that failed mostly due to overall low quality. Duplication rate ~8.4%. All looks good!

- KMBP004F: 97.6% of reads pass filters. Failures again mostly due to low read quality overall. 6.9% duplication. This is the sample with a shit-ton of reads, there are 933 million total reads after filtering.

Unfortunately fastp doesn't read out how many single reads vs. paired end reads there are, so I ended up running fastQC on the data as well just to get those counts and get some other information on the read quality.
I'll also save the read counts here: `~/HellsCanyon/dataEdited/metagenomes/2019_intensive_ancillary_info/metagenome_read_count.tsv`, then download it to `dataEdited/metagenomes/reports/metagenome_read_count.tsv`.


**Count reads in metagenome pre-trimming**

I counted the reads in each pre-trimmed metagenome.
I used the linux function zgrep to count the lines that contained the "@" sign in each file ("@" should only be found at the start of each header line for a sequence in a fastq file), then read out that information to a report file.


**Count reads in metagenome**

Next I counted the number of reads in each of the trimmed files.
I did this using the same method as with the pre-trimmed files, just now on the forward, reverse, single, and merged reads.


**Coverage post-trimming**

I also wanted to determine the overall coverage of each of the metagenomes for normalization purposes later.
For this, I used the kseq_fastq_base function in [readfq program](https://github.com/billzt/readfq), which reads the number of bases in a nucleotide file.
I did this for each set of fastq files within a metagenome (forward, reverse, single, and merged), since we use all of those for mapping.
I downloaded all these files to my local computer (`dataEdited/metagenomes/reports/metagenome_coverage.tsv`).

In the R script `code/scripts/coverage_normalization.R`, I generated an R object containing a vector that has normalization values for each metagenome, saved out to `dataEdited/metagenomes/reports/metagenome_normalization_vector.rds`.
For normalization, we assumed 100 million 150bp PE sequencing reads.


**Cluster metagenomes using Mash**

Next I wanted to cluster the metagenomes within a year using the program Mash as a way to determine clusters for assemblies.
We'll also compare all of them (across years) just to get a sense of the similarity.

*Generate sketches for each metagenome*

First I generated an individual sketch for each metagenome.
I used 20 threads, set a seed value of 50, and flagged with -r to indicate we are using reads here.
I'll require a minimum of 2 copies of each kmer (flag: -m 2).
Finally, I'll use a kmer size of 21 and the sketch size will be 100,000.
When running it with p and r flags, it warns that "The option p will be ignored with r."

*Paste sketches*

Just "paste" them all together.

*Run distance calculation*

Then we calculated the Jaccard distance between each of the metagenomes using `mash dist`.
This returns stdout that we saved into a tsv file.
The first two columns are the names of the metagenome files.
The third column is the Jaccard distance.
The fourth column is the p-value and the last column is the number of hashes shared between the reference and the query files.


**Metagenome assembly**

*Grouping metagenomes*

Next, I needed to assign the metagenomes to clusters that I could assemble together, since I did not want to either run all individual assemblies or do it as a large co-assembly.

*Assemble the metagenomes*

Once the metagenomes were trimmed and ready to go, I assembled them using metaSPADes.
I assembled the metagenomes in the clusters that I defined above, except for in 2019 where each metagenome was assembled individually.
I also performed coassemblies of the metagenomes from 2017 and 2018.
I included the paired reads, the merged reads, and the single reads.
For kmers, I included 21,33,55,77,99, and 127.

*Clean the assembly*

I then cleaned the assembly to only include scaffolds over 1000bp long, and simplified the names of the scaffolds.
I saved out a report file so we can tie the scaffolds back to their original names, if we need to use the assembly graph.

*Get assembly stats*

Then I wanted to get the assembly stats for each assembly.
For this, I used the abyss-fac perl script from the ABySS assembler (`code/assembly_processing/abyss-fac.pl`).
This returns the usual metagenome stats, such as total length, number of scaffolds, n50, things like that.
Then I aggregated the assembly stats into one file (`all_assemblies_stats.txt`) and download that to my local computer (`dataEdited/metagenomes/reports/`).

**Predict proteins**

Next I wanted to predict proteins so that we could search for the hgcA sequences, as well as other metabolic proteins.
I used Prodigal for this (version 2.6.3), on the metagenome mode (`-p meta`).
I saved the amino acid and nucleic acid sequences, as well as the GFF files.
Then, I cleaned both the sequence files with cleanFASTA.
Finally, I counted up all the ORFs we had in the different assemblies, and stored it here: `dataEdited/metagenomes/reports/`.



**Next steps**

From this script, proceed to the `mapping_assembly` workflow, under `code/analysis_workflows`.
