
#### Assemble needed sequences using EMIRGE


**Set up**

I set up a tab-separated text file with mean insert size and standard deviation deviations.
I pulled the data from the library prep sheets here: `dataEdited/dnaSequencing/libraryPrep_seqOutput/`
We'll save that here: `dataEdited/16S_from_MG/metagenome_inserts.txt`.
I uploaded it here: `~/HellsCanyon/dataEdited/16S_from_MG`.


**Run EMIRGE on metagenomes**

I individually ran EMIRGE on each of the metagenomes that we have (from 2017 and 2019).
I unzipped the fastq files temporarily for this analysis, then deleted them once it was done.
The analysis was based on the SILVA database (SILVA_138_SSURef_NR99_tax_silva_trunc).
The insert size and standard deviation was called from the Excel sheet that I had previously set up (see above).


**Pull out all FASTA sequences**

Once EMIRGE was done, I pulled out all the fasta sequences that it had assembled and renamed them using the `emirge_rename_fasta.py` python script in EMIRGE.
This add the prior abundances to the names.
I then added the metagenome ID to each name to ensure they were unique and moved them to a single folder.


**Cluster 16S sequences across metagenomes**

I only wanted one set of 16S sequences, so I used CD-HIT to cluster the fasta files.
I used a cut-off of 0.97 to cluster them, then generated a key file, where the 16S sequence name and the representative sequences are stored in two separate columns.


**Extract abundance of 16S sequences**

I based this section on the workflow from [Vincent Denef's group](https://github.com/DenefLab/EMIRGE).
I extracted the total number of reads mapped to each metagenome's set of 16S sequences using samtools.
I also pulled out the fractional abundance of each sequence from the name of the 16S sequence (stored as "NormPrior"), and then cleaned up the count file.


**Rename seqs in abundance file with representative ID**

I then took the relative abundance file I generated in the previous step and renamed each sequence with the representative sequence generated using CD-HIT.


#### Classify 16S sequences using mothur

I then classified all of my 16S sequences using [TaxAss](https://github.com/McMahonLab/TaxAss).
This is also outlined on the Denef github workflow, cited above.
I used the FreshTrain dataset for my primary database (accessed March 10th, 2019) and Silva for the secondary one (also accessed March 10th, 2019).
I won't outline all the details here since I've mostly just followed the usual workflow for this.
