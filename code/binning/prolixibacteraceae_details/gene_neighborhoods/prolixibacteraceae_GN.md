### Notes on the gene neighborhood analysis for the hgcA+ Prolixibacteraceae

First I just looked at the gene neighborhood of the *hgcA* sequence from anvio_hgcA_0130, the bin from this study.
The ORF ID for the gene was HC18ME02_000000002532_3.
I just manually inspected this in Geneious, pulling out the amino acid sequences to either end and blasting them on the online NCBI server or through the MOTIF program (also online).

This showed some very interesting data, interesting enough to warrant further investigation into the other hgcA+ Prolixibacteraceae bins.

**Extract gene neighborhood of *hgcA* genes from reference genomes**

Next I wanted to pull out the gene neighborhoods of the seven *hgcA* sequences I identified in the reference genomes.
My solution in here is pretty wonky.
For some reason, the gene IDs of the reference genomes don't match up with the scaffolds as I'd expect.
I really don't like the way GFF files assign an ID that's different than the number that is assigned the ORF in the fasta file.
