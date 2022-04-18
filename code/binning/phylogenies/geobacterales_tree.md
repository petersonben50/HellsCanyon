### Protocol for identifying phylogeny of hgcA+ Geobacterales bin



**Look at GTDB tree**

First we'll investigate the GTDB tree to identify nearby references.
Going 6 levels back from the bin, we found 37 bins, divided into two distinct clusters.
For the tree, let's pull out the branch that includes anvio_hgcA_0210, then two from the other cluster for an outgroup.
For the outgroup, we used two that were also used in the Mendota paper: GCF_000016745.1 and GCF_000007985.2.
Save this out to a list: `GTDB_ref_accession_numbers.txt`, which then gets saved to an Excel file: `GTDB_ref_accession_numbers.xlsx`.

I removed anvio_hgcA_0210 from this list, then used NCBI's E-Utilities to pull out the assemblies for the reference bins.
I had to do this separately for the RefSeq vs. GenBank sequences.
For this, I first retrieved the entries using epost, then summarized and extracted the FTP paths.
Then I used `wget` to pull down the needed genomes.
Then, I ran them through the usual ORF prediction and rp16 searches to find the rp16 genes.
These were aligned, then concatenated and masked in Geneious.
I generated a ML tree with FastTree, checked it out in Geneious, then went ahead with a tree from RAxML.
They were nearly identical in this case.

I also identified *hgcA* genes from these reference genomes and predicted their taxonomy with GTDB.

Finally, I visualized this data together in a tree.
Looks like our bin is definitely Pelobacteraceae, possibly of an unknown genus and definitely closely related to the GEO_0030 from the Mendota study.
Unlike the Bacteroidetes bin, there is no shortage of hgcA+ genomes in this tree.
It seems to be widely distributed throughout this family.
