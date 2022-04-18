### Protocol for identifying phylogeny of hgcA+ Geobacterales bin



**Look at GTDB tree**

First we'll investigate the GTDB tree to identify nearby references.
Going 6 levels back from the bin, we found 37 bins, divided into two distinct clusters.
For the tree, let's pull out the branch that includes anvio_hgcA_0210, then two from the other cluster for an outgroup.
For the outgroup, we used two that were also used in the Mendota paper: GCF_000016745.1 and GCF_000007985.2.
Save this out to a list: `GTDB_ref_accession_numbers.txt`, which then gets saved to an Excel file: `GTDB_ref_accession_numbers.xlsx`.

I removed anvio_hgcA_0210 from this list, then used NCBI's E-Utilities to pull out the assemblies for the reference bins.
For this, I first retrieved the entries using epost.
