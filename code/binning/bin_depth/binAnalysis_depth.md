#### Bin analysis - depth

This protocol file outlines how we calcuated the depth of the bins, both the hgcA+ bins and the auto-generated high quality bins that will be used for Nick Scheel's project.


**Pull out coverage data for hgcA+ bins**

First I pulled out the coverage data for the hgcA+ bins from the mapping files.
I used a submission file for this on GLBRC.
I generated individual files for each metagenome, each containing all the mapping info for the bins that were mapped to.
Like my other analyses, I trimmed off 150 bp on both ends of each contig.


**Pull out coverage data for HQ bins**

Here, I pulled out the coverage data for all the automatically generated bins that Nick is going to use.
This was also done by submission, and uses the same executable as the above analysis.
