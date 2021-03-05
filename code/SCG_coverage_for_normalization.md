### Protocol to identify and quantify SCGs in Hells Canyon metagenomes

**Get set up**

First I needed to prep the assembly and metagenome list files from the keys that I made.

**Submit jobs**

I wrote these analyses as submission files.
Each job goes after a different gene.
First it pulls the copies of the gene out of all the different assemblies.
It then dereplicates those genes (97% identity).
It then calculates the coverage of all the identified genes in all the different metagenomes.



**Clean up**

Remove the assembly/metagenome lists that I generated in this step.
