#!/bin/sh

######################
# code/binning/binAnalysis_metabolism.sh
# Benjamin D. Peterson

# This set of scripts will contain the bash
# code we need to look for metabolic genes
# in our hgcA+ bins.
######################

cd ~/HellsCanyon/dataEdited
mkdir binAnalysis
cp binning/bins_hgcA_keepers/taxonomy_summary.txt binAnalysis
cp binning/bins_hgcA_keepers/checkM_stats.csv binAnalysis
cp binning/bins_hgcA_keepers/ORFs.faa binAnalysis
cp binning/bins_hgcA_keepers/ORFs_G2B.tsv binAnalysis


####################################################
####################################################
# Run METABOLIC
####################################################
####################################################
cd ~/HellsCanyon/dataEdited
mkdir binAnalysis/metabolism
mkdir binAnalysis/metabolism/METABOLIC_bins

cat ~/HellsCanyon/dataEdited/binning/bins_hgcA_keepers/bins_hgcA_keepers_list.txt | while read binID
do
  echo "Moving" $binID".fna to" $binID".fasta"
  cp binning/bins_hgcA_keepers/DNA/$binID.fna binAnalysis/metabolism/METABOLIC_bins/$binID.fasta
done
ls -d /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/binAnalysis/metabolism/METABOLIC_bins/*.fasta
# Edit the METABOlIC_processing.sub file to include all these bins
condor_submit ~/HellsCanyon/code/submission/METABOLIC_processing.sub


####################################################
####################################################
# Custom set of metabolic HMMs
####################################################
####################################################

screen -S HCC_metabolic_HMMs
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate batch_HMMs
PYTHONPATH=''
PERL5LIB=''
binAnalysis=~/HellsCanyon/dataEdited/binAnalysis
scripts=~/HellsCanyon/code/generalUse
metabolic_HMMs=~/HellsCanyon/references/metabolic_HMMs
cd $binAnalysis

#chmod +x $scripts/batch_HMMs.py
python $scripts/batch_HMMs.py --orf_file $binAnalysis/ORFs.faa \
                              --g2b $binAnalysis/ORFs_G2B.tsv \
                              --hmm_folder $metabolic_HMMs\
                              --hmm_csv $metabolic_HMMs.csv \
                              --output $binAnalysis/metabolism/batch_HMMs
conda deactivate




#########################
# Confirm dsrA phylogeny
#########################
screen -S HCC_dsrA_tree
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''
scripts=~/HellsCanyon/code/generalUse
binAnalysis=~/HellsCanyon/dataEdited/binAnalysis

# Set up directory
cd $binAnalysis/metabolism/batch_HMMs
mkdir dsrA

# Copy over my sequences and align them
sed 's/*//' alignments/hmm_hits/dsrA.faa > dsrA/dsrA.faa
cd dsrA
grep '>' dsrA.faa | \
  sed 's/>//' \
  > dsrA_list.txt
muscle -in dsrA.faa \
        -out dsrA.afa

# Copy in reference sequences
cp ~/references/metabolicProteins/sulfur/dsrA/dsrA_karthik_clean.afa \
    dsrA_karthik_clean.afa

# Align seqs
muscle -profile \
        -in1 dsrA_karthik_clean.afa \
        -in2 dsrA.afa \
        -out dsrA_phylogeny.afa
# Trim alignment
trimal -in dsrA_phylogeny.afa \
        -out dsrA_phylogeny_trimmed.afa \
        -gt 0.5

# Generate ML tree
#screen -S EG_dsrA_tree
python $scripts/cleanFASTA.py dsrA_phylogeny_trimmed.afa
mv -f dsrA_phylogeny_trimmed.afa_temp.fasta dsrA_phylogeny_trimmed.afa
FastTree dsrA_phylogeny_trimmed.afa \
    > dsrA_phylogeny_trimmed.tree


####################################################
####################################################
# Search for MHCs
####################################################
####################################################

binAnalysis=~/HellsCanyon/dataEdited/binAnalysis
scripts=~/HellsCanyon/code/generalUse
mkdir $binAnalysis/metabolism/MHCs
cd $binAnalysis/metabolism/MHCs
$scripts/Find_multiheme_protein.py $binAnalysis/ORFs.faa 3
mv $binAnalysis/ORFs_3_heme* .

echo -e "binID\tgeneID\themeCount" > heme_count_bins.tsv
tail -n +2 ORFs_3_heme_count.txt | awk -F '\t' '{ print $1 }' | while read geneID
do
  binID=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' $binAnalysis/ORFs_G2B.tsv`
  hemeCount=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' ORFs_3_heme_count.txt`
  echo -e $binID"\t"$geneID"\t"$hemeCount
  echo -e $binID"\t"$geneID"\t"$hemeCount >> heme_count_bins.tsv
done

##########################
# Search for MHCs with 10 heme-binding sites
##########################
mkdir $binsGood/metabolism/MHCs_10
cd $binsGood/metabolism/MHCs_10
$scripts/Find_multiheme_protein.py $binsGood/ORFs.faa 10
mv $binsGood/ORFs_10_heme* .

echo -e "binID\tgeneID\themeCount" > heme_count_bins_10.tsv
tail -n +2 ORFs_10_heme_count.txt | awk -F '\t' '{ print $1 }' | while read geneID
do
  binID=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' $binsGood/binsGood_G2B.tsv`
  hemeCount=`awk -F '\t' -v geneID="$geneID" '{ if ($1 == geneID) print $2 }' ORFs_10_heme_count.txt`
  echo -e $binID"\t"$geneID"\t"$hemeCount
  echo -e $binID"\t"$geneID"\t"$hemeCount >> heme_count_bins.tsv
done


####################################################
####################################################
# Search for BB-OMPs
####################################################
####################################################

##########################
# Pull out genes adjacent to MHCs
##########################

# Pull out names of adjacent genes
binAnalysis=~/HellsCanyon/dataEdited/binAnalysis
cd $binAnalysis/metabolism
mkdir PCC
mkdir PCC/list
rm -f $binAnalysis/metabolism/PCC/adjacent_genes_all_list.txt

cat MHCs/ORFs_3_heme_list.txt | while read gene
do
  echo "Working on" $gene
  scaffold=$(echo $gene | rev | cut -d"_" -f2- | rev)
  assembly=$(echo $gene | rev | cut -d"_" -f3- | rev)
  ORFnumber=$(echo $gene | rev | cut -d"_" -f1 | rev)
  preceedingORFnumber=$(expr $ORFnumber - 1)
  followingORFnumber=$(expr $ORFnumber + 1)
  echo $scaffold"_"$preceedingORFnumber >> PCC/list/adjacent_genes_all_list.txt
  echo $scaffold"_"$followingORFnumber >> PCC/list/adjacent_genes_all_list.txt
done

# Find unique gene names
cd PCC
wc -l list/adjacent_genes_all_list.txt
sort list/adjacent_genes_all_list.txt | \
  uniq \
  > list/adjacent_genes_unique_list.txt

# Pull out adjacent genes
rm -f adjacent_genes.faa
cat list/adjacent_genes_unique_list.txt | while read geneID
do
  assembly=$(echo $geneID | rev | cut -d"_" -f3- | rev)
  echo "Looking for" $geneID "in" $assembly
  grep -A 1 -m 1 $geneID$ $binAnalysis/ORFs.faa >> adjacent_genes.faa
done


#########################
# Search adjacent genes for BBOMPs
#########################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""
cd $binAnalysis/metabolism/PCC

# Run the PCC HMM
pcc_omp_HMM=~/HellsCanyon/references/metabolicProteins/EET/pcc_omp.HMM
hmmsearch --tblout pcc_omp_custom.out \
          -T 40 \
          $pcc_omp_HMM \
          adjacent_genes.faa \
          > pcc_omp_custom_output.txt
# Pull out the sequences of interest
scripts=~/HellsCanyon/code/generalUse/
python $scripts/extract_protein_hitting_HMM.py \
        pcc_omp_custom.out \
        adjacent_genes.faa \
        pcc_omp_custom.faa


#########################
# Generate a phylogeny
#########################

cd $binAnalysis/metabolism/PCC
mkdir references

# Find references from RefSeq
blastp -query pcc_omp_custom.faa \
        -db ~/references/ncbi_db/refseq/refseq_protein \
        -evalue 0.005 \
        -outfmt '6 sseqid sseq' \
        -max_target_seqs 5 \
        -num_threads 3 \
        -out references/refseq_pcc_omp.txt
# Pull out amino acid sequences and dereplicate them

cd references
awk -F '\t' '{ print ">"$1"\n"$2 }' refseq_pcc_omp.txt > refseq_pcc_omp.faa
sed -i 's/ref|//' refseq_pcc_omp.faa
sed -i 's/|//' refseq_pcc_omp.faa
sed -i 's/-//g' refseq_pcc_omp.faa
# Dereplicate refseq sequences against reference dataset
cd-hit-2d -i ~/HellsCanyon/references/metabolicProteins/EET/pcc_omp.faa \
          -i2 refseq_pcc_omp.faa \
          -o blast_pcc_omp_uniq.faa \
          -c 0.97 \
          -n 5
# Get final list of refseq references
grep '>' blast_pcc_omp_uniq.faa | \
  sed 's/>//' \
  > blast_pcc_omp_uniq_list.txt
# Concatenate sequences
cd $binAnalysis/metabolism/PCC
mkdir phylogeny
cat references/blast_pcc_omp_uniq.faa \
    ~/HellsCanyon/references/metabolicProteins/EET/pcc_omp.faa \
    pcc_omp_custom.faa \
    > phylogeny/bbomp_phylogeny.faa

# Generate alignment
cd phylogeny
muscle -in bbomp_phylogeny.faa \
        -out bbomp_phylogeny.afa
# Trim the alignment
trimal -in bbomp_phylogeny.afa \
        -out bbomp_phylogeny_trimmed.afa \
        -gt 0.5
FastTree bbomp_phylogeny_trimmed.afa > bbomp.tree


# Retrieve reference information on local computer
HCC
cd dataEdited/bins/binAnalysis/metabolism/PCC
epost -db protein -input blast_pcc_omp_uniq_list.txt | \
    esummary | \
    xtract -pattern DocumentSummary -element AccessionVersion,Organism > refseq_bbomp_metadata.tsv



####################################################
####################################################
# Check PCC_porin hits from batch HMMs
####################################################
####################################################
binAnalysis=~/HellsCanyon/dataEdited/binAnalysis
cd $binAnalysis/metabolism
mkdir PCC_batchHMMs
mkdir PCC_batchHMMs/PCC

cp batch_HMMs/alignments/hmm_hits/pcc_porin.faa PCC_batchHMMs/PCC/pcc_porin.faa
cd PCC_batchHMMs/PCC/

grep '>' pcc_porin.faa | sed 's/>//' | awk '{ print $1 }' | while read gene
do
  #echo "Working on" $gene
  scaffold=$(echo $gene | rev | cut -d"_" -f2- | rev)
  assembly=$(echo $gene | rev | cut -d"_" -f3- | rev)
  ORFnumber=$(echo $gene | rev | cut -d"_" -f1 | rev)
  preceedingORFnumber=$(expr $ORFnumber - 1)
  followingORFnumber=$(expr $ORFnumber + 1)
  grep $scaffold"_"$preceedingORFnumber $binAnalysis/metabolism/MHCs/ORFs_3_heme_count.txt
  grep $scaffold"_"$followingORFnumber $binAnalysis/metabolism/MHCs/ORFs_3_heme_count.txt
done



####################################################
####################################################
# Identify CAZymes in each bin
####################################################
####################################################
# On local computer
cd /Users/benjaminpeterson/Documents/research/HellsCanyon/dataEdited/bins/binAnalysis/bin_ORFs
split -l 60000 ORFs.faa ORFs_split_



####################################################
####################################################
# Characterize MoORs in bins
####################################################
####################################################
#mkdir ~/HellsCanyon/references/custom_hmms
#cp ~/5M/references/MoORs/MoOR.HMM ~/HellsCanyon/references/custom_hmms
#cp ~/5M/references/MoORs/MoOR_final_reference.afa ~/HellsCanyon/references/custom_hmms/MoOR_reference.afa

# Search bins for MoORs
screen -S HCC_MoORs
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
scripts=~/HellsCanyon/code/generalUse/
binAnalysis=~/HellsCanyon/dataEdited/binAnalysis

cd $binAnalysis
mkdir metabolism/MoORs
hmmsearch --tblout metabolism/MoORs/MoOR.out \
          -T 50 \
          ~/HellsCanyon/references/custom_hmms/MoOR.HMM \
          ORFs.faa \
          > metabolism/MoORs/MoOR_log.txt

# Pull out the gene names, use this to find correct portion of G2B file.
cd metabolism/MoORs/
rm -f MoOR_G2B.tsv
rm -f putative_MoORs.faa
grep -v '#' MoOR.out | \
  awk '{ print $1 }' | while read gene
  do
    awk -F '\t' -v gene="$gene" '$1 == gene { print $0 }' ../../ORFs_G2B.tsv >> MoOR_G2B.tsv
    grep -A 1 $gene\$ ../../ORFs.faa >> putative_MoORs.faa
    tail -n 1 MoOR_G2B.tsv
  done
sed -i 's/*//g' putative_MoORs.faa
grep ">" putative_MoORs.faa | \
  sed 's/>//' > putative_MoORs_list.txt

# Generate alignment
muscle -in putative_MoORs.faa \
        -out putative_MoORs.afa

# Consensus alignment between study proteins and references
muscle -profile \
        -in1 putative_MoORs.afa \
        -in2 ~/HellsCanyon/references/custom_hmms/MoOR_reference.afa \
        -out putative_MoORs_ref_1.afa

# Clean it up
python $scripts/cleanFASTA.py putative_MoORs_ref_1.afa
mv -f putative_MoORs_ref_1.afa_temp.fasta putative_MoORs_ref_1.afa

# Mask the alignment at 50% gaps
trimal -in putative_MoORs_ref_1.afa \
        -out putative_MoORs_ref_1_masked.afa \
        -gt 0.5
FastTree putative_MoORs_ref_1_masked.afa > putative_MoORs_1.tree


#########################
# MoOR tree take 2
#########################

cd ~/HellsCanyon/dataEdited/binAnalysis/metabolism/MoORs
mkdir tree_take_2
rm -f tree_take_2/MoORs.faa
cat MoORs_list.txt | while read MoOR
do
  grep -A 1 $MoOR$ putative_MoORs.faa >> tree_take_2/MoORs.faa
done
cd tree_take_2
# Generate alignment of my sequences
muscle -in MoORs.faa \
        -out MoORs.afa
# Consensus alignment between study proteins and references
muscle -profile \
        -in1 MoORs.afa \
        -in2 ~/HellsCanyon/references/custom_hmms/MoOR_reference.afa \
        -out MoORs_ref.afa

# Clean it up
python $scripts/cleanFASTA.py MoORs_ref.afa
mv -f MoORs_ref.afa_temp.fasta MoORs_ref.afa

# Mask the alignment at 50% gaps
trimal -in MoORs_ref.afa \
        -out MoORs_ref_masked.afa \
        -gt 0.5

# Generate tree using RAxML
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -m PROTGAMMAAUTO \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s MoORs_ref_masked.afa \
        -n MoORs
