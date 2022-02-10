#!/bin/sh

#########################
# code/metabolic_analyses/metabolic_hmms.sh
# Benjamin D. Peterson
#########################


##################################################
##################################################
# Identification of metabolic proteins
##################################################
##################################################

#########################
# Run HMMs
#########################

screen -S HCC_metabolic_HMMs

cd ~/HellsCanyon/dataEdited/
mkdir metabolic_analyses
mkdir metabolic_analyses/identification

scripts=~/HellsCanyon/code/generalUse
metabolic_HMMs=~/HellsCanyon/references/metabolic_HMMs
ORFs=~/HellsCanyon/dataEdited/assemblies/ORFs

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

tail -n +2 $metabolic_HMMs.csv | awk -F ',' '{ print $2 }' | while read HMM
do
  cd ~/HellsCanyon/dataEdited/metabolic_analyses
  geneName=$(awk -F ',' -v HMM="$HMM" '$2 == HMM { print $1 }' $metabolic_HMMs.csv)

  if [ ! -e identification/$geneName.afa ]; then
    echo "Searching for" $geneName
    mkdir identification/$geneName

    cat ~/HellsCanyon/metadata/lists/assembly_key.csv | while read -r line
    do
      assembly=`echo $line | awk -F ',' '{ print $1 }'`
      year=`echo $line | awk -F ',' '{ print $2 }'`
      if [ ! -d identification/$geneName/$year ]; then
        mkdir identification/$geneName/$year
      fi

      if [ ! -e identification/$geneName/$year/$assembly.out ]; then
        echo "Searching for" $geneName "in" $assembly
        hmmsearch --tblout identification/$geneName/$year/$assembly\_$geneName.out \
                  --cpu 4 \
                  --cut_tc \
                  $metabolic_HMMs/$HMM \
                  $ORFs/$assembly.faa \
                  > identification/$geneName/$year/$assembly\_$geneName.txt
        lineCount=`wc -l < identification/$geneName/$year/$assembly\_$geneName.out`
        if [ $lineCount -eq 13 ]; then
          echo "No" $gene "hits in" $geneName
        else
          echo "Pulling" $geneName "sequences out of" $assembly
          python $scripts/extract_protein_hitting_HMM.py \
                  identification/$geneName/$year/$assembly\_$geneName.out \
                  $ORFs/$assembly.faa \
                  identification/$geneName/$year/$assembly\_$geneName.faa
        fi
      else
        echo "Search for" $geneName "is already done in" $assembly
      fi
    done

    # Aggregate all sequences and align to HMM
    cd identification
    awk -F ',' '{ print $2 }' ~/HellsCanyon/metadata/lists/assembly_key.csv | \
      sort | uniq | \
      while read year
      do
        shopt -s nullglob
        for i in $geneName/$year/*$geneName.faa; do FOUND=$i;break;done

        if [ ! -z $FOUND ]; then
          echo "Concatenating and aligning" $geneName
          cat $geneName/$year/*$geneName.faa \
              > $geneName\_from_$year\_all.faa
          hmmalign -o $geneName\_from_$year.sto \
                      $metabolic_HMMs/$HMM \
                      $geneName\_from_$year\_all.faa
          $scripts/convert_stockhold_to_fasta.py $geneName\_from_$year.sto
          grep '>' $geneName\_from_$year\_all.faa | \
            sed 's/>//' \
            > $geneName\_from_$year\_all_list.txt
        else
          echo "No" $geneName "sequences found at all :("
        fi
    done
    FOUND=""
  else
    echo "Already pulled out" $geneName "sequences"
  fi
done

conda deactivate

exit




#########################
# Extract depths of all scaffolds
#########################

# Generate list of all scaffolds that have potential metabolic genes on them
cd ~/HellsCanyon/dataEdited/metabolic_analyses
mkdir depth
grep '>' -h identification/201*/*.afa | \
  sed 's/>//' | \
  cut -d"_" -f1,2 | \
  sort | \
  uniq \
  > depth/scaffold_all_list.txt

# Pull out all depths
chmod +x /home/GLBRCORG/bpeterson26/HellsCanyon/code/executables/aggregate_depth_proteins.sh
condor_submit /home/GLBRCORG/bpeterson26/HellsCanyon/code/submission/aggregate_depth_proteins.sub


#########################
# Dereplicate sequences
#########################

cd ~/HellsCanyon/dataEdited/metabolic_analyses
mkdir dereplication
cdhit=~/programs/cdhit-master

echo "geneID,geneName" > metabolic_gene_key.csv

tail -n +2 $metabolic_HMMs.csv | awk -F ',' '{ print $1 }' | while read geneName
do
  cd ~/HellsCanyon/dataEdited/metabolic_analyses
  awk -F ',' '{ print $2 }' ~/HellsCanyon/metadata/lists/assembly_key.csv | \
    sort | uniq | \
    while read year
    do
      $cdhit/cd-hit -g 1 \
                    -i identification/$geneName\_from_$year\_all.faa \
                    -o dereplication/$geneName\_from_$year.faa \
                    -c 0.97 \
                    -n 5 \
                    -d 0
      $cdhit/clstr2txt.pl dereplication/$geneName\_from_$year.faa.clstr \
        > dereplication/$geneName\_from_$year.tsv
    done
  cd dereplication
  cat $geneName\_from_*.faa > $geneName\_derep.faa
  grep '>' $geneName\_derep.faa | \
    sed 's/>//' \
    > $geneName\_derep_list.txt
  awk -v geneName="$geneName" '{ print $0","geneName }' $geneName\_derep_list.txt >> ../metabolic_gene_key.csv
done


############################################
############################################
# Verify and classify narG gene
############################################
############################################

#########################
# Prep reference dataset
#########################
# Locally
cd ~/Documents/research/HellsCanyon/references/narG
grep '>' narG_luke_database.faa | \
  awk '{ print $1"\t"$2 }' | \
  sed 's/>//' | \
  sed 's/\[organism=//' | \
  sed 's/\]//' > narG_luke_database.tsv
grep '>' narG_vs_nxrA_1.faa narG_vs_nxrA_2.faa | \
  awk -F '_ ' '{ print $1"\t"$2 }' | \
  sed 's/nxrA_[1:2].faa:>//' > narG_vs_nxrA.tsv


# On GLBRC
screen -S HCC_narG_ref_tree
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""

cd ~/HellsCanyon/references/narG
python ~/HellsCanyon/code/generalUse/cleanFASTA.py narG_luke_database.faa
mv -f narG_luke_database.faa_temp.fasta narG_luke_database_clean.faa

# Generate alignment
muscle -in narG_luke_database_clean.faa \
        -out narG_luke_database_clean.afa
FastTree narG_luke_database_clean.afa > narG_luke_database_fastTree.tree


#########################
# Compare hits against narG HMM to reference dataset
#########################
cd ~/HellsCanyon/dataEdited/metabolic_analyses
mkdir N
mkdir N/narG
cd N/narG

# Generate FastTree
cp ~/HellsCanyon/dataEdited/metabolic_analyses/dereplication/narG_derep.faa .
cat ~/HellsCanyon/references/narG/narG_luke_database_clean.faa \
    narG_derep.faa > narG_for_tree.faa
muscle -in narG_for_tree.faa \
        -out narG_for_tree.afa
FastTree narG_for_tree.afa > narG_for_tree.tree

# Generate RAxML ML tree
trimal -in narG_for_tree.afa \
        -out narG_for_tree_trimmed.afa \
        -gt 0.5

# RAxML tree
# Used the reduced fasta file
raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -N autoMRE \
        -m PROTGAMMALG \
        -x 2381 \
        -T 20 \
        -s narG_for_tree_trimmed.afa.reduced \
        -n narG


############################################
############################################
# Classify dsrA genes
############################################
############################################

screen -S HCC_dsrA_tree
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=''
scripts=~/HellsCanyon/code/generalUse
# Set up directory
cd ~/HellsCanyon/dataEdited/metabolic_analyses
mkdir sulfur
mkdir sulfur/dsrA

# Copy over my sequences and align them
sed 's/*//' dereplication/dsrA_derep.faa > sulfur/dsrA/dsrA.faa
grep '>' sulfur/dsrA/dsrA.faa | \
  sed 's/>//' \
  > sulfur/dsrA/dsrA_list.txt
metabolic_HMMs=~/HellsCanyon/references/metabolic_HMMs
cd sulfur/dsrA
hmmalign -o dsrA.sto \
            $metabolic_HMMs/TIGR02064.HMM \
            dsrA.faa
$scripts/convert_stockhold_to_fasta.py dsrA.sto

# Copy in reference sequences
cp ~/references/metabolicProteins/sulfur/dsrA/dsrA_karthik_clean.afa \
    dsrA_karthik_clean.afa

# Align seqs
muscle -profile \
        -in1 dsrA_karthik_clean.afa \
        -in2 dsrA.afa \
        -out dsrA_phylogeny.afa
# Trim alignment

# Generate ML tree
#screen -S EG_dsrA_tree
cd ~/HellsCanyon/dataEdited/metabolic_analyses/sulfur/dsrA/
python $scripts/cleanFASTA.py dsrA_phylogeny_masked.afa
mv -f dsrA_phylogeny_masked.afa_temp.fasta dsrA_phylogeny_masked.afa
FastTree dsrA_phylogeny_masked.afa \
    > dsrA_phylogeny_masked.tree



raxml=/opt/bifxapps/raxml-8.2.11/raxmlHPC-PTHREADS
$raxml -f a \
        -p 283976 \
        -N autoMRE \
        -x 2381 \
        -T 20 \
        -s dsrA_phylogeny_masked.afa \
        -n dsrA


############################################
############################################
# Identify potential PCCs
############################################
############################################

screen -S HCC_MHCs
cd ~/HellsCanyon/dataEdited/metabolic_analyses
mkdir MHCs
cd MHCs
scripts=~/HellsCanyon/code/generalUse
ORFs=~/HellsCanyon/dataEdited/assemblies/ORFs

# Let's look for proteins with at least 3 haem-binding sites
cat ~/HellsCanyon/metadata/lists/assembly_key.csv | while read -r line
do
  assembly=`echo $line | awk -F ',' '{ print $1 }'`
  year=`echo $line | awk -F ',' '{ print $2 }'`
  if [ ! -d $year ]; then
    mkdir $year
  fi
  echo "Searching for MHCs in" $assembly
  $scripts/Find_multiheme_protein.py $ORFs/$assembly.faa 3
  # Move the scripts to the correct directory
  mv $ORFs/$assembly\_3_heme* $year/
done
cat */*3_heme_list* > MHC_list.txt

# Pull out names of adjacent genes
cd ~/HellsCanyon/dataEdited/metabolic_analyses
mkdir BBOMP
rm -f BBOMP/adjacent_genes_all_list.txt
ORFs=~/HellsCanyon/dataEdited/assemblies/ORFs
cat MHCs/MHC_list.txt | while read gene
do
  echo "Working on" $gene
  scaffold=$(echo $gene | rev | cut -d"_" -f2- | rev)
  assembly=$(echo $gene | rev | cut -d"_" -f3- | rev)
  ORFnumber=$(echo $gene | rev | cut -d"_" -f1 | rev)
  preceedingORFnumber=$(expr $ORFnumber - 1)
  followingORFnumber=$(expr $ORFnumber + 1)
  echo $scaffold"_"$preceedingORFnumber >> BBOMP/adjacent_genes_all_list.txt
  echo $scaffold"_"$followingORFnumber >> BBOMP/adjacent_genes_all_list.txt
done

# Find unique gene names
wc -l BBOMP/adjacent_genes_all_list.txt
sort BBOMP/adjacent_genes_all_list.txt | \
  uniq \
  > BBOMP/adjacent_genes_unique_list.txt

# Pull out adjacent genes
rm -f BBOMP/adjacent_genes.faa
cat BBOMP/adjacent_genes_unique_list.txt | while read geneID
do
  assembly=$(echo $geneID | rev | cut -d"_" -f3- | rev)
  echo "Looking for" $geneID "in" $assembly
  grep -A 1 -m 1 $geneID$ $ORFs/$assembly.faa >> BBOMP/adjacent_genes.faa
done


#########################
# Search adjacent genes for BBOMP
#########################
cd ~/HellsCanyon/dataEdited/metabolic_analyses/BBOMP
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=""
PERL5LIB=""

# Run the BBOMP HMM
pcc_omp_HMM=~/HellsCanyon/references/metabolicProteins/EET/pcc_omp.HMM
hmmsearch --tblout pcc_omp_custom.out \
          -T 50 \
          $pcc_omp_HMM \
          adjacent_genes.faa \
          > pcc_omp_custom_output.txt
# Pull out the sequences of interest
scripts=~/HellsCanyon/code/generalUse/
python $scripts/extract_protein_hitting_HMM.py \
        pcc_omp_custom.out \
        adjacent_genes.faa \
        pcc_omp_custom.faa
# Dereplicate sequences
cd-hit -i pcc_omp_custom.faa \
        -o pcc_omp_custom_derep.faa \
        -g 1 \
        -c 0.97 \
        -n 5 \
        -d 0
clstr2txt.pl pcc_omp_custom_derep.faa.clstr > pcc_omp_custom_derep.tsv
grep '>' pcc_omp_custom_derep.faa | sed 's/>//' > pcc_omp_custom_derep_list.txt

# Align sequences to HMM
hmmalign -o pcc_omp_custom_derep.sto \
            $pcc_omp_HMM \
            pcc_omp_custom_derep.faa
# Convert alignment to fasta format
scripts=~/HellsCanyon/code/generalUse
$scripts/convert_stockhold_to_fasta.py pcc_omp_custom_derep.sto
python $scripts/cleanFASTA.py pcc_omp_custom_derep.afa
mv -f pcc_omp_custom_derep.afa_temp.fasta pcc_omp_custom_derep.afa


#########################
# Prepare reference sequences
#########################
cd ~/HellsCanyon/dataEdited/metabolic_analyses/BBOMP
mkdir references
# Find references from RefSeq
blastp -query pcc_omp_custom_derep.faa \
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
cd-hit -g 1 \
        -i refseq_pcc_omp.faa \
        -o refseq_pcc_omp_derep.faa \
        -c 0.97 \
        -n 5 \
        -d 0
rm -f refseq_pcc_omp_clean.faa
grep '>' refseq_pcc_omp_derep.faa | sort | uniq | while read fastaHeader
do
  grep -A 1 -m 1 $fastaHeader refseq_pcc_omp_derep.faa >> refseq_pcc_omp_clean.faa
done
# Dereplicate refseq sequences against reference dataset
cd-hit-2d -i ~/HellsCanyon/references/metabolicProteins/EET/pcc_omp.faa \
          -i2 refseq_pcc_omp_clean.faa \
          -o blast_pcc_omp_uniq.faa \
          -c 0.97 \
          -n 5
# Get final list of refseq references
grep '>' blast_pcc_omp_uniq.faa | \
  sed 's/>//' \
  > blast_pcc_omp_uniq_list.txt
# Retrieve reference information on local computer
HCC
cd dataEdited/metabolic_analyses/BBOMP
epost -db protein -input blast_pcc_omp_uniq_list.txt | \
    esummary | \
    xtract -pattern DocumentSummary -element AccessionVersion,Organism > refseq_bbomp_metadata.tsv


#########################
# Generate a phylogeny
#########################
# Concatenate sequences
cd ~/HellsCanyon/dataEdited/metabolic_analyses/BBOMP
mkdir phylogeny
cat references/blast_pcc_omp_uniq.faa \
    ~/HellsCanyon/references/metabolicProteins/EET/pcc_omp.faa \
    pcc_omp_custom_derep.faa \
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


#########################
# Quantify PCC depth
#########################
source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
cd /home/GLBRCORG/bpeterson26/HellsCanyon/dataEdited/metabolic_analyses/BBOMP
mkdir depth
awk -F '_' '{ print $1"_"$2 }' bbomp_for_abundance.txt | sort | uniq > depth/scaffold_all_list.txt
# Run submission files
chmod +x /home/GLBRCORG/bpeterson26/HellsCanyon/code/executables/aggregate_depth_proteins.sh
condor_submit /home/GLBRCORG/bpeterson26/HellsCanyon/code/submission/aggregate_depth_proteins_bbomp.sub



############################################
############################################
# Run FeGenie HMMs on the assemblies
############################################
############################################

#########################
# Get set up
#########################
mkdir ~/HellsCanyon/references/FeGenie_hmms
cd ~/HellsCanyon/references/FeGenie_hmms
cp /home/GLBRCORG/cnolmsted/FeGenie_EET_Stuff/FeGenie/hmms/iron/iron_reduction/*hmm .
#ls *hmm | sed 's/.hmm//' > ../FeGenie_hmms.txt
echo -e 'ImcH\nCbcA' > ../FeGenie_hmms.txt

#########################
# Identification of reductive EET genes
#########################

screen -S HCC_FeGenie_HMMs
cd ~/HellsCanyon/dataEdited/metabolic_analyses
mkdir FeGenie
mkdir FeGenie/identification

scripts=~/HellsCanyon/code/generalUse
metabolic_HMMs=~/HellsCanyon/references/FeGenie_hmms
ORFs=~/HellsCanyon/dataEdited/assemblies/ORFs

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''

cat $metabolic_HMMs.txt | while read geneName
do
  cd ~/HellsCanyon/dataEdited/metabolic_analyses/FeGenie

  if [ ! -e identification/$geneName.afa ]; then
    echo "Searching for" $geneName
    mkdir identification/$geneName

    cat ~/HellsCanyon/metadata/lists/assembly_key.csv | while read -r line
    do
      assembly=`echo $line | awk -F ',' '{ print $1 }'`
      year=`echo $line | awk -F ',' '{ print $2 }'`
      if [ ! -d identification/$geneName/$year ]; then
        mkdir identification/$geneName/$year
      fi

      if [ ! -e identification/$geneName/$year/$assembly.out ]; then
        echo "Searching for" $geneName "in" $assembly
        hmmsearch --tblout identification/$geneName/$year/$assembly\_$geneName.out \
                  --cpu 4 \
                  -E 1e-30 \
                  $metabolic_HMMs/$geneName.hmm \
                  $ORFs/$assembly.faa \
                  > identification/$geneName/$year/$assembly\_$geneName.txt
        lineCount=`wc -l < identification/$geneName/$year/$assembly\_$geneName.out`
        if [ $lineCount -eq 13 ]; then
          echo "No" $geneName "hits in" $assembly
        else
          echo "Pulling" $geneName "sequences out of" $assembly
          python $scripts/extract_protein_hitting_HMM.py \
                  identification/$geneName/$year/$assembly\_$geneName.out \
                  $ORFs/$assembly.faa \
                  identification/$geneName/$year/$assembly\_$geneName.faa
        fi
      else
        echo "Search for" $geneName "is already done in" $assembly
      fi
    done

    # Aggregate all sequences and align to HMM
    cd identification
    awk -F ',' '{ print $2 }' ~/HellsCanyon/metadata/lists/assembly_key.csv | \
      sort | uniq | \
      while read year
      do
        shopt -s nullglob
        for i in $geneName/$year/*$geneName.faa; do FOUND=$i;break;done

        if [ ! -z $FOUND ]; then
          echo "Concatenating and aligning" $geneName
          cat $geneName/$year/*$geneName.faa \
              > $geneName\_from_$year\_all.faa
          hmmalign -o $geneName\_from_$year.sto \
                      $metabolic_HMMs/$geneName.hmm \
                      $geneName\_from_$year\_all.faa
          $scripts/convert_stockhold_to_fasta.py $geneName\_from_$year.sto
          grep '>' $geneName\_from_$year\_all.faa | \
            sed 's/>//' \
            > $geneName\_from_$year\_all_list.txt
        else
          echo "No" $geneName "sequences found at all :("
        fi
        FOUND=""
      done

    grep -h -v "#" $geneName/*/*$geneName.out | sort -r -k 6 > $geneName\_scores.txt

  else
    echo "Already pulled out" $geneName "sequences"
  fi
done



#########################
# Extract depths of all scaffolds
#########################

# Generate list of all scaffolds that have potential metabolic genes on them
cd ~/HellsCanyon/dataEdited/metabolic_analyses/FeGenie
mkdir depth
cat identification/*_from_201*_all_list.txt | \
  cut -d"_" -f1,2 | \
  sort | \
  uniq \
  > depth/scaffold_all_list.txt

# Pull out all depths
condor_submit /home/GLBRCORG/bpeterson26/HellsCanyon/code/submission/aggregate_depth_proteins_FeGenie.sub



#########################
# Dereplicate sequences
#########################
cd ~/HellsCanyon/dataEdited/metabolic_analyses/FeGenie
mkdir dereplication
cdhit=~/programs/cdhit-master

echo "geneID,geneName" > FeGenie_gene_key.csv

cat $metabolic_HMMs.txt | while read geneName
do
  cd ~/HellsCanyon/dataEdited/metabolic_analyses/FeGenie
  awk -F ',' '{ print $2 }' ~/HellsCanyon/metadata/lists/assembly_key.csv | \
    sort | uniq | \
    while read year
    do
      $cdhit/cd-hit -g 1 \
                    -i identification/$geneName\_from_$year\_all.faa \
                    -o dereplication/$geneName\_from_$year.faa \
                    -c 0.97 \
                    -n 5 \
                    -d 0
      $cdhit/clstr2txt.pl dereplication/$geneName\_from_$year.faa.clstr \
        > dereplication/$geneName\_from_$year.tsv
    done
  cd dereplication
  cat $geneName\_from_*.faa > $geneName\_derep.faa
  grep '>' $geneName\_derep.faa | \
    sed 's/>//' \
    > $geneName\_derep_list.txt
  awk -v geneName="$geneName" '{ print $0","geneName }' $geneName\_derep_list.txt >> ../FeGenie_gene_key.csv
done






############################################
############################################
# Inspection of nirS gene
############################################
############################################

cd ~/HellsCanyon/dataEdited/metabolic_analyses/identification
mkdir ../nirS
grep -h -v "#" nirS/201*/*_nirS.out | sort -r -k 6 > ../nirS/nirS_scores.txt

source /home/GLBRCORG/bpeterson26/miniconda3/etc/profile.d/conda.sh
conda activate bioinformatics
PYTHONPATH=''
PERL5LIB=''
