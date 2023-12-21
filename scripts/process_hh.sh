###################################
##### Process HH Mtb seq data #####
###################################

# Request node: 
salloc -N 1 -n 1 -t 2:00:00 --mem=48G

#### Generate preliminary data for within-host Mtb study ####
PROJ_DIR=/uufs/chpc.utah.edu/common/home/walter-group2/tb/hh/
SCRIPTS_DIR=${PROJ_DIR}scripts/
DATA_DIR=${PROJ_DIR}data/
META_DIR=${PROJ_DIR}metadata/
# Define Singularity cache.
NXF_SINGULARITY_CACHEDIR=${PROJ_DIR}cache/
export NXF_SINGULARITY_CACHEDIR

cd ${DATA_DIR}

module load gatk/4.1
module load htslib
module load iqtree/2.2.2.4
module load nextflow; module load singularity

## 0 (updated). Download data from capture studies with Nextflow fetchngs. 
# Test data
nextflow run nf-core/fetchngs -profile test,singularity --outdir test --force_sratools_download # works

# Download single capture study
nextflow run nf-core/fetchngs --input input_PRJEB37609.csv --outdir PRJEB37609 -profile singularity --nf_core_pipeline viralrecon --force_sratools_download #done
nextflow run nf-core/fetchngs --input input_PRJEB21685.csv --outdir PRJEB21685 -profile singularity --nf_core_pipeline viralrecon --force_sratools_download # done 
nextflow run nf-core/fetchngs --input input_PRJNA339209.csv --outdir PRJNA339209 -profile singularity --nf_core_pipeline viralrecon --force_sratools_download -resume -r 1.10.0
nextflow run nf-core/fetchngs --input input_PRJEB9206.csv --outdir PRJEB9206 -profile singularity --nf_core_pipeline viralrecon --force_sratools_download -resume -r 1.10.0
nextflow run nf-core/fetchngs --input input_goig_culture.csv --outdir goig_culture -profile singularity --nf_core_pipeline viralrecon --force_sratools_download -r 1.11.0 # not complete

# Check which files are missing
files=$(awk '{print $1}' PRJEB9206/metadata/PRJEB9206.runinfo_ftp.tsv  | sed '1d')
shopt -s failglob
for file in $files; do 
  
  if echo PRJEB9206/fastq/${file}* &>/dev/null
  then
    echo 'exists'
   else
    echo ${file} missing
  fi
done




## 0. Create conda environment for parallel-fastq-dump
#module load anaconda
#mamba install -c bioconda parallel-fastq-dump -n parallel-fastq-dump -m
#mamba install -c bioconda sra-tools=2.10 -n parallel-fastq-dump

## 1. Download Colangeli et al. 2020 data. SRA:PRJNA475130. (Pair identifier is Library name.)
DATA_DIR=data/vitoria/
while read sra; do
  if [ -e $DATA_DIR${sra}_1.fastq.gz ]; then echo '';
  else  echo -e 'submit job' $sra
  sbatch ${SCRIPTS_DIR}download_sra.sh ${sra} ${DATA_DIR}
  fi
done < ${META_DIR}PRJNA475130_accessions.txt

## 1b. Download BC data. SRA:PRJNA413593
DATA_DIR=data/bc/
while read sra; do
  if [ -e $DATA_DIR${sra}_1.fastq.gz ]; then echo '';
  else  echo -e 'submit job' $sra
  sbatch ${SCRIPTS_DIR}download_sra.sh ${sra} ${DATA_DIR}
  fi
done < ${META_DIR}PRJNA413593_accessions.txt

## 1c. Download Lee data. (Pair identifier is Library name.)
DATA_DIR=data/lee/
while read sra; do
  if [ -e $DATA_DIR${sra}_1.fastq.gz ]; then echo '';
  else  echo -e 'submit job' $sra
  sbatch ${SCRIPTS_DIR}download_sra.sh ${sra} ${DATA_DIR}
  fi
done < ${META_DIR}PRJNA549270_accessions.txt
# Note: PRJNA549270, SRR9588640 is the PacBio sequence data, don't need for ref genome generation, don't include
# Run fastq-dump
while read sra; do
  # Delete carriage return
  sra_mod=$(echo $sra | tr -d '\r')
  #echo $sra_mod
  if [ -e $DATA_DIR${sra}_1.fastq.gz ]; then echo $sra_mod;
 #  fasterq-dump ${sra_mod} --outdir ${DATA_DIR} --split-files --progress
  fi
done < ${META_DIR}PRJNA549270_accessions.txt

## 1d. Download Walker data.
# ssh to the DTN
ssh u6045141@dtn05.hpc.utah.edu
source activate parallel-fastq-dump

while read line; do
  echo $line
  $line &
done < /uufs/chpc.utah.edu/common/home/walter-group1/tb/hh/scripts/ena-file-download-20230706-2237.sh

## 1e. Download Goig data. 
# Use ENA FTP downloader
DATA_DIR=data/goig/
accessions_list1=$(awk -F, 'NR!=1 {print $1}' metadata/Goig2020_sputum_matched_cultures.csv | paste -sd,)
accessions_list2=$(awk -F, 'NR!=1 {print $43}' metadata/Goig2020_sputum_matched_cultures.csv | grep -v NA | grep -v -e '^$' |  paste -sd,)
accessions_list=${accessions_list1},${accessions_list2}

java -jar scripts/ena-file-downloader.jar --accessions=$accessions_list --format=READS_FASTQ \
--location=${PROJ_DIR}${DATA_DIR} --protocol=FTP --asperaLocation=null \
--email=katharine.walter@hsc.utah.edu

## 1f. Download Doyle - done
DATA_DIR=data/doyle/
java -jar scripts/ena-file-downloader.jar --accessions=PRJEB21685 --format=READS_FASTQ \
--location=${PROJ_DIR}${DATA_DIR} --protocol=FTP --asperaLocation=null \
--email=katharine.walter@hsc.utah.edu


## 1g. Download Brown - running
DATA_DIR=data/brown/
java -jar scripts/ena-file-downloader.jar --accessions=PRJEB9206 --format=READS_FASTQ \
--location=${PROJ_DIR}${DATA_DIR} --protocol=FTP --asperaLocation=null \
--email=katharine.walter@hsc.utah.edu

## 1h. Download Votintseva - running
DATA_DIR=data/votintseva
java -jar scripts/ena-file-downloader.jar --accessions=PRJNA339209 --format=READS_FASTQ \
--location=${PROJ_DIR}${DATA_DIR} --protocol=FTP --asperaLocation=null \
--email=katharine.walter@hsc.utah.edu

## 2. Run snp calling pipeline.

# Create an input list of files to process in pipeline: Vitoria
input_list=/uufs/chpc.utah.edu/common/home/walter-group1/tb/mtb-call2/resources/input/vitoria.tsv
echo -e "sample\tfastq_1\tfastq_2\tbatch" > $input_list
for f1 in $(ls -d ${PROJ_DIR}${DATA_DIR}*_1.fastq.gz); do
  f2=${f1/_1.fastq.gz/_2.fastq.gz}
  samp=$(basename ${f1/_1.fastq.gz})
  bat=vitoria
  echo -e "$samp\t$f1\t$f2\t$bat"
done >> $input_list

# Create an input list of files to process in pipeline: BC
DATA_DIR=data/bc/
input_list=/uufs/chpc.utah.edu/common/home/walter-group1/tb/mtb-call2/resources/input/bc.tsv
echo -e "sample\tfastq_1\tfastq_2\tbatch" > $input_list
for f1 in $(ls -d ${PROJ_DIR}${DATA_DIR}*_1.fastq.gz); do
  f2=${f1/_1.fastq.gz/_2.fastq.gz}
  samp=$(basename ${f1/_1.fastq.gz})
  bat=bc
  echo -e "$samp\t$f1\t$f2\t$bat"
done >> $input_list

# Process shorter list of input files to start--now running shorter list of 108 samples.
hh_list="/uufs/chpc.utah.edu/common/home/walter-group1/tb/hh/metadata/PRJNA413593_possible_hh_accessions.tsv"
input_list=/uufs/chpc.utah.edu/common/home/walter-group1/tb/mtb-call2/resources/input/bc_hh_members.tsv
echo -e "sample\tfastq_1\tfastq_2\tbatch" > $input_list
while read samp; do
  #echo $samp
  f1=/uufs/chpc.utah.edu/common/home/walter-group1/tb/hh/data/bc/${samp}_1.fastq.gz
  f2=${f1/_1.fastq.gz/_2.fastq.gz}
  bat=bc
  echo -e "$samp\t$f1\t$f2\t$bat"
done < "$hh_list" >> $input_list

# Update nextflow.config; run.
cd /uufs/chpc.utah.edu/common/home/walter-group1/tb/mtb-call2
module load singularity # or docker (not necessary on Stanford SCG)
module load java nextflow

sbatch scripts/submit_mtb_pipeline.sh

# Get mean cov for BC samples (need to redo the nextflow run, excluding the missing sample because the summary file isn't generated if one samples fails)
samps=$(cat resources/input/bc_hh_members.tsv | tail -n26 | awk '{print $1}')
for samp in $samps; do
  #echo $samp
  cat results/bc/${samp}/stats/${samp}_coverage_stats.txt  | grep -v '#' | head -n3 | awk '{print $2'} | tail -n1
done > mean_cov.out

# Sample list for Walker files
DATA_DIR=data/walker/
input_list=/uufs/chpc.utah.edu/common/home/walter-group1/tb/mtb-call2/resources/input/walker.tsv
echo -e "sample\tfastq_1\tfastq_2\tbatch" > $input_list
for f1 in $(ls -d ${PROJ_DIR}${DATA_DIR}*_1.fastq.gz); do
  f2=${f1/_1.fastq.gz/_2.fastq.gz}
  samp=$(basename ${f1/_1.fastq.gz})
  bat=walker
  echo -e "$samp\t$f1\t$f2\t$bat"
done >> $input_list

# Create an input list of files to process in pipeline: Goig
input_list=/uufs/chpc.utah.edu/common/home/walter-group2/tb/mtb-call2/resources/input/goig.tsv
DATA_DIR=data/goig/reads_fastq/
echo -e "sample\tfastq_1\tfastq_2\tbatch" > $input_list
for f1 in $(ls -d ${PROJ_DIR}${DATA_DIR}*/*_1.fastq.gz); do
  f2=${f1/_1.fastq.gz/_2.fastq.gz}
  samp=$(basename ${f1/_1.fastq.gz})
  bat=goig
  echo -e "$samp\t$f1\t$f2\t$bat"
done >> $input_list

# Create an input list of files to process in pipeline: capture studies
input_list=/uufs/chpc.utah.edu/common/home/walter-group2/tb/mtb-call2/resources/input/capture.tsv

echo -e "sample\tfastq_1\tfastq_2\tbatch" > $input_list
for project in PRJNA339209 PRJEB9206 PRJEB21685 PRJEB37609 ; do

for f1 in $(ls -d ${PROJ_DIR}data/${project}/fastq/*_1.fastq.gz); do
  f2=${f1/_1.fastq.gz/_2.fastq.gz}
  samp=$(basename ${f1/_1.fastq.gz})
  bat=${project}
  echo -e "$samp\t$f1\t$f2\t$bat"
done >> $input_list

done

# Create an input list of files to process in pipeline: Goig culture studies
input_list=/uufs/chpc.utah.edu/common/home/walter-group2/tb/mtb-call2/resources/input/goig_culture.tsv

echo -e "sample\tfastq_1\tfastq_2\tbatch" > $input_list

for f1 in $(ls -d ${PROJ_DIR}data/goig_culture/fastq/*_1.fastq.gz); do
  f2=${f1/_1.fastq.gz/_2.fastq.gz}
  samp=$(basename ${f1/_1.fastq.gz})
  bat=${project}
  echo -e "$samp\t$f1\t$f2\t$bat"
done >> $input_list

# Run pipeline by modifying nextflow.config.
sbatch scripts/submit_mtb_pipeline.sh 

## 3. Fasta analysis
## 3.0 Create environment phylo
module load iqtree/2.2.0
#mamba create -c bioconda snp-sites -n phylo
source activate phylo

# Create concatenated msa for vitoria
cat ../mtb-call2/results/vitoria/SRR7276*/fasta/*PPEmask.fa > results/fasta/vitoria_PPEmask.fa
# Snps only
snp-sites -m results/fasta/vitoria_PPEmask.fa > results/fasta/vitoria_PPEmask_snps.fa

# Create concatenated msa for bc
cat ../mtb-call2/results/bc/SRR*/fasta/*PPEmask.fa > results/fasta/bc_PPEmask.fa
# Snps only
snp-sites -m results/fasta/bc_PPEmask.fa > results/fasta/bc_PPEmask_snps.fa

# Create concatenated msa for lee
cat ../mtb-call2/results/lee/SRR*/fasta/*PPEmask.fa > results/fasta/lee_PPEmask.fa
# Snps only
snp-sites -m results/fasta/lee_PPEmask.fa > results/fasta/lee_PPEmask_snps.fa

# Create concatenated msa for walker
cat ../mtb-call2/results/walker/ERR*/fasta/*PPEmask.fa > results/fasta/walker_PPEmask.fa
# Snps only
snp-sites -m results/fasta/walker_PPEmask.fa > results/fasta/walker_PPEmask_snps.fa

## 4. IQTree phylogeny building
iqtree2 -s  results/fasta/vitoria_PPEmask_snps.fa --prefix results/tree/vitoria
iqtree2 -s  results/fasta/bc_PPEmask_snps.fa --prefix results/tree/bc
iqtree2 -s  results/fasta/lee_PPEmask_snps.fa --prefix results/tree/lee
iqtree2 -s  results/fasta/walker_PPEmask_snps.fa --prefix results/tree/walker

# Redo IQTree fitting all with the same GTR model.
model=GTR+G+ASC
sbatch scripts/iqtree.sh results/fasta/vitoria_PPEmask_snps.fa results/tree/vitoria_GTR ${model}
sbatch scripts/iqtree.sh results/fasta/bc_PPEmask_snps.fa results/tree/bc_GTR ${model}
sbatch scripts/iqtree.sh results/fasta/lee_PPEmask_snps.fa results/tree/lee_GTR ${model}
sbatch scripts/iqtree.sh results/fasta/walker_PPEmask_snps.fa results/tree/walker_GTR ${model}

## 5. Create VCF tables.
## Vitoria ##
# Convert variants to table.
vcf_list=results/vars/vcfs_vitoria.list
ls ../mtb-call2/results/vitoria/SRR*/vars/*gatk_ann.vcf.gz > $vcf_list

while read vcf; do
  echo $vcf
  tabix $vcf
  output=${vcf/_ann.vcf.gz/.table}
  gatk VariantsToTable -V $vcf -F CHROM -F POS -F REF -F ALT -F QUAL -F ANN -F FILTER -F TYPE -GF PPE -GF GT -GF AD -GF DP -GF GQ -O ${output} --show-filtered
done < $vcf_list

# List tables
vars_tables_list=results/vars/var_table_vitoria.list
sed 's/_ann.vcf.gz/.table/' $vcf_list > $vars_tables_list

# Convert LoFreq variants to table.
lofreq_list=results/vars/vcfs_lofreq_vitoria.list
ls ../mtb-call2/results/vitoria/SRR*/vars/*_lofreq_ann.vcf.gz > $lofreq_list

while read vcf; do
  echo $vcf
  tabix $vcf
  output=${vcf/_ann.vcf.gz/.table}
  gatk VariantsToTable -V $vcf -F CHROM -F POS -F REF -F ALT -F QUAL -F DP -F DP4 -F ANN -F FILTER -F TYPE -GF PPE -GF GT -GF AD -GF DP -GF GQ -O ${output} --show-filtered
done < $lofreq_list

## BC ##
# Convert variants to table.
vcf_list=results/vars/vcfs_bc.list
ls ../mtb-call2/results/bc/SRR*/vars/*gatk_ann.vcf.gz > $vcf_list

while read vcf; do
  echo $vcf
  tabix $vcf
  output=${vcf/_ann.vcf.gz/.table}
  gatk VariantsToTable -V $vcf -F CHROM -F POS -F REF -F ALT -F QUAL -F ANN -F FILTER -F TYPE -GF PPE -GF GT -GF AD -GF DP -GF GQ -O ${output} --show-filtered
done < $vcf_list

# List tables
vars_tables_list=results/vars/var_table_bc.list
sed 's/_ann.vcf.gz/.table/' $vcf_list > $vars_tables_list

# Convert LoFreq variants to table.
lofreq_list=results/vars/vcfs_lofreq_bc.list
ls ../mtb-call2/results/bc/SRR*/vars/*_lofreq_ann.vcf.gz > $lofreq_list

while read vcf; do
  echo $vcf
  tabix $vcf
  output=${vcf/_ann.vcf.gz/.table}
  gatk VariantsToTable -V $vcf -F CHROM -F POS -F REF -F ALT -F QUAL -F DP -F DP4 -F ANN -F FILTER -F TYPE -GF PPE -GF GT -GF AD -GF DP -GF GQ -O ${output} --show-filtered
done < $lofreq_list

## Walker ##
# Convert variants to table.
vcf_list=results/vars/vcfs_walker.list
ls ../mtb-call2/results/walker/ERR*/vars/*gatk_ann.vcf.gz > $vcf_list

while read vcf; do
  echo $vcf
  tabix $vcf
  output=${vcf/_ann.vcf.gz/.table}
  gatk VariantsToTable -V $vcf -F CHROM -F POS -F REF -F ALT -F QUAL -F ANN -F FILTER -F TYPE -GF PPE -GF GT -GF AD -GF DP -GF GQ -O ${output} --show-filtered
done < $vcf_list

# List tables
vars_tables_list=results/vars/var_table_walker.list
sed 's/_ann.vcf.gz/.table/' $vcf_list > $vars_tables_list

# Convert LoFreq variants to table.
lofreq_list=results/vars/vcfs_lofreq_bc.list
ls ../mtb-call2/results/bc/SRR*/vars/*_lofreq_ann.vcf.gz > $lofreq_list

while read vcf; do
  echo $vcf
  tabix $vcf
  output=${vcf/_ann.vcf.gz/.table}
  gatk VariantsToTable -V $vcf -F CHROM -F POS -F REF -F ALT -F QUAL -F DP -F DP4 -F ANN -F FILTER -F TYPE -GF PPE -GF GT -GF AD -GF DP -GF GQ -O ${output} --show-filtered
done < $lofreq_list

## Lee ##
# Convert variants to table.
vcf_list=results/vars/vcfs_lee.list
ls ../mtb-call2/results/lee/SRR*/vars/*gatk_ann.vcf.gz > $vcf_list

while read vcf; do
  echo $vcf
  tabix $vcf
  output=${vcf/_ann.vcf.gz/.table}
  gatk VariantsToTable -V $vcf -F CHROM -F POS -F REF -F ALT -F QUAL -F ANN -F FILTER -F TYPE -GF PPE -GF GT -GF AD -GF DP -GF GQ -O ${output} --show-filtered
done < $vcf_list

# List tables
vars_tables_list=results/vars/var_table_lee.list
sed 's/_ann.vcf.gz/.table/' $vcf_list > $vars_tables_list

# Convert LoFreq variants to table.
lofreq_list=results/vars/vcfs_lofreq_lee.list
ls ../mtb-call2/results/lee/SRR*/vars/*_lofreq_ann.vcf.gz > $lofreq_list

while read vcf; do
  echo $vcf
  tabix $vcf
  output=${vcf/_ann.vcf.gz/.table}
  gatk VariantsToTable -V $vcf -F CHROM -F POS -F REF -F ALT -F QUAL -F DP -F DP4 -F ANN -F FILTER -F TYPE -GF PPE -GF GT -GF AD -GF DP -GF GQ -O ${output} --show-filtered
done < $lofreq_list

## Goig ##
# Convert variants to table.
vcf_list=results/vars/vcfs_goig.list
ls ../mtb-call2/results/goig/ERR*/vars/*gatk_ann.vcf.gz > $vcf_list

while read vcf; do
  echo $vcf
  tabix $vcf
  output=${vcf/_ann.vcf.gz/.table}
  gatk VariantsToTable -V $vcf -F CHROM -F POS -F REF -F ALT -F QUAL -F ANN -F FILTER -F TYPE -GF PPE -GF GT -GF AD -GF DP -GF GQ -O ${output} --show-filtered
done < $vcf_list

## Capture ##
# Convert variants to table.
vcf_list=results/vars/vcfs_capture.list
ls ../mtb-call2/results/PRJ*/*/vars/*gatk_ann.vcf.gz > $vcf_list
ls ../mtb-call2/results/goig_culture*/*/vars/*gatk_ann.vcf.gz >> $vcf_list

while read vcf; do
  echo $vcf
  tabix $vcf
  output=${vcf/_ann.vcf.gz/.table}
  gatk VariantsToTable -V ${vcf} -F CHROM -F POS -F REF -F ALT -F QUAL -F ANN -F FILTER -F TYPE -GF PPE -GF GT -GF AD -GF DP -GF GQ -O ${output} --show-filtered
done < $vcf_list

# List tables
vars_tables_list=results/vars/var_table_capture.list
sed 's/_ann.vcf.gz/.table/' $vcf_list > $vars_tables_list

## 5. Run process vars R script for each dataset.
module load R/4.2.2

Rscript scripts/0-process_vcfs.R results/vars/var_table_vitoria.list vitoria
Rscript scripts/0-process_vcfs.R results/vars/var_table_bc.list bc
Rscript scripts/0-process_vcfs.R results/vars/var_table_lee.list lee
Rscript scripts/0-process_vcfs.R results/vars/var_table_walker.list walker # requires more memory than others
Rscript scripts/0-process_vcfs.R results/vars/var_table_capture.list capture
#Rscript scripts/0-process_vcfs.R results/vars/var_table_goig.list goig

## Troubleshooting: are hard filters being applied? 
fasta1=/uufs/chpc.utah.edu/common/home/walter-group2/tb/mtb-call2/results/bc/SRR6398092/fasta/SRR6398092_gatk_PPEmask.fa
fasta2=/uufs/chpc.utah.edu/common/home/walter-group2/tb/mtb-call2/results/bc/SRR6397994/fasta/SRR6397994_gatk_PPEmask.fa
vcf1=/uufs/chpc.utah.edu/common/home/walter-group2/tb/mtb-call2/results/bc/SRR6398092/vars/SRR6398092_gatk_filt.vcf.gz 
vcf2=/uufs/chpc.utah.edu/common/home/walter-group2/tb/mtb-call2/results/bc/SRR6397994/vars/SRR6397994_gatk_filt.vcf.gz 

grep -o N $fasta | wc -l # 0 Ns.
zcat $vcf | grep -v '#' | grep lowDepth 

# many positions with lowDepth filter 
# As currently, written the bcftools consensus is not setting these lowDepth sites to N. 
# Proposed a solution/ pipeline update to Marco. 

# Look at fastas, sites of differenence w/ snp-sites, create a vcf
# Look at coverage at these sites
# is this issue because we are under-filtering or because they are over filtering? 
source activate phylo
cat $fasta1 $fasta2 > tmp.fa
snp-sites tmp.fa
less $fasta