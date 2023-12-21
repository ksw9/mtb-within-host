###################################
##### Process HH Mtb seq data #####
###################################

# Request node: 
salloc -N 1 -n 1 -t 2:00:00 --mem=48G

#### Generate preliminary data for within-host Mtb study ####
VARCALL_DIR=/uufs/chpc.utah.edu/common/home/walter-group2/
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

# Download SRA data with nf-core
nextflow run nf-core/fetchngs --input input_PRJNA475130.csv --outdir PRJNA475130 -profile singularity --nf_core_pipeline viralrecon --force_sratools_download #done
nextflow run nf-core/fetchngs --input input_PRJNA413593.csv --outdir PRJNA413593 -profile singularity --nf_core_pipeline viralrecon --force_sratools_download #done
nextflow run nf-core/fetchngs --input input_PRJNA549270.csv --outdir PRJNA549270 -profile singularity --nf_core_pipeline viralrecon --force_sratools_download #done

# Capture studies
nextflow run nf-core/fetchngs --input input_PRJEB37609.csv --outdir PRJEB37609 -profile singularity --nf_core_pipeline viralrecon --force_sratools_download #done
nextflow run nf-core/fetchngs --input input_PRJEB21685.csv --outdir PRJEB21685 -profile singularity --nf_core_pipeline viralrecon --force_sratools_download # done 
nextflow run nf-core/fetchngs --input input_PRJNA339209.csv --outdir PRJNA339209 -profile singularity --nf_core_pipeline viralrecon --force_sratools_download -resume -r 1.10.0
nextflow run nf-core/fetchngs --input input_PRJEB9206.csv --outdir PRJEB9206 -profile singularity --nf_core_pipeline viralrecon --force_sratools_download -resume -r 1.10.0
nextflow run nf-core/fetchngs --input input_goig_culture.csv --outdir goig_culture -profile singularity --nf_core_pipeline viralrecon --force_sratools_download -r 1.11.0 # not complete

## 2. Run snp calling pipeline.

# Create an input list of files to process in pipeline: Vitoria
input_list=${VARCALL_DIR}/resources/input/vitoria.tsv
echo -e "sample\tfastq_1\tfastq_2\tbatch" > $input_list
for f1 in $(ls -d ${PROJ_DIR}${DATA_DIR}*_1.fastq.gz); do
  f2=${f1/_1.fastq.gz/_2.fastq.gz}
  samp=$(basename ${f1/_1.fastq.gz})
  bat=vitoria
  echo -e "$samp\t$f1\t$f2\t$bat"
done >> $input_list

# Create an input list of files to process in pipeline: BC
DATA_DIR=data/bc/
input_list=/${VARCALL_DIR}/resources/input/bc.tsv
echo -e "sample\tfastq_1\tfastq_2\tbatch" > $input_list
for f1 in $(ls -d ${PROJ_DIR}${DATA_DIR}*_1.fastq.gz); do
  f2=${f1/_1.fastq.gz/_2.fastq.gz}
  samp=$(basename ${f1/_1.fastq.gz})
  bat=bc
  echo -e "$samp\t$f1\t$f2\t$bat"
done >> $input_list

# Process shorter list of input files to start--now running shorter list of 108 samples.
hh_list="{PROJ_DIR}/metadata/PRJNA413593_possible_hh_accessions.tsv"
input_list=/${VARCALL_DIR}/resources/input/bc_hh_members.tsv
echo -e "sample\tfastq_1\tfastq_2\tbatch" > $input_list
while read samp; do
  #echo $samp
  f1=/uufs/chpc.utah.edu/common/home/walter-group1/tb/hh/data/bc/${samp}_1.fastq.gz
  f2=${f1/_1.fastq.gz/_2.fastq.gz}
  bat=bc
  echo -e "$samp\t$f1\t$f2\t$bat"
done < "$hh_list" >> $input_list

# Update nextflow.config; run.
cd ${VARCALL_DIR}
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
input_list=${VARCALL_DIR}/resources/input/walker.tsv
echo -e "sample\tfastq_1\tfastq_2\tbatch" > $input_list
for f1 in $(ls -d ${PROJ_DIR}${DATA_DIR}*_1.fastq.gz); do
  f2=${f1/_1.fastq.gz/_2.fastq.gz}
  samp=$(basename ${f1/_1.fastq.gz})
  bat=walker
  echo -e "$samp\t$f1\t$f2\t$bat"
done >> $input_list

# Create an input list of files to process in pipeline: Goig
input_list=${VARCALL_DIR}/resources/input/goig.tsv
DATA_DIR=data/goig/reads_fastq/
echo -e "sample\tfastq_1\tfastq_2\tbatch" > $input_list
for f1 in $(ls -d ${PROJ_DIR}${DATA_DIR}*/*_1.fastq.gz); do
  f2=${f1/_1.fastq.gz/_2.fastq.gz}
  samp=$(basename ${f1/_1.fastq.gz})
  bat=goig
  echo -e "$samp\t$f1\t$f2\t$bat"
done >> $input_list

# Create an input list of files to process in pipeline: capture studies
input_list=${VARCALL_DIR}/resources/input/capture.tsv

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
input_list=${VARCALL_DIR}/resources/input/goig_culture.tsv

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

# Create concatenated msa for walker
cat ../mtb-call2/results/walker/ERR*/fasta/*PPEmask.fa > results/fasta/walker_PPEmask.fa
# Snps only
snp-sites -m results/fasta/walker_PPEmask.fa > results/fasta/walker_PPEmask_snps.fa

## 4. IQTree phylogeny building
iqtree2 -s  results/fasta/vitoria_PPEmask_snps.fa --prefix results/tree/vitoria
iqtree2 -s  results/fasta/bc_PPEmask_snps.fa --prefix results/tree/bc
iqtree2 -s  results/fasta/walker_PPEmask_snps.fa --prefix results/tree/walker

# Redo IQTree fitting all with the same GTR model.
model=GTR+G+ASC
sbatch scripts/iqtree.sh results/fasta/vitoria_PPEmask_snps.fa results/tree/vitoria_GTR ${model}
sbatch scripts/iqtree.sh results/fasta/bc_PPEmask_snps.fa results/tree/bc_GTR ${model}
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
