#!/usr/bin/env Rscript

# Input:
# arg1: vcf_list
# arg2: output prefix

# Main steps:
# Read in vcf tables
# Filter vcf tables
# Save filtered vcf table
# Pairwise shared snps
# Save pairwise shared snps
# Summarize shared snps by threshold
# Save summarized shared snps

# RUN:
# Rscript scripts/0-process_vcfs.R <input_directory> <output_prefix>

args <- commandArgs(trailingOnly = TRUE)
print(args)

input_directory = args[1]
output_prefix = args[2]

# Source config file
source('scripts/config.R')
source('scripts/process_vcf_functions.R')

# Define output files
snps_table_file = paste0(results_dir,'processed/',output_prefix,'_snps.csv')
shared_snps_table_file = paste0(results_dir,'processed/',output_prefix,'_shared_snps.csv')
shared_snps_summary_file = paste0(results_dir,'processed/',output_prefix,'_shared_snps_summary.csv')

##### Create vars table ####
vars_table = read_vcf_tables(input_directory)

##### Filter vars table #####
snps_table = filter_vcf_table(vars_table)

##### Save table to directory #####
print(snps_table_file)
print(output_prefix)
write_csv(snps_table, file = snps_table_file)

##### Pairwise shared iSNVs #####
sample_names = unique(snps_table$Sample)
shared_snps_table <- map(sample_names, variants = snps_table, .f = left_join_snps, .progress = TRUE) %>%
  map_dfr(~ as.data.frame(.)) %>%
  mutate(dp_flag = case_when(dp_flag.x == 'Expected' & dp_flag.y == 'Expected' ~ 'Expected', 
                                                          TRUE ~ 'Low/high depth')) %>%
  rename(ppe_flag = ppe_flag.x) %>%
  select(-ppe_flag.y)

##### Save table to directory #####
write_csv(shared_snps_table, file = shared_snps_table_file)
