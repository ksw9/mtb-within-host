###############################
#### Process VCF functions ####
###############################

# 1. Read in VCF tables in an input directory.
# List GATK VCF files to read in. 
read_vcf_tables <- function(vars_list) {
  var_files = readLines(vars_list)
  #var_files = dir(path = Sys.glob(input_directory), pattern = '*gatk.table', full.names = TRUE)
  names(var_files) = str_remove(basename(var_files), "_.+")

  # Read in all tables. 
  vars_table = purrr::map_dfr(var_files, 
                            ~read_tsv(.x) %>% 
                              rename_with(~str_remove(., ".+\\.")) , .id = 'Sample') # Remove sample name appended to column name

  vars_table = vars_table %>% 
    separate_wider_delim(cols = AD, names = c("AD1","AD2"), delim = ",") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(AD1 = as.numeric(AD1), AD2 = as.numeric(AD2), 
                maf = min(AD1,AD2)/DP) %>% 
    mutate(Sample = str_remove(Sample, '_.+'))
  
  return(vars_table)
}

# 2. Filter VCF table. 

# Select snps for analysis
filter_vcf_table <- function(vars_table) {
  snps_table = vars_table %>% filter(TYPE == 'SNP')

  # Add column so that pairwise function works and column for dp_flag.
  snps_table <- snps_table %>% group_by(Sample) %>%
  dplyr::mutate(mean_dp = mean(DP), sd_dp = (sd(DP)),
                dp_flag = case_when((DP >= (mean_dp - 2*sd_dp) & DP <= (mean_dp + 2*sd_dp)) & (AD1 >= min_AD & AD2 >= min_AD) ~ 'Expected',
                                     TRUE ~ 'Low depth'), 
                ppe_flag = case_when(!is.na(PPE) ~ 'PE/PPE',
                                     TRUE ~ 'Outside')) 
  return(snps_table)
  
}

# 3. Function to left join single sample df with full variants df. (This will include all duplicates.)
left_join_snps <- function(sample_name,variants) {
  print(sample_name)
  variants %>% filter(Sample == sample_name) %>% left_join(variants, by = c('POS' = 'POS', 'ALT' = 'ALT'), multiple = "all") %>%
    rowwise() %>%
    mutate(pair_name = paste(sort(c(Sample.x,Sample.y)), collapse = '_')) %>%
    distinct()
}

# 4. Summarize shared iSNVs
thresh_fun <- function(df, threshold){
  df %>%
    distinct_at(vars(pair_name,POS), .keep_all = TRUE) %>% # Need to filter out duplicates
    ungroup() %>%
    group_by(pair_name, ppe_flag, dp_flag, pair_id) %>% 
    filter(maf.x >= threshold, maf.y >= threshold) %>% 
    count()
}
