#!/usr/bin/env Rscript

# Main steps:
# Combine shared iSNVs with metadata
# Summarize

# RUN:
# Rscript scripts/1-summarize_shared_isnvs.R 

# Source config file
source('scripts/config.R')

#### Vitoria ####
# Earlier error wrt to duplicates by pair_name & POS.

meta_vitoria_file = here('metadata/vitoria.csv')
shared_snps_table_vitoria = 'results/processed/vitoria_shared_snps.csv'
meta_vitoria = read_csv(meta_vitoria_file)
shared_snps_table_vitoria = vroom(shared_snps_table_vitoria)
shared_snps_table_vitoria_dt = lazy_dt(shared_snps_table_vitoria)

shared_snps_table_vitoria_ann <- shared_snps_table_vitoria_dt %>% 
  left_join(meta_vitoria[,c("Index case","HHC","months","pair_name")], by = "pair_name") %>%
  mutate(pair_id = case_when(Sample.x == Sample.y ~ 'Sample',
                            !is.na(months) ~ 'Household',
                             TRUE ~ 'Outside household'))
# Write annotated shared snps
write_csv(shared_snps_table_vitoria_ann %>% as_tibble(), file = here('results/processed/vitoria_shared_snps_ann.csv'))

# Summarize flags
shared_snps_summary_vitoria <- map_df(alt_freq_thresholds, ~thresh_fun(as_tibble(shared_snps_table_vitoria_ann), .), .id = 'threshold')

shared_snps_summary_vitoria <- shared_snps_summary_vitoria %>% 
  ungroup() %>%
  complete(nesting(pair_name,pair_id), threshold, ppe_flag, dp_flag, fill = list(n = 0)) %>% 
  mutate_at(vars(threshold), .funs = as.numeric)

# Check correct dim
length(unique(shared_snps_summary_vitoria$pair_name))*length(alt_freq_thresholds)*4
dim(shared_snps_summary_vitoria)
write_csv(shared_snps_summary_vitoria, file = 'results/processed/shared_snps_summary_vitoria.csv')

#### BC ####
meta_bc_file = 'metadata/bc.csv'
shared_snps_table_bc = 'results/processed/bc_shared_snps.csv'
meta_bc = read_csv(meta_bc_file)
shared_snps_table_bc = vroom(shared_snps_table_bc)
shared_snps_table_bc_dt = lazy_dt(shared_snps_table_bc)

shared_snps_table_bc_ann = shared_snps_table_bc_dt %>% 
  left_join(meta_bc[,c('Run','library','hh')], by = c('Sample.x' = 'Run')) %>%
  left_join(meta_bc[,c('Run','library','hh')], by = c('Sample.y' = 'Run')) %>%
  dplyr::rename(sample_name.x = library.x, sample_name.y= library.y,house.x = hh.x, house.y = hh.y) %>%
  mutate(pair_id = case_when(sample_name.x == sample_name.y ~ 'Sample',
                             house.x == house.y ~ 'Household',
                             TRUE ~ 'Outside household'))

# Write annotated shared snps
write_csv(shared_snps_table_bc_ann %>% as_tibble(), file = here('results/processed/bc_shared_snps_ann.csv'))

# Summarize flags
shared_snps_summary_bc <- map_df(alt_freq_thresholds, ~thresh_fun(as_tibble(shared_snps_table_bc_ann), .), .id = 'threshold')
shared_snps_summary_bc <- shared_snps_summary_bc %>% 
  ungroup() %>%
  complete(nesting(pair_name,pair_id), threshold, ppe_flag, dp_flag, fill = list(n = 0)) %>% 
  mutate_at(vars(threshold), .funs = as.numeric)

# Check correct dim
length(unique(shared_snps_summary_bc$pair_name))*length(alt_freq_thresholds)*4
dim(shared_snps_summary_bc)
write_csv(shared_snps_summary_bc, file = 'results/processed/shared_snps_summary_bc.csv')

#### Lee ####
meta_lee_file = 'metadata/lee.csv'
shared_snps_table_lee = 'results/processed/lee_shared_snps.csv'
meta_lee = read_csv(meta_lee_file)
shared_snps_table_lee = vroom(shared_snps_table_lee)
shared_snps_table_lee_dt = lazy_dt(shared_snps_table_lee)

shared_snps_table_lee_ann = shared_snps_table_lee_dt %>% 
  left_join(meta_lee[,c('Run','SampleName','sublineage','predicted_index')], by = c('Sample.x' = 'Run')) %>%
  left_join(meta_lee[,c('Run','SampleName','sublineage','predicted_index')], by = c('Sample.y' = 'Run')) %>%
  dplyr::rename(sample_name.x = SampleName.x, sample_name.y = SampleName.y) 

# Add information. exclude site: 2996085 (only 2% pass the DP flag anyway)
shared_snps_summary_lee = shared_snps_table_lee_ann %>% 
  mutate(pair_id = case_when(sample_name.x == sample_name.y ~ 'Sample',
                      sublineage.x == sublineage.y ~ 'Sublineage',
                       TRUE ~ 'Outside sublineage'))
    # pair_id = case_when(sample_name.x == sample_name.y ~ 'Sample',
    #                     #sublineage.x == 'IIIB' & sublineage.y == 'IIIB' & (sample_name.x == 'MT-504' | sample_name.y == 'MT-504') ~ 'IIIB: index-contact',
    #                     #sublineage.x == 'IIIB' & sublineage.y == 'IIIB' ~ 'IIIB: other',
    #                     sublineage.x == sublineage.y ~ 'Other sublineage',
    #                     TRUE ~ 'Outside sublineage'))


# Summarize flags
shared_snps_summary_lee <- map_df(alt_freq_thresholds, ~thresh_fun(as_tibble(shared_snps_summary_lee), .), .id = 'threshold')
shared_snps_summary_lee <- shared_snps_summary_lee %>% 
  ungroup() %>%
  complete(nesting(pair_name,pair_id), threshold, ppe_flag, dp_flag, fill = list(n = 0)) %>% 
  mutate_at(vars(threshold), .funs = as.numeric)

# Check correct dim
length(unique(shared_snps_summary_lee$pair_name))*length(alt_freq_thresholds)*4
dim(shared_snps_summary_lee)
write_csv(shared_snps_summary_lee, file = 'results/processed/shared_snps_summary_lee.csv')

#### Walker ####
meta_walker_file = 'metadata/walker.csv'
shared_snps_table_walker = 'results/processed/walker_shared_snps.csv'
meta_walker = read_csv(meta_walker_file)
shared_snps_table_walker = vroom(shared_snps_table_walker, delim = ',') # faster read in
shared_snps_table_walker_dt = lazy_dt(shared_snps_table_walker)

shared_snps_table_walker_ann = shared_snps_table_walker_dt %>%
  left_join(meta_walker %>% 
              dplyr::select(pair_name, epi_link,epi_connection), by = 'pair_name') %>%
  dplyr::mutate(pair_id = case_when(Sample.x == Sample.y ~ 'Sample',
                                    epi_connection == 'Household / family' ~ 'Household',
                                    epi_link == 'Yes' ~ 'Linked',     
                                    TRUE ~ 'Unlinked'))

# Write annotated shared snps
write_csv(shared_snps_table_walker_ann %>% as_tibble(), file = here('results/processed/walker_shared_snps_ann.csv'))

# Summarize flags
shared_snps_summary_walker <- map_df(alt_freq_thresholds, ~thresh_fun(as_tibble(shared_snps_table_walker_ann), .), .id = 'threshold')
shared_snps_summary_walker <- shared_snps_summary_walker %>% 
  ungroup() %>%
  complete(nesting(pair_name,pair_id), threshold, ppe_flag, dp_flag, fill = list(n = 0)) %>% 
  mutate_at(vars(threshold), .funs = as.numeric)

# Check correct dim
length(unique(shared_snps_summary_walker$pair_name))*length(alt_freq_thresholds)*4
dim(shared_snps_summary_walker)
write_csv(shared_snps_summary_walker, file = 'results/processed/shared_snps_summary_walker.csv')


#### Goig ####
metadata_file_goig = 'metadata/Goig2020_sputum_matched_cultures.csv'
shared_snps_table_goig = 'results/processed/goig_shared_snps.csv'
meta_goig = read_csv(metadata_file_goig)
shared_snps_table_goig = vroom(shared_snps_table_goig, delim = ',') # faster read in
shared_snps_table_goig_dt = lazy_dt(shared_snps_table_goig)

shared_snps_table_goig_ann = shared_snps_table_goig %>% 
  left_join(meta_goig[,c('Run','Matched_Culture_Accession','sequencing')], by = c('Sample.x' = 'Run')) %>%
  left_join(meta_goig[,c('Run','Matched_Culture_Accession','sequencing')], by = c('Sample.y' = 'Run')) %>%
  dplyr::rename(sample_name.x = Sample.x, sample_name.y = Sample.y) %>%
  mutate(pair_id = case_when((sample_name.y == Matched_Culture_Accession.x) | (sample_name.x == Matched_Culture_Accession.y) ~ 'Matched pair',
                              sample_name.x == sample_name.y & sequencing.x == 'Bait capture' ~ 'Bait capture',
                             sample_name.x == sample_name.y & sequencing.x == 'Direct sequencing'~ 'Direct sequencing',
                             sample_name.x == sample_name.y ~ 'Culture',
                             TRUE ~ 'Random pair'))

# Summarize flags
shared_snps_summary_goig <- map_df(alt_freq_thresholds, ~thresh_fun(as_tibble(shared_snps_table_goig_ann), .), .id = 'threshold')
shared_snps_summary_goig <- shared_snps_summary_goig %>% 
  ungroup() %>%
  complete(nesting(pair_name,pair_id), threshold, ppe_flag, dp_flag, fill = list(n = 0)) %>% 
  mutate_at(vars(threshold), .funs = as.numeric)

# Check correct dim
length(unique(shared_snps_summary_goig$pair_name))*length(alt_freq_thresholds)*4
dim(shared_snps_summary_goig)
write_csv(shared_snps_summary_goig, file = 'results/processed/shared_snps_summary_goig.csv', col_names = TRUE)
