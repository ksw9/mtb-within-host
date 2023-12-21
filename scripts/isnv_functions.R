### Functions to source ####
##### Function to read in VCF, convert to tibble and bind columns for the gt and fix slots.####
vcf_to_tbl <- function(input_vcf){
  tmp = vcfR2tidy(read.vcfR(input_vcf))
  tmp_tbl = cbind(tmp$fix,tmp$gt)
  return(tmp_tbl)
}
#### Summarize data. ####
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  # library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


#### Get pairwise overlaps of iSNVs with shared ALT allele. ####
calc_pairwise_overlaps <- function(alt_freq_threshold, input_vars) {
  #pb$tick()
  print(alt_freq_threshold)
  n_sets <- length(unique(input_vars$Sample))
  set_names <- unique(input_vars$Sample)
  
  # Filter by alt_freq_threshold. Find matches on POS + ALT ALLELE
  df_filt = input_vars %>% group_by(Sample) %>% arrange(Sample) %>% filter(maf >= alt_freq_threshold) %>% 
    distinct() #%>% select(POS,ALT,Sample)%>% distinct()
  
  
  if (dim(df_filt)[1] > 0 ){
    vec_name1 <- character()
    vec_name2 <- character()
    vec_num_shared <- integer()
    
    for (i in seq_len(n_sets)) {
      #print(i)
      name1 <- set_names[i]
      set1 <- df_filt %>% filter(Sample == name1) 
      for (j in seq(i, n_sets)) { # updated so we get self-comparison.
        name2 <- set_names[j]
        set2 <- df_filt %>% filter(Sample == name2)
        set_intersection <- inner_join(set1,set2, by = c('POS','ALT'))
        num_shared <- dim(set_intersection)[1]
        
        vec_name1 <- c(vec_name1, name1)
        vec_name2 <- c(vec_name2, name2)
        vec_num_shared <- c(vec_num_shared, num_shared)
      }
    }
    
    result <- data.frame(name1 = vec_name1,
                         name2 = vec_name2,
                         num_shared = vec_num_shared,
                         threshold = alt_freq_threshold,
                         stringsAsFactors = FALSE) } 
  else{ combs <- data.frame(t(combn(set_names, 2)))
  result <- data.frame(name1 = combs$X1,
                       name2 = combs$X2,
                       num_shared = 0,
                       threshold = alt_freq_threshold,
                       stringsAsFactors = FALSE) 
  }
  return(result)
 # pb$tick()$print()
}

#### Function to count number of shared iSNVs.Doesn't match based on ALT allele.####
calc_pairwise_overlaps_no_matching <- function(alt_freq_threshold) {
  
  # Create sets - need to retain set even when empty
  df_filt = vars_all %>% select(!GFF_FEATURE) %>% group_by(Sample) %>% filter(ALT_FREQ > alt_freq_threshold) %>% distinct()
  sets = split(df_filt$POS,df_filt$Sample)
  n_sets <- length(sets)
  set_names <- names(sets)
  
  vec_name1 <- character()
  vec_name2 <- character()
  vec_num_shared <- integer()
  
  for (i in seq_len(n_sets - 1)) {
    name1 <- set_names[i]
    set1 <- sets[[i]]
    for (j in seq(i, n_sets)) { # updated so we get self-comparison.
      name2 <- set_names[j]
      set2 <- sets[[j]]
      set_intersection <- intersect(set1, set2)
      num_shared <- length(set_intersection)
      
      vec_name1 <- c(vec_name1, name1)
      vec_name2 <- c(vec_name2, name2)
      vec_num_shared <- c(vec_num_shared, num_shared)
    }
  }
  
  result <- data.frame(name1 = vec_name1,
                       name2 = vec_name2,
                       num_shared = vec_num_shared,
                       threshold = alt_freq_threshold,
                       stringsAsFactors = FALSE)
  return(result)
}
# # Code from : https://blog.jdblischak.com/posts/pairwise-overlaps/


#### Function to output all shared iSNVs with no filtering; match based on ALT allele. ####

# Function to get intersection iSNVs by POS and ALT allele.
# identify_shared_isnvs <-function(name1,name2, vars_df,p){
#   set1 <- vars_df %>% filter(Sample == name1) 
#   set2 <- vars_df %>% filter(Sample == name2) 
#   set_intersection <- inner_join(set1,set2, by = c('POS','ALT'))
#   return(set_intersection)
# }

# Use datatable to find intersection to find intersection iSNVS to speed things. Use p for progress bar.
identify_shared_isnvs <- function(name1,name2, vars_df, p = p){
  p() # update progressr
  #print(paste(name1,name2))
  set1 <- vars_df[vars_df$Sample == name1,]
  set2 <- vars_df[vars_df$Sample == name2,]
  set_intersection <- merge(set1,set2, by = c('POS','ALT'))
  #set_intersection <- set_intersection[set_intersection$ALT == set_intersection$i.ALT,]
  return(set_intersection)
}

# Use datatable to find intersection to find transition sites.  Use p for progress bar.
identify_fixation_sites <- function(name1,name2, vars_df, p = p){
  p() # update progressr
  #print(paste(name1,name2))
  set1 <- vars_df[vars_df$Sample == name1,]
  set2 <- vars_df[vars_df$Sample == name2,]
  set_intersection <- merge(set1,set2, by = c('POS','ALT'))
  set_intersection <- set_intersection[(set_intersection$ALT_FREQ.x <= .5 & set_intersection$ALT_FREQ.y >= .7) | 
                                         (set_intersection$ALT_FREQ.y <= .5 & set_intersection$ALT_FREQ.x >= .7),]
  return(set_intersection)
}
#name1 = 'ST038_01_D5_S37';name2='ST038_02_D1_S4'

# Function to filter by maf and then summarize. 
stats_by_threshold = function(alt_threshold, input, all_pairs = all_pairs){
  # Summarize iSNVs by pair.
  isnvs_passing_threshold = input %>%
    ungroup() %>% rowwise() %>%
    filter(maf.x >= alt_threshold & maf.y >= alt_threshold) %>%
    mutate(sample_names = paste(sort(c(Sample.x,Sample.y)),collapse = ';')) %>%
    group_by(sample_names) %>%
    dplyr::summarize(sum_maf.x = sum(maf.x), 
                     sum_maf.y = sum(maf.y),
                     gmean_maf = sqrt(sum_maf.x*sum_maf.y),
                     n_isnvs = n())
  # Add back pairs with 0 iSNVs passing filter.
  out = all_pairs %>%
    left_join(isnvs_passing_threshold,by = c('sample_names')) %>%
    replace_na(list(sum_maf.x = 0, sum_maf.y = 0, gmean_maf = 0, n_isnvs = 0))
  
  # Add warning
  if(dim(out)[1] != dim(all_pairs)[1]) warning('Incomplete output!')
  return(out)
}

# Get number of fixed sites by threshold .
fixation_by_threshold = function(alt_threshold, input, all_pairs = all_pairs){
  # Summarize iSNVs by pair.
  isnvs_passing_threshold = input %>%
    ungroup() %>% rowwise() %>%
    filter(maf.x >= alt_threshold & maf.y >= alt_threshold) %>%
    group_by(sample_names) %>%
    dplyr::summarize(n_fixed= n())
  # Add back pairs with 0 iSNVs passing filter.
  out = all_pairs %>%
    left_join(isnvs_passing_threshold,by = c('sample_names')) %>%
    replace_na(list(n_fixed = 0))
  
  # Add warning
  if(dim(out)[1] != dim(all_pairs)[1]) warning('Incomplete output!')
  return(out)
}

#### Calculate pairwise sum shared MAF ####
calc_pairwise_diversity <- function(alt_freq_threshold) {
  print(alt_freq_threshold)
  # Find matches on POS + ALT ALLELE
  df_filt = vars_all %>% select(!GFF_FEATURE) %>% group_by(Sample) %>% arrange(Sample) %>% filter(maf > alt_freq_threshold) %>% distinct() #%>% select(POS,ALT,Sample)%>% distinct()
  
  n_sets <- length(unique(df_filt$Sample))
  set_names <- unique(df_filt$Sample)
  
  vec_name1 <- character()
  vec_name2 <- character()
  vec_num_shared <- integer()
  vec_sum_maf.x <- double()
  vec_sum_maf.y <- double()
  vec_shannon.x <- double()
  vec_shannon.y <- double()
  vec_simpson.x <- double()
  vec_simpson.y <- double()
  
  for (i in seq_len(n_sets - 1)) {
    #print(i)
    name1 <- set_names[i]
    set1 <- df_filt %>% filter(Sample == name1) 
    for (j in seq(i, n_sets)) { # updated so we get self-comparison.
      name2 <- set_names[j]
      set2 <- df_filt %>% filter(Sample == name2)
      set_intersection <- inner_join(set1,set2, by = c('POS','ALT'))
      num_shared <- dim(set_intersection)[1]
      
      # Additional metrics
      sum_maf.x = sum(set_intersection$maf.x)
      sum_maf.y = sum(set_intersection$maf.y)
      
      shannon.x = -sum(set_intersection$maf.x/sum_maf.x*log(set_intersection$maf.x/sum_maf.x))
      shannon.y = -sum(set_intersection$maf.y/sum_maf.y*log(set_intersection$maf.y/sum_maf.y))
      simpson.x = 1/sum((set_intersection$maf.x/sum_maf.x)^2)
      simpson.y = 1/sum((set_intersection$maf.y/sum_maf.y)^2)
      
      vec_name1 <- c(vec_name1, name1)
      vec_name2 <- c(vec_name2, name2)
      vec_num_shared <- c(vec_num_shared, num_shared)
      vec_sum_maf.x <- c(vec_sum_maf.x,sum_maf.x)
      vec_sum_maf.y <- c(vec_sum_maf.y,sum_maf.y)
      vec_shannon.x <- c(vec_shannon.x,shannon.x)
      vec_shannon.y <- c(vec_shannon.y,shannon.y)
      vec_simpson.x <- c(vec_simpson.x,simpson.x)
      vec_simpson.y <- c(vec_simpson.y,simpson.y)
      
      
    }
  }
  
  result <- data.frame(name1 = vec_name1,
                       name2 = vec_name2,
                       num_shared = vec_num_shared,
                       sum_maf.x = vec_sum_maf.x,
                       sum_maf.y = vec_sum_maf.y,
                       shannon.x = vec_shannon.x,
                       shannon.y = vec_shannon.y,
                       simpson.x = vec_simpson.x,
                       simpson.y = vec_simpson.y,
                       threshold = alt_freq_threshold,
                       stringsAsFactors = FALSE)
  return(result)
}
#### Calculate probability of a household pair for more than X+ iSNVs ####
calc_prob_hh_pair <- function(num_isnvs, pair_df){
  tt <- pair_df %>% 
    #filter(dist <=1 ) %>% 
    filter(!pair_id %in% c('sample','participant', 'Sample','Within-host')) %>%
    filter(n_isnvs >= num_isnvs) %>%
    ungroup() %>% group_by(threshold) %>%
    dplyr::summarize(n = n(),
                     n_hh = length(which(pair_id %in% c('Household','Transmission cluster','Transmission pair'))),
                     prob_hh = n_hh/n, 
                     isnv_threshold = num_isnvs)
  return(tt)
}
calc_prob_hh_pair_maf <- function(maf, pair_df){
  tt <- pair_df %>% 
    #filter(dist <=1 ) %>% 
    filter(!pair_id %in% c('sample','participant', 'Sample','Within-host')) %>%
    filter(gmean_maf >= maf) %>%
    ungroup() %>% group_by(threshold) %>%
    dplyr::summarize(n = n(),
                     n_hh = length(which(pair_id %in% c('Household','Transmission cluster','Transmission pair'))),
                     prob_hh = n_hh/n, 
                     maf_threshold = maf)
  return(tt)
}
#### Calculate sensitivity of a household pair for more than X+ iSNVs ####
calc_sens_hh_pair <- function(num_isnvs){
  tt <- pair_model_df %>% filter(hh ==1) %>%
    ungroup() %>% group_by(threshold) %>%
    dplyr::summarize(n = n(),
                     n_pass_threshold = length(which(num_shared >= num_isnvs)),
                     sens =n_pass_threshold/n, 
                     isnv_threshold = num_isnvs)
  return(tt)
}

#### Calculate specificity of a household pair for more than X+ iSNVs ####
calc_spec_hh_pair <- function(num_isnvs){
  tt <- pair_model_df %>% filter(hh ==0) %>%
    ungroup() %>% group_by(threshold) %>%
    dplyr::summarize(n = n(),
                     n_pass_threshold = length(which(num_shared >= num_isnvs)),
                     prob_pass_threshold = n_pass_threshold/n, 
                     spec = 1-prob_pass_threshold,
                     isnv_threshold = num_isnvs)
  return(tt)
}



#### Jaccard concordance ####
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}
