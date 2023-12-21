#####################
#### Config file ####
#####################

library(tidyverse)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(vroom)
library(here)

project_dir = '/uufs/chpc.utah.edu/common/home/walter-group2/tb/hh/'
plot_dir = '/uufs/chpc.utah.edu/common/home/walter-group2/tb/hh/plots/'
scripts_dir = '/uufs/chpc.utah.edu/common/home/walter-group2/tb/hh/scripts/'
meta_dir = '/uufs/chpc.utah.edu/common/home/walter-group2/tb/hh/metadata/'
results_dir = '/uufs/chpc.utah.edu/common/home/walter-group2/tb/hh/results/'
process_dir = '/uufs/chpc.utah.edu/common/home/walter-group2/tb/hh/processed/'
sequencing_dir = '/uufs/chpc.utah.edu/common/home/walter-group2/tb/mtb-call2/results/'

# Source functions
source(paste0(scripts_dir,'process_vcf_functions.R'))

# Set variables
min_AD = 5
alt_freq_thresholds <- c(.005,.01,.02,.05,.1, .2,.5)
names(alt_freq_thresholds) = c(.005,.01,.02,.05,.1, .2,.5)
