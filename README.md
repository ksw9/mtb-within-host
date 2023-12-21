# Signatures of transmission in within-host *M. tuberculosis* variation

Scripts accomanying our manuscript examining shared within-host *M. tuberculosis* variation in epidemiologically linked individuals. 

# Directory

```
├── scripts
│   ├── process_hh.sh # Driver script
│   ├── 0-process_vcfs.R # read and filter VCFs
│   ├── 1-summarize_shared_isnvs.R # summarize pairwise shared iSNVs
│   ├── 2-organize_metadata.Rmd # organize metadata accomapnying SRA data
│   ├── 3-plots.Rmd # all manuscript plots
│   ├── iqtree.sh # fit IQ-Tree maximum likelihood trees
│   ├── config.R # define variables and paths
│   ├── process_vcf_functions.R # functions for VCF input parsing

```

# Abstract
Because M. tuberculosis evolves slowly, transmission clusters often contain multiple individuals with identical consensus genomes, making it difficult to reconstruct transmission chains. Finding additional sources of shared M. tuberculosis variation could help overcome this problem. Previous studies have reported M. tuberculosis diversity within infected individuals, however, whether within-host variation improves transmission inferences remains unclear.

## Methods.
To evaluate the transmission information present in within-host M. tuberculosis variation, we re-analyzed publicly available sequence data from three household transmission studies, using household membership as a proxy for transmission linkage between donor-recipient pairs. 

## Findings. 
We found moderate levels of minority variation present in M. tuberculosis sequence data from cultured isolates that varied significantly across studies (mean: 7-170 minority variants >1% minor allele frequency, outside of PE/PPE genes). Isolates from household members shared more minority variants than did isolates from unlinked individuals in the three studies (median 95 shared minority variants vs. 1; 3 vs. 1, and 0 vs. 0, respectively). Shared within-host variation was significantly associated with household membership in each study (OR: 3.00 [2.40-3.89]; OR: 1.69 [1.21-2.37], and OR: 1.37 [1.04-1.60], for one standard deviation increase in shared minority variants), and shared minority variants predicted household membership with an AUC from 0.68-0.97. 

## Interpretation. 
Within-host M. tuberculosis variation could enhance the resolution of transmission inferences. Pathogen enrichment approaches improve the recovery of within-host variation, highlighting the need to optimize approaches to recover and incorporate within-host variation into automated phylogenetic and transmission inference. 
