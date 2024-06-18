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
│   ├── config.R # define variables and paths. **Needs to be updated for specific users/clusters.**
│   ├── process_vcf_functions.R # functions for VCF input parsing

```

# Abstract
Because the Mycobacterium tuberculosis complex (MTBC) evolves slowly, isolates from individuals linked in transmission often have identical or nearly identical genomes, making it difficult to reconstruct transmission chains. Finding additional sources of shared MTBC variation could help overcome this problem. Previous studies have reported MTBC diversity within infected individuals; however, whether within-host variation improves transmission inferences remains unclear.

## Methods.
To evaluate the transmission information present in within-host MTBC variation, we conducted a retrospective genomic epidemiology study in which we re-analyzed publicly available sequence data from household transmission studies for which both genomic and epidemiological contact data were available, using household membership as a proxy for transmission linkage. We quantified minority variants, positions with two or more alleles each supported by at least 5-fold coverage and above a 1% minor allele frequency, outside of PE/PPE genes, within individual samples and shared across samples. We compared the performance of a general linear model for household membership that included shared minority variants and one that included only fixed genetic differences. 

## Findings. 
We identified three MTBC household transmission studies with publicly available whole genome sequencing data and epidemiological linkages: a household transmission study in Vitória, Brazil (Colangeli et al.)1, a retrospective population-based study of pediatric tuberculosis in British Columbia, Canada (Guthrie et al.),2and a retrospective population-based study in Oxfordshire, England (Walker et al.).3 We found moderate levels of minority variation present in MTBC sequence data from cultured isolates that varied significantly across studies: mean 168.6 (95% CI: 151.4-185.9) minority variants for the Colangeli et al. dataset, 5.8 (95% CI: 1.5-10.2) for Guthrie et al., and 7.1 (95% CI: 2.4,11.9) for Walker et al.). Isolates from household pairs shared more minority variants than did randomly selected pairs of isolates: mean 97.7 (95% CI: 79.1-116.3) shared minority variants vs. 9.8 (95% CI: 8.6-11.0) in Colangeli et al.; 0.8 (95% CI: 0.1-1.5) vs. 0.2 (95% CI: 0.1-0.2) in Guthrie et al.; and 0.7 (95% CI: 0.1-1.3) vs. 0.2 (0.2-0.2) in Walker et al.; all p<0.0001, Wilcoxon rank sum test) (Table 2; Fig. 3). Shared within-host variation was significantly associated with household membership (OR: 1.51, 95% CI: 1.30,1.71, p < 0.0001), for one standard deviation increase in shared minority variants. Models that included shared within-host variation improved the accuracy of predicting household membership in all three studies as compared to models without within-host variation: AUC: 0.95 versus 0.92, 0.99 versus 0.95, and 0.93 versus 0.91, for the Colangeli et al., Guthrie et al., and Walker et al. studies, respectively. 

## Interpretation. 
Within-host MTBC variation persists through culture and could enhance the resolution of transmission inferences. The substantial differences in minority variation recovered across studies highlights the need to optimize approaches to recover and incorporate within-host variation into automated phylogenetic and transmission inference. 
