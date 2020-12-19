# HM450K Analysis

The folder **fun/** contains functions dedicated to HM450K data analysis.  
The folder **data/** contains data that could be useful for HM450K analysis.


## Functions  
* `probeseq_match.R` - Counts HM450K probe sequence matches on hg19 genome assembly scaffolds.  

## Data
`HM450K_probes_annotation.csv` - A CSV file listing HM450K probes purposes.
* The first column "ID" contains the ID of probes.
* The second column "Target" contains probes' target types. The target can be either:
   * A type of methylation: CpG, a CpH (any other C-dinucleotide different from a CpG) or a SNP captured by genotyping probes.
   * A preparation step involved in the Human Methylation 450K method.
* The third column "Purpose" defines the purpose of probes: The intensity of a probe can either be used for methylation analysis (Analysis), or for quality control (QC).
`HM450K_probes_filters.csv` - A CSV file listing HM450K probes that are either cross-reactive, or polymorphic, or both cross-reactive and polymorphic.
* The first column provide the ID of CG probes.  
* The second column specify wether the probe is part of the cross-reactive probes using BLAT results (TRUE) as defined in the paper [*Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray* - Chen Y. et al. Epigenetics 2013](https://pubmed.ncbi.nlm.nih.gov/23314698/) or not (FALSE).  
* The third column specify wether the probe is part of the polymorphic probes using BLAT results (TRUE) as defined in the paper [*Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray* - Chen Y. et al. Epigenetics 2013](https://pubmed.ncbi.nlm.nih.gov/23314698/) or not (FALSE).  
* the fourth column specify wether the probe is part of the cross-reactive probes using BOWTIE2 results (TRUE) as defined in the Github repository [*illumina450k_filtering* - Miles Benton. Github 2015](https://github.com/sirselim/illumina450k_filtering) or not (FALSE).  
