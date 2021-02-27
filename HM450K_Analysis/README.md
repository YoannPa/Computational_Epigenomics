# HM450K Analysis

The folder **fun/** contains functions dedicated to HM450K data analysis.  
The folder **data/** contains data that could be useful for HM450K analysis.


## Functions  
* `probeseq_match.R` - Counts HM450K probe sequence matches on hg19 genome assembly scaffolds.  

## Data
`HM450K_probes_annotation.csv` - A CSV file listing HM450K probes purposes.  
* The first column "ID" contains the ID of probes.  
* The second column "Target" contains probes' target types. The target can be either:  
   * A type of methylation: "CpG", "CpH" (any other C-dinucleotide different from a CpG) or "SNP" (when captured by genotyping probes).  
   * A preparation step involved in the Infinium Human Methylation 450K method.  
* The third column "Purpose" defines the purpose of probes: The intensity of a probe can either be used for methylation analysis (Analysis), or for quality control (QC).  
  
  
`HM450K_probes_filters.csv` - A CSV file listing HM450K probes that are either cross-reactive, or polymorphic, or both cross-reactive and polymorphic.
* The first column provides the ID of CG probes.  
* The second column specifies whether the probe is part of the cross-reactive probes using BLAT results (TRUE) as defined in the paper [*Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray* - Chen Y. et al. Epigenetics 2013](https://pubmed.ncbi.nlm.nih.gov/23314698/) or not (FALSE).  
* The third column specifies whether the probe is part of the polymorphic probes using BLAT results (TRUE) as defined in the paper [*Discovery of cross-reactive probes and polymorphic CpGs in the Illumina Infinium HumanMethylation450 microarray* - Chen Y. et al. Epigenetics 2013](https://pubmed.ncbi.nlm.nih.gov/23314698/) or not (FALSE).  
* the fourth column specifies whether the probe is part of the cross-reactive probes using BOWTIE2 results (TRUE) as defined in the Github repository [*illumina450k_filtering* - Miles Benton. Github 2015](https://github.com/sirselim/illumina450k_filtering) or not (FALSE).

`HM450K_genotyping_probes_hg19.bed` - A BED file containing the coordinates of the HM450K genotyping 'rs' probes on the hg19 genome assembly version.
* The first column provides the ID of rs probes, which also match Human dbSNP ID.
* The second column states on which chromosome the probe is located.
* The third column matches the 1-based position of the SNP.  

A more detailed introduction to the 'rs' genotyping probes used in HM450K can be found in this [*paper*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5984806/):  
_"There are 65 probes placed on the 450K chip querying high-frequency SNPs (with 59 of these on the EPIC chip; their probe identifiers start with “rs”). Just as for CpG sites, a β-value is calculated for each SNP locus, based on fluorescence intensities from two probes targeting either the wild type or the common mutant variant. These β-values usually fall into one of three disjunct clusters, corresponding to the heterozygous and the two homozygous genotypes (AB, AA, or BB). The specific combination of SNPs across these 65 probes serves as a genetic fingerprint: fingerprints of samples from the same donor match but differ between individuals – with the exception of monozygotic twins – thereby enabling one to check for discrepancies with the metadata."_  
(_N.B.: The SNP probe ID **rs13369115** has been merged with another SNP (rs60784560) under the new ID **rs10155413**. More detailed [here](https://www.ncbi.nlm.nih.gov/snp/rs10155413#history)._)  
