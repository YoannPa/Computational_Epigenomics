# Subtelomeric regions annotation

_"Subtelomeres are segments of DNA between telomeric caps and chromatin. Each chromosome has two subtelomeres immediately adjacent to the long (TTAGGG)n repeats. Subtelomeres are considered to be the most distal (farthest from the centromere) region of unique DNA on a chromosome, and they are unusually dynamic and variable mosaics of multichromosomal blocks of sequence."_ - [**Wikipedia**](https://en.wikipedia.org/wiki/Subtelomere)  

While this definition seems rather simple and consensual, span of regions referred as "subtelomeric" regions on human chromosomes are extremely variable from one paper to another.  
Consequently, there is no consensual genomic annotation of subtelomeric regions provided by oftenly used annotation repositories such as the UCSC Genome Browser, the NCBI, the ENCODE project or ENSEMBL.  

Here, I attempt to provide an annotation of subtelomeric regions in the hg19 human genome assembly in a BED format.  

The annotation has been generated based on 500kb sequences provided in a [supplementary FASTA file](https://genome.cshlp.org/content/suppl/2014/04/16/gr.166983.113.DC1/Supplemental_FileS1.txt) from the paper [*Subtelomeric CTCF and cohesin binding site organization using improved subtelomere assemblies and a novel annotation pipeline* - Nicholas S. et al. Genome Research 2014](https://genome.cshlp.org/content/24/6/1039.full).  

50bp have been extracted on both extremities of each subtelomere sequences.  
These 50 bp sequences have been then BLASTed to the hg19 human genome assembly (GRCh37.p13). Sequences submission has been done through the NCBI BLAST API using the R package [**NCBI.BLAST2DT**](https://github.com/YoannPa/NCBI.BLAST2DT). BLAST hits have been collected using the same R package.  

I selected the hits of interest based on the following conditions:  
* The hit must be located on the GRCh37.p13 primary assembly.  
* The hit must be located on the same chromosome as the one from which the BLASTed 50bp sequences are coming from.  
* The hit must have the maximum score of alignment (score = 100) where it happens.  

I then designed a [function](function_subtelomeres_bed_annotation.R) able to extract positions of these perfect alignments on distal (most distant hits from chromosomes centromere) and proximal (closest hits from chromosomes centromere) parts of both p and q arms of each chromosome.  

As there was not always a perfect hit available I completed the missing coordinates using the **Gaps** annotation provided by the UCSC genome browser: The Gaps annotation contain coordinates of telomeres and centromeres on the hg19 human genome assembly. Missing coordinates where so replaced by telomere ends on p arms, and telomere starts on q arms.  

Another issue to overcome was the coexistence of multiple hits on the same chromosome, all with a perfect alignment score of 100.  
In such case I chosed to select the hit to be the closest to 500kb from either the distal or the proximal coordinate (knowing that we expect final subtelomeric regions to have similar length to sequences provided in the paper).  

The resulting BED annotation after these steps is available in the file [hg19_subtelomeres.bed](hg19_subtelomeres.bed).  

As one can spot in this version of the subtelomeres annotation, some gaps remain between telomeres and the newly defined subtelomeric regions. One reason that could explain these gaps is the fact that the hg19 assembly wasn't the same back when the paper defined the subtelomere they described. To fill gaps between telomeres and subtelomeres, I extended the subtelomeric regions, on the p arm where the telomeres end, and on the q arm where the telomeres start.  

The resulting BED annotation after subtelomeric regions extension is available in the file [hg19_extended_subtelomeres.bed](hg19_extended_subtelomeres.bed).  

If you wish to reproduce this work, the selected BLAST hits used here are available in [subtelomere_selected_BLAST_hits.csv](subtelomere_selected_BLAST_hits.csv) and the complete script to reproduce the subtelomeres annotation on the hg19 assembly is available in [make_subtelomeres_bed_annotation.R](make_subtelomeres_bed_annotation.R).  

