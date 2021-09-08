# Pericentromeric regions annotation
## Introduction
Pericentromeric regions are often described as segments of DNA situated close to a chromosome centromere or on each chromosome arm alongside the centromere.  
While this definition seems rather simple and consensual, span of regions referred as "pericentromeric" regions on human chromosomes are extremely variable from one paper to another.  
Consequently, there is no consensual genomic annotation of pericentromeric regions provided by oftenly used annotation repositories such as the UCSC Genome Browser, the NCBI, the ENCODE project or ENSEMBL.  

Here, I attempt to provide an annotation of pericentromeric regions in the hg19 human genome assembly in a BED format based on an analysis of chromosomes content in specific repetitive elements.  
I focused mainly on 2 classes of repeats:  
* Satellite sequences: Alpha satellites, Beta satellites and HSATII.  
* Transposable elements: LINEs, SINEs and LTRs.  

Some satellite sequences, such as alpha and beta satellites, have already largely been described in the litterature to be more densely present around centromeres, rather than throughout the rest of the human genome where they are less frequent **<sup>(1,2,3,4)</sup>**.  
Their location, close to the centromere seems to make them more prone to follow specific regulations and having specific roles related to chromosomes organisation during mitosis and chromatin stability around the centromeres.
These satellites are also known to have an expansion and an expression activity from the centromeric areas which can be disrupted in some cancers like in the case of HSATII **<sup>(2,5)</sup>**.

Satellites are not the only class of repeat with the ability to "expand" within the human genome and to be associated with pericentromeric regions.
Through transposition mecanisms, transposable elements such as LINEs, SINEs and some ERVs can duplicate and insert within genomes. One specific type of sequence described as Pericentromeric Interspersed Repeats (PIRs) have been found to remain more conserved around centromeres in a subset of human chromosomes **<sup>(7)</sup>**. From an evolution point of view, they are associated with a high duplication activity event which has been estimated to have occured around 5 to 7 millions of years ago before the separation between great apes and human species **<sup>(6,7)</sup>**.  

Other types of repeats have also been described to be associated with centromeric and pericentromeric regions **<sup>(8)</sup>**.  

In the production of this pericentromeric annotation I only considered the 2 classes mentionned before: satellites and transposable elements.  

## Selecting transposable elements of interest.
One sequence provided as a starting point by [J. E. Horvath et al.](https://academic.oup.com/mbe/article/20/9/1463/976868) to study association of interspersed elements and pericentromeric is [AC073318](https://www.ncbi.nlm.nih.gov/nuccore/AC073318.8/): starting from 71,401 bp and ending at 120,576 bp.  
I extracted the sequence, splitted it into 7kb smaller sequences and submitted these for a BLAST against hg19 genome assembly on NCBI servers using the R package [NCBI.BLAST2DT](https://github.com/YoannPa/NCBI.BLAST2DT).  
I collected the resulting BLAST hits using the same tool, and considered only hits with an alignment length strictly superior to 700, an E-value strictly inferior to 0.01 and an overall score strictly superior to 600 (Figures [1](figures/Fig1.pdf), [2](figures/Fig2.pdf) and [3](figures/Fig3.pdf)).  
Then I extracted alignements from the different BLASTED smaller sequences in order to keep only those which were following each other (the same way the smaller sequence 2, should follow the smaller sequence 1, and is followed by smaller sequence 3, 4, 5, ...) in both directions (forward and reverse alignments). None of the alignments were perfect alignment, and I did not obtain more than 2 smaller sequences following each other at best.  
After looking at query sequences parts aligned, I realized that 4 sub-sequences within the submissions had a higher probabily of alignment on hg19. Thus, I extracted these 4 sub-sequences, and submitted them to NCBI servers for another BLAST (You can find a FASTA of these 4 sub-sequences in the file [PIR4_subsequences_A_B_C_D.fasta](data/PIR4_subsequences_A_B_C_D.fasta)).  
These sub-sequences do not seem to have much in common in terms of similarities: You can check the Clustal Omega multiple alignment result available in the file [MView_ClustalO_multiple_alignment_PIR4_A_B_C_D.html](https://htmlpreview.github.io/?https://github.com/YoannPa/Computational_Epigenomics/blob/master/Pericentromeric_regions_annotation/data/MView_ClustalO_multiple_alignment_PIR4_A_B_C_D.html) and the resulting phylogenetic tree based on the multiple alignment in [Figure 4](figures/Fig4.png).
In order to select the best hits from each sub-sequences BLASTed, I decided to select a cut-off on a minimum overall score of the alignments returned. The cut-off is specific to each sub-sequence and for each chromosome: it is the minimum score necessary so that all the alignments above this score fall within both arms proximal halves of a chromosome for a given sub-sequence.  
This cut-off strategy is based on 2 facts:  
1. Whatever the span that pericentromeric regions might have, by definition, they should not expand beyond the proximal halves of both arms of a chromosomes. Otherwise these regions could also be considered as subtelomeric with the same reasoning.  
2. From the analysis of alignments distribution, and previous knowledge reported about the origin of these sub-sequence, closer alignments are to centromeres, higher the alignment score is supposed to be.  
Indeed, using this strategy, I have been able to select alignments with the best scores for each sub-sequences and on most of the chromosomes. We can clearly see that closer alignments are to centromeres, higher their scores are ([Figure 5](figures/Fig5.png)).   
Furthermore, sub-sequences alignments selected were not ubiquitous of all chromosomes ([Figure 6](figures/Fig6.png)). For example:  
* Chromosome 9 contains exclusively sub-sequences PIR4A and PIR4B.  
* Chromosome 20 contains exclusively sub-sequences PIR4D and PIR4B.  
* Chromosome 13 contains exclusively sub-sequences PIR4B.  
* Chromosome X doesn't contain any sub-sequences PIR4B.  
Unexpectedly, even before setting a cut-off no alignments were returned from BLAST of these sub-sequences on chromosome 4.  
Regarding chromosomes 11, 14 and 19, despite the fact that some hits were available, the cut-off strategy applied did not keep any of these. Possible reasons are that alignment scores of these hits were simply too low and/or that the alignment scores of hits closer to these chromosomes' centromeres were not higher than the rest of the hits.  
Alignement frames ratio on the genome is also pretty balanced: almost 50% in forward, and 50% in reverse ([Figure 7](figures/Fig7.png)).  
You can find an annotation of these selected alignments in the file [hg19_PIR4.bed](data/hg19_PIR4.bed).  
Since these sub-sequences cannot be aligned quite well to each other, and that they are supposed to be overlapping with interspersed elements I also checked the repeat composition of the selected alignments.
By looking at overlapping repeats in these regions, I found that each sub-sequences was composed of a specific chaining of repetitive elements.  
All sub-sequences have in common to overlap on their entire length with a transposable element from the LINES L1 family.  
The vast majority of other identified repeats overlapping are transposables elements from the SINEs or ERVs families.  
For each sub-sequence, I have been able to associate a specific chaining pattern of repetitive elements overlapping their alignments, independent from the frame of alignments considered.
Basically, long pieces of L1 overlapping with the selected alignements are interrupted at some specific locations by SINEs or ERVs. The fact that we are able to find this pattern to appear in multiple alignments, and in different directions, makes me believe that, not only the L1 elements overlapping might have transposed throughtout primate genomes millions of years ago, but they might also have carried prior insertions of other transposable elements (namely here SINEs and ERVs) which occured before these transposition events took place.  
With the alignment score decreasing, this chaining patterns of repetitive elements tend to disappear to only keep an overlapping with an L1 element or an Alu sequence. Again, the highest alignement scores and most conserved repeats patterns can be observed closer to centromeres. The probable explanation for the disappearance of these patterns in alignments more distant from centromere might be the increased occurence of mutations. As these sequenced are in non-coding regions, and do not need to be as conserved as centromeric regions (necessary for kinetochore scaffolds attachment during Mitosis) they might have accumulated more mutation than alignments found closer to centromeres. This accumulation of mutations can make RepeatMasker fail to identify highly mutated repeat sequences, and by consequence make it fail to annotate them.  
You can find the chaining patterns associated to each sub-sequences for all selected alignements in the file [PIR4_repeats_chains.tsv](data/PIR4_repeats_chains.tsv): Alignments are ordered by decreasing alignment scores (top scores in top rows).  
Based on their location and composition I have called these selected alignements **PICS: Pericentromeric Interspersed Composite Sequences**.  

## Selecting satellite sequences of interest.
After defining PICS I also selected some satellite sequences in order to define those who could be pericentromeric. Based on the litterature introduced before, I selected alpha, beta and HSATII satellites. As these repeats are already annotated on hg19, I retrieved there coordinates from the UCSC genome browser databases. After assessing their distribution throughout hg19, and seeing they were more densely annotated when closer to centromeres (see [Figure 8](figures/Fig8.png)), I decided to set a cut-off on their density in the genome. For each chromosome I considered as pericentromeric, satellites from these 3 families:  
1. when there was more than 3 satellites annotated within a window of N Kb.  
2. when all the N Kb windows matching previous conditions were located within the proximal halves of both chromosome arms.  
The window size N considered depended on the chromosome: Some chromosome containing less of these satellites needed larger windows to be considered. However, regardless of the size of windows considered for each chromosome, all windows had to fit within the proximal haves of both chromosome arms anyway, for the sattelites contained to be considered pericentromeric.  
Using these density cut-off, I have been able to retain most of the sattelite considered initially (see [Figure 9](figures/Fig9.png)): unexpectedly contrary to what was described in papers I selected prior to the analysis, a larger portion of the beta satellite was kept with this approach. The portion of alpha satellite, the family of satellites often presented as the pericentromeric, was smaller (in terms of percentage of the total alpha satellites annotated before applying the cut-off). Even more surprising, almost all HSATII satellite annotations were kept after applying the cut-off. Therefore, here HSATII satellites are considered more pericentromeric than beta satellites, these latter being more pericentromeric than alpha satellites.  
Using this selection strategy worked on most chromosomes in hg19 but notably failed for chromosomes 4 and 5 where no satellites were retained, because their density on these chromosomes was too low. A peculiarity of chromosome 4 is this high satellite density region located on the distal part of its p arm, outside the proximal halves of both arms. It consist of some alpha satellites close to each others, making the selection based on satellites density fail on this chromosome.

## Defining pericentromeric regions from PICS and selected satellites
Based on all the selected PICS and satellites (see [Figure 10](figures/Fig10.png)) I defined pericentromeric regions as following for each chromosome:  
* On p arm: the pericentromeric region starts with the start position of the first element selected, and ends with the start position of the centromere.  
* On q arm: the pericentromeric region starts with the end position of the centromere, and ends with the end position of the last element selected.  
The pericentromeric regions defined this way are sumarized on [Figure 11](figures/Fig11.png).  
Pericentromeric regions have been defined on almost all chromosome, except one: chromosome 4.  
On chromosome 4:  
* No PICS was aligned, and consequently none alignment has been selected.
* regarding satellites, only few alpha satellite sequences are annotated and their distribution doesn't make it possible to select any, based on their density, since their density is not higher around the centromere.
In terms of repeats chromosome 4 must have a specific composition to fail both attempts to identify pericentromeric sequences. Besides, the definition of pericentromeric regions failed on the q arm of chromosome 19, and on the p arm of acrocentric chromosomes 13, 14, 15 and 22. Chromosome 21 is the only acrocentric chromosome for which a pericentromeric has been defined on its p arm.  

Albaits this annotation of pericentromeric regions is probably not perfect it provides a new approach toward defining pericentromeric regions based on chromosome composition in satellite sequences and specific transposable elements. The fact that during the analysis each chromosome is processed separately also emphatize the fact that, whatever expansion mecanisms are, or have been, in play in these areas, the activity of these expansions might not be, or have been, the same for all chromosomes.
We can also note similarities with the schematic map provided in [J.M. Mudge and M.S. Jackson 2005](https://www.karger.com/Article/Pdf/80801) especially about acrocentric chromosomes p arms, and the uncommon composition of pericentromeric regions of chromosome 4.  

The outlook for this annotation is to be completed, maybe with more repeats families, and probably with other types of data (epigenomic and transcriptomic data for example). With respect to this, **outside contributions and feedbacks are all very welcome.**  

## Reproducing the annotation
The work shared here can be reproduced using scripts available in the **/src** directory.  
All BED annotations used here will be available in the hg19_annotation folder of this repository.  
In case you encounter some issues executing the code, please fill in an issue in this repository describing the error message you have or the issue you are facing.  
Also if some files are missing, it is possible that I forgot to upload some, so feel free to ask.  

## References
1. [_J.M. Mudge and M.S. Jackson. Evolutionary implications of pericentromeric gene expression in humans._](https://www.karger.com/Article/Pdf/80801)
2. [_Ksenia Smurova and Peter De Wulf. Centromere and Pericentromere Transcription: Roles and Regulation in Sickness and in Health._](https://www.frontiersin.org/articles/10.3389/fgene.2018.00674/full)
3. [_M. A. Dobrynin et al. Human pericentromeric tandemly repeated DNA is transcribed at the end of oocyte maturation and is associated with membraneless mitochondria-associated structures._](https://www.nature.com/articles/s41598-020-76628-8)
4. [_V.A. Shepelev et al. Annotation of suprachromosomal families reveals uncommon types of alpha satellite organization in pericentromeric regions of hg38 human genome assembly._](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4496801/pdf/main.pdf)
5. [_Francesca Bersani et al. Pericentromeric satellite repeat expansions through RNA-derived DNA._](https://www.pnas.org/content/112/49/15148)
6. [_M. Pita et al. A Highly Conserved Pericentromeric Domain in Human and Gorilla Chromosomes._](https://www.karger.com/Article/PDF/251962)
7. [_J. E. Horvath et al. Using a Pericentromeric Interspersed Repeat to Recapitulate the Phylogeny and Expansion of Human Centromeric Segmental Duplications._](https://academic.oup.com/mbe/article/20/9/1463/976868)
8. [_Juliann E. Horvath et al. The Mosaic Structure of Human Pericentromeric DNA: A Strategy for Characterizing Complex Regions of the Human Genome._](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC310890/)
9. [_Lakshay Anand. chromoMap-An R package for Interactive Genomic Visualization of Multi-Omics Data_](https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html)
