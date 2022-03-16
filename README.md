# UnlockingFunctionalPotentialofPolyploidYeasts
Repository contains code used in the work "Unlocking the functional potential of polyploid yeasts"

**Abstract:** 
Breeding and domestication have generated widely exploited crops, animals and microbes. However, many Saccharomyces cerevisiae industrial strains have complex polyploid genomes and are sterile, preventing genetic improvement strategies based on breeding. Here, we present a novel strain improvement approach based on the budding yeasts’ property to promote genetic recombination when meiosis is interrupted and cells return-to-mitotic-growth (RTG). We demonstrated that two unrelated sterile industrial strains with complex triploid and tetraploid genomes were RTG-competent and developed a visual screening for easy and high-throughput identification of recombined RTG clones based on colony phenotypes. Sequencing of the evolved clones revealed unprecedented levels of RTG-induced genome-wide recombination. We generated and extensively phenotyped a RTG library and identified clones with superior biotechnological traits. Thus, we propose the RTG-framework as a fully non-GMO workflow to rapidly improve industrial yeasts that can be easily brought to the market.


The folders contain the scripts used for analysis of RTG derived from sterile polyploid strains. The scripts usually contained
either the id of some of the RTG analysed as example, or a commented entry in which it is specified where the id of the samples
need to be indicated.
The scripts need to be launched in order following the number prefix in the script (e.g. 1.xxx, then 1.1xxx, then 2...).
Scripts without numbering use either the output of others scripts of are standalone and are made to be launched at the end
of the numbered scripts.
Feel free to contact me for any questions regarding the uploaded material.

OS1364 folder
Script prefix / explanation

1./1.1/ -> mapping short reads sequencing of RTG samples on S288C reference, duplicate removal and coverage calculation

2./ -> extraction of coverage data from the files produced above and creation of .txt file with coverage by position

2.1./ -> Plot of coverage profile across all 16 chromosomes based on S288C coordinates

3.1.2./3.1. -> indexing of vcf of OS1364 obtained by mapping short reads from https://doi.org/10.1038/s41586-018-0030-5 mapped against S288C reference
              p3 is for ploidy3 chromosome, p4 is for ploidy 4 chromosomes (chromosome III)
              
4.2./4.4. -> Variant calling in the RTG experiments to identify the same set of markers of the parental strain and eventual allele frequency (AF) variation

5./5.1 -> Filter of variants called based on quality and type of variant

6./ -> Intersect the set of variants called in each RTG dataset with the variants genotyped for the parental strain to use the same
       set of heterozygous markers for AF shift detection
       
7./ -> Evaluation of stretches of markers with allele frequency variation which implies recombination upon RTG
       before any filtering step (CNV filter, shared LOH filter).
       This file create intermediate file used as control (e.g. table of "LOH, AF shift regions", Table of markers
       with detected AF shift) and used for plotting regions with AF shift.


