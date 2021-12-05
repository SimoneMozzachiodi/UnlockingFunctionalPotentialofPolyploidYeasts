# UnlockingFunctionalPotentialofPolyploidYeasts
Repository contains code used in the work "Unlocking the functional potential of polyploid yeasts"

This folder contains the scripts used for loss of heterozygosity regions detections in RTG derived from OS1364 strain.
The scripts need to be launched in order following the number prefix in the script (e.g. 1.xxx, then 1.1xxx, then 2...).

Script prefix / explanation

1./1.1/ -> mapping short reads sequencing of RTG samples on S288C reference, duplicate removal and coverage calculation

2./ -> extraction of coverage data from the files produced above and creation of .txt file with coverage by position

3.1.1./ -> Plot of coverage profile across all 16 chromosomes based on S288C coordinates

3.1.2./3.1. -> indexing of vcf of OS1364 obtained by mapping short reads from https://doi.org/10.1038/s41586-018-0030-5 mapped against S288C reference
              p3 is for ploidy3 chromosome, p4 is for ploidy 4 chromosomes (chromosome III)
              
4.2./4.4. -> Variant calling in the RTG experiments to identify the same set of markers of the parental strain and eventual allele frequency (AF) variation

5./5.1 -> Filter of variants called based on quality and type of variant

6./ -> Intersect the set of variants called in each RTG dataset with the variants genotyped for the parental strain to use the same
       set of heterozygous markers for AF shift detection
       
7./ -> Evaluation of stretches of markers with allele frequency variation which implies recombination upon RTG.
       This file create intermediate file used as control (e.g. table of "LOH" based on the length or table
       of markers with tabulated the respective allele frequency shift or genotype)


