# HHV8-non-B-DNA
R scripts used in publication: Genomic changes in Kaposi Sarcoma-associated Herpesvirus and their clinical correlates

To replicate, one needs SIST (Stress-Induced Structural Transitions) (Zhabinskaya, et al, 2015) for calculation of cruciform, denaturation and Z-DNA formation probabilities, which is run in Bash. R packages triplex (Hon, et al, 2013) and pqsfinder (Hon, et al, 2017) is used for calculating triplex and G-quadruplex scores, respectively. Refer to their respective manuals.

SIST - https://bitbucket.org/benhamlab/sist_codes/src/master/

pqsfinder - https://bioconductor.org/packages/release/bioc/html/pqsfinder.html

triplex - https://bioconductor.org/packages/release/bioc/html/triplex.html

common_functions.R
initializes many functions used in the study

plot_nB-DNA.R
For plotting predicted non-B DNA probabilities or scores along KSHV sequences and combining plots together. Codes to make Figure 4A - D are in here.

controls_&sims.R
Used in making Table 3, and for 1000 simulations for testing significance of association of observed breakpoints with G4

odds.R
Calculate odds that a random 17.5 kb slice of the KSHV genome GK18 will include the 'minimum region' frequently over-covered in whole viral genome sequencing of many KS tumors.

motifs_test.R
Counting Motifs and their association with breakpoints. Tested to see if motifs identified in (Brown, 2014; see reference list in manuscript) are present and associated with the observed breakpoints. This code was not used in the publication because sample size (31) is not big enough to give a robust answer.
