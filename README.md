This code implements a Snakemake pipeline for performing Mendelian randomization analyses on Genome-Wide Association Studies (GWAS).
The aim is to assess whether there is a causal relationship between a genotype known to be associated with a particular disease and a specific phenotype within the GWAS population.

The GWAS used in this study was conducted in a population of patients with calcific aortic stenosis. These patients underwent blood testing, during which their lipoprotein concentrations were measured. Mendelian randomization analyses were then performed to assess the causal effect of lipoproteins on calcific aortic stenosis.

INSTALLATION :
To run it, you need :

Snakemake: [Install Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

Plink: [Install Plink](https://www.cog-genomics.org/plink/2.0/)

R:  [Install R](https://cran.r-project.org/)
And the following R packages : 
- GRAPPLEofNEW To install this package, you need to build it :
```bash
R CMD build your/path/to/folder_of/GRAPPLEofNEW```
Then ```install.packages("your/path/to/GRAPPLEofNEW_0.2.2.tar.gz")``` 
- ieugwasr
- dplyr
- TwoSampleMR
- RadialMR
- data.table
- ggplot2
- remotes
- MRCIEU
- MRBEE
- RMVMR
- MVMR
