This code implements a Snakemake pipeline for performing Mendelian randomization analyses on Genome-Wide Association Studies (GWAS).
The aim is to assess whether there is a causal relationship between a genotype known to be associated with a particular disease and a specific phenotype within the GWAS population.

The GWAS used in this study was conducted in a population of patients with calcific aortic stenosis. These patients underwent blood testing, during which their lipoprotein concentrations were measured. Mendelian randomization analyses were then performed to assess the causal effect of lipoproteins on calcific aortic stenosis.

INSTALLATION :
To run it, you need :

Snakemake: [Install Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

Plink: [Install Plink](https://www.cog-genomics.org/plink/2.0/)

R:  [Install R](https://cran.r-project.org/)
And the following R packages : 
- GRAPPLEofNEW
  ###### To install this package, you need to build it
```bash
R CMD build your/path/to/folder_of/GRAPPLEofNEW
```
  ###### Then install it on R like the others
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

To set your parameters, you have two config files : 

##### A snakemale.yaml : 

Specify each of your folders paths (folder data for selection, exposure, outcome, covariates if you want to do multivariate mendelian randomisation, and the folder where you want the results to be)


##### A R.yaml : 

Here you can set up all the parameters of your analysis : Your p value treshold choosing SNPs, the reference panel for clumping, parameters for multivariate / univariate analysis, kind of metrics you want for results...



then just run the code with : 

`snakemake --snakefile snakefile --use-conda -np ( a dry run to just create the result folders)`

`snakemake --snakefile snakefile --use-conda -c 1 -j 1` 


###### MVMR Horse is also implemented for multivariate analysis 
###### Before using it, you need to install the application JAGS, and the library(R2jags) in R
###### The method take a lot of time to run, so you have to accept to run it in the R config file
###### Results are stored in another file than the rest of the methods, and no sensitivity analysis is made

