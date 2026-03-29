
- ### Interpretation for each result  


name : name of the exposure     

Egger_snp : Number of SNPs selected for Egger and IVW estimation       

GRAPPLE_snp : Number of SNPs selected for GRAPPLE estimation     

GRAPPLE_modes : Number of modes of the profile likelihood of GRAPPLE   

GRAPPLE_est : Beta estimation of GRAPPLE     

GRAPPLE_pval : Pvalue for the beta estimation of GRAPPLE    

Egger_est : Beta estimation of Radial Egger        

Egger_pval : Pvalue for the beta estimation of Radial Egger      

Egger_int : Intercept estimation of Radial Egger       

Egger_int_pval : Pvalue for the intercept estimation of Radial Egger  

Egger_Qstat : Heterogeneity Qstatistic of Radial Egger     

Egger_Qpval : Pvalue of the heterogeneity Qstat of Radial Egger    

Egger_I2GX : I2GX dilution term for Radial Egger     

IVW_est : Beta estimation of Radial IVW 

IVW_pval : Pvalue for the beta estimation of Radial IVW        

IVW_Qstat : Heterogeneity Qstatistic of Radial IVW

IVW_Qpval : Pvalue of the heterogeneity Qstat of Radial IVW      

GRAPPLE_se : Standard error of GRAPPLE beta estimation     

Egger_se : Standard error of Radial Egger beta estimation        

IVW_se : Standard error of Radial IVW beta estimation 

GRAPPLE_OR : Odd ratio of GRAPPLE     

Egger_OR : Odd ration of Radial Egger        

IVW_OR : Odd ratio of Radial IVW  

IVW_LL_OR : Lower limit odd ratio of Radial IVW       

IVW_UL_OR : Upper limit odd ratio of Radial IVW      

Egger_UL_OR : Upper limit odd ratio of Radial Egger    

Egger_LL_OR : Lower limit odd ratio of Radial Egger    

GRAPPLE_LL_OR : Lower limit odd ratio of GRAPPLE  

GRAPPLE_UL_OR : Upper limit odd ratio of GRAPPLE


############################################################## UNIV SENSITIVITY 

name : Name of the exposure 

Egger_snp : Number of SNPs for Egger and IVW estimation 

out_snp : Number of SNPs for Egger and IVW estimation after removing outliers at Qstat < 0.01

assoc_out_snp : Number of SNPs for Egger and IVW estimation after removing SNPs associated with outcome at p < 5e-6

conf_snp : Number of SNPs for Egger and IVW estimation after removing SNPs associated with the confounding factor at p < 5e-8 (Remnant_C and LDL_C)

***_est : Estimation of the method ***

***_est_out : Estimation of the method *** after removing outliers at Qstat < 0.01

***_est_confound : Estimation of the method *** after removing SNPs associated with the confounding factor at p < 5e-8 (Remnant_C and LDL_C)

***_est_assoc_out : Estimation of the method *** after removing SNPs associated with outcome at p < 5e-6

***_pval : Pvalue of the estimation of the method ***

***_pval_out : Pvalue of the estimation of the method *** after removing outliers at Qstat < 0.01

***_pval_confound : Pvalue of the estimation of the method *** after removing SNPs associated with the confounding factor at p < 5e-8 (Remnant_C and LDL_C)

***_pval_assoc_out : Pvalue of the estimation of the method *** after removing SNPs associated with outcome at p < 5e-6


################################################################ MULTI_RESULTS


name...1 : name of the exposure      

n_snps : Number of SNPs selected for thr stringent P threshold  

Horse_est : Beta estimation of MVMR horse       

Horse_sd : Standard deviation of the beta estimation of MVMR horse        

Horse_pval : Pvalue of the beta estimation of MVMR horse      

Horse_converge : Convergence term for the MVMR horse method   

name...7 : name of the exposure        

GRAPPLE_snp : Number of SNPs selected at the GRAPPLE pvalue threshold     

IVW_snp : Number of SNPs selected for the Radial MVMR  

GRAPPLE_beta : Beta estimation of GRAPPLE   

GRAPPLE_pval : Pvalue of the beta estimation of GRAPPLE   

GRAPPLE_Qstat : INSTRUMENT STRENGTH MODIFIED Qstat for GRAPPLE estimation (a low pvalue indicate good instrument strength)

GRAPPLE_Qpval : Pvalue of the INSTRUMENT STRENGTH MODIFIED QSTAT Qstatistic for GRAPPLE estimation   

G_converge : Convergence of GRAPPLE method. FALSE indicate CONVERGENCE and TRUE indicate NON-CONVERGENCE     

IVW_beta : Beta estimation of Radial MVMR         

IVW_pval : Pvalue of the beta estimation of Radial MVMR  

IVW_Qstat : Heterogenity Qstat adapted for Multivariate, for the Radial MVMR estimation        

IVW_Qpval : Pvalue of the heterogenity Qstat adapted for Multivariate, for the Radial MVMR estimation       

IVW_str_Q1 : INSTRUMENT STRENGTH MODIFIED Qstat for Radial MVMR estimation for the 1st covariate (a low pvalue indicate good instrument strength)    

IVW_str_Q2 : INSTRUMENT STRENGTH MODIFIED Qstat for Radial MVMR estimation for the 2nd covariate (a low pvalue indicate good instrument strength)    

IVW_str_Q3 : INSTRUMENT STRENGTH MODIFIED Qstat for Radial MVMR estimation for the 3rd covariate (a low pvalue indicate good instrument strength)



######################################################################## MULTI SENSITIVITY 


Same nomenclature as UNIV SENSITIVITY




