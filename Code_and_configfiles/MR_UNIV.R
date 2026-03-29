

##### DATA RECUP

###Loading the required libraries 

library(data.table)
library(ggplot2)
library(GRAPPLEofNEW)
library(ieugwasr)



### Recup variables from the config file 

working_dir <- getwd()

source(paste0(working_dir, "/config_R.yaml"))

## Recup arguments from the snakefile 

args = commandArgs(trailingOnly=TRUE)





### Assign the file pathways for selection, exposure and outcome 

sel_file = args[1]
exp_file = args[2]
out_file = args[3]
res_path = args[5]


name = args[4]
n_conf = name_confounding


sel_cov_file_M = args[6:(5+n_cov)]
exp_cov_file_M = args[(6+n_cov):(5+2*n_cov)]
sel_file_M = c(args[1], sel_cov_file_M)
exp_file_M = c(args[2], exp_cov_file_M)
out_file_M = out_file

print(sel_cov_file_M)

print(exp_cov_file_M)

print(sel_file_M)

print(exp_file_M)

print(out_file_M)


###### SNPs FILTERING 

### Filter SNPs according to GRAPPLE p-threshold and clump

Input <- getInput(sel.files = sel_file,
  exp.files = exp_file,
  out.files = out_file,
  plink_refdat = plink_ref,
  max.p.thres = GRAPPLE_P_thresh,
  cal.cor = T,
  p.thres.cor = 0.5,
  get.marker.candidates = F,
  marker.p.thres = GRAPPLE_P_thresh,
  clump_r2 = 0.001,
  plink_exe = plink_bin
)

library(dplyr)
library(TwoSampleMR)
library(RadialMR)

H_data <- Input[[1]]


## Save the selected and pruned SNPs, if something went wrong, or for further analysis

write.table(H_data, paste0(res_path,name,"/selected_and_pruned_SNPs_", GRAPPLE_P_thresh), quote = F, row.names = F, sep = "\t") 


##### UNIVARIATE ANALYSIS : PRIMARY ANALYSIS 

## Find the number of modes of the GRAPPLE function if we are in the univariate setting

GR_modes <- findModes(H_data, p.thres = GRAPPLE_P_thresh, map.marker = F)


## Estimation results of GRAPPLE at the GRAPPLE p-threshold

GR_est <- grappleRobustEst(
  H_data,
  p.thres = GRAPPLE_P_thresh,
  cor.mat = Input[[3]],
  niter = 20,
  diagnosis = T,
  plot.it = T
)  


## Estimation of Radial_IVW and Radial_Egger at the stringent p-threshold


#Filter SNPs at stringent p-threshold
R_data <- H_data[H_data$selection_pval < MR_P_thresh,]


# Format for radial estimation 
rad_dat <- format_radial(BXG = R_data$gamma_exp1, BYG = R_data$gamma_out1, seBXG = R_data$se_exp1, seBYG = R_data$se_out1, RSID = R_data$SNP)


#Perform Radial_IVW analysis
ivw_res <- ivw_radial(rad_dat, alpha = 0.01, weights = 3, tol = 0.0001)


#Perform Radial_Egger analysis
egg_res <- egger_radial(rad_dat,alpha = 0.01, weights = 3)


## I2GX dilution parameter for Egger (~1 showing that the Egger estimate is also a good prediction compared to IVW)

IGX_egg <- Isq(rad_dat$beta.exposure/rad_dat$se.outcome, rad_dat$se.exposure/rad_dat$se.outcome)


###### UNIVARIATE ANALYSIS : PLOTS 


#Comboplot of Radial_IVW, Radial_Egger and GRAPPLE 
Comboplot <- plot_radial(c(ivw_res, egg_res), T, F, F)

Comboplot <- Comboplot + geom_abline(intercept = 0, slope = GR_est$beta.hat, col = "green")


#QQplot of GRAPPLE
GR_plot <- GR_est$p

###### UNIVARIATE ANALYSIS : SENSITIVITY ANALYSIS - OUTLIERS  


## Sensitivity analysis without outliers for GRAPPLE
if(nrow(GR_est$outliers) != 0){
	
	#Remove outlier SNPs
        H_data_out <- H_data[!(H_data$SNP %in% rownames(GR_est$outliers)),]

	#Perform GRAPPLE analysis without outliers 
        GR_est_out <- grappleRobustEst(
        H_data_out,
        p.thres = GRAPPLE_P_thresh,
        cor.mat = Input[[3]],
        niter = 40,
        diagnosis = T,
        plot.it = T
        )
	
	#QQplot of GRAPPLE without outliers
	GR_out_plot <- GR_est_out$p
	
	#Save the outlier SNPs for further analysis
        write.table(rownames(GR_est$outliers), paste0(res_path,name,"/UNIV_GRAPPLE_outliers"), quote = F, row.names = F, col.names = F, sep = "\t")
}


## Markers SNPs of GRAPPLE for each mode, interesting if there is multiple modes
if(nrow(GR_modes$markers) != 0){
        write.table(rownames(GR_modes$markers), paste0(res_path,name,"/UNIV_GRAPPLE_markers"), quote = F, row.names = F, col.names = F, sep = "\t")
}



## Sensitivity analysis without Egger outliers for IVW and Egger
if(class(egg_res$outliers) != "character"){
	
	#Remove outliers SNPs
        rad_dat_out <- rad_dat[!(rad_dat$SNP %in% egg_res$outliers[,1]),]

	#Perform Radial_IVW analysis without outliers 
        ivw_res_o <- ivw_radial(rad_dat_out, alpha = 0.05, weights = 3, tol = 0.0001)

	#Perform Radial_Egger analysis without outliers 
        egg_res_o <- egger_radial(rad_dat_out,alpha = 0.05, weights = 3)

	#Save the outliers SNPs for further analysis 
        write.table(egg_res$outliers$SNP, paste0(res_path,name,"/UNIV_egger_outliers"), quote = F, row.names = F, col.names = F, sep = "\t")
}




###### UNIVARIATE ANALYSIS : SENSITIVITY ANALYSIS - OVER ASSOCIATION OF FILTERED_SNPs/OUTCOME (removing SNPs associated with the outcome)

# Load the outcome file 
high_out <- data.frame(fread(out_file))


# Remove SNPs associated with the outcome
high_out <- high_out[high_out$V8 < p_thresh_assoc,]

H_high <- H_data[!(H_data$SNP %in% high_out$V3),]


## Perform GRAPPLE analysis without SNPs associated with the outcome 
GR_high <- grappleRobustEst(
  H_high,
  p.thres = GRAPPLE_P_thresh,
  cor.mat = Input[[3]],
  niter = 20,
  diagnosis = T,
  plot.it = T
)


# Filter SNPs with a stringent p-threshold 
H_high <- H_high[H_high$selection_pval < 5e-8,]


# Format for radial estimation
rad_high <- format_radial(BXG = H_high$gamma_exp1, BYG = H_high$gamma_out1, seBXG = H_high$se_exp1, seBYG = H_high$se_out1, RSID = H_high$SNP)

# Perform Radial_IVW without SNPs associated with the outcome
ivw_high <- ivw_radial(rad_high, alpha = 0.01, weights = 3, tol = 0.0001)


# Perform Radial_Egger without SNPs associated with the outcome 
egg_high <- egger_radial(rad_high,alpha = 0.01, weights = 3)




###### UNIVARIATE ANALYSIS : SENSITIVITY ANALYSIS - REMOVING SNPs OF ANOTHER EXPOSURE (to ensure that if there is an effect, it doesn't come from this exposure)


if(confounding == T){
	# Load the confounding exposure file 
        confound <- data.frame(fread(confound_1))
	
	confound <- confound[confound$pval < 5e-8,]
	
	# Remove SNPs associated with the confounding exposure
        H_conf <- H_data[!(H_data$SNP %in% confound$SNP),]
        rad_conf <- rad_dat[!(rad_dat$SNP %in% confound$SNP),]

	# Perform GRAPPLE estimation
        GR_est_conf <- grappleRobustEst(
         H_conf,
         p.thres = GRAPPLE_P_thresh,
         cor.mat = Input[[3]],
         niter = 20,
         diagnosis = T,
         plot.it = T
         )
	
	# Perform Radial_IVW estimation
	ivw_res_conf <- ivw_radial(rad_conf, alpha = 0.05, weights = 3, tol = 0.0001)

	# Perform Radial_Egger estimation
        egg_res_conf <- egger_radial(rad_conf,alpha = 0.05, weights = 3)
}




###### UNIVARIATE ANALYSIS : WRITING RESULTS



## Save analysis results in a dataframe
# Number of SNPs for each threshold used 
# Beta estimation , se and pvalue 
# Heterogeneity tests (Qstats) and I2GX stat
# Intercept of Radial_Egger 
# Number of modes for GRAPPLE

results <- data.frame(name = name, Egger_snp = nrow(R_data), GRAPPLE_snp = nrow(H_data), GRAPPLE_modes = length(GR_modes$raw.modes),
 GRAPPLE_est = GR_est$beta.hat, GRAPPLE_pval = GR_est$beta.p.value,
 Egger_est = egg_res$coef[2,1], Egger_pval = egg_res$coef[2,4],
 Egger_int = egg_res$coef[1,1], Egger_int_pval = egg_res$coef[1,4],
 Egger_Qstat = egg_res$qstatistic, Egger_Qpval = 1-pchisq(egg_res$qstatistic, egg_res$df), Egger_I2GX = IGX_egg,
 IVW_est =ivw_res$coef[4,1], IVW_pval =ivw_res$coef[4,4], 
 IVW_Qstat = ivw_res$qstatistic, IVW_Qpval = 1-pchisq(ivw_res$qstatistic, ivw_res$df))

write.table(results, paste0(res_path, name, "/UNIV_results"), quote = F, row.names = F, sep = "\t")




## If some sensitivity analysis is not required, ensure to run the code without giving the result

noout <- function(var_name){
	tryCatch(
  	{
    	var_name
  	},
  	error = function(e) {
    	NA
  	})
}


## Save sensitivity analysis in a dataframe 
# Number of SNPs for each threshold used
# Beta estimation, se and pvalue for each sensitivity analysis

sensitivity <- data.frame(name = name, Egger_snp = nrow(R_data), out_snp = noout(nrow(rad_dat_out)), conf_snp = noout(nrow(rad_conf)), assoc_out_SNP = noout(nrow(H_high)),
 GRAPPLE_est = GR_est$beta.hat, GRAPPLE_est_out =noout(GR_est_out$beta.hat),GRAPPLE_est_confound = noout(GR_est_conf$beta.hat), GRAPPLE_est_assoc_out = noout(GR_high$beta.hat),
 GRAPPLE_pval = GR_est$beta.p.value, GRAPPLE_pval_out=noout(GR_est_out$beta.p.value), GRAPPLE_pval_confound = noout(GR_est_conf$beta.p.value), GRAPPLE_pval_assoc_out = noout(GR_high$beta.p.value),
 Egger_est = egg_res$coef[2,1], Egger_est_out = noout(egg_res_o$coef[2,1]), Egger_est_confound = noout(egg_res_conf$coef[2,1]), Egger_est_assoc_out = egg_high$coef[2,1],
 Egger_pval = egg_res$coef[2,4], Egger_pval_out = noout(egg_res_o$coef[2,4]), Egger_pval_confound = noout(egg_res_conf$coef[2,4]), Egger_pval_assoc_out = egg_high$coef[2,4],
 IVW_est = ivw_res$coef[4,1], IVW_est_out = noout(ivw_res_o$coef[4,1]), IVW_est_confound = noout(ivw_res_conf$coef[4,1]), IVW_assoc_out = ivw_high$coef[4,1],
 IVW_pval = ivw_res$coef[4,4], IVW_pval_out = noout(ivw_res_o$coef[4,4]), IVW_pval_confound = noout(ivw_res_conf$coef[4,4]), IVW_pval_assoc_out = ivw_high$coef[4,4])  

write.table(sensitivity, paste0(res_path, name, "/UNIV_sensitivity_", n_conf), quote = F, row.names = F, sep = "\t")



## Save plots in a pdf 

pdf(paste0(res_path, name, "/UNIV_plots.pdf"))

# QQplot of GRAPPLE
print(GR_plot)

## Comboplot of radial IVW, radial Egger and GRAPPLE estimation
print(Comboplot)
dev.off()


#### MULTIVARIATE ANALYSIS 

if(make_MULTI){
	## Load the needed libraries 
	
	library(MRBEE)
	library(RMVMR)
	library(MVMR)
	
	
	### Assign the file pathways for selection, exposure and outcome
	
	sel_cov_file_M = args[6:(5+n_cov)]
	exp_cov_file_M = args[(6+n_cov):(5+2*n_cov)]
	sel_file_M = c(args[1], sel_cov_file_M)
	exp_file_M = c(args[2], exp_cov_file_M)
	out_file_M = out_file
	
	name_exp = name
	cov_name = name_covariate
	name_conf <- paste(cov_name, collapse="_")
	
	###### SNPs FILTERING
	
	### Filter SNPs associated with at least one of the exposure, according to GRAPPLE p-threshold and clump
	
	
	Input_M <- getInput(sel.files = sel_file_M,
  	exp.files = exp_file_M,
  	out.files = out_file_M,
  	plink_refdat = plink_ref,
  	max.p.thres = GRAPPLE_P_thresh,
  	cal.cor = T,
  	p.thres.cor = 0.5,
  	get.marker.candidates = T,
  	marker.p.thres = 1e-05,
  	clump_r2 = 0.001,
  	plink_exe = plink_bin
	)
	
	
	
	H_data_M <- Input_M[[1]]
	
	## Models are not really happy with these SNPs, so better to remove them to have a lot of model hapiness
	
	H_data_M <- H_data_M[!(H_data_M$se_out1 == 0),]
	H_data_M <- H_data_M[!(H_data_M$gamma_out1 == 0),]
	
	
	## Save the selected and pruned SNPs, if something went wrong, or for further analysis 
	
	write.table(H_data_M, paste0(res_path,name,"/MULTI_selected_and_pruned_SNPs_", GRAPPLE_P_thresh), quote = F, row.names = F, sep = "\t")
	
	#### MULTIVARIATE ANALYSIS : PRIMARY ANALYSIS
	
	## Estimation results of GRAPPLE at the GRAPPLE p-threshold
	
	GR_est_M <- grappleRobustEst(
  	H_data_M,
  	p.thres = GRAPPLE_P_thresh,
  	cor.mat = Input_M[[3]],
  	niter = 30,
  	diagnosis = T,
  	plot.it = T
	)
	
	## Estimation of Radial IVW MVMR at the stringent p-threshold
	
	
	#Filter SNPs at stringent p-threshold
	H_rad_M <- H_data_M[H_data_M$selection_pval < 5e-8,]
	
	
	#Format for radial estimation
	RM_dat_M <- format_rmvmr(BXGs = H_rad_M %>% select(matches("gamma_exp")),
 	BYG = H_rad_M %>% select(matches("gamma_out")),
 	seBXGs = H_rad_M %>% select(matches("se_exp")),
 	seBYG = H_rad_M %>% select(matches("se_out")), RSID = H_rad_M$SNP)
	
	
	#Perform Radial IVW MVMR analysis
	RM_ivw_M <- ivw_rmvmr(RM_dat_M)
	
	
	##### MULTIVARIATE ANALYSIS : SENSITIVITY ANALYSIS - INSTRUMENT STRENGTH TEST 
	
	
	## Perform the modified cochran Q stat for instrument strength (Are overall SNPs good predictor for an exposure, conditionnaly on each covariate exposure in turn ?)
	
	# Take the covariance matrix for each exposure from the GRAPPLE covariance matrix estimation
	phe_cov <- Input_M[[3]][-nrow(Input_M[[3]]), -nrow(Input_M[[3]])]
	
	
	#Approximate the covariance matrix for each SNPs with this covariance matrix (just an approximation compared to the estimation of the covariance matrix with individual data)
	phe_cov_snp <- list()
	
	for(x in 1:nrow(RM_dat_M)){
        	phe_cov_snp[[x]] <- phe_cov
	}
	
	
	#Perform the modified cochran Q analysis for each exposure
	RM_str <- noout(strength_rmvmr(RM_dat_M, phe_cov_snp))
	
	if(class(RM_str) != "logical"){
		
		
		###### MULTIVARIATE ANALYSIS : SENSITIVITY ANALYSIS - HETEROGENEITY TESTS
		
		
		## Perform the Cochran Q stat adapted to the multivariate setting for GRAPPLE 
		
		# Assign the function 
		
		computeQ <- function(dat.list, p.thres = NULL) {
		
    		data <- dat.list$data
    		if (!is.null(p.thres))
        		data <- data[data$selection_pvals < p.thres, ]
		
    		nsnps <- nrow(data)
    		cor.mat <- dat.list$cor.mat
    		nn <- colnames(data)
    		b_exp <- as.matrix(data[, grep("gamma_exp", nn), drop = F])
    		se_exp <-  as.matrix(data[, grep("se_exp", nn), drop = F])
    		k <- ncol(b_exp)
		
    		Q.stats <- sapply(1:k, function(i) {
        		print(paste("Computing the conditional Q-stat for exposure ", i))
        		temp.dat <- cbind(b_exp[, -i], b_exp[, i], se_exp[, -i], se_exp[, i])
        		colnames(temp.dat) <- c(paste0("gamma_exp", 1:(k-1)), "gamma_out1",
                                		paste0("se_exp", 1:(k-1)), "se_out1")
        		ss <- 1:k
        		temp.cor <- cor.mat[c(ss[-i], i), c(ss[-i], i)]
		
        		delta <- grappleRobustEst(temp.dat, cor.mat = temp.cor,
                                  		diagnosis = F)$beta.hat
		
        		upper <- (b_exp[, i] - b_exp[, -i, drop = F] %*% delta)^2
        		temp <- t(t(cbind(se_exp[, -i], se_exp[, i])) * c(-delta, 1))
        		lower <- rowSums((temp %*% temp.cor) * temp)
        		Q.stat <- sum(upper/lower)
        		return(Q.stat)
    		})
    		names(Q.stats) <- paste0("exposure", 1:k)
		
    		df <- nrow(data) - 1
    		Q.stats.pval <- pchisq(Q.stats, df, lower.tail = F)
    		Q.stats.log.pval <- pchisq(Q.stats, df, lower.tail = F, log.p = T)
    		return(list(Q.stats = Q.stats,
                		df = df,
                		Q.stats.pval = Q.stats.pval,
                		Q.stats.log.pval = Q.stats.log.pval))
		}
		
		#Perform the Cochran Q stat
		GR_Q <- noout(computeQ(Input_M))
		
		## Perform the Cochran Q stat adapted to the multivariate setting
		RM_pleio <- noout(pleiotropy_rmvmr(RM_dat_M, RM_ivw_M))
		
		
		
		###### MULTIVARIATE ANALYSIS : SENSITIVITY ANALYSIS - OUTLIERS 
		
		## Sensitivity analysis without outliers for GRAPPLE 
		
		if(nrow(GR_est_M$outliers) != 0){
		
			#Remove outlier SNPs
        		H_data_out_M <- H_data_M[!(H_data_M$SNP %in% rownames(GR_est_M$outliers)),]
		
			#Perform GRAPPLE analysis without outliers
        		GR_est_out_M <- grappleRobustEst(
        		H_data_out_M,
        		p.thres = GRAPPLE_P_thresh,
        		cor.mat = Input_M[[3]],
        		niter = 40,
        		diagnosis = T,
        		plot.it = T
        		)
		
			#Save the outlier SNPs for further analysis
        		write.table(rownames(GR_est_M$outliers), paste0(res_path,name_exp,"/MULTI_GRAPPLE_outliers_", name_conf), quote = F, row.names = F, col.names = F, sep = "\t")
		}
		
		
		## Sensitivity analysis without outliers for Radial IVW
		
		if(length(RM_pleio$qdat[which(RM_pleio$qdat$qj_p < 1e-10),"snp"]) == 0){
		
        		rad_dat_out_M <- RM_dat_M[!(RM_dat_M$SNP %in% RM_pleio$qdat[which(RM_pleio$qdat$qj_p < 1e-2),"snp"]),]
		
        		ivw_res_o_M <- ivw_rmvmr(rad_dat_out_M)
		
        		write.table(RM_ivw_M$outliers$SNP, paste0(res_path,name_exp,"/MULTI_ivw_outliers", name_conf), quote = F, row.names = F, col.names = F, sep = "\t")
		}
		
		
		
		###### MULTIVARIATE ANALYSIS : SENSITIVITY ANALYSIS - OVER ASSOCIATION OF FILTERED_SNPs/OUTCOME (removing SNPs associated with the outcome)
		
		
		# Load the outcome file
		high_out_M <- data.frame(fread(out_file_M))
		
		# Remove SNPs associated with the outcome
		high_out_M <- high_out_M[high_out_M$V8 < 5e-6,]
		
		H_high_M <- H_data_M[!(H_data_M$SNP %in% high_out_M$V3),]
		
		
		## Perform GRAPPLE analysis without SNPs associated with the outcome
		GR_high_M <- grappleRobustEst(
  		H_high_M,
  		p.thres = GRAPPLE_P_thresh,
  		cor.mat = Input_M[[3]],
  		niter = 20,
  		diagnosis = T,
  		plot.it = T
		)
		
		
		# Filter SNPs with a stringent p-threshold
		H_high_M <- H_high_M[H_high_M$selection_pval < 5e-8,]
		
		
		# Format for radial estimation
		rad_high_M <- format_rmvmr(BXGs = H_high_M %>% select(matches("gamma_exp")),
 		BYG = H_high_M %>% select(matches("gamma_out")),
 		seBXGs = H_high_M %>% select(matches("se_exp")),
 		seBYG = H_high_M %>% select(matches("se_out")), RSID = H_high_M$SNP)
		
		
		# Perform Radial_IVW MVMR without SNPs associated with the outcome
		ivw_high_M <- ivw_rmvmr(rad_high_M)
		
		
		
		
		###### MULTIVARIATE ANALYSIS : PLOTS 
		
		RM_plot <- plot_rmvmr(RM_dat_M, RM_ivw_M)
		
		plot1 <- RM_plot$p1
		
		plot2 <- RM_plot$p2
		
		#QQplot of GRAPPLE without outliers
        		GR_out_plot_M <- GR_est_out_M$p
		
		################################################# MR HORSE IF YOU HAVE TIME 

		if(horse_try == T){

			## Set up needed libraries
			library(R2jags)
			source(paste0(working_dir, "/mr_horse.R"))
		

			## Format data for Horse 

			Horse_dat <- H_rad_M[,!(names(H_rad_M) %in% c("SNP","effect_allele","other_allele", "selection_pvals"))]

			beta_horse <- paste0("betaX", 1:(n_cov+1))
			se_horse <- paste0("betaX", 1:(n_cov+1), "se")

				
			colnames(Horse_dat) <- c(beta_horse, se_horse, "betaY", "betaYse")

			## Perform the MVMR Horse algorithm estimation 
			MVMREx <- mvmr_horse(Horse_dat)

			
			## Save the analysis of MR horse in a separated dataframe
			## Results contain beta, se, pval and the degree of convergence of the algorithm
				
			sensitivity <- data.frame(name = c(name_exp, cov_name), n_snps = nrow(Horse_dat),
 			 Horse_est = MVMREx[[1]][,2],
 			 Horse_sd = MVMREx[[1]][,3],
 			 Horse_pval = 2*(1-pt(abs(MVMREx[[1]][,2]/MVMREx[[1]][,3]), nrow(Horse_dat)-(length(cov_name)+1))),
 			 Horse_converge = MVMREx[[1]][,6])

				
			write.table(sensitivity, paste0(res_path, name_exp, "/MULTI_horse_", name_conf), quote = F, row.names = F, sep = "\t")		
		}		


		###### MULTIVARIATE ANALYSIS : WRITING RESULTS
		
		
		
		## Save analysis results in a dataframe
		# Number of SNPs for each threshold used
		# Beta estimation , se and pvalue
		# Heterogeneity test
		# Instrument strength test

		
		## Recup some list items 
		rec_list <- function(lst, r, c) {
 	 		selected_elements <- lapply(lst, function(df) df[[c]][r])
  			return(selected_elements)
		}

		r = 1:n_cov

		list_str <- rec_list(RM_str[[2]], r = r, c = "p_value")
		str_values <-  do.call(rbind, list_str)
		colnames(str_values) <- paste0("IVW_strength_", r)
		
		
		results <- data.frame(name = c(name_exp, cov_name), GRAPPLE_snp = nrow(H_data_M), IVW_snp = nrow(H_rad_M),
 		GRAPPLE_beta = GR_est_M$beta.hat, GRAPPLE_pval = GR_est_M$beta.p.value,
 		GRAPPLE_Qstat = GR_Q$Q.stats, GRAPPLE_Qpval = GR_Q$Q.stats.pval,
 		IVW_beta = RM_ivw_M$coef[,1], IVW_pval = RM_ivw_M$coef[,4],
 		IVW_Qstat = RM_pleio$gq[c(1:(n_cov+1)),1],
 		IVW_Qpval = RM_pleio$gq[c(1:(n_cov+1)),2],
 		str_values)
		
		write.table(results, paste0(res_path, name_exp, "/MULTI_results_", name_conf), quote = F, row.names = F, sep = "\t")
		
		
		## Save sensitivity analysis in a dataframe 
		# Number of SNPs for each threshold used
		# Beta estimation, se and pvalue for each sensitivity analysis
		
		
		sensitivity <- data.frame(name = c(name_exp, cov_name), ivw_snp = nrow(RM_dat_M), out_snp = noout(nrow(rad_dat_out_M)), assoc_out_SNP = noout(nrow(H_high_M)),
 		GRAPPLE_est = GR_est_M$beta.hat, GRAPPLE_est_out =noout(GR_est_out_M$beta.hat), GRAPPLE_est_assoc_out = noout(GR_high_M$beta.hat),
 		GRAPPLE_pval = GR_est_M$beta.p.value, GRAPPLE_pval_out=noout(GR_est_out_M$beta.p.value), GRAPPLE_pval_assoc_out = noout(GR_high_M$beta.p.value),
 		IVW_est = RM_ivw_M$coef[,1], IVW_est_out = noout(ivw_res_o_M$coef[,1]), IVW_assoc_out = ivw_high_M$coef[,1],
 		IVW_pval = RM_ivw_M$coef[,4], IVW_pval_out = noout(ivw_res_o_M$coef[,4]), IVW_pval_assoc_out = ivw_high_M$coef[,4])
		
		write.table(sensitivity, paste0(res_path, name_exp, "/MULTI_sensitivity_", name_conf), quote = F, row.names = F, sep = "\t")
		
		
		
		
		
		## Save plots in a pdf
		
		pdf(paste0(res_path, name_exp,"/MULTI_","plots_",name_conf,".pdf"))
		
		#QQplot of GRAPPLE 
		print(GR_out_plot_M)
		
		## First plot of radial IVW MVMR without correction
		print(plot1 + geom_abline(intercept = 0, slope = GR_est_M$beta.hat[1], col = "green") + geom_abline(intercept = 0, slope = GR_est_M$beta.hat[2], col = "purple"))
		
		## 2nd plot of radial IVW MVMR with correction (correcting with heterogenity to better assign each SNP to it's corresponding exposure)
		print(plot2 + geom_abline(intercept = 0, slope = GR_est_M$beta.hat[1], col = "green") + geom_abline(intercept = 0, slope = GR_est_M$beta.hat[2], col = "purple"))
		
		dev.off()
	}	
}	
	


