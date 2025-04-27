library(TwoSampleMR)

exposureQC <- function(exposure_data_path, exposure_name, exposure_clump_path, sampleSize=NULL,snp_col='SNP',beta_col='BETA',se_col='SE',
    ea_col='ALLELE1',ra_col='ALLELE0',eaf_col='A1FREQ',pval_col='P_BOLT_LMM_INF'){

    #read gwas file
    if (!file.exists(exposure_data_path)){
        print(paste(exposure_data_path,'not find!',sep=""))
        q()
    } else {
        exposure_data <- read.table(exposure_data_path,sep="\t",header=T)
        print(paste('gwas',nrow(exposure_data),sep=": "))
        if (nrow(exposure_data)>0){
            exposure_data$PHENOTYPE <- exposure_name
            if (!is.null(sampleSize)){
                exposure_data$N <- sampleSize 
            }
        }
    }
    #clumped
    if (!file.exists(exp_clump_path)){
        print(paste(exp_clump_path,'not find!',sep=""))
        print(paste('clumped failed, quit!', exposure_name, sep=""))
        q()
    } else {
        clump_data <- read.table(exp_clump_path,sep=" ",header=T)
        exposure_data_clumped <- exposure_data[exposure_data[[snp_col]] %in% clump_data$SNP,]
        print(paste('after clump',nrow(exposure_data_clumped),sep=": "))
    }
    # Ftest
    if (nrow(exposure_data_clumped)>0){
        exposure_data_F <- exposure_data_clumped
        exposure_data_F$R2 <- 2*(exposure_data_clumped[[beta_col]]^2)*exposure_data_clumped[[eaf_col]]*(1-exposure_data_clumped[[eaf_col]])
        exposure_data_F <- transform(exposure_data_F,F=(N-2)*R2/(1-R2))
        exposure_data_F <- exposure_data_F[exposure_data_F$F > Ffilter,]
        print(paste('after F',nrow(exposure_data_F),sep=": "))
    } else {
        print(paste('F filter failed, quit!', exposure_name, sep=""))
        q()
    }
    #convert to exposure data format
    if (nrow(exposure_data_F)>0){
        exp_dat <- format_data (
            dat = exposure_data_F,
            type = "exposure",
            phenotype_col = "PHENOTYPE",
            snp_col = snp_col,
            beta_col =  beta_col,
            se_col = se_col,
            effect_allele_col = ea_col,
            other_allele_col = ra_col,
            eaf_col = eaf_col,
            pval_col = pval_col,
            samplesize_col = "N",
        ) 
        return (exp_dat)
    } else {
        print(paste('No snp passed F filter, quit!', exposure_name, sep=""))
        q()
    }
}

outcomeAdd <- function(outcome_data_path,outcome_name,exp_dat,snp_col='rsids',beta_col='beta',se_col='sebeta',
    ea_col='alt',ra_col='ref',eaf_col='af_alt',pval_col='pvalue',samplesize_col='af_alt_cases',sampleSize=NULL){
    #convert to outcome data format
    outcome_dat <- read_outcome_data (
            snps = exp_dat$SNP,
            filename = outcome_data_path,
            sep = "\t",
            snp_col = snp_col,
            beta_col = beta_col,
            se_col = se_col,
            effect_allele_col = ea_col,
            other_allele_col = ra_col,
            eaf_col = eaf_col,
            pval_col = pval_col,
            samplesize_col = samplesize_col)

    if (nrow(outcome_dat)>0){
        outcome_dat$outcome <- outcome_name
        if (!is.null(sampleSize)){
            outcome_dat$samplesize.outcome <- sampleSize 
        }
        return (outcome_dat)
    } else {
        print('No common snp in outcome data, quit!')
        q()
    }
}

run_mr <- function(data, mr_method, label=''){

    #run multiple mr method in once
    res <- data.frame()
    for (m in mr_method){
        tryCatch({
            restmp <- mr(data,method_list=c(m))
            res <- rbind(res,restmp)
        },error=function(e){
            cat("MR ERROR :",conditionMessage(e),"\n")
        })
    }

    #MR-PRESSO test
    presso <- 'Failed'
    tryCatch({
        presso <- run_mr_presso(data)
    },error = function(e){
        cat("PRESSO ERROR :",conditionMessage(e),"\n") 
    }) 

    #Heterogeneity test
    het <- 'Failed'
    tryCatch({
        het <- mr_heterogeneity(data)
    },error = function(e){
        cat("Heterogeneity ERROR :",conditionMessage(e),"\n") 
    }) 

    #pleiotropy test
    pleio <- 'Failed'
    tryCatch({
        pleio <- mr_pleiotropy_test(data)
    },error = function(e){
        cat("Pleiotropy ERROR :",conditionMessage(e),"\n") 
    }) 

    #direction test
    steiger <- 'Failed'
    tryCatch({
        steiger<-directionality_test(data)
    },error=function(e){
        cat("Steiger ERROR :",conditionMessage(e),"\n") 
    })


    #single snp test
    sin <- 'Failed'
    tryCatch({
        sin <- mr_singlesnp(data)
    },error = function(e){
        cat("Single SNP ERROR :",conditionMessage(e),"\n") 
    }) 

    #leave one out test
    leaveoneout <- 'Failed'
    tryCatch({
        leaveoneout <- mr_leaveoneout(data)
    },error = function(e){
        cat("Leave-one-out ERROR :",conditionMessage(e),"\n") 
    }) 

    save(data,res,presso,het,pleio,steiger,sin,leaveoneout,file=paste(exp_name,'/',exp_name,'_',outcome_name,label,'.Rdata',sep=""))

    if (class(presso) != "character"){
        presso_indice=presso[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
        if (!is.null(presso_indice)&(!'No significant outliers' %in% presso_indice)){
            data <- data[-presso_indice,]
            run_mr(data, mr_method, label=paste('_PRESSO_',nrow(data),sep=''))
        }  
    }
}



pvalue <- '5e-8'
Ffilter <- 10

exposure_gwas <- "~/ecgGWAS/"
exposure_clump <- "~/ecgGWAS/5e-8_snplist" #clumping: use plink clump function
outcome_gwas <- "~/finngenGWAS/"

mr_method <- c('mr_ivw','mr_ivw_mre','mr_wald_ratio','mr_weighted_median','mr_weighted_mode','mr_egger_regression')

args <- commandArgs(trailingOnly=TRUE)

exp_name <- args[1]
if (!dir.exists(exp_name)){
    dir.create(exp_name)
}

outcome_name <- args[2]  

exp_data_path <- paste(exposure_gwas, paste(exp_name,'.txt',sep=""),sep = "/")
exp_clump_path <- paste(exposure_clump, paste(exp_name,'.txt.clumped.snplist',sep=""),sep = "/")
exp_dat <- exposureQC(exp_data_path, exp_name, exp_clump_path,snp_col='rsid',beta_col='beta',se_col='se',
    ea_col='a2',ra_col='a1',eaf_col='af',pval_col='P')

outcome_data_path <- paste(outcome_gwas,outcome_name,sep="")
finngenSampleSize <- args[3]
outcome_dat <- outcomeAdd(outcome_data_path,outcome_name,exp_dat,snp_col='rsids',beta_col='beta',se_col='sebeta',
    ea_col='alt',ra_col='ref',eaf_col='af_alt',pval_col='pvalue',samplesize_col='af_alt_cases',sampleSize=finngenSampleSize)

#harmonise
data <- harmonise_data(exposure_dat=exp_dat, outcome_dat=outcome_dat)

#run 
run_mr(data, mr_method)

