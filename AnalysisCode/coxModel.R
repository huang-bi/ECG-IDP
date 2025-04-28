library('survival')
library('survminer')
library(cutoff)


data<-read.table('~/cox/idp25001AF.txt',sep="\t",header=TRUE) #file contians at least 5 columns: sample ID, target trait value, disease state, survival time, sample keep state amd multiple columns corresponding covariates
data<-data[data$keep=='True',] #exclude sample with target disease event prior to attending assessment 
print(dim(data[data$diseaseState==1,]))
print(dim(data[data$diseaseState==0,]))

#scale continous variates
vars_to_scale <- c("height", "weight", "waist", "hip", "whr", "bmi", "age", "pulse", 
                   "diastolic", "systolic", "activityDuration","targetTrait")

data[paste0(vars_to_scale, ".scaled")] <- lapply(data[vars_to_scale], scale)


library(car)

resultFun <- function(model){
  sum.cox<-summary(model)
  #print(sum.cox)

  conf <- as.data.frame(sum.cox$conf.int)
  conf <- conf[,-which(names(conf)%in%c("exp(coef)","exp(-coef)"))]
  cox_result<-cbind(as.data.frame(sum.cox$coefficients),conf)
  cox_result$label<- rownames(cox_result)

  #schoenfeld residuals 
  res.ph <- cox.zph(model)
  zph_result<-as.data.frame(res.ph$table)
  #print(zph_result)
  zph_result$label<- rownames(zph_result)
  cox_result<-merge(cox_result,zph_result,by='label',all=TRUE,sort=FALSE)

  #Calculating vif
  vif_values <- vif(model)
  vif_result <- as.data.frame(vif_values)
  vif_result$label<- rownames(vif_result)
  #print(vif_result)
  cox_result<-merge(cox_result,vif_result,by='label',all=TRUE,sort=FALSE)

  #Cindex
  cIndex <- as.data.frame(sum.cox$concordance)
  cox_result$cIndex<-cIndex['C',]
  cox_result$'cIndex(SE)'<-cIndex['se(C)',]

  #Significance test
  logtest <- as.data.frame(sum.cox$logtest)
  cox_result$logtest<-logtest['test',]
  cox_result$logtestDf<-logtest['df',]
  cox_result$logtestP<-logtest['pvalue',]

  #sample size
  cox_result$nCase <- sum.cox$nevent
  cox_result$n <- sum.cox$n
  return(cox_result)
}



#Cox model
testTrait <- c('targetTrait.scaled','bmi.scaled','pulse.scaled','diastolic.scaled','systolic.scaled',
              'activityDuration.scaled','alcoholFreq','smokeFreq','insomnia')

#base model
univ_formulas <- sapply(testTrait,function(x) as.formula(paste('Surv(time,diseaseState)~',x,'+age.scaled+sex')))  
univ_models <- lapply(univ_formulas,function(x) coxph(x,data=data))
univ_results <- lapply(univ_models,function(x) resultFun(x))

# also check covariate
for (trait in testTrait){
  resTmp<-univ_results[[trait]]
  if (resTmp[resTmp$label==trait,'Pr(>|z|)']>=0.05){
    print(trait)
  } 
  write.csv(resTmp,paste0('~/cox/AF_',trait,'baseModel.csv'))
}

#additional model
trait = '25001.scaled' #target idp trait
finFormula<-as.formula(paste('Surv(time,diseaseState)~',trait,'+age.scaled+sex+bmi.scaled+pulse.scaled+alcoholFreq+insomnia'))
res.cox<-coxph(finFormula,data=data)
result<-resultFun(res.cox)
write.csv(result,paste0('~/cox/AF_',trait,'additionallyModel.csv'))


#find best cut off
res.cut <- logrank(data = data,
                           time = "time",
                           y = "diseaseState", 
                           x = trait,
                           cut.numb = 1, 
                           n.per = 0.1, 
                           y.per = 0.005, 
                           p.cut = 0.05,
                           strict=TRUE,
                           round=20) 

cutoffValue <- res.cut[order(res.cut$pvalue,decreasing  = F),][1,1])
print(cutoffValue)

data$Group = ifelse(data[[trait]] <= cutoffValue,"low","high")

fit <- survfit(Surv(time, diseaseState) ~ Group, data = data)
p<-ggsurvplot(fit, data = data,  
           conf.int = TRUE,  pval = TRUE, 
           
           legend.labs=c("IDP High     ", "IDP Low"),
           risk.table = TRUE, palette = c('#ff6672','#00d2f9')) 

print(p)

p$plot<-p$plot+coord_cartesian(ylim=c(0.85,1))+labs(x='Time (Days)',y='Non-Disease Probability')+
        annotate('text',x=1000,y=0.9,size=5,label = paste0("C-Index: ",round(as.data.frame(summary(res.cox)$concordance)['C',],3),'\n P-value: ',sprintf('%.2e',as.data.frame(summary(res.cox)$logtest)['pvalue',])))+
        theme(axis.title = element_text(size = 20),axis.text = element_text(size = 18),legend.text=element_text(size=14),legend.title=element_blank())

print(p$plot)
