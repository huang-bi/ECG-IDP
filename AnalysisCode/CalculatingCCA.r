library('CCA')
library('CCP')


ecgDf <- read.csv('~/traitGroup/ecgRS.txt',sep="\t")   #input file contains at least two columns: a Sample ID column and multiple columns corresponding to all traits within the same group.
idpDf <- read.csv('~/traitGroup/idpWhiteMatter.txt',sep="\t")    #input file contains at least two columns: a Sample ID column and multiple columns corresponding to all traits within the same group.

sampleID <- intersect(ecgDf$'sample ID',idpDf$'sample ID')

ecg <- ecgDf[ecgDf$'sample ID' %in% sampleID,-1]
idp <- idpDf[idpDf$'sample ID' %in% sampleID,-1]

#scaled
ecg <- data.frame(scale(ecg))
idp <- data.frame(scale(idp))

#CCA
ccaRes<-cc(ecg,idp)
print(ccaRes)
rho <- ccaRes$cor
n <- dim(ecg)[1]
p <- length(ecg)
q <- length(idp)
#significant test
sigTest<-p.asym(rho,n,p,q,tstat="Wilks")
print(sigTest)

