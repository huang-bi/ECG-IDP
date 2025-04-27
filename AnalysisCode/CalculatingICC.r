library(irr)

#calculating test-retest reliability
inputFile = '~/revisitedData/ecgV1RS.txt' #input file contians three columns: sample ID, trait value in first visit, trait value in revisit

traitDf <- read.csv(inputFile,sep="\t")
icc_result <- icc(traitDf[,-1], model="twoway",type="agreement",unit="single")
print(icc_result)


#calculating degree of agreement within the same grou
inputFile = '~/traitGroup/ecgRS.txt' #input file contains at least two columns: a Sample ID column and multiple columns corresponding to all traits within the same group.

traitDf <- read.csv(inputFile,sep="\t")

icc_result <- icc(traitDf[,c(-1)], model="twoway",type="consistency",unit="average")
print(icc_result)
