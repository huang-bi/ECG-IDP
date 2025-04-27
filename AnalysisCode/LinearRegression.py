import pandas as pd
import numpy as np
from scipy.stats import median_abs_deviation
import os
import statsmodels.api as sm


def filterMAD5(df,colName):
	tmp = df.dropna()
    mad5_upper = tmp[colName].median()+5*median_abs_deviation(tmp[colName])
    mad5_lower = tmp[colName].median()-5*median_abs_deviation(tmp[colName])
    tmp = tmp[(tmp[colName]>=mad5_lower)&(tmp[colName]<=mad5_upper)]
    return tmp


ecgDf = pd.read_csv('~/ecgTrait/ecgV1RS.txt') #input file contians two columns: sample ID, trait value 
idpDf = pd.read.csv('~/idpTrait/idp25001.txt') #input file contians two columns: sample ID, trait value 
covDf = pd.read.csv('~/covariate/ukbCovariate.txt') #input file contians multiple columns: sample ID and multiple columns corresponding to covariates
covList = list(covDf.columns[1:])

#exclude outliers
ecgCol = 'V1RS'
ecgDfMAD5 = filterMAD5(ecgDf, ecgCol)
idpCol = '25001'
idpDfMAD5 = filterMAD5(idpDf, idpCol)

#match sample
dataMatrix = pd.merge(ecgDfMAD5,idpDfMAD5,on='sample ID')
dataMatrix = pd.merge(dataMatrix,covDf,on='sample ID')
dataMartix = dataMatrix.drop(columns=['sample ID'])

#normalized
dataMatrix = dataMatrix.apply(lambda x: (x - np.mean(x)) / np.std(x),axis=0)
dataMatrix = dataMatrix.fillna(0)

#linear regression
Y = dataMatrix[[idpCol]]
X = dataMatrix[[ecgCol]+covList]
X = sm.add_constant(X)
model = sm.OLS(Y,X)
results = model.fit()
print(results)

