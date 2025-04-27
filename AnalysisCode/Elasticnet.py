import warnings
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.datasets import fetch_california_housing
from sklearn.linear_model import ElasticNet
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split,GridSearchCV
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr
from joblib import dump,load

warnings.filterwarnings('ignore')

ecgDf = pd.read_csv('~/ecgTrait/ecg95.txt') #input file contians 96 columns: sample ID, and 95 columns corresponding 95 replicable ECG trait value 
ecgList = list(ecgDf.columns)[1:]
idpDf = pd.read.csv('~/idpTrait/idp2095.txt') #input file contians 2096 columns: sample ID, and 2095 columns corresponding 2096 replicable IDP value 
ecgIDP = pd.merge(ecgDf, idpDf, on='sample ID')

#Split: Train set: 80% ~20145, Test set: 20%~5000
train = ecgIDP.sample(20145, random_state=888)
print(len(train))
test = ecgIDP[~ecgIDP['sample ID'].isin(train['sample ID'].tolist())]
print(len(test))

#add cov
cov0Df = pd.read.csv('~/covariate/Cox0Covariate.txt') #input file contians multiple columns: sample ID, sex, age
covDf = pd.read.csv('~/covariate/CoxCovariate.txt') #input file contians multiple columns: sample ID and multiple columns corresponding to covariates (exclude sex and age)


train = pd.merge(train, cov0Df, on='sample ID')
train = pd.merge(train, covDf, on='sample ID')

test = pd.merge(test, cov0Df, on='sample ID')
test = pd.merge(test, covDf, on='sample ID')

#scaled
idp=['25001']
cov1=['age','bmi','whr','pulse','activityDuration'] # Continuous variables
cov2=['sex','alcoholFreq','smokeFreq','insomnia'] # Categorical variables

train_scaled = train[['sample ID']+idp+ecgList+cov1+cov2]
scaler = StandardScaler()
train_scaled.iloc[:,1:-4] = scaler.fit_transform(train_scaled.iloc[:,1:-4])
train_X = np.array(train_scaled[ecg+otherCov1+otherCov2])
train_Y = np.array(train_scaled[idp])

test_scaled = test[['sample ID']+idp+ecgList+cov1+cov2]
test_scaled.iloc[:,1:-4] = scaler.fit_transform(test_scaled.iloc[:,1:-4])
test_X = np.array(test_scaled[ecg+otherCov1+otherCov2])
test_Y = np.array(test_scaled[idp])

#elasticnet model
model=ElasticNet()
param_grid={
    'alpha':[0.001,0.01,0.1,1,10],
    'l1_ratio':[0.2,0.5,0.8],
}

#Turing hyperparameters of the model through five-fold cross-validation on the training dataset
grid_search=GridSearchCV(model, param_grid,cv=5,return_train_score=True)
grid_search.fit(train_X,train_Y)

print('best_params:',grid_search.best_params_)
best_model=grid_search.best_estimator_
dump(best_model,"~/model/elasticnet_idp_%s.pkl"%idp[0])

#best model performance
train_predictions=best_model.predict(train_X)
test_predictions=best_model.predict(test_X)
train_rmse=np.sqrt(mean_squared_error(train_Y,train_predictions))
test_rmse=np.sqrt(mean_squared_error(test_Y,test_predictions))
train_r2=r2_score(train_Y,train_predictions)
test_r2=r2_score(validate_Y,test_predictions)
train_pcc=pearsonr(train_Y,train_predictions)
test_pcc=pearsonr(validate_Y,test_predictions)
record_train={'idp':idp[0],'type':'train','RMSE':train_rmse,'R2':train_r2,'PCC':train_pcc[0],'PCC_pvalue':train_pcc[1]}
record_test={'idp':idp[0],'type':'test','RMSE':test_rmse,'R2':test_r2,'PCC':test_pcc[0],'PCC_pvalue':test_pcc[1]}
print(f'train RMSE: {train_rmse:.4f}, R2: {train_r2:.4f}, PCC: {train_pcc[0]:.4f},{train_pcc[1]}')
print(f'test RMSE: {test_rmse:.4f}, R2: {test_r2:.4f}, PCC: {test_pcc[0]:.4f},{test_pcc[1]}')

#plot 
cv_results=grid_search.cv_results_
mean_train_score=cv_results['mean_train_score']
mean_test_score=cv_results['mean_test_score']

plt.figure(figsize=(10,6))
plt.plot(mean_train_score,label='Train score')
plt.plot(mean_test_score,label='Test score')
plt.xlabel('params')
plt.ylabel('score')
plt.legend()
#plt.show()
plt.savefig('~/performance_png/idp%s_performace.png'%idp[0])
plt.figure(figsize=(10,6))
plt.scatter(validate_Y,test_predictions)
plt.xlabel('True')
plt.ylabel('Predict')
plt.plot([min(train_Y.min(),validate_Y.min()),max(train_Y.max(),validate_Y.max())],[min(train_Y.min(),validate_Y.min()),max(train_Y.max(),validate_Y.max())],'k--',lw=4)
#plt.show()
plt.savefig('~/true_vs_predict_png/idp%s_True.vs.Predict.png'%idp[0])

errors=validate_Y-test_predictions
plt.figure(figsize=(10,6))
sns.histplot(errors,kde=True)
plt.xlabel('predicted error')
plt.ylabel('Frequency')
#plt.show()
plt.savefig('~/error_distribution_png/idp%s_Error_distribution.png'%idp[0])

