# AutoLasso: A Fast and Automatic Variable Selection Method for Multiple Linear Regression

Linear regression is one of the most widely used statistical methods and with the ease of modern data collection, often many predictors are associated with a response variable. However, most of the available variable selection methodologies (e.g., elastic net and its special case lasso) involve careful selection of the tuning parameter. This can cause numerous problems in practice and may lead to biased estimates. We present an automatic variable selection method using anew risk criteria which is simple and compares favorably with many of the currently available variable selection methods for linear models. Numerical illustrations based on simulated and real data are presented based on output from an user friendly R package `AutoLasso` that is made available for easy implementation.

Required R packages: `glmnet` and `MASS`

Author: Sujit Ghosh, Kaiyuan(Carl) Duan, and Guangjie Yu

## Proposed Methodology

Consider the familiar linear (multiple) regression model:  
<img src="https://render.githubusercontent.com/render/math?math=\LARGE{y_i = \beta_0 + \sum_{j=1}^p x_{ij}\beta_j + \sigma\epsilon_i}">,   
where <img src="https://render.githubusercontent.com/render/math?math=\epsilon_i">’s satisfy the usual Gauss-Markov assumptions. 


- Step 1: Standardize the response and predictor variables:  
<img src="https://render.githubusercontent.com/render/math?math=\sum_{i=1}^n x_{ij}=0"> and <img src="https://render.githubusercontent.com/render/math?math=\sum_{i=1}^n x_{ij}^2=n"> for all <img src="https://render.githubusercontent.com/render/math?math=j=1,\ldots,p"> and <img src="https://render.githubusercontent.com/render/math?math=\sum_{i=1}^n y_i=0\\">   
Let <img src="https://render.githubusercontent.com/render/math?math=Y=(y_1, y_2,\ldots,y_n)^t"> denote the centered response vector and <img src="https://render.githubusercontent.com/render/math?math=X=((x_{ij}))"> denotes the <img src="https://render.githubusercontent.com/render/math?math=n\times p"> centered and scaled design matrix   

- Step 2: Get an initial (consistent) estimator of <img src="https://render.githubusercontent.com/render/math?math=\beta"> (e.g., using the `glmnet`} R package)
and denote it by <img src="https://render.githubusercontent.com/render/math?math=\hat{\beta}^{(0)}">        

- Step 3: Choose a subset of variables minimizing the risks:  
(i) Order the absolute values of <img src="https://render.githubusercontent.com/render/math?math=\hat{\beta}^{(0)}">'s and choose first <img src="https://render.githubusercontent.com/render/math?math=k"> variables   
(ii) Compute the risks: <img src="https://render.githubusercontent.com/render/math?math=R_k = k||Y-X_k\hat{\beta}_k||^2/(n-k)"> where <img src="https://render.githubusercontent.com/render/math?math=\hat{\beta}_k"> is the least square estimate based on top <img src="https://render.githubusercontent.com/render/math?math=k"> variables <img src="https://render.githubusercontent.com/render/math?math=X_k"> selected in step (i) above   
(iii) Chose <img src="https://render.githubusercontent.com/render/math?math=k\in \{1,2,\ldots,r\}"> that minimizes the 
<img src="https://render.githubusercontent.com/render/math?math=R_k">, where <img src="https://render.githubusercontent.com/render/math?math=r<rank(X)">  

- Step 4: Output the <img src="https://render.githubusercontent.com/render/math?math=\hat{\beta}"> with <img src="https://render.githubusercontent.com/render/math?math=\hat{k}"> non-zero entries obtained from Step 3(iii) and rest entries set to zeros. Also output <img src="https://render.githubusercontent.com/render/math?math=\hat{\sigma}^2=||Y-X_{\hat{k}}\hat{\beta}_{\hat{k}}||^2/(n-\hat{k})">  

## Simulated Data Design

Generate data using following scenarios:  
<img src="https://render.githubusercontent.com/render/math?math=\LARGE{y_i =\beta_0 + \sum_{j=1}^p x_{ij}\beta_j + \sigma\epsilon_i}">, <img src="https://render.githubusercontent.com/render/math?math=\LARGE{\epsilon_i\stackrel{i.i.d}{\sim} N(0, 1)}">  

- Two Choices for correlation matrix <img src="https://render.githubusercontent.com/render/math?math=\Sigma_0">:  
 (i)  AutoRegressive of order 1 (AR(1)): <img src="https://render.githubusercontent.com/render/math?math=\Sigma_0=((\rho^{|i-j|}))">  
 (ii) Compound symmetry: <img src="https://render.githubusercontent.com/render/math?math=\Sigma_0=(1-\rho)I_p+\rho 1_p 1_p^t"> where <img src="https://render.githubusercontent.com/render/math?math=\rho\in [0, 1)">     

- Error variance <img src="https://render.githubusercontent.com/render/math?math=\sigma^2"> is chosen to have desired signal-to-noise ratio (SNR) :  
<img src="https://render.githubusercontent.com/render/math?math=\sigma^2 =||X\beta||_2^2/(n(snr)^2)">(e.g., <img src="https://render.githubusercontent.com/render/math?math=snr=2, 5, 10"> etc.)   

We explore various combinations of (<img src="https://render.githubusercontent.com/render/math?math=n,p,\rho, snr">) for above scenarios each based on 1000 simulated data replicates

## Performance Metrics for Variable Selection

We concerned Overall Accuracy Rates and Biases:

- Overall Accuracy Rates:  
  - True Positive Rate
  - False Positive Rate
  - False Discovery Rate

- Biases (Performance in reference)
  - Biase of active variables: <img src="https://render.githubusercontent.com/render/math?math=\{\hat{\beta_j}-\beta_j^{true}: \beta_j^{true} \neq 0, j=1,...,p\}">     
  - Relative error bias of <img src="https://render.githubusercontent.com/render/math?math=\sigma^2">: <img src="https://render.githubusercontent.com/render/math?math=\frac{\hat{\sigma}^2}{\sigma_{true}^2}-1">
  - All Biases: <img src="https://render.githubusercontent.com/render/math?math=\{\hat{\beta_j}-\beta_j^{true}: j=1,...,p\}">
 
## Accuracy Measures for two methods
- Case 1: n=100, p.true=5, p=2n, cor.type=Compound.Symmetry, rho=0.85, snr=5

| Methods | TPR  |   FPR  | FDR  |
|---------|------| -------|------|
|True     |  1   | 0.000  | 0.000|
|AutoLASSO|  1   | 0.000  | 0.002|
|    LASSO|  1   | 0.044  | 0.631|

- Case 2: Case 2: n=100, p.true=5, p=2n, cor.type=Auto.Regressive, rho=0.85, snr=5

| Methods | TPR  |   FPR  | FDR  |
|---------|------| -------|------|
|True     |  1   | 0.000  | 0.000|
|AutoLASSO|  1   | 0.000  | 0.000|
|    LASSO|  1   | 0.064  | 0.715|


## Application

To study the application of Automated LASSO approach, we compare with linear regression `lm` and Classic LASSO method. Data was taken from the R package `COVID19` available to the public for educational and scientific use. Due to the outbreak of corona-virus in the United States, exploding in March and stable in the beginning of June, we selected the data based on the date between March 1st and June 10th, 2020. The dataset presented in this paper is four predictors (selected from a preliminary analysis of a set of five variables) modeled on a selected response variable. There were 5202 observations. Conclusions about the data were not desired, the data was used only to observe how linear regression `lm`, Automated LASSO, and classic LASSO models perform on highly correlated data sets. All computation was done with R software using the package `glmnet` and `AutoLasso`.

- Number of variables selected by two methods:

| day.ahead | AutoLASSO|   LASSO  | lm      |  
|---------  |--------- | -------  |------   |
|1          |  2(0.99) | 3(0.99)  | 57(0.97)|
|5          |  4(0.98) | 46(0.97) | 57(0.95)|
|    7      |  8(0.98) | 47(0.96) | 57(0.94)|
|    10     |  21(0.95)| 46(0.95) | 57(0.93)|   

Kendall’s Tau measures the association predicted and hold observations      
**AutoLasso achieves almost same accuracy of prediction using much smaller subset of predictor variables**

## Reference

USPROC 2020 Electronic Undergraduate Statistics Research Conference: https://www.causeweb.org/usproc/eusrc/2020/virtual-posters/6


## Examples:
  > x=model.matrix(log1p(confirmed) ~ log1p(confirmed.lag)+day+cancel+internal+state-1,data=covid.data)  
  > y=log1p(confirmed)
# LASSO
  > cv.lasso.fit=cv.glmnet(x,y)  
  > lasso.fit.glmnet=glmnet(x,y)  
  > log.confirmed.pred=predict(cv.lasso.fit,newx=x[test,],type='response',s="lambda.1se")  
  > p.cor=cor(log1p(confirmed[test]),log.confirmed.pred)
# AutoLasso
  > autolasso.fit=auto.lasso(x,y)  
  > log.confirmed.pred=cbind(rep(1,length(test)),x[test,])%*%autolasso.fit$coef  
  > p.cor=cor(log1p(confirmed[test]),log.confirmed.pred)
# Model Compare
  > model.compare(n=100,ratio=0.5,p.true=5,rho=0,cor.type=1,snr=10,N.sim=100)
# Simulation Study  
  > sim.study(n=100,N.sim=100,ratio=2,p.true=5,rho=0.25,cor.type=1,snr=10,method=1)
