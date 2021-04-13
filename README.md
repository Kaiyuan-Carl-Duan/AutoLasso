# AutoLasso: A Fast and Automatic Variable Selection Method for Multiple Linear Regression

Linear regression is one of the most widely used statistical methods and with the ease of modern data collection, often many predictors are associated with a response variable. However, most of the available variable selection methodologies (e.g., elastic net and its special case lasso) involve careful selection of the tuning parameter. This can cause numerous problems in practice and may lead to biased estimates. We present an automatic variable selection method using anew risk criteria which is simple and compares favorably with many of the currently available variable selection methods for linear models. Numerical illustrations based on simulated and real data are presented based on output from an user friendly R package `AutoLasso` that is made available for easy implementation.

Methodology and Performance of Simulations: https://www.causeweb.org/usproc/eusrc/2020/virtual-posters/6

Required R packages: `glmnet` and `MASS`

Author: Sujit Ghosh, Kaiyuan(Carl) Duan, and Guangjie Yu

## Application

To study the application of Automated LASSO approach, we compare with linear regression `lm` and Classic LASSO method. Data was taken from the R package `COVID19` available to the public for educational and scientific use. Due to the outbreak of corona-virus in the United States, exploding in March and stable in the beginning of June, we selected the data based on the date between March 1st and June 10th, 2020. The dataset presented in this paper is four predictors (selected from a preliminary analysis of a set of five variables) modeled on a selected response variable. There were 5202 observations. Conclusions about the data were not desired, the data was used only to observe how linear regression `lm`, Automated LASSO, and classic LASSO models perform on highly correlated data sets. All computation was done with R software using the package `glmnet` and `AutoLasso`.

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
