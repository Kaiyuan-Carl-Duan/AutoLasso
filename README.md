# AutoLasso

Linear regression is one of the most widely used statistical methods and with the ease of modern data collection, often many predictors are associated with a response variable. However, most of the available variable selection methodologies (e.g., elastic net and its special case lasso) involve careful selection of the tuning parameter. This can cause numerous problems in practice and may lead to biased estimates. We present an automatic variable selection method which is simple and compares favorably with many of the currently available variable selection methods for linear models. Numerical illustrations based on simulated and real data are presented based on output from an user friendly R package `AutoLasso` that is made available for easy implementation.

Required R packages: `glmnet` and `MASS`

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
