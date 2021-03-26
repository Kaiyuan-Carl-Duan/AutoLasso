##############################################################
#Simulation Study: Variable Selection via Best Subset Method #
#Copyright (2019): Sujit K. Ghosh, NC State University       #
##############################################################
#Last modified on: 06/27/2020 (by Carl Duan)
############################################

#########################################
#Variable selection by automatic LASSO  #
#########################################
lasso1=function(x,y){
  library(glmnet)
  fit1=cv.glmnet(x,y)
  beta.est=as.numeric(coef(fit1,s="lambda.1se"))
  #beta.est=as.numeric(coef(fit1,s="lambda.min"))
  x=cbind(rep(1,nrow(x)),x)
  sigma2.est=mean((y-as.vector(x%*%beta.est))^2)
  return(list(coef=beta.est,sigma2=sigma2.est))
}

lasso2=function(x,y){
  library(glmnet)
  fit1=cv.glmnet(x,y,relax=TRUE)
  beta.est=as.numeric(coef(fit1,s="lambda.1se",gamma="gamma.1se"))
  #beta.est=as.numeric(coef(fit1,s="lambda.min"))
  x=cbind(rep(1,nrow(x)),x)
  sigma2.est=mean((y-as.vector(x%*%beta.est))^2)
  return(list(coef=beta.est,sigma2=sigma2.est))
}

auto.lasso=function(x,y,ratio=0.5){
  risk=c(); p=ncol(x); n=nrow(x)
  ratio=min(ratio,0.95)
  m.max=ceiling(ratio*qr(x)$rank)

  beta.lasso=lasso1(x=x,y=y)$coef

  for(m in 1:m.max){
    #find m largest beta's in absolute values
    active.indx=order(-abs(beta.lasso[-1]))[1:m]
    #obtain LS estimates with active m variables;
    fit.m=lm(y~x[,active.indx])
    risk[m]=(m+1)*((summary(fit.m)$sigma)^2)
  }

  #select m.opt
  m.opt=which.min(risk)

  #find m from the minimum of risk
  active.indx.opt=order(-abs(beta.lasso[-1]))[1:m.opt]
  #obtain LS estimates with active m variables;
  fit.opt=lm(y~x[,active.indx.opt]); beta.est=rep(0,p)
  beta.est[active.indx.opt]=fit.opt$coef[-1]
  beta.est=c(fit.opt$coef[1],beta.est)
  #MLE of sigma2:
  sigma2.est=mean(fit.opt$residuals^2)
  #UMVUE of sigma2:
  #sigma2.est=sum(fit.opt$residuals^2)/(n-m.opt-1)

  #return estimate beta and estimate sigma
  return(list(coef=beta.est,sigma2=sigma2.est,best.fit=fit.opt))
}

auto.lasso.2=function(x,y,ratio=0.5){
  risk=c(); p=ncol(x); n=nrow(x)
  ratio=min(ratio,0.95)
  m.max=ceiling(ratio*qr(x)$rank)

  beta.lasso=lasso2(x=x,y=y)$coef

  for(m in 1:m.max){
    #find m largest beta's in absolute values
    active.indx=order(-abs(beta.lasso[-1]))[1:m]
    #obtain LS estimates with active m variables;
    fit.m=lm(y~x[,active.indx])
    risk[m]=(m+1)*((summary(fit.m)$sigma)^2)
  }

  #select m.opt
  m.opt=which.min(risk)

  #find m from the minimum of risk
  active.indx.opt=order(-abs(beta.lasso[-1]))[1:m.opt]
  #obtain LS estimates with active m variables;
  fit.opt=lm(y~x[,active.indx.opt]); beta.est=rep(0,p)
  beta.est[active.indx.opt]=fit.opt$coef[-1]
  beta.est=c(fit.opt$coef[1],beta.est)
  #MLE of sigma2:
  sigma2.est=mean(fit.opt$residuals^2)
  #UMVUE of sigma2:
  #sigma2.est=sum(fit.opt$residuals^2)/(n-m.opt-1)

  #return estimate beta and estimate sigma
  return(list(coef=beta.est,sigma2=sigma2.est,best.fit=fit.opt))
}

###########################
#Perform a data generation#
###########################
gen.data=function(n=100,ratio=0.5,p.true=5,rho=0,cor.type=1,snr=10){
  library(MASS)
  p=floor(ratio*n); p.true=min(max(p.true,3),p)

  if(cor.type==1){Sigma0=toeplitz(rho^c(0:(p-1)))}#auto regressive
  else{Sigma0=(1-rho)*diag(p)+rho*matrix(1,p,p)} #compound symmetry

  #Set active and inactive beta values:
  #nzero.indx=1:p.true
  nzero.indx=sample(p,p.true)
  beta=rep(0,p)
  beta[nzero.indx]=c(-1,-1,rep(1,p.true-2))
  beta0=0.5; beta.true=c(beta0,beta)
  sigma.true=sqrt(beta0^2+as.vector(t(beta)%*%Sigma0%*%beta))/snr

  #set.seed(27695)
  x=mvrnorm(n,mu=rep(0,p),Sigma=Sigma0); X=cbind(rep(1,n),x)
  y=mvrnorm(mu=as.vector(X%*%beta.true),Sigma=(sigma.true^2)*diag(n))
  #collect the output
  output=list(y=y,x=x,beta.true=beta.true,sigma.true=sigma.true)
  return(output)
}

#try the code:
#mydata=gen.data(n=100,ratio=1.5,rho=0.25)

#test the code
#Case1: easy
#mydata=gen.data(n=100,ratio=0.5,rho=0.25,snr=10,cor.type=2)
#Case2: moderate
#mydata=gen.data(n=100,ratio=1,rho=0.5,snr=5,cor.type=2)
#Case2: tough
#mydata=gen.data(n=100,ratio=2,rho=0.85,snr=1,cor.type=2)

#out1=lasso1(x=mydata$x,y=mydata$y)
#out2=auto.lasso(x=mydata$x,y=mydata$y)
#rbind(mydata$beta.true,out1$coef,out2$coef)
#c(mydata$sigma.true^2,out1$sigma2,out2$sigma2)

###################
#Compare the Model#
###################
model.compare=function(n=100,ratio=0.5,p.true=5,rho=0,cor.type=1,snr=10,N.sim=100,method=1){

  #provide the data generation
  mydata=gen.data(n=n,ratio=ratio,p.true=p.true,rho=rho,cor.type=cor.type,snr=snr)

  #collect the data from data generation function
  x=mydata$x;sigma.true=mydata$sigma.true;beta.true=mydata$beta.true
  p=dim(x)[2];p.true=sum(mydata$beta.true==1)+sum(mydata$beta.true==-1)
  beta0=beta.true[1];beta=beta.true[-1];true.indx=1-(beta==0)

  #begin the stimulation program
  start.time=proc.time(); start.date=date()

  p.est.auto.lasso=c(); bias.auto.lasso=c(); bias.sigma2.auto.lasso=c()
  selected.indx.auto.lasso=c(); miss.class.auto.lasso=c()

  p.est.lasso1=c(); bias.beta.lasso1=c(); bias.sigma2.lasso1=c()
  selected.indx.lasso1=c(); miss.class.lasso1=c()

  for(i in 1:N.sim){

    #generate response variables:
    y=beta0+as.vector(x%*%beta)+rnorm(n,0,sd=sigma.true)

    #Estimate coef's using the default lasso method:

    fit.auto.lasso=auto.lasso(x=x,y=y)
    beta.hat.auto.lasso=fit.auto.lasso$coef
    sigma2.hat.auto.lasso=fit.auto.lasso$sigma2

    if(method==1){fit.lasso1=lasso1(x=x,y=y)}
    else{fit.lasso1=auto.lasso.2(x=x,y=y)}
    beta.hat.lasso1=fit.lasso1$coef
    sigma2.hat.lasso1=fit.lasso1$sigma2

    #Summary of variables selected by the lasso method:
    TPR.true=sum(true.indx[which(true.indx!=0)])/p.true
    FPR.true=sum(true.indx[which(true.indx==0)])/(p-p.true)
    FDR.true=sum(true.indx[which(true.indx==0)])/max(p.true,1)

    selected.auto.lasso=as.numeric(1-(beta.hat.auto.lasso[-1]==0))
    selected.indx.auto.lasso=rbind(selected.indx.auto.lasso,selected.auto.lasso)

    p.est.auto.lasso=c(p.est.auto.lasso,sum(selected.auto.lasso))
    bias.auto.lasso=rbind(bias.auto.lasso,beta.hat.auto.lasso-beta.true)
    bias.sigma2.auto.lasso=c(bias.sigma2.auto.lasso,(sigma2.hat.auto.lasso/sigma.true^2)-1)
    miss.class.auto.lasso=c(miss.class.auto.lasso,mean(abs(selected.auto.lasso-true.indx)))

    selected.lasso1=as.numeric(1-(beta.hat.lasso1[-1]==0))
    selected.indx.lasso1=rbind(selected.indx.lasso1,selected.lasso1)

    p.est.lasso1=c(p.est.lasso1,sum(selected.lasso1))
    miss.class.lasso1=c(miss.class.lasso1,mean(abs(selected.lasso1-true.indx)))
    bias.beta.lasso1=rbind(bias.beta.lasso1,beta.hat.lasso1-beta.true)
    bias.sigma2.lasso1=c(bias.sigma2.lasso1,(sigma2.hat.lasso1/sigma.true^2)-1)
  }

  TPR.auto.lasso=sum(colMeans(selected.indx.auto.lasso)[which(true.indx!=0)])/p.true
  FPR.auto.lasso=sum(colMeans(selected.indx.auto.lasso)[which(true.indx==0)])/(p-p.true)
  FDR.auto.lasso=sum(colMeans(selected.indx.auto.lasso)[which(true.indx==0)])/max(sum(colMeans(selected.indx.auto.lasso)),1)


  TPR.lasso1=sum(colMeans(selected.indx.lasso1)[which(true.indx!=0)])/p.true
  FPR.lasso1=sum(colMeans(selected.indx.lasso1)[which(true.indx==0)])/(p-p.true)
  FDR.lasso1=sum(colMeans(selected.indx.lasso1)[which(true.indx==0)])/max(sum(colMeans(selected.indx.lasso1)),1)

  #####end of simulation study
  time.elapsed=proc.time()-start.time
  end.date=date()

  #######################
  #Graphical summaries  #
  #######################
  par(mfrow=c(2,4),cex=0.4)
  bar.col=rep("grey",p+1)
  bar.col[which(beta.true!=0)]=rep("green",p.true+1)

  #graphical summary for bestmodel
  colnames(bias.auto.lasso)=0:p
  boxplot(bias.auto.lasso, col=bar.col); abline(h=0,col="red")
  title(paste("Biases of All Coefficients auto.lasso(n=",n,", p=",p,")"))
  barplot(colMeans(selected.indx.auto.lasso),col=bar.col, names.arg=1:p)
  title(paste("Propn of correct selection auto.lasso(p.true=",p.true,", N.sim=",N.sim,")"))
  boxplot(bias.auto.lasso[,which(beta.true!=0)], col="green"); abline(h=0,col="red")
  title("Biases of active variables (auto.lasso)")
  #boxplot(bias.auto.lasso[,which(beta.true==0)], col="grey"); abline(h=0,col="red")
  #title("Biases of zero coefficients (auto.lasso)")
  boxplot(bias.sigma2.auto.lasso,col="lightblue"); abline(h=0,col="red")
  title("Rel. Bias of sigma^2 (auto.lasso)")

  #graphical summary for lasso1
  colnames(bias.beta.lasso1)=0:p
  boxplot(bias.beta.lasso1, col=bar.col); abline(h=0,col="red")
  title(paste("Biases of All Coefficients lasso1 (n=",n,", p=",p,")"))
  barplot(colMeans(selected.indx.lasso1),col=bar.col, names.arg=1:p)
  title(paste("Propn of correct selection lasso1(p.true=",p.true,", N.sim=",N.sim,")"))
  boxplot(bias.beta.lasso1[,which(beta.true!=0)], col="green"); abline(h=0,col="red")
  title("Biases of active variables (lasso1)")
  #boxplot(bias.beta.lasso1[,which(beta.true==0)], col="grey"); abline(h=0,col="red")
  #title("Biases of zero coefficients (lasso1)")
  boxplot(bias.sigma2.lasso1,col="lightblue"); abline(h=0,col="red")
  title("Rel. Bias of sigma^2 (lasso1)")

  #######################
  #Numerical summaries  #
  #######################
  #Set up a function to find mode
  getmode=function(v) {
    uniqv=unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }

  #percentage of times each variables selected:
  colnames(selected.indx.auto.lasso)=1:p ;ordr=order(-true.indx)
  est.propn=round(colMeans(selected.indx.auto.lasso),3)[ordr]
  est.propn1=round(colMeans(selected.indx.lasso1),3)[ordr]
  selected.propn=rbind(true.indx[ordr],est.propn,est.propn1)
  rownames(selected.propn)=c("true","auto.lasso","lasso1")
  cat("Proportion of correct selection:",fill=T)
  print(selected.propn)
  cat("--------------------------",fill=T)

  #TPR&FPR of different methods:
  TPR.disp=rbind(TPR.true,TPR.auto.lasso,TPR.lasso1)
  FPR.disp=rbind(FPR.true,FPR.auto.lasso,FPR.lasso1)
  FDR.disp=rbind(FDR.true,FDR.auto.lasso,FDR.lasso1)
  TFPR.disp=round(cbind(TPR.disp,FPR.disp,FDR.disp),3)
  rownames(TFPR.disp)=c("true","auto.lasso","lasso1")
  colnames(TFPR.disp)=c("TPR","FPR","FDR")
  cat("TPR & FPR for different methods:",fill=T)
  print(TFPR.disp)
  cat("--------------------------",fill=T)

  #Distibution of number of selected variables:
  cat("Sampling distribution of number of selected variables:",fill=T)
  cat("auto.lasso method",fill=T)
  print(round(table(p.est.auto.lasso)/N.sim,3))
  cat("lasso1 method",fill=T)
  print(round(table(p.est.lasso1)/N.sim,3))
  mean.mode=rbind(cbind(mean(p.est.auto.lasso),getmode(p.est.auto.lasso)),cbind(mean(p.est.lasso1),getmode(p.est.lasso1)))
  rownames(mean.mode)=c("auto.lasso","lasso1")
  colnames(mean.mode)=c("mean","mode")
  print(mean.mode)
  cat("--------------------------",fill=T)

  #Comparasion between methods of number of active variables found
  cat("active variables in different methods")
  act.num=sum(beta.true!=0)-1
  act.var.auto.lasso.feq=round(table(p.est.auto.lasso)/N.sim, 3)[which((as.data.frame(round(table(p.est.auto.lasso)/N.sim, 3)))$p.est.auto.lasso==act.num)]
  act.var.lasso1.feq=round(table(p.est.lasso1)/N.sim, 3)[which((as.data.frame(round(table(p.est.lasso1)/N.sim, 3)))$p.est.lasso1==act.num)]
  act.var.summary=cbind(table(act.num), act.var.auto.lasso.feq, act.var.lasso1.feq)
  print(act.var.summary)
  cat("--------------------------",fill=T)

  #Biases of active variables:
  cat("bias of active cofficients:",fill=T)
  cat("auto.lasso method",fill=T)
  print(round(apply(bias.auto.lasso[,which(beta.true!=0)],2,summary),4))
  cat("lasso1 method",fill=T)
  print(round(apply(bias.beta.lasso1[,which(beta.true!=0)],2,summary),4))
  cat("--------------------------",fill=T)

  #Relative Bias of sigma2:
  cat("relative bias of sigma^2:",fill=T)
  bias.sigma2.summary=rbind(round(summary(bias.sigma2.auto.lasso),4), round(summary(bias.sigma2.lasso1),4))
  rownames(bias.sigma2.summary)=c("auto.lasso", "lasso1")
  print(bias.sigma2.summary)
  cat("--------------------------",fill=T)

  #summary of overall miss classification:
  cat("Overall missclassification rate:",fill=T)
  miss.class.summary=rbind(round(summary(miss.class.auto.lasso),4), round(summary(miss.class.lasso1),4))
  rownames(miss.class.summary)=c("auto.lasso","lasso1")
  print(miss.class.summary)
  cat("--------------------------",fill=T)

  #Real time used (subtract start time from completed):
  cat(c("started at: ",start.date),fill=T)
  cat(c("completed at: ",end.date),fill=T)
  cat("--------------------------",fill=T)
}

##################
#Simulation Study#
##################
sim.study=function(n=100,N.sim=100,ratio=2,p.true=5,rho=0.25,cor.type=1,snr=10,method=1){
  #Perform a data generation
  mydata=gen.data(n=n,ratio=ratio,p.true=p.true,rho=rho,snr=snr,cor.type=cor.type)

  #Perform a simulation process
  model.compare(n=n,ratio=ratio,rho=rho,cor.type=cor.type,snr=snr,N.sim=N.sim,
                method=method)
}
