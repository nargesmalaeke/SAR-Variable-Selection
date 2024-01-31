library(Metrics)
library(MASS)
library(stats4)
library(sn)
library(geoR)
library(nlme)
library(mvtnorm)
library(pracma)
library(Matrix)
library(statmod)
library(HyperbolicDist)
library(spatialprobit)
library(spatstat)
library(coda)
gibbsSAR=function(p,niter,burn){
  n=100
  #### Data Generating
  X=matrix(rep(1,p),nrow=n,ncol=p)
  for(j in 2:p){
    X[,j]=rnorm(n)}
  beta1=matrix(rnorm(5,0,5),ncol=1)
  beta2=matrix(rnorm(5,0,1),ncol =1)
  zero=matrix(rep(0,p-10),ncol=1)
  betat=rbind(beta1,beta2,zero)
  sigma2t=1
  rhot=rnorm(1,0,3)
  Omegat=sigma2t*diag(n)
  epsilon0=matrix(rmvnorm(1,rep(0,n),Omegat))
  ## generate random points within a unit square
  pp <- runifpoint(n)
  ## row-stochastic matrix W based on the five-nearest neighborhood
  W=kNearestNeighbors(pp$x,pp$y,k=5)
  ## generate y
  y0=(exp(-rhot*W))%*%(X%*%betat+epsilon0)
  ## Standardization of X and Y
  for(j in 2:p){
    for(i in 1:n){X[i,j]=(X[i,j]-mean(X[,j]))/sd(X[,j])}}
  y=(y0-mean(y0))/sd(y0)
  ## Metroplis Algorithm function for Generate Rho
  Omg=diag(n)
  drho=function(rho){
    e1=exp(rho*W)%*%y-X%*%beta
prho=dnorm(rho,0,10)
ans=exp(-(t(e1)%*%(solve(Omg))%*%e1)/2)*prho
return(as.numeric(ans))}
  metropolisRho=function(niterrho,burnin,startval,si,thin){
    rhochain = matrix(1,nrow=niterrho)
    rhochain[1,]=startval
    for(i in 2:niterrho){
      currentrho=rhochain[i-1,]
      proposedrho=rnorm(1,currentrho,si) #rho*
      A=(drho(proposedrho))/(drho(currentrho))
      R=min(1,A)
      u=runif(1)
      if(u<R){
        rhochain[i,] = proposedrho
      } else {rhochain[i,] = currentrho
      } }
    train=niterrho-burnin
    trho=as.matrix(rhochain[(burnin+1):niterrho,])
    Mrho=matrix(rep(1,train/thin),ncol=1)
    for(i in 1:(train/thin)){
      Mrho[i,]=trho[(thin*i),]}
    rho=mean(Mrho)
    return(rho)
  }
  #### estimate with train and burn-in
  train=niter-burn
  trainfun=function(burn){
    trainbeta=mat[(burn+1):niter,1:p]
    trainrho=mat[(burn+1):niter,p+1]
    trainrho=as.matrix(trainrho,ncol=1)
    trainmat=cbind(trainbeta,trainrho)
    return(trainmat)
  }
  ### input initial values for parameters
  beta0=matrix(runif(p,-5,5),ncol=1)
  for(j in 1:p){
    if(beta0[j,]==0)
      beta0[j,]=runif(1,1,2)}
  rho0=runif(1,-2,2)
  if(rho0==0){
    rho0=runif(1,1,3)
  }
  ## Parameters for NG Prior
  psi0=matrix(runif(p,1,10),ncol=1)
  lambda20=runif(1,1,10)
  ## Parameters for DL Prior
  phi0=matrix(runif(p,1,10),ncol=1)
  varphi0=matrix(runif(p,1,10),ncol=1)
  tu0=runif(1,1,10)
  ### RMSE for Point Estimates
  result=matrix(rep(0,6),nrow=3,ncol=2)
  ### MCMC Algorithm / None / Without Shrinkage
  beta=beta0
  rho=rho0
  mat <- matrix(0, ncol = p+1 , nrow = niter)
  mat[1,]=cbind(t(beta),rho)
  startnone=Sys.time()
  for(i in 2:niter){
    ### beta
    Sigma_=1000*diag(p)
    Sigmabar=solve(solve(Sigma_)+t(X)%*%solve(Omg)%*%X)
    betabar=Sigmabar%*%(t(X)%*%solve(Omg)%*%exp(rho*W)%*%y)
    beta=matrix(rmvnorm(1,betabar,Sigmabar,checkSymmetry = FALSE))
    for(j in 1:p){
      if(0 < beta[j,] && beta[j,] < 0.0001){
        beta[j,]=0.001}
      if(-0.0001 < beta[j,] && beta[j,] < 0){
        beta[j,]=-0.001}}
    ### rho
    rho1=metropolisRho(200,100,5,30,10)
    rho=metropolisRho(1000,500,rho1,30,10)
    mat[i,]=cbind(t(beta),rho)
  }
  endnone=Sys.time()
  TimeNone=endnone-startnone
  trainmat=trainfun(burn)
  estimate=function(thin){
    trainbeta=trainmat[,1:p]
    betahat=matrix(1,nrow=p,ncol=1)
    nr=round(train/thin)
    betathin=matrix(0,nrow=nr,ncol=p)
    for(j in 1:p){
      for(i in 2:nr){
        betathin[i-1,j]=trainbeta[(thin*(i-1)),j]}
      betathin[nr,j]=trainbeta[train,j]
      betahat[j,]=mean(betathin[,j])}
    trainrho=as.matrix(trainmat[,p+1],ncol=1)
    rhothin=matrix(0,nrow=nr,ncol=1)
    for(i in 2:nr){
      rhothin[i,]=trainrho[(thin*(i-1)),]}
    rhothin[nr,]=trainrho[train,]
    rhohat=mean(rhothin)
    est=rbind(betahat,rhohat)
    return(est)}
  beta=as.matrix(estimate(thin)[1:p,],ncol=1)
  rho=as.matrix(estimate(thin)[p+1,])
  betaa=betathin[,p]
  result[1,1]=rmse(rhot,rho)
  result[1,2]=rmse(betat,beta)
  ### MCMC Algorithm / NG Prior
  beta=beta0
  rho=rho0
  lambda2=lambda20
  psi=psi0
  ### hyperparameters
  theta=0.1
  d0=0.01
  d1=0.01
  ### input initial values in matrix
  mat <- matrix(0,ncol =2*p+2 , nrow = niter)
  mat[1,]=cbind(t(beta),rho,lambda2,t(psi))
  startNG=Sys.time()
  for(i in 2:niter){
    ### lambda2
    S1=0
    for(j in 1:p){
      S1=S1+psi[j,]}
    lambda2=rgamma(1,d0+(theta*p),d1+(2^(-1)*theta*S1))
    ### psi
    for(j in 1:p){
      psi[j,]=rgig(1,Theta=c(theta-1/2,beta[j,]^2,theta*lambda2))
    }
    ### beta
    Sigma_=Diagonal(p,psi)
    Sigmabar=solve(solve(Sigma_)+t(X)%*%solve(Omg)%*%X)
    betabar=Sigmabar%*%(t(X)%*%solve(Omg)%*%exp(rho*W)%*%y)
    beta=matrix(rmvnorm(1,betabar,Sigmabar,checkSymmetry = FALSE))
    for(j in 1:p){
      if(0 < beta[j,] && beta[j,] < 0.0001){
        beta[j,]=0.001}
      
      if(-0.0001 < beta[j,] && beta[j,] < 0){
        beta[j,]=-0.001}}
    ### rho
    rho1=metropolisRho(200,100,5,30,10)
    rho=metropolisRho(1000,500,rho1,30,10)
    mat[i,]=cbind(t(beta),rho,lambda2,t(psi))
  }
  endNG=Sys.time()
  TimeNG=endNG-startNG
  trainmat=trainfun(burn)
  beta=as.matrix(estimate(thin)[1:p,],ncol=1)
  rho=as.matrix(estimate(thin)[p+1,])
  result[2,1]=rmse(rhot,rho)
  result[2,2]=rmse(betat,beta)
  #### MCMC Algorithm / DL
  beta=beta0
  rho=rho0
  tu=tu0
  phi=phi0
  varphi=varphi0
  ## hyperparameter
  a=1/p
  mat <- matrix(0,ncol =3*p+2 , nrow = niter)
  mat[1,]=cbind(t(beta),rho,tu,t(phi),t(varphi))
  ### gibbs sampling for DL
  startDL=Sys.time()
  for(i in 2:niter){
    ### varphi
    mu=matrix(rep(0,p),ncol=1)
    for (j in 1:p) {
      mu[j,]=phi[j,]*tu/abs(beta[j,])
      Evarphi=rgig(1,Theta = c(-1/2,1,mu[j,]^(-2)))
      varphi[j,]=1/Evarphi }
    ### tu
    S1=0
    for(j in 1:p){
      S1=S1+abs(beta[j,])/phi[j,]}
    tu=rgig(1,Theta = c(1-p,2*S1,1))
    ### phi
    t=matrix(rep(0,p),ncol=1)
    S2=0
    for(j in 1:p){
      t[j,]=rgig(1,Theta = c(a-1,2*abs(beta[j,]),1))
      S2=S2+t[j,]}
    for(j in 1:p){
      phi[j,]=t[j,]/S2}
    ### beta
    d=matrix(rep(0,p),ncol=1)
    for (j in 1:p) {
      d[j,]=varphi[j,]*phi[j,]^2*tu^2}
    Sigma_=Diagonal(p,d)
    Sigmabar=solve(solve(Sigma_)+t(X)%*%solve(Omg)%*%X)
    betabar=Sigmabar%*%(t(X)%*%solve(Omg)%*%exp(rho*W)%*%y)
    beta=matrix(rmvnorm(1,betabar,Sigmabar,checkSymmetry = FALSE))
    for(j in 1:p){
      if(0 < beta[j,] && beta[j,] < 0.0001){
        beta[j,]=0.001}
      if(-0.0001 < beta[j,] && beta[j,] < 0){
        beta[j,]=-0.001}}
    ### rho
    rho1=metropolisRho(200,100,5,30,10)
    rho=metropolisRho(1000,500,rho1,30,10)
    mat[i,]=cbind(t(beta),rho,tu,t(phi),t(varphi))
  }
  endDL=Sys.time()
  TimeDL=endDL-startDL
  beta=as.matrix(estimate(thin)[1:p,],ncol=1)
  rho=as.matrix(estimate(thin)[p+1,])
  result[3,1]=rmse(rhot,rho)
  result[3,2]=rmse(betat,beta)
  Time=rbind(TimeNone,TimeNG,TimeDL)
  Method=matrix(c("None","NG","DL"),ncol=1)
  Answer=cbind(Method,result,Time)
  return(Answer)
}
EmpRMSE=function(iter,p,niter,burn){
  Time=matrix(0,nrow=3,ncol=iter)
  EmpRMSE=matrix(0,nrow=3,ncol=iter)
  for (k in 1:iter) {
    EmpRMSE[,k]=matrix(as.numeric(gibbsSAR(p,niter,burn)[,2]),ncol=1)
    Time[,k]=matrix(as.numeric(gibbsSAR(p,niter,burn)[,3]),ncol=1)
  }
  EstTime=as.matrix(rep(0,3))
  EstRMSE=as.matrix(rep(0,3))
  for(i in 1:3){
    EstRMSE[i,]=mean(EmpRMSE[i,])
    EstTime[i,]=mean(Time[i,])}
  Method=matrix(c("None","NG","DL"),ncol=1)
  Answer=cbind(Method,EstRMSE,EstTime)
  return(Answer)}