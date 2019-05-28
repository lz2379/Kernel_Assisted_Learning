#########################R script for Kernel Assisted Learning: Part A, source functions###########################
#########################Authors: Liangyu Zhu, Wenbin Lu, Michael Kosorok, Rui Song################################
#########################Contact: lzhu12@ncsu.edu##################################################################


####Load dependent packages#########
library(pdist)
library(optimr)

##################################################
####1. Supplementary functions for calculation####
##################################################

#g function
g=function(x){
  return(1/(1+exp(-x)))
}
#first derivative of g()
g1=function(x){
  return(exp(x)/((1+exp(x))^2))
}
#second derivative of g()
g2=function(x){
  return((exp(x)-exp(2*x))/((1+exp(x))^3))
}
#Gaussian kernel function
K=function(x){
  return(exp(-x^2/2)/(sqrt(2*pi)))
}
#first derivative of the kernel function
K1=function(x){
  return((-x*exp(-x^2/2))/(sqrt(2*pi)))
}
#second derivative of the kernel function
K2=function(x){
  return(((x^2-1)*exp(-x^2/2))/sqrt(2*pi))
}

#nonlinear baseline function
mu2=function(x){
  return(1+0.5*cos(2*pi*x))
}



##################################################
### 2. Generating data
##################################################

#Randomized trials, settings 1-4
generate_sample_random1=function(C=-10, N=500, beta=c(0,1), sd0=0.5, 
                                 mu0=function(x) 0 , 
                                 constX, constA, expo){
  p=length(beta)-1
  A=runif(N)
  X=matrix(rnorm(N*p,0,1),nrow=N)
  Q=C*(A-g(cbind(1,X) %*% beta))^2
  mu=mu0(X)+Q
  Y=rnorm(rep(1,N),mean=mu,sd=rep(sd0,N))
  hhx=constX*sd(X)*sqrt(p)*N^(-1/expo)
  hha=constA*sd(A)*N^(-1/expo)
  return(finddose(Y,X,A,hx=hhx,ha=hha,sampling=3000,
                  g=g, g1=g1,g2=g2,
                  init=diag(rep(1,p+1)) ))
}

#Randomized trial, setting 5 
generate_sample_random5=function(C=-25, N=500, beta=c(1,0.5,0.5,0),
                                 sd0=1, mbeta=c(8,4,-2,-2), 
                                 constX, constA, expo){
  p=length(beta)-1
  A=runif(N,min=0,max=2)
  X=matrix(runif(N*p,min=-1,max=1),nrow=N)
  Q=C*(A-cbind(1,X)%*% beta)^2
  mu=Q+cbind(1,X)%*% mbeta
  Y=rnorm(rep(1,N),mean=mu,sd=rep(sd0,N))
  hhx=constX*apply(as.matrix(X),2,sd)*N^(-1/expo)
  hha=constA*sd(A)*N^(-1/expo)
  return(finddose(Y,X,A=A,hx=hhx,ha=hha,sampling=3000,xlim=c(-1,1),
                  bootstrap = 100, 
                  g=function(x) x, g1=function(x) 1, g2=function(x) 0,
                  init=diag(rep(1,p+1)), seed=100 ))
}


#Observational study, settings 1-4  
generate_sample_obs1=function(C=-10,N=500,beta=c(0,1),sd0=0.5,mu0=function(x) 0 ,constX,constA,expo=5){
  p=length(beta)-1
  X=matrix(rnorm(N*p,0,1),nrow=N)
  fopt=g(cbind(1,X) %*% beta)
  A=rbeta(n=N,shape1=2*exp(cbind(1,X) %*% beta),shape2=2)
  Q=C*(A-fopt)^2
  mu=mu0(X)+Q
  Y=rnorm(rep(1,N),mean=mu,sd=rep(sd0,N))
  hhx=constX*sd(X)*sqrt(p)*N^(-1/expo)
  hha=constA*sd(A)*N^(-1/expo)
  return(finddose(Y,X,A,hx=hhx,ha=hha,sampling=3000,
                  g=g, g1=g1,g2=g2))
}



##################################################
### 3. calculate the objective function to maximize
##################################################
#'beta' are the parameters
#'A' is a vector of the observed doses
#'Y' is a vector of the observed outcomes
#'KX' is the calculated kernel function for covariates, 'logKX' is the log of KX, only one of these two are needed
#'CX' is the kernel estimated marginal density function for covariates
#'Xj1' is the sampled grids for calculating integration in the covariate space
#'ha' is the bandwidth for doses
#'N' is the number of observations
#'sampling' is number of sampled grids for calculating integrations
#'g' is the link function
conc=function(beta,Y,A,logKX=NA,KX=NA,CX,Xj1,ha,N,sampling,link,track=FALSE){
  gg=link(Xj1%*%beta)
  dA=matrix(rep(A,sampling)-rep(gg,each=N),ncol=sampling)
  if (sum(!is.na(logKX))>0){
    ee=exp(-((dA/ha)^2/2)-logKX)
  } else{
    ee=exp(-((dA/ha)^2/2))*KX
  }
  loss=-mean(t(Y)%*%ee/colSums(ee)*CX)
  if (track) {
    cat("beta:", beta, "loss:", loss,"\n")
  }
  return(loss)
}

##################################################
###4.  Dose finding.  ??????????????????????combine thse two
##################################################
###'Y' is a vector of the observed outcomes, 'X' is a matrix of covariates, 
###''A' is a vector of the observed treatments.
###'hx' is the bandwidth for X, 'ha' is the bandwidth for A
###sampling is the number of grid points used to estimate the integration
###xlim is the range of X for integration
###g is the link function, g1 is the first derivative of g, g2 is the second derivative of g
###init is the initial points
###seed is the random seed

finddose=function(Y,X,A,hx=0.1,ha=0.03,sampling=10000,xlim=c(-3,3),
                   bootstrap=0,g,g1,g2,init,seed=NA){
  N=dim(X)[1]
  p=dim(X)[2]
  X1=cbind(1,X)
  
  #sample grid points for estimating the integration
  if (!is.na(seed)){
    set.seed(seed)
  }
  Xj=matrix(runif(p*sampling,min=xlim[1],max=xlim[2]),ncol=p)
  Xj1=cbind(1,Xj)
  
  if (p>1){
    X0=as.matrix(X)%*%diag(as.numeric(1/hx))
    Xj0=as.matrix(Xj)%*% diag(as.numeric(1/hx))
  } else{
    X0=as.matrix(X)/hx
    Xj0=as.matrix(Xj)/hx
  }
  
  dX=as.matrix(pdist(X0,Xj0))
  logKX=dX^2/2
  CX=(colMeans(exp(-logKX))/sqrt(2*pi))/exp(sum(log(hx)))
  
  op=multistart(parmat=init, fn=conc,  Y=Y, A=A,  logKX=logKX, CX=CX, Xj1=Xj1, ha=ha, N=N, sampling=sampling, link=g)
  argmin=which.min(op$value)
  solution=as.numeric(op[argmin,1:(p+1)])
  
  #calculate sd using bootstrap
  if (bootstrap>0){
    bootstrap_sol=matrix(rep(NA,4*bootstrap),nrow=bootstrap)
    for (i in 1:bootstrap){
      cat("bootstrap",i,"\n")
      index=sample(N,N,replace=TRUE)
      logKXb=logKX[index,]
      CXb=(colMeans(exp(-logKXb))/sqrt(2*pi))/exp(sum(log(hx)))
      op=optim(solution,fn=conc,A=A[index],logKX=logKX[index,], CX=CXb,Xj1=Xj1,Y=Y[index],ha=ha,N=N,sampling=sampling,link=g)
      bootstrap_sol[i,]=as.numeric(op$par)
    }
    bootstrap_sd=apply(bootstrap_sol,2,sd)
  } else{
    bootstrap_sd=NA
  }
  
  ##Calculate the estimated cov using formula
  beta_X=Xj1%*%solution
  gg=g(beta_X)
  gg1=g1(beta_X)
  gg2=g2(beta_X)
  dA=matrix(rep(gg,each=N)-rep(A,sampling),ncol=sampling)
  KX=K(dX)/exp(sum(log(hx)))
  KA=K(dA/ha)/ha
  KA1=K1(dA/ha)/(ha^2)
  KA2=K2(dA/ha)/(ha^3)
  mY=matrix(rep(Y,sampling),ncol=sampling)
  A_=colMeans(mY*KX*KA)
  B_=colMeans(KX*KA)
  tA =colMeans(mY*KX*KA1)*gg1 #first derivative with respect to beta
  tB=colMeans(KX*KA1)*gg1
  ttA=colMeans(mY*KX*KA2)*(gg1^2)+ colMeans(mY*KX*KA1)*gg2 #second derivative with respect to beta
  ttB=colMeans(KX*KA2)*(gg1^2)+ colMeans(KX*KA1)*gg2 #second derivative with respect to beta
  CX=colMeans(KX)
  Dn=t(Xj1)%*%diag(as.numeric(  ((ttA/B_)-2* (tA*tB/(B_^2))-A_*ttB/(B_^2)+
                                   2*(A_*(tB^2)/(B_^3)))*CX   ))%*%(Xj1)/sampling*(xlim[2]-xlim[1])^p
  Dn_1=solve(Dn) #inverse of Dn
  
  Phi=((mY* matrix(rep(gg1*CX/B_,each=N),nrow=N)-    
        matrix(rep(gg1*CX*A_/(B_)^2,each=N),nrow=N))*KX*KA1) %*% Xj1/sampling*(xlim[2]-xlim[1])^p
  covS=cov(Phi)/N
  cov=Dn_1%*%covS%*%Dn_1
  
  return(list(beta=solution,cov=cov,N=N,hx=hx,ha=ha,bootstrap_sd=bootstrap_sd))
  
}


#Indicator function of I(x==0)
I=function(x,h=1){
  return((x==0)*h+(x!=0)*(1-h))
}
#sample from discrete dataset
sample_discrete=function(Xd,sampling,seed=100){
  pd=length(Xd[1,])
  base=2^(c(1:pd)-1)
  combo=unique(as.matrix(Xd)%*% base)
  combo=combo[order(combo)]
  sp=sample(combo,size=sampling,replace=TRUE)
  out=matrix(NA,nrow=sampling,ncol=pd)
  for (i in pd:1){
    out[,i]=floor(sp/base[i])
    sp=sp%%base[i]
  }
  return(out)
}

#find dose with data containing both continuous and discrete values
finddose2=function(Y,A,Xc,Xd,hc,hd=1,ha=0.03,sampling=2000,xlim,resample=200,
                   init, #parscale #can add this to control step size when optimizing
                   seed=100,g,g1,g2, calc_cov=TRUE){
  X=cbind(Xc,Xd)
  N=dim(X)[1]
  p=dim(X)[2]
  pc=dim(Xc)[2]
  pd=dim(Xd)[2]
  X1=cbind(1,Xc,Xd)
  
  set.seed(seed)
  #sample grid for calculating integrals
  Xjc=matrix(runif(pc*sampling,min=xlim[,1],max=xlim[,2]),ncol=pc,byrow=TRUE)
  Xjd=sample_discrete(Xd,sampling)
  Xj=cbind(Xjc,Xjd)
  Xj1=cbind(1,Xj)
  
  if (pc>1){
    Xc0=as.matrix(Xc)%*%diag(as.numeric(1/hc))
    Xjc0=as.matrix(Xjc)%*% diag(as.numeric(1/hc))
  } else{ 
    Xc0=as.matrix(Xc)/hc
    Xjc0=as.matrix(Xjc)/hc
  }
  
  dXc=as.matrix(pdist(as.matrix(Xc0),as.matrix(Xjc0)))
  
  KXc=K(dXc)/exp(sum(log(hc)))
  IXd=I(as.matrix(pdist(Xd,Xjd)))
  KX=KXc*IXd
  CX=colMeans(KX)
  
  
  op=multistart(parmat=init,method="BFGS",fn=conc,A=A,KX=KX, CX=CX,Xj1=Xj1,
                Y=Y,ha=ha,N=N,sampling=sampling,link=g, track=TRUE#, control=list(parscale=parscale)
                )
  argmin=which.min(op$value)
  solution=as.numeric(op[argmin,1:(p+1)])
  
  if (calc_cov==TRUE){
    
    ##Calculate the estimated cov
    beta_X=Xj1%*%solution
    gg=g(beta_X)
    gg1=g1(beta_X)
    gg2=g2(beta_X)
    dA=matrix(rep(gg,each=N)-rep(A,sampling),ncol=sampling)
    #KX=K(dX,d=p)/exp(sum(log(hx)))
    KA=K(dA/ha)/ha
    KA1=K1(dA/ha)/(ha^2)
    KA2=K2(dA/ha)/(ha^3)
    mY=matrix(rep(Y,sampling),ncol=sampling)
    
    A_=colMeans(mY*KX*KA)
    B_=colMeans(KX*KA)
    tA =colMeans(mY*KX*KA1)*gg1 #first derivative with respect to beta
    tB=colMeans(KX*KA1)*gg1
    ttA=colMeans(mY*KX*KA2)*(gg1^2)+ colMeans(mY*KX*KA1)*gg2 #second derivative with respect to beta
    ttB=colMeans(KX*KA2)*(gg1^2)+ colMeans(KX*KA1)*gg2
    
    s_A_=sign(A_)
    s_B_=sign(B_)
    s_tA=sign(tA)
    s_tB=sign(tB)
    s_ttA=sign(ttA)
    s_ttB=sign(ttB)
    
    l_A_=log(abs(A_))
    l_B_=log(abs(B_))
    l_tA =log(abs(tA)) #first derivative with respect to beta
    l_tB=log(abs(tB))
    l_ttA=log(abs(ttA)) #second derivative with respect to beta
    l_ttB=log(abs(ttB)) #second derivative with respect to beta
    #CX=colMeans(KX)
    sign1=s_ttA*s_B_
    sign2=s_tA*s_tB
    sign3=s_A_*s_ttB
    sign4=s_A_*s_B_
    Dn=t(Xj1)%*%diag(as.numeric(  (  sign1*exp(l_ttA-l_B_)-  sign2*2* exp(l_tA+l_tB-2*l_B_)  -  sign3*exp(l_A_+l_ttB-2*(l_B_))  +
                                       sign4*2*exp(l_A_+2*l_tB-3*l_B_)                                                   )*CX   ))%*%(Xj1)/sampling*(exp(sum(log(xlim[,2]-xlim[,1]))))
    Dn_1=solve(Dn)
    
    s_C=sign(CX)
    l_C=log(abs(CX))
    sign5=s_C*s_B_
    sign6=s_C*s_A_
    Phi=((mY* matrix(rep(gg1*sign5*exp(l_C-l_B_),each=N),nrow=N)-    matrix(rep(gg1*sign6*exp(l_C+l_A_-2*l_B_),each=N),nrow=N))*KX*KA1) %*% Xj1/sampling*(exp(sum(log(xlim[,2]-xlim[,1]))))
    covS=cov(Phi)/N
    cov=Dn_1%*%covS%*%Dn_1
  } else{
    cov=NA
  }
  return(list(beta=solution,cov=cov,N=N,hc=hc,ha=ha,allresult=op))
}

##################################################
#### 5. cross validation for finding the h
##################################################
#for continuous data
cross_val1=function(par,X,A,Y, expo, fold=5,track=FALSE){
  NN=dim(X)[1]
  p=dim(X)[2]
  testsize=floor(NN/fold)
  trainsize=NN-testsize
  constX=par[1:p]
  constA=par[(p+1)]
  
  if (track) {
    cat(par,"\n")
  }
  mse=0
  for (i in 1:fold){
    set.seed(i)
    train_index=sample(NN,size=trainsize)
    
    trainX=X[train_index,]
    trainY=Y[train_index]
    trainA=A[train_index]
    
    testX=X[-train_index,]
    testY=Y[-train_index]
    testA=A[-train_index]
    
    testX1=cbind(1,testX)
    trainX1=cbind(1,trainX)
    
    hx=constX*apply(as.matrix(trainX),2,sd)*trainsize^(-1/expo)
    ha=constA*sd(trainA)*trainsize^(-1/expo)
    if (p>1){
      trainX0=as.matrix(trainX)%*%diag(as.numeric(1/hx)) 
      testX0=as.matrix(testX)%*%diag(as.numeric(1/hx))
    } else {
      trainX0=as.matrix(trainX)/hx 
      testX0=as.matrix(testX)/hx
    }
    dX=as.matrix(pdist(as.matrix(trainX0),as.matrix(testX0)))
    KX=K(dX)/exp(sum(log(hx)))
    KA=K(as.matrix(pdist(as.matrix(trainA),as.matrix(testA)))/ha)
    predictY=t(KX*KA)%*%trainY/trainsize/colMeans(KX*KA)
    mse=mse+mean((predictY-testY)^2)
  }
  return(sqrt(mse/fold/testsize))
}

# for dataset with continuous variables and discrete variables
cross_val2=function(par,X,Xc=NA,Xd,A,Y, expo, fold=5){
  NN=dim(X)[1]
  pc=dim(Xc)[2]
  
  testsize=floor(NN/fold)
  trainsize=NN-testsize
  constX=par[1:pc]
  constA=par[pc+1]
  
  cat(par,"\n")
  
  mse=0
  for (i in 1:fold){
    set.seed(i)
    train_index=sample(NN,size=trainsize)
    
    trainX=X[train_index,]
    trainY=Y[train_index]
    trainA=A[train_index]
    trainXc=Xc[train_index,]
    trainXd=Xd[train_index,]
    
    testX=X[-train_index,]
    testY=Y[-train_index]
    testA=A[-train_index]
    testXc=Xc[-train_index,]
    testXd=Xd[-train_index,]
    
    testX1=cbind(1,testX)
    trainX1=cbind(1,trainX)
    
    hc=constX*apply(as.matrix(trainXc),2,sd)*trainsize^(-1/expo)
    ha=constA*sd(trainA)*trainsize^(-1/expo)
    if (pc>1){
      trainXc0=as.matrix(trainXc)%*%diag(as.numeric(1/hc)) 
      testXc0=as.matrix(testXc)%*%diag(as.numeric(1/hc))
    } else {
      trainXc0=as.matrix(trainXc)/hc 
      testXc0=as.matrix(testXc)/hc
    }
    
    dXc=as.matrix(pdist(as.matrix(trainXc0),as.matrix(testXc0)))
    KXc=K(dXc,d=pc)/exp(sum(log(hc)))
    IXd=I(as.matrix(pdist(trainXd,testXd)))
    KX=KXc*IXd
    KA=K(as.matrix(pdist(as.matrix(trainA),as.matrix(testA)))/ha)
    predictY=t(KX*KA)%*%trainY/trainsize/colMeans(KX*KA)
    mse=mse+sum((testY-predictY)^2)
  }
  print(sqrt(mse/rep/testsize))
  return(sqrt(mse/rep/testsize))
}


##################################################
#### 6. estimate the value function of the suggested dose on a test dataset
##################################################


testvalue=function(betahat,outdose=NA, testY,testA,testXc,testXd,testsize,hc,ha){
  #estimate value function using kernel with test dataset
  if (sum(!is.na(outdose))>0){
    predictA=outdose
  } else{
    testX1=cbind(1,testXc,testXd)
    predictA=g(as.matrix(testX1)%*%betahat)
  }
  Aj=matrix(rep(testA,each=testsize),nrow=testsize)
  
  predictAi=matrix(rep(predictA,testsize),nrow=testsize)
  Ka=as.matrix(K((Aj-predictAi)/ha))/ha
  if (length(hc)>1){
    dXc=as.matrix(dist(as.matrix(testXc)%*%diag(as.numeric(1/hc))))
  } else {
    dXc=as.matrix(dist(as.matrix(testXc)/hc))
  }
  Kxc=K(dXc)/exp(sum(log(hc)))
  Ixd=I(as.matrix(dist(testXd)))
  Kx=Kxc*Ixd
  
  Yj=matrix(rep(testY,each=testsize),nrow=testsize)
  return(mean(rowMeans(Yj*Ka*Kx)/rowMeans(Ka*Kx)))
}

