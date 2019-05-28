##################################################
##### simulation Examples ##############
##################################################

##################################################
###setting 1
##################################################
source('KAl-functions.r')
expo=4.5
constX=1.25
constA=1.75

for (samplesize in c(400,800)){
  set.seed(2000)
  numberofsets=500
  result_01=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  cov_01=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  results_01=list()
  for (i in 1:numberofsets){
    print(i)
    result=generate_sample_random1(beta=c(0,0.5),N=samplesize,constX=constX,constA=constA,expo=expo)
    result_01[i,]=result$beta
    cov_01[i,1]=result$cov[1,1]
    cov_01[i,2]=result$cov[2,2]
    results_01[[i]]=result
    write.csv(result_01, paste("beta_01", samplesize, ".csv"))
    write.csv(cov_01,paste("cov_01", samplesize, ".csv"))
    save.image(paste("result_01", samplesize, ".RData"))
  }
  
  colMeans(result_01) #mean of the estimated parameter
  apply(result_01,2,sd) # sd of the estimated parameters
  colMeans(sqrt(cov_01))  # mean of estimated sd
  #Confidence Intervals
  lowerCI=result_01-sqrt(cov_01)*1.96
  upperCI=result_01+sqrt(cov_01)*1.96
  #coverage of confidence intervals
  mean(lowerCI[,1]*upperCI[,1]<0)
  mean((lowerCI[,2]-0.5)*(upperCI[,2]-0.5)<0)
}

##################################################
####Scenario 2
##################################################
for (samplesize in c(400,800)){
  set.seed(2000)
  numberofsets=500
  result_02=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  cov_02=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  results_02=list()
  for (i in 1:numberofsets){
    print(i)
    result=generate_sample_random1(beta=c(0,1),N=samplesize,constX=constX,constA=constA,expo=expo)
    result_02[i,]=result$beta
    cov_02[i,1]=result$cov[1,1]
    cov_02[i,2]=result$cov[2,2]
    results_02[[i]]=result
    write.csv(result_02, paste("beta_02", samplesize, ".csv"))
    write.csv(cov_02,paste("cov_02", samplesize, ".csv"))
    save.image(paste("result_02", samplesize, ".RData"))
  }
  
  colMeans(result_02) #mean of the estimated parameter
  apply(result_02,2,sd) # sd of the estimated parameters
  colMeans(sqrt(cov_02))  # mean of estimated sd
  #Confidence Intervals
  lowerCI=result_02-sqrt(cov_02)*1.96
  upperCI=result_02+sqrt(cov_02)*1.96
  #coverage of confidence intervals
  mean(lowerCI[,1]*upperCI[,1]<0)
  mean((lowerCI[,2]-1)*(upperCI[,2]-1)<0)
}

##################################################
####Scenario 3
##################################################

for (samplesize in c(400,800)){
  set.seed(2000)
  numberofsets=500
  result_03=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  cov_03=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  results_03=list()
  for (i in 1:numberofsets){
    print(i)
    result=generate_sample_random1(beta=c(0,0.5),mu0=mu2, N=samplesize,constX=constX,constA=constA,expo=expo)
    result_03[i,]=result$beta
    cov_03[i,1]=result$cov[1,1]
    cov_03[i,2]=result$cov[2,2]
    results_03[[i]]=result
    write.csv(result_03, paste("beta_03", samplesize, ".csv"))
    write.csv(cov_03,paste("cov_03", samplesize, ".csv"))
    save.image(paste("result_03", samplesize, ".RData"))
  }
  
  colMeans(result_03) #mean of the estimated parameter
  apply(result_03,2,sd) # sd of the estimated parameters
  colMeans(sqrt(cov_03))  # mean of estimated sd
  #Confidence Intervals
  lowerCI=result_03-sqrt(cov_03)*1.96
  upperCI=result_03+sqrt(cov_03)*1.96
  #coverage of confidence intervals
  mean(lowerCI[,1]*upperCI[,1]<0)
  mean((lowerCI[,2]-0.5)*(upperCI[,2]-0.5)<0)
}




##################################################
##########Scenario 4
##################################################

for (samplesize in c(400,800)){
  set.seed(2000)
  numberofsets=500
  result_04=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  cov_04=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  results_04=list()
  for (i in 1:numberofsets){
    print(i)
    result=generate_sample_random1(beta=c(0,1), mu0=mu2, N=samplesize, constX=constX,constA=constA,expo=expo)
    result_04[i,]=result$beta
    cov_04[i,1]=result$cov[1,1]
    cov_04[i,2]=result$cov[2,2]
    results_04[[i]]=result
    write.csv(result_04, paste("beta_04", samplesize, ".csv"))
    write.csv(cov_04,paste("cov_04", samplesize, ".csv"))
    save.image(paste("result_04", samplesize, ".RData"))
  }
  
  colMeans(result_04) #mean of the estimated parameter
  apply(result_04,2,sd) # sd of the estimated parameters
  colMeans(sqrt(cov_04))  # mean of estimated sd
  #Confidence Intervals
  lowerCI=result_04-sqrt(cov_04)*1.96
  upperCI=result_04+sqrt(cov_04)*1.96
  #coverage of confidence intervals
  mean(lowerCI[,1]*upperCI[,1]<0)
  mean((lowerCI[,2]-1)*(upperCI[,2]-1)<0)
}


##################################################
#############Observational studies
###################################################

##################################################
####Scenario 1
##################################################
expo=4.5
constX=0.8
constA=3.2

for (samplesize in c(400,800)){
  set.seed(2000)
  numberofsets=500
  result_01=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  #cov_01=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  cov_01=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  results_01=list()
  for (i in 1:numberofsets){
    print(i)
    result=generate_sample_obs1(beta=c(0,0.5),N=samplesize,constX=constX,constA=constA,expo=expo)
    result_01[i,]=result$beta
    cov_01[i,1]=result$cov[1,1]
    cov_01[i,2]=result$cov[2,2]
    results_01[[i]]=result
    write.csv(result_01,paste("beta_01_obs",samplesize,".csv"))
    write.csv(cov_01,paste("cov_01_obs", samplesize,".csv"))
    save.image(paste("result_01_obs", samplesize, ".RData"))
  }
  colMeans(result_01)
  apply(result_01,2,sd)
  sqrt(colMeans(cov_01))
  lowerCI=result_01-sqrt(cov_01)*1.96
  upperCI=result_01+sqrt(cov_01)*1.96
  mean(lowerCI[,1]*upperCI[,1]<0)
  mean((lowerCI[,2]-0.5)*(upperCI[,2]-0.5)<0)
}

##################################################
####Scenario 2
##################################################

expo=4.5
constX=0.9
constA=2.35


for (samplesize in c(400,800)){
  set.seed(2000)
  numberofsets=500
  result_02=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  #cov_02=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  cov_02=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  results_02=list()
  for (i in 1:numberofsets){
    print(i)
    result=generate_sample_obs1(beta=c(0,1),N=samplesize,constX=constX,constA=constA,expo=expo)
    result_02[i,]=result$beta
    cov_02[i,1]=result$cov[1,1]
    cov_02[i,2]=result$cov[2,2]
    results_02[[i]]=result
    write.csv(result_02,paste("beta_02_obs",samplesize,".csv"))
    write.csv(cov_02,paste("cov_02_obs", samplesize,".csv"))
    save.image(paste("result_02_obs", samplesize, ".RData"))
  }
  colMeans(result_02)
  apply(result_02,2,sd)
  sqrt(colMeans(cov_02))
  lowerCI=result_02-sqrt(cov_02)*1.96
  upperCI=result_02+sqrt(cov_02)*1.96
  mean(lowerCI[,1]*upperCI[,1]<0)
  mean((lowerCI[,2]-1)*(upperCI[,2]-1)<0)
}

##################################################
####Scenario 3
##################################################
expo=4.5
constX=0.9
constA=3.05

for (samplesize in c(400,800)){
  set.seed(2000)
  numberofsets=500
  result_03=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  #cov_03=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  cov_03=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  results_03=list()
  for (i in 1:numberofsets){
    print(i)
    result=generate_sample_obs1(beta=c(0,0.5),mu0=mu2,N=samplesize,constX=constX,constA=constA,expo=expo)
    result_03[i,]=result$beta
    cov_03[i,1]=result$cov[1,1]
    cov_03[i,2]=result$cov[2,2]
    results_03[[i]]=result
    write.csv(result_03,paste("beta_03_obs",samplesize,".csv"))
    write.csv(cov_03,paste("cov_03_obs", samplesize,".csv"))
    save.image(paste("result_03_obs", samplesize, ".RData"))
  }
  colMeans(result_03)
  apply(result_03,2,sd)
  sqrt(colMeans(cov_03))
  lowerCI=result_03-sqrt(cov_03)*1.96
  upperCI=result_03+sqrt(cov_03)*1.96
  mean(lowerCI[,1]*upperCI[,1]<0)
  mean((lowerCI[,2]-0.5)*(upperCI[,2]-0.5)<0)
}


##################################################
####Scenario 4
##################################################

expo=4.5
constX=0.9
constA=2.35

for (samplesize in c(400,800)){
  set.seed(2000)
  numberofsets=500
  result_04=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  #cov_04=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  cov_04=matrix(rep(NA,2*numberofsets),nrow=numberofsets)
  results_04=list()
  for (i in 1:numberofsets){
    print(i)
    result=generate_sample_obs1(beta=c(0,1),mu0=mu2,N=samplesize,constX=constX,constA=constA,expo=expo)
    result_04[i,]=result$beta
    cov_04[i,1]=result$cov[1,1]
    cov_04[i,2]=result$cov[2,2]
    results_04[[i]]=result
    write.csv(result_04,paste("beta_04_obs",samplesize,".csv"))
    write.csv(cov_04,paste("cov_04_obs", samplesize,".csv"))
    save.image(paste("result_04_obs", samplesize, ".RData"))
  }
  colMeans(result_04)
  apply(result_04,2,sd)
  sqrt(colMeans(cov_04))
  lowerCI=result_04-sqrt(cov_04)*1.96
  upperCI=result_04+sqrt(cov_04)*1.96
  mean(lowerCI[,1]*upperCI[,1]<0)
  mean((lowerCI[,2]-1)*(upperCI[,2]-1)<0)
}

##################################################
##########Scenario 5, sample size 400
##################################################


#cross validation for finding the optimal constants for bandwidths
#generate a dataset for cross-validation

C=-25
N=800
beta=c(1,0.5,0.5,0)
mbeta=c(8,4,-2,-2)
sd0=1
p=length(beta)-1
A=runif(N,min=0,max=2)
X=matrix(runif(N*p,min=-1,max=1),nrow=N)
Q=C*(A-cbind(1,X)%*% beta)^2
mu=Q+cbind(1,X)%*% mbeta
Y=rnorm(rep(1,N),mean=mu,sd=rep(sd0,N))

par=optim(rep(1,4),cross_val1, X=X, A=A, Y=Y, expo=8.5, fold=5, track=TRUE)$par
#constX=par[1:p]
#constA=par[p+1]
#the constants we got from cross validation
expo=8.5
constX=c(0.5674214, 0.5818246, 3.4937701)
constA= 0.3179026

##start simulation

for (samplesize in c(400,800)){
  numberofsets=500
  result_05=matrix(rep(NA,4*numberofsets),nrow=numberofsets)
  cov_05=matrix(rep(NA,4*numberofsets),nrow=numberofsets)
  cov_05=matrix(rep(NA,4*numberofsets),nrow=numberofsets)
  bootsd_05=matrix(rep(NA,4*numberofsets),nrow=numberofsets)
  results_05=list()
  for (i in 1:numberofsets){
    print(i)
    set.seed(1000+i)
    result=generate_sample_random5(N=samplesize,constX=constX,constA=constA,expo=expo)
    result_05[i,]=result$beta
    cov_05[i,1]=result$cov[1,1]
    cov_05[i,2]=result$cov[2,2]
    cov_05[i,3]=result$cov[3,3]
    cov_05[i,4]=result$cov[4,4]
    bootsd_05[i,]=result$bootstrap_sd
    results_05[[i]]=result
    write.csv(result_05,paste("beta_05",samplesize,".csv"))
    write.csv(cov_05,paste("cov_05",samplesize,".csv"))
    write.csv(bootsd_05,paste("bootsd_05",samplesize,".csv"))
    save.image(paste("result_05",samplesize,".RData"))
  }
  
  colMeans(result_05)
  apply(result_05,2,sd)
  sqrt(colMeans(cov_05)) #mean of estimated sd
  colMeans(bootsd_05) # mean of bootstrap sd 
  #confidence intervals
  lowerCI=result_05-sqrt(cov_05)*1.96
  upperCI=result_05+sqrt(cov_05)*1.96
  mean((lowerCI[,1]-1)*(upperCI[,1]-1)<0)
  mean((lowerCI[,2]-0.5)*(upperCI[,2]-0.5)<0)
  mean((lowerCI[,3]-0.5)*(upperCI[,3]-0.5)<0)
  mean((lowerCI[,4])*(upperCI[,4])<0)
  #bootstrap confidence intervals
  lowerCI=result_05-bootsd_05*1.96
  upperCI=result_05+bootsd_05*1.96
  mean((lowerCI[,1]-1)*(upperCI[,1]-1)<0)
  mean((lowerCI[,2]-0.5)*(upperCI[,2]-0.5)<0)
  mean((lowerCI[,3]-0.5)*(upperCI[,3]-0.5)<0)
  mean((lowerCI[,4])*(upperCI[,4])<0)
  
  
  #calculate value function
  C=-25
  testN=1000
  beta=c(1,0.5,0.5,0)
  mbeta=c(8,4,-2,-2)
  V=rep(NA,numberofsets)
  for (i in 1:numberofsets){
    print(i)
    set.seed(5000+i)
    #generate a separate dataset for evaluting value
    p=length(beta)-1
    A=runif(testN,min=0,max=2)
    X=matrix(runif(testN*p,min=-1,max=1),nrow=testN)
    V[i]=mean(C*(cbind(1,X)%*%t(result_05[i,])-cbind(1,X)%*% beta)^2+cbind(1,X)%*% mbeta)
  }
  print(mean(V))
  
}





