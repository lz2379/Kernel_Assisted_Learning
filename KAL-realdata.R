#########Kernel Assisted Learning applied to the warfarin dataset#############################################


##################################################################
#######1. dose finding using all real data
####################################################################
source('KAl-functions.r')
Data <- read.csv("data_processed.csv")
NN=dim(Data)[1]

c_index=c(3)    #column index for continuous variables
d_index=c(6,9)  #column index for discrete variables
var_index=c(c_index,d_index)
Xc=as.matrix(Data[,c_index])
Xd=as.matrix(Data[,d_index])
X=Data[,var_index]
A=Data$Dose_standardized
Y=Data$Y

p=length(var_index)  
pc=length(c_index)   #number of continuous variables
pd=length(d_index)   #number of discrete variables
#setting bandwidths
par=c(1.25,1.75)
constX=par[1:pc]     
constA=par[pc+1]
expo=4.5
hc=constX*apply(as.matrix(Xc),2,sd)*NN^(-1/expo)
ha=constA*sd(A)*NN^(-1/expo)
#find the range of the continuous variable
xlim=as.matrix(cbind(apply(Xc,2,min),apply(Xc,2,max)))

#inital points
beta0=rbind(rep(0,p+1),diag(rep(0.5,p+1)),diag(rep(-0.5,p+1)))
scale=c(1,1/apply(as.matrix(X),2,sd))

results=finddose2(Y=Y,A=A,Xc=as.matrix(Xc),
                  Xd=as.matrix(Xd),hc=hc,ha=ha,init=beta0,xlim=xlim,sampling=4000,
                  g=g,g1=g1,g2=g2)
beta_hat=as.numeric(results$beta)
#here we got the solution
#c(0.052298357, -0.002358858,  0.109252696, -0.672216357)
dose_hat=g(as.matrix(cbind(1,X))%*%beta_hat)
#save.image("KAL-real.Rdata")
write.csv(dose_hat,"results/KAL-real-dose.csv",row.names=FALSE)

########################################################
###########2. Bootstrap for estimating sd of the solution
##############################################################
beta0=t( c(0.052298357, -0.002358858,  0.109252696, -0.672216357))
for (i in (1:200)){
  print(i)
  set.seed(10000+i)
  train_index=sample(NN,size=NN, replace=TRUE)
  train=Data[train_index,]
  trainX=train[,var_index]
  trainY=train$Y
  trainA=train$Dose_standardized
  trainXc=train[,c_index]
  trainXd=train[,d_index]
  hc=constX*apply(as.matrix(trainXc),2,sd)*NN^(-1/expo)
  h2=constA*sd(trainA)*NN^(-1/expo)
  
  results[[i]]=finddose2(Y=trainY,A=trainA,Xc=as.matrix(trainXc),
                         Xd=as.matrix(trainXd),hc=hc,ha=ha,init=beta0,xlim=xlim,sampling=8000,calc_cov=FALSE,
                         g=g,g1=g1,g2=g2)
  beta_hat[i,]=as.numeric(results[[i]]$beta)
  write.csv(beta_hat,"results/KAL-real-bootstrap.csv", row.names=FALSE)
}    

########################################################
###########3. Estimate value function
########################################################
beta0=t( c(0.052298357, -0.002358858,  0.109252696, -0.672216357))
###Using one third as testing dataset
testsize=floor(NN/3)
trainsize=NN-testsize
beta_hat=matrix(nrow=200,ncol=(p+1))
value_hat=rep(NA,200)
results=list()
for (i in (1:200)){
  print(i)
  set.seed(10000+i)
  train_index=sample(NN,size=trainsize)
  train=Data[train_index,]
  test=Data[-train_index,]
  
  trainX=train[,var_index]
  trainY=train$Y
  trainA=train$Dose_standardized
  trainXc=train[,c_index]
  trainXd=train[,d_index]
  
  testX=test[,var_index]
  testY=test$Y
  testA=test$Dose_standardized
  testXc=test[,c_index]
  testXd=test[,d_index]
  
  hc=constX*apply(as.matrix(trainXc),2,sd)*trainsize^(-1/expo)
  ha=constA*sd(trainA)*trainsize^(-1/expo)
  results[[i]]=finddose2(Y=trainY,A=trainA,Xc=as.matrix(trainXc),
                         Xd=as.matrix(trainXd),hc=hc,ha=ha,xlim=xlim,init=beta0,sampling=4000,
                         g=g,g1=g1,g2=g2)
  beta_hat[i,]=as.numeric(results[[i]]$beta)
  
  hc=constX*apply(as.matrix(testXc),2,sd)*testsize^(-1/expo)
  ha=constA*sd(testA)*testsize^(-1/expo)
  
  value_hat[i]=testvalue(beta_hat[i,],testY,testA,testXc,testXd,testsize,hc=hc,ha=ha)
  write.csv(value_hat,"results/KAL-real-value.csv",row.names=FALSE)
}



#########################################
#########################################

