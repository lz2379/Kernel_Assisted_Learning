###############################################################
#######Outcome weighted learning for real data
source('O-learning-functions.r')
source('KAL-functions.R')
#### Load dependent packages ############

Data <- read.csv("data-processed.csv")
#Data=Data[complete.cases(Data),]
NN=dim(Data)[1]
c_index=c(3)
d_index=c(6,9)
var_index=c(c_index,d_index)

Xc=Data[,c_index]
Xd=Data[,d_index]
X=Data[,var_index]
A=Data$Dose_standardized
Y=Data$Y
p=length(var_index)
pc=length(c_index)
pd=length(d_index)

par=c(1.25, 1.75)
constX=par[1:pc]
constA=par[pc+1]
expo=4.5
####################################################
#########1. Dose finding with all data
####################################################
ind=5
allloop =0
datagen = Scenario1
covsize=30
epsilon = 0.1

results_lol_s1 <- list()
  
constant = min(quantile(Y,0.6),0)
trainweight = Y - constant
index = which(Y > quantile(Y,0.6))
#estimate propensity score
### We use the penalized quantile regression to enhance the model fitting.
#rqmodel = rq(trainA[index] ~ trainX[index,],.5,method="lasso",weights=trainweight[index],lambda = 1)
#predict dose
############### propensity score part
mydata = data.frame(T = A,X = X)
model.num = lm(T~1,data = mydata)
ps.num= dnorm((mydata$T-model.num$fitted)/(summary(model.num))$sigma,0,1)
model.den = gbm(T~.,data = mydata, shrinkage = 0.0005,
                  interaction.depth = 4, distribution = "gaussian",n.trees = 20000)
opt = optimize(F.aac.iter,interval = c(1,20000), data = mydata, ps.model =
                   model.den,
                 ps.num = ps.num,rep = 50,criterion = "pearson")
best.aac.iter = opt$minimum
best.aac = opt$objective
# Calculate the inverse probability weights
model.den$fitted = predict(model.den,newdata = mydata,
                             n.trees = floor(best.aac.iter), type = "response")
ps.den = dnorm((mydata$T-model.den$fitted)/sd(mydata$T-model.den
                                                $fitted),0,1)
weight.gbm = ps.num/ps.den
  
### model with estimated propensity score times reward/outcome as weigths
  
rqmodel = rq(A[index] ~ X[index,],.5,method="lasso", w= (trainweight*weight.gbm)[index],lambda=1)
tmpvalue = rqmodel$coefficients[1] + X %*% rqmodel$coefficients[-1]
  
lol_dose_hat=as.numeric(tmpvalue)
### model with estimated propensity score times reward/outcome as weigths
  
model_pen = svm(x = trainX[index,], y = trainA[index], w= (trainweight*weight.gbm)[index], type="eps-regression",
                  epsilon = 0.15, scale=FALSE)
tmpvalue2 = predict(model_pen,testX)
tmpvalue2 = pmin(pmax(tmpvalue2,0),1)
  
kol_dose_hat=as.numeric(tmpvalue2)
#plot(density(lol_dose_hat),col=2,main="Dose_hat using all data",ylim=c(0,15),xlim=c(0,1))
#lines(density(kol_dose_hat),col=3)
#legend(x=0.4,y=10,legend=c("lol","kol"),lty=c(1,1),col=c(2,3))
write.csv(lol_dose_hat,"results/lol-real-dose.csv",row.names=FALSE)
write.csv(kol_dose_hat,"results/kol-real-dose.csv",row.names=FALSE)
#save.image("O-learning-real.Rdata")


####################################################
#########2. Estimate value function
####################################################

ind=5
allloop =0
datagen = Scenario1
covsize=30
epsilon = 0.1

results_lol_s1 <- list()

testsize=floor(NN/3)
trainsize=NN-testsize
lol_beta_hat=matrix(nrow=200,ncol=p+1)
lol_value_hat=rep(NA,200)
lol_dose_hat=matrix(nrow=200,ncol=testsize)
kol_value_hat=rep(NA,200)
kol_dose_hat=matrix(nrow=200,ncol=testsize)
for (i in 1:200){
  print(i)
  set.seed(10000+i)
  train_index=sample(NN,size=trainsize)
  train=Data[train_index,]
  test=Data[-train_index,]
  
  trainX=as.matrix(train[,var_index])
  trainY=train$Y
  trainA=train$Dose_standardized
  trainXc=as.matrix(train[,c_index])
  trainXd=as.matrix(train[,d_index])
  
  testX=as.matrix(test[,var_index])
  testY=test$Y
  testA=test$Dose_standardized 
  testXc=as.matrix(test[,c_index])
  testXd=as.matrix(test[,d_index])
  
  constant = min(quantile(trainY,0.6),0)
  trainweight = trainY - constant
  index = which(trainY > quantile(trainY,0.6))
  #estimate propensity score
  ### We use the penalized quantile regression to enhance the model fitting.
  #rqmodel = rq(trainA[index] ~ trainX[index,],.5,method="lasso",weights=trainweight[index],lambda = 1)
  #predict dose
  
  ############### propensity score part
  mydata = data.frame(T = trainA,X = trainX)
  model.num = lm(T~1,data = mydata)
  ps.num= dnorm((mydata$T-model.num$fitted)/(summary(model.num))$sigma,0,1)
  model.den = gbm(T~.,data = mydata, shrinkage = 0.0005,
                  interaction.depth = 4, distribution = "gaussian",n.trees = 20000)
  opt = optimize(F.aac.iter,interval = c(1,20000), data = mydata, ps.model =
                   model.den,
                 ps.num = ps.num,rep = 50,criterion = "pearson")
  best.aac.iter = opt$minimum
  best.aac = opt$objective
  # Calculate the inverse probability weights
  model.den$fitted = predict(model.den,newdata = mydata,
                             n.trees = floor(best.aac.iter), type = "response")
  ps.den = dnorm((mydata$T-model.den$fitted)/sd(mydata$T-model.den
                                                $fitted),0,1)
  weight.gbm = ps.num/ps.den
  
  ### model with estimated propensity score times reward/outcome as weights
  
  rqmodel = rq(trainA[index] ~ trainX[index,],.5,method="lasso", w= (trainweight*weight.gbm)[index],lambda=1)
  tmpvalue = rqmodel$coefficients[1] + testX %*% rqmodel$coefficients[-1]
  
  
  #estimate value
  hc=constX*apply(as.matrix(testXc),2,sd)*testsize^(-1/expo)
  h2=constA*sd(testA)*testsize^(-1/expo)
  
  Aj=matrix(rep(testA,each=testsize),nrow=testsize)
  predictAi=matrix(rep(tmpvalue,testsize),nrow=testsize)
  Ka=as.matrix(K((Aj-predictAi)/h2))
  
  if (length(hc)>1){
    dXc=as.matrix(dist(as.matrix(testXc)%*%diag(as.numeric(1/hc))))
  } else {
    dXc=as.matrix(dist(as.matrix(testXc)/hc))
  }
  Kxc=K(dXc)/exp(sum(log(hc)))
  Ixd=I(as.matrix(dist(testXd)))
  Kx=Kxc*Ixd
  
  Yj=matrix(rep(testY,each=testsize),nrow=testsize)
  lol_beta_hat[i,]=as.numeric(rqmodel$coefficients)
  lol_value_hat[i]=mean(rowMeans(Yj*Ka*Kx)/rowMeans(Ka*Kx))
  lol_dose_hat[i,]=as.numeric(tmpvalue)
  write.csv(lol_value_hat,"results/lol-real-value.csv",row.names=FALSE)
  
  ### model with estimated propensity score times reward/outcome as weigths
  
  model_pen = svm(x = trainX[index,], y = trainA[index], w= (trainweight*weight.gbm)[index], type="eps-regression",
                  epsilon = 0.15, scale=FALSE)
  tmpvalue2 = predict(model_pen,testX)
  tmpvalue2 = pmin(pmax(tmpvalue2,0),1)
  

  predictAi=matrix(rep(tmpvalue2,testsize),nrow=testsize)
  Ka=as.matrix(K((Aj-predictAi)/h2))

  kol_value_hat[i]=mean(rowMeans(Yj*Ka*Kx)/rowMeans(Ka*Kx))
  kol_dose_hat[i,]=as.numeric(tmpvalue2)
  write.csv(kol_value_hat,"results/kol-real-value.csv",row.names=FALSE)
}




