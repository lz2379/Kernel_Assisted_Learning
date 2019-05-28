########################################################################
######################O-learning for setting5###########################
########################################################################

source('O-learning-functions.r')
ind=5
allloop =0
datagen = Scenario1
covsize=30
epsilon = 0.1


solution0=c(1,0.5,0.5,0)

numberofsets=500



for (samplesize in c(400,800)){
  results_lol_s1 <- list()
  N=samplesize
  lol_beta_hat=matrix(nrow=numberofsets,ncol=4)
  lol_value_hat=rep(NA,numberofsets)
  kol_value_hat=rep(NA,numberofsets)
  C=-25
  beta=c(1,0.5,0.5,0)
  mbeta=c(8,4,-2,-2)
  sd0=1
  for (i in 1:numberofsets){
    #generate data
    print(i)
    set.seed(1000+i)
    p=length(beta)-1
    A=runif(N,min=0,max=2)
    X=matrix(runif(N*p,min=-1,max=1),nrow=N)
    Q=C*(A-cbind(1,X)%*% beta)^2
    mu=Q+cbind(1,X)%*% mbeta
    Y=rnorm(rep(1,N),mean=mu,sd=rep(sd0,N))
    
    #O-learning
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
    
    #generate a new dataset for evaluating value function 
    set.seed(5000+i)
    C_new=-25
    N_new=1000
    beta_new=c(1,0.5,0.5,0)
    mbeta_new=c(8,4,-2,-2)
    p_new=length(beta_new)-1
    A_new=runif(N_new,min=0,max=2)
    X_new=matrix(runif(N_new*p_new,min=-1,max=1),nrow=N_new)
    
    tmpvalue = rqmodel$coefficients[1] + X_new %*% rqmodel$coefficients[-1]
    #estimate value
    Yhat=C*(tmpvalue-cbind(1,X_new)%*% beta)^2+cbind(1,X_new)%*% mbeta
    
    lol_value_hat[i]=mean(Yhat)
    write.csv(lol_value_hat,paste("results/lol-sim05-value",samplesize, ".csv"),row.names=FALSE)
    
    ### model with estimated propensity score times reward/outcome as weigths
    
    model_pen = svm(x = X[index,], y = A[index], w= (trainweight*weight.gbm)[index], type="eps-regression",
                    epsilon = 0.15, scale=FALSE)
    tmpvalue2 = predict(model_pen,X_new)
    tmpvalue2 = pmin(pmax(tmpvalue2,0),1)
    Yhat2=C_new*(tmpvalue2-cbind(1,X_new)%*% beta_new)^2+cbind(1,X_new)%*% mbeta_new
    kol_value_hat[i]=mean(Yhat2)
    write.csv(kol_value_hat,paste("results/kol-sim05-value",samplesize,".csv"),row.names=FALSE)
    
  }
}
