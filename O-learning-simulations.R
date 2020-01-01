

source('O-learning-functions.r')
source('KAL-functions.r')
ind = 5
allloop  = 0
datagen  =  Scenario1
epsilon  =  0.1

########################################################################
######################O-learning for setting1-4 ###########################
########################################################################

##############################
#LOL SETTING 1-4
##############################
#Randomized trials


beta0 = c(0,0.5) #beta0 = c(0,1) for setting 2 4

numberofsets = 500

for (samplesize in c(400,800)){
  results_lol_s1 <- list()
  N = samplesize
  lol_beta_hat = matrix(nrow = numberofsets,ncol = length(beta0)) 
  cat("\n")
  for (i in 1:numberofsets){
    #generate data
    cat(i," ")
    set.seed(1000+i)
    dataset = generate_sample_random1(beta = beta0, N = samplesize) #mu0 = mu2 for setting 3 4
    Y = dataset$Y
    X = as.matrix(dataset$X)
    A = dataset$A
    N = dataset$N
    p = dataset$p
    #O-learning
    constant  =  min(quantile(Y, 0.6), 0)
    trainweight  =  Y - constant
    index  =  which(Y > quantile(Y, 0.6))
    #estimate propensity score
    ### We use the penalized quantile regression to enhance the model fitting.
    #rqmodel  =  rq(trainA[index] ~ trainX[index,],.5,method = "lasso",weights = trainweight[index],lambda  =  1)
    #predict dose
    
    ############### propensity score part
    mydata  =  data.frame(T  =  A,X  =  X)
    model.num  =  lm(T~1,data  =  mydata)
    ps.num =  dnorm((mydata$T - model.num$fitted) / (summary(model.num))$sigma,0,1)
    model.den  =  gbm(T~.,data  =  mydata, shrinkage  =  0.0005,
                    interaction.depth  =  4, distribution  =  "gaussian",n.trees  =  20000)
    opt  =  optimize(F.aac.iter,interval  =  c(1,20000), data  =  mydata, ps.model  = 
                     model.den,
                   ps.num  =  ps.num,rep  =  50,criterion  =  "pearson")
    best.aac.iter  =  opt$minimum
    best.aac  =  opt$objective
    # Calculate the inverse probability weights
    model.den$fitted  =  predict(model.den,newdata  =  mydata,
                               n.trees  =  floor(best.aac.iter), type  =  "response")
    ps.den  =  dnorm((mydata$T - model.den$fitted) / sd(mydata$T-model.den
                                                  $fitted),0,1)
    weight.gbm  =  ps.num / ps.den
    
    ### model with estimated propensity score times reward/outcome as weigths
    
    rqmodel  =  rq(A[index] ~ X[index,],.5,method = "lasso", w =  (trainweight*weight.gbm)[index],lambda = 1)
    
    lol_beta_hat[i,] = rqmodel$coefficients

  }
  cat("\n Sample Size:",samplesize,"\n")
  
  #calculate value
  v = sim_value1(betahat = lol_beta_hat,
                 dose_function =  dose_lol,
                 beta = beta0,
                 testN = 1000) #mu0 = mu2 for setting 3 4
}

#########
#Observational Studies

beta0 = c(0,1) #beta0 = c(0,1) for setting 2 4

numberofsets = 500

for (samplesize in c(400,800)){
  results_lol_s1 <- list()
  N = samplesize
  lol_beta_hat = matrix(nrow = numberofsets,ncol = length(beta0)) 
  #lol_value_hat = rep(NA,numberofsets)
  #kol_value_hat = rep(NA,numberofsets)
  cat("\n")
  for (i in 1:numberofsets){
    #generate data
    cat(i," ")
    set.seed(1000+i)
    dataset = generate_sample_obs1(beta = beta0, 
                                   N = samplesize,
                                   mu0 = mu2) #mu0 = mu2 for setting 3 4
    Y = dataset$Y
    X = as.matrix(dataset$X)
    A = dataset$A
    N = dataset$N
    p = dataset$p
    #O-learning
    constant  =  min(quantile(Y,0.6),0)
    trainweight  =  Y - constant
    index  =  which(Y > quantile(Y,0.6))
    #estimate propensity score
    ### We use the penalized quantile regression to enhance the model fitting.
    #rqmodel  =  rq(trainA[index] ~ trainX[index,],.5,method = "lasso",weights = trainweight[index],lambda  =  1)
    #predict dose
    
    ############### propensity score part
    mydata  =  data.frame(T  =  A,X  =  X)
    model.num  =  lm(T~1,data  =  mydata)
    ps.num =  dnorm((mydata$T-model.num$fitted)/(summary(model.num))$sigma,0,1)
    model.den  =  gbm(T~.,data  =  mydata, shrinkage  =  0.0005,
                    interaction.depth  =  4, distribution  =  "gaussian",n.trees  =  20000)
    opt  =  optimize(F.aac.iter,interval  =  c(1,20000), data  =  mydata, ps.model  = 
                     model.den,
                   ps.num  =  ps.num,rep  =  50,criterion  =  "pearson")
    best.aac.iter  =  opt$minimum
    best.aac  =  opt$objective
    # Calculate the inverse probability weights
    model.den$fitted  =  predict(model.den,newdata  =  mydata,
                               n.trees  =  floor(best.aac.iter), type  =  "response")
    ps.den  =  dnorm((mydata$T-model.den$fitted)/sd(mydata$T-model.den
                                                  $fitted),0,1)
    weight.gbm  =  ps.num/ps.den
    
    ### model with estimated propensity score times reward/outcome as weigths
    
    rqmodel  =  rq(A[index] ~ X[index,],.5,method = "lasso", w =  (trainweight*weight.gbm)[index],lambda = 1)
    
    lol_beta_hat[i,] = rqmodel$coefficients

  }
  cat("\n Sample Size:",samplesize,"\n")
  v = sim_value1(betahat = lol_beta_hat,
                 dose_function =  dose_lol,
                 beta = beta0,testN = 1000,
                 mu0 = mu2) #mu0 = mu2 for setting 3 4
}



##############################
#KOL SETTING 1-4
##############################
#Randomized trials

beta0 = c(0,1) #beta0 = c(0,1) for setting 2 4

numberofsets = 500
meanfunction = mu2
for (samplesize in c(800)){
  results_lol_s1 <- list()
  N = samplesize
  #lol_beta_hat = matrix(nrow = numberofsets,ncol = length(beta0)) 
  #lol_value_hat = rep(NA,numberofsets)
  kol_value_hat = rep(NA,numberofsets)
  cat("\n")
  for (i in 1:numberofsets){
    #generate data
    cat(i," ")
    set.seed(1000+i)
    dataset = generate_sample_random1(beta = beta0, 
                                      N = samplesize,
                                      mu0 = meanfunction) #mu0 = mu2 for setting 3 4
    Y = dataset$Y
    X = as.matrix(dataset$X)
    A = dataset$A
    N = dataset$N
    p = dataset$p
    #O-learning
    constant  =  min(quantile(Y,0.6),0)
    trainweight  =  Y - constant
    index  =  which(Y > quantile(Y,0.6))
    #estimate propensity score
    ### We use the penalized quantile regression to enhance the model fitting.
    #rqmodel  =  rq(trainA[index] ~ trainX[index,],.5,method = "lasso",weights = trainweight[index],lambda  =  1)
    #predict dose
    
    ############### propensity score part
    mydata  =  data.frame(T  =  A,X  =  X)
    model.num  =  lm(T~1,data  =  mydata)
    ps.num =  dnorm((mydata$T-model.num$fitted)/(summary(model.num))$sigma,0,1)
    model.den  =  gbm(T~.,data  =  mydata, shrinkage  =  0.0005,
                    interaction.depth  =  4, distribution  =  "gaussian",n.trees  =  20000)
    opt  =  optimize(F.aac.iter,interval  =  c(1,20000), data  =  mydata, ps.model  = 
                     model.den,
                   ps.num  =  ps.num,rep  =  50,criterion  =  "pearson")
    best.aac.iter  =  opt$minimum
    best.aac  =  opt$objective
    # Calculate the inverse probability weights
    model.den$fitted  =  predict(model.den,newdata  =  mydata,
                               n.trees  =  floor(best.aac.iter), type  =  "response")
    ps.den  =  dnorm((mydata$T-model.den$fitted)/sd(mydata$T-model.den
                                                  $fitted),0,1)
    weight.gbm  =  ps.num/ps.den
    
    ### model with estimated propensity score times reward/outcome as weigths
    

    #generate a new dataset for evaluating value function 
    set.seed(5000+i)
    C_new = -10
    N_new = 1000
    beta_new = beta0
    mbeta_new = meanfunction
    p_new = length(beta_new)-1
    X_new = matrix(rnorm(N_new*p,0,1),nrow = N_new)
    dose_optim = g(cbind(1,X_new)%*% beta_new)
    ### model with estimated propensity score times reward/outcome as weigths
    
    model_pen  =  svm(x  =  X[index,], y  =  A[index], w =  (trainweight*weight.gbm)[index], type = "eps-regression",
                    epsilon  =  0.15, scale = FALSE)
    tmpvalue2  =  predict(model_pen,X_new)
    tmpvalue2  =  pmin(pmax(tmpvalue2,0),1)
    Yhat2 = mbeta_new(X_new)+C_new*(tmpvalue2-dose_optim)^2
    kol_value_hat[i] = mean(Yhat2)
    write.csv(kol_value_hat,paste("results/simulation/kol-sim04-value",samplesize,".csv"),row.names = FALSE)
    
  }
  cat("\n Sample Size:",samplesize,"\n")
  cat("Mean of V:", mean(kol_value_hat),"\n")
  cat("Sd of V:", sd(kol_value_hat),"\n")
}

########
#Observational Studies

beta0 = c(0,1) #beta0 = c(0,1) for setting 2 4

numberofsets = 500
meanfunction = mu2
for (samplesize in c(800)){
  results_lol_s1 <- list()
  N = samplesize
  #lol_beta_hat = matrix(nrow = numberofsets,ncol = length(beta0)) 
  #lol_value_hat = rep(NA,numberofsets)
  kol_value_hat = rep(NA,numberofsets)
  cat("\n")
  for (i in 1:numberofsets){
    #generate data
    cat(i," ")
    set.seed(1000+i)
    dataset = generate_sample_obs1(beta = beta0, 
                                      N = samplesize,
                                      mu0 = meanfunction) #mu0 = mu2 for setting 3 4
    Y = dataset$Y
    X = as.matrix(dataset$X)
    A = dataset$A
    N = dataset$N
    p = dataset$p
    #O-learning
    constant  =  min(quantile(Y, 0.6), 0)
    trainweight  =  Y - constant
    index  =  which(Y > quantile(Y, 0.6))
    #estimate propensity score
    ### We use the penalized quantile regression to enhance the model fitting.
    #rqmodel  =  rq(trainA[index] ~ trainX[index,],.5,method = "lasso",weights = trainweight[index],lambda  =  1)
    #predict dose
    
    ############### propensity score part
    mydata  =  data.frame(T  =  A,X  =  X)
    model.num  =  lm(T~1,data  =  mydata)
    ps.num =  dnorm((mydata$T-model.num$fitted) / (summary(model.num))$sigma,0,1)
    model.den  =  gbm(T~.,data  =  mydata, shrinkage  =  0.0005,
                      interaction.depth  =  4, distribution  =  "gaussian",n.trees  =  20000)
    opt  =  optimize(F.aac.iter,interval  =  c(1,20000), data  =  mydata, ps.model  = 
                       model.den,
                     ps.num  =  ps.num,rep  =  50,criterion  =  "pearson")
    best.aac.iter  =  opt$minimum
    best.aac  =  opt$objective
    # Calculate the inverse probability weights
    model.den$fitted  =  predict(model.den,newdata  =  mydata,
                                 n.trees  =  floor(best.aac.iter), type  =  "response")
    ps.den  =  dnorm((mydata$T-model.den$fitted)/sd(mydata$T-model.den
                                                    $fitted),0,1)
    weight.gbm  =  ps.num/ps.den
    
    ### model with estimated propensity score times reward/outcome as weigths
    
    
    #generate a new dataset for evaluating value function 
    set.seed(5000 + i)
    C_new = -10
    N_new = 1000
    beta_new = beta0
    mbeta_new = meanfunction
    p_new = length(beta_new) - 1
    X_new = matrix(rnorm(N_new * p,0,1),nrow = N_new)
    dose_optim = g(cbind(1,X_new) %*% beta_new)
    ### model with estimated propensity score times reward/outcome as weigths
    
    model_pen  =  svm(x  =  X[index,], y  =  A[index], w =  (trainweight * weight.gbm)[index], type = "eps-regression",
                      epsilon  =  0.15, scale = FALSE)
    tmpvalue2  =  predict(model_pen,X_new)
    tmpvalue2  =  pmin(pmax(tmpvalue2,0),1)
    Yhat2 = mbeta_new(X_new) + C_new*(tmpvalue2-dose_optim)^2
    kol_value_hat[i] = mean(Yhat2)
    write.csv(kol_value_hat, paste("results/simulation/kol-sim04-obs-value", samplesize, ".csv"),row.names = FALSE)
    
  }
  cat("\n Sample Size:",samplesize,"\n")
  cat("Mean of V:", mean(kol_value_hat),"\n")
  cat("Sd of V:", sd(kol_value_hat),"\n")
}