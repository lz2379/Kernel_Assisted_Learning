###################### R script for Personalized dose finding: Part A, source funcions  ####
###################### using Outcome weighted learning         #############################
###################### Authors: Guanhua Chen, Donglin Zeng, Michael R. Kosorok #############
###################### Contact: g.chen@vanderbilt.edu          #############################


#### Load dependent packages ############
library(gbm)
library(polycor)
library(kernlab)
## The 'SVMW' package is for weighted supporting vector machines/regression, Dr. Holloway (NCSU) kindly 
## implement the modification of 'e1071' package and share the package with us.
library(SVMW)  
library(glmnet)
library(quantreg)
library(truncnorm)


### 1. Data Generating for the simulation data.

### 'size' is the sample size, 'ncov' is the covariate dimention, 'seed' is for random seed.
### 'R' is the reward, 'A' is the received dose. 
### Scenario 1, the optimal rule is linear in X.
Scenario1 <- function(size,ncov,seed){
  set.seed(seed)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  A = runif(size,0,2)
  D_opt = 1 + 0.5*X[,2] + 0.5*X[,1]
  mu =   8 + 4*X[,1] - 2*X[,2] - 2*X[,3] - 25*((D_opt-A)^2)
  R = rnorm(length(mu),mu,1)
  datainfo = list(X=X,A=A,R=R,D_opt=D_opt,mu=mu)
  return(datainfo)
}

### Scenario 2, the optimal rule is nonlinear in X.

Scenario2 <- function(size,ncov,seed){
  set.seed(seed)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  A = runif(size,0,2)
  D_opt = I(X[,1] > -0.5)*I(X[,1] < 0.5)*0.6 + 1.2*I(X[,1] > 0.5) + 1.2*I(X[,1] < -0.5) + 
    X[,4]^2 + 0.5*log(abs(X[,7])+1) - 0.6
  #D_opt = 1.1/(abs(X[,1])+1) + 0.8*X[,4]^2 + 0.9*log(abs(X[,7])+1) - 0.5 
  mu =   8 + 4*cos(2*pi*X[,2]) - 2*X[,4] - 8*X[,5]^3 - 15*abs(D_opt-A)
  #mu =   1 + X[,1] + 0.5*X[,2] + 8*exp(-abs(D_opt-A))
  #R = runif(length(mu),mu-1,mu+1)
  R = rnorm(length(mu),mu,1)
  datainfo = list(X=X,A=A,R=R,D_opt=D_opt,mu=mu)
  return(datainfo)
}


### Scenario 4 similar to Scenario 2, but it is a observational study
### that is A is dependent on X, and we do not know the exact relationship of A and X
### when train the data.


Scenario4 <- function(size,ncov,seed){
  set.seed(seed)
  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
  D_opt = I(X[,1] > -0.5)*I(X[,1] < 0.5)*0.6 + 1.2*I(X[,1] > 0.5) + 1.2*I(X[,1] < -0.5) + 
    X[,4]^2 + 0.5*log(abs(X[,7])+1) - 0.6
  A = rtruncnorm(size,a=0,b=2,mean=D_opt,sd=0.5)
  mu =   8 + 4*cos(2*pi*X[,2]) - 2*X[,4] - 8*X[,5]^3 - 15*abs(D_opt-A)
  R = rnorm(length(mu),mu,1)
  datainfo = list(X=X,A=A,R=R,D_opt=D_opt,mu=mu)
  return(datainfo)
}


### 2. For evaluating the value of estimated rule using testind data   ####
### Due to the different output model, L-O-Learning and K-O-leanring     ##
### use different prediction function                                    ##
### Scenario3 and Scenario4 can use similar prediction functions     ######

pred_s1 <- function(model,test){
  if(sum(grep('rq',model$call)) ==0){tmpvalue = predict(model,test$X)}
  if(sum(grep('rq',model$call)) ==1){tmpvalue = model$coefficients[1] + test$X %*% model$coefficients[-1]}
  ### to make sure the predicted dose is in the reasonable range
  pred = pmin(pmax(tmpvalue,0),2)
  pred_value = mean(8 + 4*test$X[,1] - 2*test$X[,2] - 2*test$X[,3] - 25*((test$D_opt-pred)^2))
  results = list(pred_dose=pred,pred_value=pred_value)
  return(results)
}

pred_s2 <- function(model,test){
  if(sum(grep('rq',model$call)) ==0){tmpvalue = predict(model,test$X)}
  if(sum(grep('rq',model$call)) ==1){tmpvalue = model$coefficients[1] + test$X %*% model$coefficients[-1]}
  pred = pmin(pmax(tmpvalue,0),2)
  pred_value = mean(8 + 4*test$X[,2] - 2*test$X[,4] - 8*test$X[,5]^3 - 15*abs(test$D_opt-pred))
  results = list(pred_dose=pred,pred_value=pred_value)
  return(results)
}


### 3. other utility functions                                      ####
###                                                                 ####

### create the design matrix for the lasso, which is main effect and pairwise
### interaction.

design <- function(train){
  mat = cbind(train$A,train$X)
  colnames(mat) = c("dose",paste("X",seq(1,ncol(train$X)),sep=""))
  int_mat=model.matrix(~.^2,data=data.frame(mat))
  sqmat = mat^2
  colnames(sqmat) = paste0(colnames(sqmat),":",colnames(sqmat),sep="")
  design_mat = cbind(int_mat[,-1],sqmat)
  return(design_mat)
}

### The terms which involve treatment and covariate interaction,
### and is used for determining the optimal dose for lasso, can 
### also use analytic solution for the lasso model in the paper.

lasso_fn <- function(d,coefs,X,ncov){
  obj = sum(coefs[1],coefs[(ncov+2):(2*ncov+1)]*X[2:(ncov+1)])*d+ coefs[2*ncov + (ncov-1)*ncov/2 + 2]*(d^2)  
  return(obj) 
}

### 4. Functions for propensity score estimation, Code is copied from the appendix of Zhu et .al (2015)

F.aac.iter = function(i,data,ps.model,ps.num,rep,criterion) {
  # i: number of iterations (trees)
  # data: dataset containing the treatment and the covariates
  # ps.model: the boosting model to estimate p(T_iX_i)
  # ps.num: the estimated p(T_i)
  # rep: number of replications in bootstrap
  # criterion: the correlation metric used as the stopping criterion
  GBM.fitted = predict(ps.model,newdata = data,n.trees = floor(i),type = "response")
  ps.den = dnorm((data$T - GBM.fitted)/sd(data$T-GBM.fitted),0,1)
  wt = ps.num/ps.den
  aac_iter = rep(NA,rep)
  for (i in 1:rep){
    bo = sample(1:dim(data)[1],replace = TRUE,prob = wt)
    newsample = data[bo,]
    j.drop = match(c("T"),names(data))
    j.drop = j.drop[!is.na(j.drop)]
    x = newsample[,-j.drop]
    if(criterion == "spearman"| criterion == "kendall"){
      ac = apply(x, MARGIN = 2, FUN = cor, y = newsample$T,
                 method = criterion)
    } else if (criterion == "distance"){
      ac = apply(x, MARGIN = 2, FUN = dcor, y = newsample$T)
    } else if (criterion == "pearson"){
      ac = matrix(NA,dim(x)[2],1)
      for (j in 1:dim(x)[2]){
        ac[j] = ifelse (!is.factor(x[,j]), cor(newsample$T, x[,j],
                                               method = criterion),polyserial(newsample$T, x[,j]))
      }
    } else print("The criterion is not correctly specified")
    aac_iter[i] = mean(abs(1/2*log((1+ac)/(1-ac))),na.rm = TRUE)
  }
  aac = mean(aac_iter)
  return(aac)
}


#### 5. Cross validation and get solution for K-O-Learning and L-O-Learning

### loss function for O-learning.
loss_O_learning <- function(model,test,epsilon){
  epsilon = 0.05
  Y = test$A
  R = test$R
  ### L-O-Learning use quantile regression, while K-O-Learning use SVM, hence
  ### the prediction is a little different.
  if(sum(grep('rq',model$call)) ==0){pred = predict(model,test$X)}
  if(sum(grep('rq',model$call)) ==1){pred = model$coefficients[1] + test$X %*% model$coefficients[-1]}
  mse = mean(R*(abs(Y - pred) -  eps_insen(Y - pred,epsilon))/epsilon)
  return(mse)
}

# epsilon_insensative loss
eps_insen <- function(z,epsilon){
  out = (z > epsilon)*(z - epsilon) + (z < -epsilon)*(-epsilon - z)
  return(out)
}


# Just for calculating the distance of two matrix, may be useful for gaussian bandwidth selection
distmat <- function(X,U){
  a = as.matrix(rowSums(X^2))     
  b = as.matrix(rowSums(U^2))
  one.a = matrix(1, ncol = nrow(b))     
  one.b = matrix(1, ncol = nrow(a))
  K1 = one.a %x% a
  K2 = X %*% t(U)
  K3 = t(one.b %x% b)
  K =  K1 - 2 * K2 + K3
  return(K)
}

### 5(a), for K-O-Learning

### "tunefunc" is a tuning function, which evaluate the performance of the methods,
### this function can be "loss_O_learning "
dc_solution_gaussian <- function(train,tune,epsilon,sigma,lambda,tunefunc){
  loop.index = 0
  curm = list(coefs=list(),epsilon=list(),sigma=list(),lambda=list(),error=list(),obj=list(),value=list())
  for(i in 1:length(epsilon)){
    for( j in 1:length(sigma)){
      #w=rnorm(ncol(train$X),0,4)
      index = which(train$R > quantile(train$R,0.65))
      tmpmodel = svm(x = train$X[index,], y = train$A[index], w= train$R[index], type="eps-regression", 
                     gamma = sigma[j], epsilon = 0.1, scale=FALSE)
      #tmpmodel = ksvm(train$X[index,],train$A[index],kernel="rbfdot",kpar=list(sigma=sigma[j]),scale=TRUE)
      for( k in 1:length(lambda)){
        loop.index = loop.index + 1
        model = dc_loop_gl(train$X,train$A,train$R,epsilon[i],sigma[j],lambda[k],tmpmodel)
        tune.error <- tunefunc(model,tune,epsilon[i])
        #value <- pred3_kernel(coefs,train,test,sigma[j],epsilon[i])
        curm$model[[loop.index]] = model
        curm$lambda[[loop.index]] =  lambda[k]
        curm$epsilon[[loop.index]] = epsilon[i]
        curm$sigma[[loop.index]] = sigma[j]
        curm$error[[loop.index]] = tune.error
        #curm$obj[[loop.index]] = obj
        print(paste("error=",tune.error,"epsilon=",epsilon[i],"sigma=",sigma[j],"lambda=",lambda[k]))
      }
    }
  }
  return(curm)
}


dose.cv.gaussian <- function(datall,epsilon,lambda,sigma=NULL,cvfold,tunefunc,optimfunc,seed){
  set.seed(seed)
  solution = list()
  result = list()
  error = 0
  max.step = 10
  nsize = length(datall$R)/cvfold
  alldist = as.vector(distmat(datall$X,datall$X))
  #sigma = 0.5*1/quantile(alldist[which(alldist>0)])
  ## empirical ways of choosing sigma, other ways are possible
  if(is.null(sigma)){
    sigma = 1/quantile(prob=c(0.25,0.5,0.75),alldist[which(alldist>0)])
  }
  for( i in 1:cvfold){
    print(paste("===== Fold ",i,"=====",sep=" "))
    index = seq(((i-1)*nsize+1),(i*nsize))
    train = list(X=datall$X[-index,],A=datall$A[-index],R=datall$R[-index],D_opt=datall$D_opt[-index])
    tune = list(X=datall$X[index,],A=datall$A[index],R=datall$R[index],D_opt=datall$D_opt[index])
    solution[[i]] <- optimfunc(train,tune,epsilon,sigma,lambda,tunefunc)
    error =  unlist(solution[[i]]$error) + error
  }
  finalindex = which.min(error)
  tmp.result <- optimfunc(datall,datall,unlist(solution[[1]]$epsilon)[finalindex],unlist(solution[[1]]$sigma)[finalindex],unlist(solution[[1]]$lambda)[finalindex],tunefunc)
  result$model = tmp.result$model[[1]]
  result$lambda = tmp.result$lambda[[1]]
  result$epsilon = tmp.result$epsilon[[1]]
  result$sigma = tmp.result$sigma[[1]]
  result$error = tmp.result$error[[1]]
  #result$obj = tmp.result$obj[[1]]
  return(result)
}


dc_loop_gaussian <- function(x,a,r,epsilon,sigma,C0,model){
  ### set up the maximum steps of iterations
  maxstep = 20
  steps = 0
  ### observations that contribute to the loss function that in the inital steps
  index_old = which(abs(a - predict(model,x))< epsilon)
  index = seq(1,length(a),1)
  # check whether the solution change between two steps, one can check the solution
  # or just the index
  while(sum(index_old) != sum(index) && steps < maxstep){
    steps = steps + 1
    index_old = which(abs(a - predict(rqmodel,x))< epsilon)
    if(length(index)<2) {break}
    rqmodel = svm(x = x[index,], y = a[index], w = r[index], type="eps-regression", gamma=sigma, cost = C0, epsilon = 0.1,
                  scale=FALSE)
    predtmp = predict(rqmodel,x)
    index = which(abs(a - predtmp) < epsilon)
  }
  pickedmodel = rqmodel
  return(pickedmodel)
}



### 5b. For L-O-Learning

dc_loop <- function(x,a,r,epsilon,C0,models){
  w_old = w_t
  b_old = b_t
  w = w_t + 1
  b = b_old + 1
  maxstep = 20
  steps = 0
  #while(sum((w-w_old)^2,(b - b_old)^2)/(length(w)+1) > 0.0001 && steps < maxstep){
  while(max((w-w_old)^2,(b - b_old)^2) > 0.0001 && steps < maxstep){
    steps = steps + 1
    index = which(abs(a - x %*% w_t - b_t) < epsilon)
    if(length(index)<2) {break}
    rqmodel = rq(a[index] ~ x[index,],.5,weights=r[index],method="lasso",lambda = C0)
    w = rqmodel$coefficients[-1]
    b = rqmodel$coefficients[1]
    w_old = w_t
    b_old = b_t
    w_t = w
    b_t = b
  }
  pickedmodel = rqmodel
  return(pickedmodel)
}


dc_solution <- function(train,tune,epsilon,lambda,tunefunc,pred){
  loop.index = 0
  curm = list(coefs=list(),epsilon=list(),lambda=list(),error=list(),obj=list(),value=list())
  for(j in 1:length(epsilon)){
    # The selection of 0.65 is arbitary, can be other quantile. 
    index = which(train$R > quantile(train$R,0.65))
    tmpmodel = cv.glmnet(train$X[index,],train$A[index])
    w=as.vector(coef(tmpmodel))[-1]
    b=as.vector(coef(tmpmodel))[1]
    for(k in 1:length(lambda)){
      loop.index = loop.index + 1
      model = dc_loop(train$X,train$A,train$R,epsilon[j],lambda[k],w_t=w,b_t=b)
      tune.error <- tunefunc(model,tune)
      value <- pred(coefs,tune)
      curm$coefs[[loop.index]] = model
      curm$lambda[[loop.index]] =  lambda[k]
      curm$epsilon[[loop.index]] = epsilon[j]
      curm$error[[loop.index]] = tune.error
      w = coefs$w
      b = coefs$b
      print(paste("error=",tune.error,"epsilon=",epsilon[j],"lambda=",lambda[k]))
    }    
  }
  return(curm)
}



dose.cv <- function(datall,epsilon,lambda,cvfold,tunefunc,optimfunc,seed){
  set.seed(seed)
  solution = list()
  result = list()
  error = 0
  max.step = 10
  nsize = length(datall$R)/cvfold
  for( i in 1:cvfold){
    print(paste("===== Fold ",i,"=====",sep=" "))
    index = seq(((i-1)*nsize+1),(i*nsize))
    train = list(X=datall$X[-index,],A=datall$A[-index],R=datall$R[-index],D_opt=datall$D_opt[-index])
    tune = list(X=datall$X[index,],A=datall$A[index],R=datall$R[index],D_opt=datall$D_opt[index])
    solution[[i]] <- optimfunc(train,tune,epsilon,lambda,tunefunc)
    error =  unlist(solution[[i]]$error) + error
  }
  finalindex = which.min(error)
  tmp.result <- optimfunc(datall,datall,unlist(solution[[1]]$epsilon)[finalindex],unlist(solution[[1]]$lambda)[finalindex],tunefunc)
  result$model = tmp.result$model[[1]]
  result$lambda = tmp.result$lambda[[1]]
  result$epsilon = tmp.result$epsilon[[1]]
  result$error = tmp.result$error[[1]]
  #result$obj = tmp.result$obj[[1]]
  return(result)
}
