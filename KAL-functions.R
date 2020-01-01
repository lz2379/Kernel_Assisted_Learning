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
g = function(x){
  return(1 / (1 + exp(-x)))
}
#first derivative of g()
g1 = function(x){
  return(exp(x) / ((1 + exp(x))^2))
}
#second derivative of g()
g2 = function(x){
  return((exp(x)-exp(2 * x)) / ((1 + exp(x))^3))
}

#Gaussian kernel function
K = function(x){
  return(exp(-x^2 / 2) / (sqrt(2 * pi)))
}
#first derivative of the kernel function
K1 = function(x){
  return((-x * exp(-x^2 / 2)) / (sqrt(2 * pi)))
}
#second derivative of the kernel function
K2 = function(x){
  return(((x^2-1) * exp(-x^2 / 2)) / sqrt(2 * pi))
}

#nonlinear baseline function
mu2 = function(x){
  return(1 + 0.5 * cos(2 * pi * x))
}



##################################################
### 2. Generating data
##################################################

#Randomized trials, settings 1-4
generate_sample_random1 = function(
                          C = -10,        
                          N = 500,        #sample size 
                          beta = c(0,1),  #parameters for optimal dosage
                          sd0 = 0.5, 
                          mu0 = function(x) 0 #baseline function
                          ){  
  p = length(beta)-1  #dimension of X
  A = runif(N)        #doses
  X = matrix(rnorm(N * p,0,1),nrow = N) #covariates
  Q = C * (A-g(cbind(1,X) %*% beta))^2
  mu = mu0(X) + Q       #E(Y|X,A)
  Y = rnorm(rep(1,N),mean = mu,sd = rep(sd0,N)) #outcome
  return(list(Y = Y,X = X,A = A,N = N,p = p))
}

#Observational study, setting 1-4 
generate_sample_obs1 = function(
                       C = -10,
                       N = 500,
                       beta = c(0,1),
                       sd0 = 0.5,
                       mu0 = function(x) 0
                       ){
  p = length(beta)-1
  X = matrix(rnorm(N * p,0,1),nrow = N)
  fopt = g(cbind(1,X) %*% beta) #optimal dosage
  #doses are related to fopt, thus related to X here
  A = rbeta(n = N,shape1 = 2 * exp(cbind(1,X) %*% beta),shape2 = 2)
  Q = C * (A-fopt)^2
  mu = mu0(X) + Q
  Y = rnorm(rep(1,N),mean = mu,sd = rep(sd0,N))
  return(list(Y = Y,X = X,A = A,N = N,p = p))
}


######################################
#Run KAL with simulated datasets
#Randomized trial, setting 1-4
KAL1 = function(
        C = -10, 
        N = 500, 
        beta = c(0,1), 
        sd0 = 0.5, 
        mu0 = function(x) 0 , 
        constX, constA, expo #parameters for kernel bandwidth
        ){
  dataset = generate_sample_random1(C = C,
                                    N = N,
                                    beta = beta,
                                    sd0 = sd0,
                                    mu0 = mu0)
  Y = dataset$Y
  X = dataset$X
  A = dataset$A
  N = dataset$N
  p = dataset$p
  #calculate bandwidth
  hhx = constX * sd(X) * sqrt(p) * N^(-1 / expo)
  hha = constA * sd(A) * N^(-1 / expo)
  return(finddose(Y, X, A, 
                  hx = hhx,
                  ha = hha,
                  sampling = 3000,
                  g = g, g1 = g1, g2 = g2,
                  init = diag(rep(1,p + 1)) #initial parameters 
                  )
         )
}


#Observational study, settings 1-4  
KAL_obs1 = function(
            C = -10,
            N = 500,
            beta = c(0,1),
            sd0 = 0.5,
            mu0 = function(x) 0 ,
            constX,constA,expo = 5
            ){
  dataset = generate_sample_obs1(C = C,
                                 N = N,
                                 beta = beta,
                                 sd0 = sd0,
                                 mu0 = mu0)
  Y = dataset$Y
  X = dataset$X
  A = dataset$A
  N = dataset$N
  p = dataset$p
  hhx = constX * sd(X) * sqrt(p) * N^(-1 / expo)
  hha = constA * sd(A) * N^(-1 / expo)
  return(finddose(Y, X, A,
                  hx = hhx,
                  ha = hha,
                  sampling = 3000,
                  g = g, g1 = g1,g2 = g2,
                  init = diag(rep(1,p + 1))
                  )
         )
}



##################################################
### 3. calculate the loss function to minimize
##################################################

calc_loss = function(
        beta, #the parameters
        Y,    #a vector of the observed outcomes
        A,    #a vector of the observed doses
        logKX = NA, #log of KX, 
        KX = NA,    #the calculated kernel function for covariates, 
        #only one of logKX and KX is needed
        CX,         #the kernel estimated marginal density function for covariates
        Xj1,  #the sampled grids for calculating integration in the covariate space
        ha,   #the bandwidth for doses
        N,    #number of observations
        sampling,   #number of sampled grids for calculating integrations
        link, #link function g()
        track = FALSE
        ){
  #calculate optimal dosage for grids Xj1
  gg = link(Xj1 %*% beta) 
  #calculate the distance between optimal dosage and the actual dose for dataset
  dA = matrix(rep(A,sampling)-rep(gg,each = N),ncol = sampling)
  
  if (sum(!is.na(logKX))>0){
    mhat = exp(-((dA / ha)^2 / 2)-logKX)
  } else{
    mhat = exp(-((dA / ha)^2 / 2)) * KX
  }
  loss = -mean(t(Y)%*%mhat / colSums(mhat) * CX)
  
  if (track) {
    cat("beta:", beta, "loss:", loss,"\n")
  }
  
  return(loss)
}

##################################################
###4.  Dose finding.  
##################################################

#Dose finding when all variables are continuous
finddose = function(Y, #vector of the observed outcomes 
                    X, #a matrix of covariates
                    A, #a vector of the observed treatments
                    hx = 0.1,  #bandwidth for X
                    ha = 0.03, #bandwidth for A
                    sampling = 10000, #the number of grid points used to estimate the integration
                    xlim = c(-3,3), #the range of X for integration
                    bootstrap = 0,  #number of bootstraps for calculating standard deviation
                    g,g1,g2, #g is the link function, g1 is the first derivative of g, g2 is the second derivative of g
                    init,  #the initial points
                    seed = NA #the random seed
                    ){
  N = dim(X)[1]
  p = dim(X)[2]
  X1 = cbind(1,X)
  
  #sample grid points for estimating the integration
  if (!is.na(seed)){
    set.seed(seed)
  }
  Xj = matrix(runif(p * sampling,min = xlim[1],max = xlim[2]),ncol = p)
  Xj1 = cbind(1,Xj)
  
  #apply bandwidthds for calculating kernel functions
  if (p>1){
    X0 = as.matrix(X)%*%diag(as.numeric(1 / hx))
    Xj0 = as.matrix(Xj)%*% diag(as.numeric(1 / hx))
  } else{
    X0 = as.matrix(X) / hx
    Xj0 = as.matrix(Xj) / hx
  }
  
  #kernel density estimation
  dX = as.matrix(pdist(X0,Xj0))
  logKX = dX^2 / 2
  CX = (colMeans(exp(-logKX)) / sqrt(2 * pi)) / exp(sum(log(hx)))
  
  #finding optimal parameters
  op = multistart(parmat = init, 
                  fn = calc_loss,  
                  Y = Y, 
                  A = A,  
                  logKX = logKX, 
                  CX = CX, 
                  Xj1 = Xj1, 
                  ha = ha, 
                  N = N, 
                  sampling = sampling, 
                  link = g)
  argmin = which.min(op$value)
  solution = as.numeric(op[argmin,1:(p + 1)])
  
  #calculate sd using bootstrap
  if (bootstrap > 0){
    bootstrap_sol = matrix(rep(NA,4 * bootstrap),nrow = bootstrap)
    for (i in 1:bootstrap){
      cat("bootstrap",i,"\n")
      #sampling bootstrap dataset
      index = sample(N, N, replace = TRUE)
      logKXb = logKX[index,]
      CXb = (colMeans(exp(-logKXb)) / sqrt(2 * pi)) / exp(sum(log(hx)))
      #finding optimal solution for bootstrap dataset
      op = optim(solution,
                 fn = calc_loss,
                 A = A[index],
                 logKX = logKX[index,], 
                 CX = CXb,
                 Xj1 = Xj1,
                 Y = Y[index],
                 ha = ha,
                 N = N,
                 sampling = sampling,
                 link = g)
      bootstrap_sol[i,] = as.numeric(op$par)
    }
    bootstrap_sd = apply(bootstrap_sol, 2, sd)
  } else{
    bootstrap_sd = NA
  }
  
  ##Calculate the estimated cov using formula
  beta_X = Xj1%*%solution
  gg = g(beta_X)
  gg1 = g1(beta_X)
  gg2 = g2(beta_X)
  dA = matrix(rep(gg, each = N) - rep(A, sampling), ncol = sampling)
  KX = K(dX) / exp(sum(log(hx)))
  KA = K(dA / ha) / ha
  KA1 = K1(dA / ha) / (ha^2)
  KA2 = K2(dA / ha) / (ha^3)
  mY = matrix(rep(Y, sampling),ncol = sampling)
  A_ = colMeans(mY * KX * KA)
  B_ = colMeans(KX * KA)
  #first derivative of A_ and B_ with respect to beta
  tA  = colMeans(mY * KX * KA1) * gg1 
  tB = colMeans(KX * KA1) * gg1  
  #second derivative
  ttA = colMeans(mY * KX * KA2) * (gg1^2) +  colMeans(mY * KX * KA1) * gg2 
  ttB = colMeans(KX * KA2) * (gg1^2) +  colMeans(KX * KA1) * gg2 
  CX = colMeans(KX)
  Dn = t(Xj1) %*% diag(as.numeric(  ((ttA / B_)-2 *  (tA * tB / (B_^2))-A_ * ttB / (B_^2) + 
                                   2 * (A_ * (tB^2) / (B_^3))) * CX   )) %*% (Xj1) / sampling * (xlim[2] - xlim[1])^p
  Dn_1 = solve(Dn) #inverse of Dn
  
  Phi = ((mY *  matrix(rep(gg1 * CX / B_,each = N),nrow = N)-    
        matrix(rep(gg1 * CX * A_ / (B_)^2,each = N),nrow = N)) * KX * KA1) %*% Xj1 / sampling * (xlim[2]-xlim[1])^p
  covS = cov(Phi) / N
  cov = Dn_1 %*% covS %*% Dn_1
  
  return(list(beta = solution,
              cov = cov,
              N = N,
              hx = hx,ha = ha,
              bootstrap_sd = bootstrap_sd))
  
}


#Indicator function of I(x == 0) for the purpose of stratifying categorial variables
I = function(x,h = 1){
  return((x == 0) * h + (x != 0) * (1-h))
}

#sample from discrete dataset
sample_discrete = function(Xd, #matrix of categorical covariates
                           sampling, #number of sampling grids
                           seed = 100){
  pd = length(Xd[1,])
  #write possible combinations of categories into single numebers
  base = 10^(c(1:pd) - 1) 
  combo = unique(as.matrix(Xd) %*% base)
  combo = combo[!is.na(combo)]
  combo = combo[order(combo)]
  #sample
  sp = sample(combo,size = sampling,replace = TRUE)
  out = matrix(NA,nrow = sampling,ncol = pd)
  #change samples into categorical vectors
  for (i in pd:1){
    ind = unique(Xd[,i])
    ind1 = matrix(rep(ind,sampling),ncol = length(ind),byrow = TRUE)
    temp = matrix(rep(sp / base[i],length(ind)),ncol = length(ind))
    temp2 =  apply( temp-ind1,1, FUN = function(x) which.min(abs(x)))
    out[,i] = ind[temp2]
    sp = sp - base[i] * out[,i]
  }
  return(out)
}

#Dose finding with data containing both continuous and discrete values
finddose2 = function(Y,A,
                     Xc, #continuous covariates
                     Xd, #discrete covariates
                     hc,hd = 1,ha = 0.03,
                     sampling = 2000,
                     xlim,
                     resample = 200,
                     init, 
                     parscale, #can add this to control step size when optimizing
                     seed = 100,
                     g,g1,g2, 
                     calc_cov = TRUE){
  X = cbind(Xc,Xd)
  N = dim(X)[1]
  p = dim(X)[2]
  pc = dim(Xc)[2] #number of continous covariates
  pd = dim(Xd)[2] #number of discrete covariates
  X1 = cbind(1,Xc,Xd)
  
  set.seed(seed)
  #sample grid for calculating integrals
  Xjc = matrix(runif(pc * sampling,min = xlim[,1],max = xlim[,2]),ncol = pc,byrow = TRUE)
  Xjd = sample_discrete(Xd,sampling)
  Xj = cbind(Xjc,Xjd)
  Xj1 = cbind(1,Xj)
  
  if (pc>1){
    Xc0 = as.matrix(Xc)%*%diag(as.numeric(1 / hc))
    Xjc0 = as.matrix(Xjc)%*% diag(as.numeric(1 / hc))
  } else{ 
    Xc0 = as.matrix(Xc) / hc
    Xjc0 = as.matrix(Xjc) / hc
  }
  
  dXc = as.matrix(pdist(as.matrix(Xc0),as.matrix(Xjc0)))
  
  KXc = K(dXc) / exp(sum(log(hc)))
  IXd = I(as.matrix(pdist(Xd,Xjd)))
  KX = KXc * IXd
  CX = colMeans(KX)
  
  #minimizing loss function
  op = multistart(parmat = init,fn = calc_loss,A = A,KX = KX, CX = CX,Xj1 = Xj1,#method = "BFGS",
                Y = Y,ha = ha,N = N,sampling = sampling,link = g, track = TRUE, control = list(parscale = parscale)
                )
  argmin = which.min(op$value)
  solution = as.numeric(op[argmin,1:(p + 1)])
  
  if (calc_cov == TRUE){
    
    ##Calculate the estimated cov
    beta_X = Xj1%*%solution
    gg = g(beta_X)
    gg1 = g1(beta_X)
    gg2 = g2(beta_X)
    dA = matrix(rep(gg,each = N)-rep(A,sampling),ncol = sampling)
    #KX = K(dX,d = p) / exp(sum(log(hx)))
    KA = K(dA / ha) / ha
    KA1 = K1(dA / ha) / (ha^2)
    KA2 = K2(dA / ha) / (ha^3)
    mY = matrix(rep(Y,sampling),ncol = sampling)
    
    A_ = colMeans(mY * KX * KA)
    B_ = colMeans(KX * KA)
    tA  = colMeans(mY * KX * KA1) * gg1 #first derivative with respect to beta
    tB = colMeans(KX * KA1) * gg1
    ttA = colMeans(mY * KX * KA2) * (gg1^2) +  colMeans(mY * KX * KA1) * gg2 #second derivative with respect to beta
    ttB = colMeans(KX * KA2) * (gg1^2) +  colMeans(KX * KA1) * gg2
    
    s_A_ = sign(A_)
    s_B_ = sign(B_)
    s_tA = sign(tA)
    s_tB = sign(tB)
    s_ttA = sign(ttA)
    s_ttB = sign(ttB)
    
    l_A_ = log(abs(A_))
    l_B_ = log(abs(B_))
    l_tA  = log(abs(tA)) #first derivative with respect to beta
    l_tB = log(abs(tB))
    l_ttA = log(abs(ttA)) #second derivative with respect to beta
    l_ttB = log(abs(ttB)) #second derivative with respect to beta
    #CX = colMeans(KX)
    sign1 = s_ttA * s_B_
    sign2 = s_tA * s_tB
    sign3 = s_A_ * s_ttB
    sign4 = s_A_ * s_B_
    Dn = t(Xj1)%*%diag(as.numeric(  (  sign1 * exp(l_ttA-l_B_)-  sign2 * 2 *  exp(l_tA + l_tB-2 * l_B_)  -  sign3 * exp(l_A_ + l_ttB-2 * (l_B_))   + 
                                       sign4 * 2 * exp(l_A_ + 2 * l_tB-3 * l_B_)                                                   ) * CX   ))%*%(Xj1) / sampling * (exp(sum(log(xlim[,2]-xlim[,1]))))
    Dn_1 = solve(Dn)
    
    s_C = sign(CX)
    l_C = log(abs(CX))
    sign5 = s_C * s_B_
    sign6 = s_C * s_A_
    Phi = ((mY *  matrix(rep(gg1 * sign5 * exp(l_C-l_B_),each = N),nrow = N)-    matrix(rep(gg1 * sign6 * exp(l_C + l_A_-2 * l_B_),each = N),nrow = N)) * KX * KA1) %*% Xj1 / sampling * (exp(sum(log(xlim[,2]-xlim[,1]))))
    covS = cov(Phi) / N
    cov = Dn_1%*%covS%*%Dn_1
  } else{
    cov = NA
  }
  return(list(beta = solution,cov = cov,N = N,hc = hc,ha = ha,allresult = op))
}


##################################################
#### 5. calculate value function of the suggested dose 
##################################################

dose_KAL = function(testX,
                    betahat
                    ){
  return( g(cbind(1,testX) %*% betahat) )
}
##for simulation settings 1-4
sim_value1 = function(betahat = NA,  
                      dose_function = dose_KAL, #dose_suggested = NA,
                      C = -10, 
                      testN = 1000, 
                      beta = c(0,1), 
                      sd0 = 0.5, 
                      mu0 = function(x) 0 ,
                      seed = 5000
                      ){
  numberofsets = dim(betahat)[1]
  V = rep(NA,numberofsets)
  p = length(beta)-1
  for (i in 1:numberofsets){
      set.seed(seed + i)
      X = matrix(rnorm(testN * p,0,1),nrow = testN)
      dose_optim = g(cbind(1,X)%*% beta)
      V[i] = mean( mu0(X) + C * (dose_optim-dose_function(X,t(betahat[i,])))^2)
  }
  cat("Mean of V:", mean(V),"\n")
  cat("Sd of V:", sd(V),"\n")
  return(V)
}

###
# for real data
# estimate the value function of the suggested dose on a test dataset
###


testvalue = function(betahat = NA, #this option is for dose with form g(beta*X)
                     outdose = NA, #suggested dose, this option is for dose without form g(beta*X)
                     testY,testA,
                     testXc,testXd,
                     testsize,
                     hc,ha){
  
  if (sum(!is.na(outdose))>0){
    predictA = outdose
  } else{
    testX1 = cbind(1,testXc,testXd)
    predictA = g(as.matrix(testX1) %*% betahat)
  }
  Aj = matrix(rep(testA,each = testsize),nrow = testsize)
  
  #Estimate E(Y|A=suggested dose, X) using kernel estimation
  predictAi = matrix(rep(predictA,testsize),nrow = testsize)
  Ka = as.matrix(K((Aj - predictAi) / ha)) / ha
  if (length(hc)>1){
    dXc = as.matrix(dist(as.matrix(testXc) %*% diag(as.numeric(1 / hc))))
  } else {
    dXc = as.matrix(dist(as.matrix(testXc) / hc))
  }
  Kxc = K(dXc) / exp(sum(log(hc)))
  Ixd = I(as.matrix(dist(testXd)))
  Kx = Kxc * Ixd
  
  Yj = matrix(rep(testY,each = testsize),nrow = testsize)
  return(mean(rowMeans(Yj * Ka * Kx) / rowMeans(Ka * Kx)))
}

