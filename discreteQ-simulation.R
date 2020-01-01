source('KAl-functions.r')
source('discreteQ-function.R')

####################################
# Simulation settings with discretized Q learning
####################################


###################################
# Randomized trial Scenario 1-4
###################################
beta0 = c(0,  0.5) #beta0 = c(0, 1) for setting 2, 4
dose_range = 1
p = 1
d = 10
for (samplesize in c(400, 800)){
  set.seed(2000)
  numberofsets = 500
  result_01 = matrix(rep(NA, d * (2 * p + 1) * numberofsets), 
                     nrow = numberofsets)
  for (i in 1:numberofsets){
    
    #generate dataset
    dataset = generate_sample_random1(beta = beta0, 
                                      N = samplesize
                                      ) #mu0 = mu2 for setting 3, 4
    Y = dataset$Y
    X = dataset$X
    A = dataset$A
    N = dataset$N
    p = dataset$p
    #dummy variables for dose intervals
    Ad = floor(A * d/(dose_range))
    #quadratic terms
    X2 = X^2
    X_all = cbind(X, X2)
    colnames(X_all) = c(paste("X", c(1:p), sep = ""), paste("X", c(1:p), "_2", sep = ""))
    #fit the Q-learning model
    lm_fit = lm(Y ~ factor(Ad) + X_all + as.factor(Ad):(X_all))
    beta_est = lm_fit$coefficients
    result_01[i, ] = beta_est
  }
  cat("\n Sample Size:", samplesize, "\n")
  #calculate value function
  v = sim_value1(betahat = result_01, 
                 dose_function =  dose_discreteQ, 
                 beta = beta0, 
                 testN = 1000)
}

###################################
# Observational Studies setting1-4
###################################
beta0 = c(0,  1) #beta0 = c(0, 0.5)
dose_range = 1
d = 10
p = length(beta0)-1
for (samplesize in c(400, 800)){
  set.seed(2000)
  numberofsets = 500
  result_01 = matrix(rep(NA, d * (2 * p + 1) * numberofsets), nrow = numberofsets)
  for (i in 1:numberofsets){
    #print(i)
    dataset = generate_sample_obs1(beta = beta0, 
                                   N = samplesize, 
                                   mu0 = mu2)
    Y = dataset$Y
    X = dataset$X
    A = dataset$A
    N = dataset$N
    p = dataset$p
    #dummy variables for dose intervals
    Ad = floor(A * d / (dose_range))
    #quadratic terms
    X2 = X^2
    X_all = cbind(X, X2)
    colnames(X_all) = c(paste("X", c(1:p), sep = ""), paste("X", c(1:p), "_2", sep = ""))
    #fit the Q-learning model
    lm_fit = lm(Y ~ factor(Ad) + X_all + as.factor(Ad):(X_all))
    beta_est = lm_fit$coefficients
    result_01[i, ] = beta_est
  }
  cat("\n Sample Size:", samplesize, "\n")
  #value function
  v = sim_value1(betahat = result_01, 
                 dose_function =  dose_discreteQ, 
                 beta = beta0, 
                 testN = 1000, 
                 mu0 = mu2)
}


