##################################################
##### Simulation Examples ##############
##################################################

##################################################
#Randomized trials #Scenario 1-4
##################################################
#-----------------
#Scenario 1
#-----------------

source('KAL-functions.r')

expo = 4.5
constX = 1.25
constA = 1.75
beta0 = c(0, 0.5)

for (samplesize in c(400,800)){
  set.seed(2000)
  numberofsets = 500
  result_01 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  cov_01 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  results_01 = list()
  for (i in 1:numberofsets){
    print(i)
    #generate dataset and calculate parameters
    result = KAL1(beta = beta0,
                  N = samplesize,
                  constX = constX,
                  constA = constA,
                  expo = expo)
    result_01[i, ] = result$beta
    cov_01[i, 1] = result$cov[1, 1]
    cov_01[i, 2] = result$cov[2, 2]
    results_01[[i]] = result
    #output
    #write.csv(result_01, paste("results/simulation/beta_01", samplesize, ".csv"))
    #write.csv(cov_01, paste("results/simulation/cov_01", samplesize, ".csv"))
    #save.image(paste("results/simulation/result_01", samplesize, ".RData"))
  }
  
  colMeans(result_01) #mean of the estimated parameter
  apply(result_01, 2, sd) # sd of the estimated parameters
  colMeans(sqrt(cov_01))  # mean of estimated sd
  
  #Confidence Intervals
  lowerCI = result_01 - sqrt(cov_01) * 1.96
  upperCI = result_01 + sqrt(cov_01) * 1.96
  
  #coverage of confidence intervals
  mean(lowerCI[, 1] * upperCI[, 1] < 0)
  mean((lowerCI[, 2] - 0.5) * (upperCI[, 2] - 0.5) < 0)
}

#calculate value function
beta_400 = read.csv("results/simulation/beta_01 400 .csv")[, 2:3]
beta_800 = read.csv("results/simulation/beta_01 800 .csv")[, 2:3]
v = sim_value1(betahat = beta_400,beta = beta0,testN = 1000)
v = sim_value1(betahat = beta_800,beta = beta0,testN = 1000)

#-----------------
#Scenario 2
#-----------------
beta0 = c(0,1)
for (samplesize in c(400,800)){
  set.seed(2000)
  numberofsets = 500
  result_02 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  cov_02 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  results_02 = list()
  for (i in 1:numberofsets){
    print(i)
    result = KAL1(beta = beta0,
                  N = samplesize,
                  constX = constX,
                  constA = constA,
                  expo = expo)
    result_02[i,] = result$beta
    cov_02[i,1] = result$cov[1,1]
    cov_02[i,2] = result$cov[2,2]
    results_02[[i]] = result
    write.csv(result_02, paste("results/simulation/beta_02", samplesize, ".csv"))
    write.csv(cov_02,paste("results/simulation/cov_02", samplesize, ".csv"))
    #save.image(paste("results/simulation/result_02", samplesize, ".RData"))
  }
  
  colMeans(result_02) #mean of the estimated parameter
  apply(result_02,2,sd) # sd of the estimated parameters
  colMeans(sqrt(cov_02))  # mean of estimated sd
  
  #Confidence Intervals
  lowerCI = result_02 - sqrt(cov_02) * 1.96
  upperCI = result_02 + sqrt(cov_02) * 1.96
  
  #coverage of confidence intervals
  mean(lowerCI[,1] * upperCI[,1] < 0)
  mean((lowerCI[,2] - 1) * (upperCI[,2] - 1) < 0)
}
#calculate value function
beta_400 = read.csv("results/simulation/beta_02 400 .csv")[,2:3]
beta_800 = read.csv("results/simulation/beta_02 800 .csv")[,2:3]
v = sim_value1(betahat = beta_400,beta = beta0,testN = 1000)
v = sim_value1(betahat = beta_800,beta = beta0,testN = 1000)


#-----------------
#Scenario 3
#-----------------
beta0 = c(0, 0.5)
for (samplesize in c(400, 800)){
  set.seed(2000)
  numberofsets = 500
  result_03 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  cov_03 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  results_03 = list()
  for (i in 1:numberofsets){
    print(i)
    result = KAL1(beta = beta0,
                  mu0 = mu2, 
                  N = samplesize,
                  constX = constX,
                  constA = constA,
                  expo = expo)
    result_03[i,] = result$beta
    cov_03[i,1] = result$cov[1,1]
    cov_03[i,2] = result$cov[2,2]
    results_03[[i]] = result
    write.csv(result_03, paste("results/simulation/beta_03", samplesize, ".csv"))
    write.csv(cov_03,paste("results/simulation/cov_03", samplesize, ".csv"))
    #save.image(paste("results/simulation/result_03", samplesize, ".RData"))
  }
  
  colMeans(result_03) #mean of the estimated parameter
  apply(result_03,2,sd) # sd of the estimated parameters
  colMeans(sqrt(cov_03))  # mean of estimated sd
  
  #Confidence Intervals
  lowerCI = result_03 - sqrt(cov_03) * 1.96
  upperCI = result_03 + sqrt(cov_03) * 1.96
  
  #coverage of confidence intervals
  mean(lowerCI[,1] * upperCI[,1] < 0)
  mean((lowerCI[,2] - 0.5) * (upperCI[,2] - 0.5) < 0)
}

beta_400 = read.csv("results/simulation/beta_03 400 .csv")[,2:3]
beta_800 = read.csv("results/simulation/beta_03 800 .csv")[,2:3]
v = sim_value1(betahat = beta_400,mu0 = mu2,beta = beta0,testN = 1000)
v = sim_value1(betahat = beta_800,mu0 = mu2,beta = beta0,testN = 1000)


#-----------------
#Scenario 4
#-----------------
beta0 = c(0,1)
for (samplesize in c(400,800)){
  set.seed(2000)
  numberofsets = 500
  result_04 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  cov_04 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  results_04 = list()
  for (i in 1:numberofsets){
    print(i)
    result = KAL1(beta = beta0, 
                  mu0 = mu2, 
                  N = samplesize, 
                  constX = constX,
                  constA = constA,
                  expo = expo)
    result_04[i,] = result$beta
    cov_04[i,1] = result$cov[1,1]
    cov_04[i,2] = result$cov[2,2]
    results_04[[i]] = result
    write.csv(result_04, paste("results/simulation/beta_04", samplesize, ".csv"))
    write.csv(cov_04,paste("results/simulation/cov_04", samplesize, ".csv"))
    #save.image(paste("results/simulation/result_04", samplesize, ".RData"))
  }
  
  colMeans(result_04) #mean of the estimated parameter
  apply(result_04,2,sd) # sd of the estimated parameters
  colMeans(sqrt(cov_04))  # mean of estimated sd
  
  #Confidence Intervals
  lowerCI = result_04 - sqrt(cov_04) * 1.96
  upperCI = result_04 + sqrt(cov_04) * 1.96
  
  #coverage of confidence intervals
  mean(lowerCI[,1] * upperCI[,1] < 0)
  mean((lowerCI[,2] - 1) * (upperCI[,2] - 1) < 0)
}

beta_400 = read.csv("results/simulation/beta_04 400 .csv")[,2:3]
beta_800 = read.csv("results/simulation/beta_04 800 .csv")[,2:3]
v = sim_value1(betahat = beta_400,mu0 = mu2,beta = beta0,testN = 1000)
v = sim_value1(betahat = beta_800,mu0 = mu2,beta = beta0,testN = 1000)



##################################################
#############Observational studies #Scenario 1-4
###################################################

#-----------------
#Scenario 1
#-----------------
expo = 4.5
constX = 0.8
constA = 3.2
beta0 = c(0,0.5)
for (samplesize in c(400,800)){
  set.seed(2000)
  numberofsets = 500
  result_01 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  cov_01 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  results_01 = list()
  
  for (i in 1:numberofsets){
    print(i)
    result = KAL_obs1(beta = c(0,0.5),
                      N = samplesize,
                      constX = constX,
                      constA = constA,
                      expo = expo)
    result_01[i,] = result$beta
    cov_01[i,1] = result$cov[1,1]
    cov_01[i,2] = result$cov[2,2]
    results_01[[i]] = result
    write.csv(result_01, paste("results/simulation/beta_01_obs",samplesize,".csv"))
    write.csv(cov_01, paste("results/simulation/cov_01_obs", samplesize,".csv"))
    #save.image(paste("results/simulation/result_01_obs", samplesize, ".RData"))
  }
  
  colMeans(result_01)
  apply(result_01,2,sd)
  colMeans(sqrt(cov_01))
  
  #confidence interval
  lowerCI = result_01 - sqrt(cov_01) * 1.96
  upperCI = result_01 + sqrt(cov_01) * 1.96
  mean(lowerCI[,1] * upperCI[,1] < 0)
  mean((lowerCI[,2] - 0.5) * (upperCI[,2] - 0.5) < 0)
}

#value function
beta_400 = read.csv("results/simulation/beta_01_obs 400 .csv")[,2:3]
beta_800 = read.csv("results/simulation/beta_01_obs 800 .csv")[,2:3]
v = sim_value1(betahat = beta_400,beta = beta0,testN = 1000)
v = sim_value1(betahat = beta_800,beta = beta0,testN = 1000)

result_01 = read.csv("results/simulation/beta_01_obs 400 .csv")[,2:3]
cov_01 = read.csv("results/simulation/cov_01_obs 400 .csv")[,2:3]

result_01 = read.csv("results/simulation/beta_01_obs 800 .csv")[,2:3]
cov_01 = read.csv("results/simulation/cov_01_obs 800 .csv")[,2:3]

#-----------------
#Scenario 2
#-----------------

expo = 4.5
constX = 0.9
constA = 2.35
beta0 = c(0,1)

for (samplesize in c(400,800)){
  set.seed(2000)
  numberofsets = 500
  result_02 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  cov_02 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  results_02 = list()
  for (i in 1:numberofsets){
    print(i)
    result = KAL_obs1(beta = c(0,1),
                      N = samplesize,
                      constX = constX,
                      constA = constA,
                      expo = expo)
    result_02[i,] = result$beta
    cov_02[i,1] = result$cov[1,1]
    cov_02[i,2] = result$cov[2,2]
    results_02[[i]] = result
    write.csv(result_02,paste("results/simulation/beta_02_obs",samplesize,".csv"))
    write.csv(cov_02,paste("results/simulation/cov_02_obs", samplesize,".csv"))
    #save.image(paste("results/simulation/result_02_obs", samplesize, ".RData"))
  }
  #summary
  colMeans(result_02)
  apply(result_02,2,sd)
  sqrt(colMeans(cov_02))
  
  #confidence interval
  lowerCI = result_02 - sqrt(cov_02) * 1.96
  upperCI = result_02 + sqrt(cov_02) * 1.96
  mean(lowerCI[,1] * upperCI[,1] < 0)
  mean((lowerCI[,2] - 1) * (upperCI[,2] - 1) < 0)
}

#value function
beta_400 = read.csv("results/simulation/beta_02_obs 400 .csv")[,2:3]
beta_800 = read.csv("results/simulation/beta_02_obs 800 .csv")[,2:3]
v = sim_value1(betahat = beta_400, beta = beta0, testN = 1000)
v = sim_value1(betahat = beta_800, beta = beta0, testN = 1000)

#-----------------
#Scenario 3
#-----------------
expo = 4.5
constX = 0.8
constA = 3.2
beta0 = c(0,0.5)
for (samplesize in c(400,800)){
  set.seed(2000)
  numberofsets = 500
  result_03 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  #cov_03 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  cov_03 = matrix(rep(NA,2 * numberofsets),nrow = numberofsets)
  results_03 = list()
  for (i in 1:numberofsets){
    print(i)
    result = KAL_obs1(beta = c(0,0.5),
                      mu0 = mu2,
                      N = samplesize,
                      constX = constX,
                      constA = constA,
                      expo = expo)
    result_03[i,] = result$beta
    cov_03[i,1] = result$cov[1,1]
    cov_03[i,2] = result$cov[2,2]
    results_03[[i]] = result
    write.csv(result_03,paste("results/simulation - test/beta_03_obs",samplesize,constX,constA,".csv"))
    write.csv(cov_03,paste("results/simulation - test/cov_03_obs", samplesize,constX,constA,".csv"))
    #save.image(paste("results/simulation - test/result_03_obs", samplesize, constX,constA,".RData"))
  }
  #summary
  colMeans(result_03)
  apply(result_03,2,sd)
  colMeans(sqrt(cov_03))
  
  #CI
  lowerCI = result_03 - sqrt(cov_03) * 1.96
  upperCI = result_03 + sqrt(cov_03) * 1.96
  mean(lowerCI[,1] * upperCI[,1] < 0)
  mean((lowerCI[,2] - 0.5) * (upperCI[,2] - 0.5) < 0)
}

#Value
beta_400 = read.csv("results/simulation - test/beta_03_obs 400 0.8 3.2 .csv")[,2:3]
beta_800 = read.csv("results/simulation - test/beta_03_obs 800 0.8 3.2 .csv")[,2:3]
v = sim_value1(betahat = beta_400,beta = beta0,mu0 = mu2,testN = 1000)
v = sim_value1(betahat = beta_800,beta = beta0,mu0 = mu2,testN = 1000)

#-----------------
#Scenario 4
#-----------------

expo = 4.5
constX = 0.9
constA = 2.35
beta0 = c(0,1)
for (samplesize in c(400, 800)){
  set.seed(2000)
  numberofsets = 500
  result_04 = matrix(rep(NA, 2 * numberofsets), nrow = numberofsets)
  cov_04 = matrix(rep(NA, 2 * numberofsets), nrow = numberofsets)
  results_04 = list()
  for (i in 1:numberofsets){
    print(i)
    result = KAL_obs1(beta = c(0, 1),
                      mu0 = mu2,
                      N = samplesize,
                      constX = constX,
                      constA = constA,
                      expo = expo)
    result_04[i,] = result$beta
    cov_04[i, 1] = result$cov[1, 1]
    cov_04[i, 2] = result$cov[2, 2]
    results_04[[i]] = result
    write.csv(result_04,paste("results/simulation/beta_04_obs",samplesize,".csv"))
    write.csv(cov_04,paste("results/simulation/cov_04_obs", samplesize,".csv"))
    #save.image(paste("results/simulation/result_04_obs", samplesize, ".RData"))
  }
  colMeans(result_04)
  apply(result_04, 2, sd)
  sqrt(colMeans(cov_04))

  #CI
  lowerCI = result_04 - sqrt(cov_04) * 1.96
  upperCI = result_04 + sqrt(cov_04) * 1.96
  mean(lowerCI[,1] * upperCI[,1] < 0)
  mean((lowerCI[,2] - 1) * (upperCI[,2] - 1) < 0)
}

#value
beta_400 = read.csv("results/simulation/beta_04_obs 400 .csv")[, 2:3]
beta_800 = read.csv("results/simulation/beta_04_obs 800 .csv")[, 2:3]
v = sim_value1(betahat = beta_400,beta = beta0,mu0 = mu2,testN = 1000)
v = sim_value1(betahat = beta_800,beta = beta0,mu0 = mu2,testN = 1000)


