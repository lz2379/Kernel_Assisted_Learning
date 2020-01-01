#####################
#1. estimate optimal treatment regime using discretized Q learning
#######################

#calculate the suggested dose given the estimated parameters 
dose_discreteQ = function(testX, 
                          beta_est, 
                          d = 10,  #number of intervals
                          dose_range = 1){
  q = dim(testX)[2]
  testX2 = testX^2
  testX_all = cbind(testX, testX2)
  testsize = dim(testX)[1]
  predY = matrix(rep(NA, d * testsize), nrow = testsize)
  predY[, 1] = beta_est[1] + testX_all %*% beta_est[(d + 1):(d + 2 * q)]
  for (j in 2:d){
    predY[, j] = beta_est[1] + beta_est[j] + 
      testX_all %*% beta_est[(d + 1):(d + 2 * q)] + 
      testX_all %*% beta_est[(d + (2 * q) * (j - 1) + 1):(d + (2 * q) * j)]
  }
  #find the interval which maximize the outcome
  outdose = (apply(predY, 1, which.max) - 0.5)/d
  
  return(outdose)
}