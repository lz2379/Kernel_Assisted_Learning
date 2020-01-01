######################Discrete Q learning with real dataset################################
############################################################################
source("KAL-functions.r")


###################################################
######1. Run with all data
###################################################
Data <- read.csv("data-processed.csv")
NN = dim(Data)[1]

#variable names
continuous_var = c("Height")
discrete_var = c("Gender", "VKORC1.AG")
dose_var = "Dose_standardized"
outcome_var = "Y2"

#variable indexes
c_index = as.numeric(sapply(continuous_var,  FUN = function(x) which(names(Data) == x)))    #column index for continuous variables
d_index = as.numeric(sapply(discrete_var,  FUN = function(x) which(names(Data) == x)))  #column index for discrete variables
var_index = c(c_index, d_index)
dose_index = which(names(Data) == dose_var)
outcome_index = which(names(Data) == outcome_var)

#extract variables
Xc = as.matrix(Data[, c_index])
Xd = as.matrix(Data[, d_index])
X = Data[, var_index]
A = Data[, dose_index]
Y = Data[, outcome_index] # use quadratic outcome

#number of variables
p = length(var_index)
pc = length(c_index)
pd = length(d_index)

X2 = as.matrix(cbind(X, Xc^2))

#number of intervals
d = 10
Ad = floor(A * d)
Ad[(Ad == d)] = d-1

#Fit data using all data
lm_fit = lm(Y~X2 + as.factor(Ad):X2)
beta_est = lm_fit$coefficients

predY = matrix(rep(NA, d * NN), nrow = NN)
predY[, 1] = beta_est[1] + X2 %*% beta_est[2:5]
for (j in 2:10){
  predY[, j] = beta_est[1] + X2 %*% beta_est[2:5] + X2 %*% beta_est[(4 * j-2):(4 * j + 1)]
}
#Estimate the optimal dose
outdose = (apply(predY, 1, which.max)-0.5) / d
write.csv(outdose, "results/realdata/discreteQ-real-dose-Y2.csv", 
          row.names = FALSE)

####################################################
#######2. Estimate value function
####################################################

value = rep(NA, 200)
par2 = c(1.25, 1.75)
constX = par2[1:pc]
constA = par2[pc + 1]
expo = 4.5
hc = constX * apply(as.matrix(Xc), 2, sd) * NN^(-1 / expo)
ha = constA * sd(A) * NN^(-1 / expo)

testsize = floor(NN / 3)
trainsize = NN - testsize

value = rep(NA, 200)
for (i in 1:200){
  print(i)
  set.seed(10000 + i)
  train_index = sample(NN, size = trainsize)
  #split data into test and training
  train = Data[train_index, ]
  test = Data[-train_index, ]
  
  trainX = train[, var_index]
  trainY = train$Y
  trainA = train$Dose_standardized
  trainXc = train[, c_index]
  trainXd = train[, d_index]
  trainX2 = X2[train_index, ]
  trainAd = Ad[train_index]
  
  testX = test[, var_index]
  testY = test$Y
  testA = test$Dose_standardized
  testXc = test[, c_index]
  testXd = test[, d_index]
  testX2 = X2[-train_index, ]
  testAd = Ad[-train_index]
  
  #fit linear models
  lm_fit = lm(trainY~trainX2 + as.factor(trainAd):trainX2)
  beta_est = lm_fit$coefficients
  
  predY = matrix(rep(NA, d * testsize), nrow = testsize)
  predY[, 1] = beta_est[1] + testX2 %*% beta_est[2:5]
  for (j in 2:10){
    predY[, j] = beta_est[1] + testX2 %*% beta_est[2:5] + testX2 %*% beta_est[(4 * j-2):(4 * j + 1)]
  }
  outdose = (apply(predY, 1, which.max) - 0.5) / d
  
  #estimate value function with test dataset
  value[i] = testvalue(outdose = outdose, 
                       testY = testY, 
                       testA = testA, 
                       testXc = testXc, 
                       testXd = testXd, 
                       testsize = testsize, 
                       hc = hc, ha = ha)
  print(value[i])
  write.csv(value, "results/realdata/discreteQ-real-value-Y2.csv", 
            row.names = FALSE)
}
