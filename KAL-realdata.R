#########Kernel Assisted Learning############################ 
#########Applied to the warfarin dataset#####################
##############################################


##################################################################
#######1. dose finding using all real data
####################################################################
source('KAl-functions.r')


Data = read.csv("data-processed.csv")

NN = dim(Data)[1]#sample size

#variable names
continuous_var = c("Height")
discrete_var = c("Gender", "VKORC1.AG")
dose_var = "Dose_standardized"
outcome_var = "Y2" #quadratic outcome

#index of variables
c_index = as.numeric(sapply(continuous_var,  FUN = function(x) which(names(Data) == x)))    #column index for continuous variables
d_index = as.numeric(sapply(discrete_var,  FUN = function(x) which(names(Data) == x)))  #column index for discrete variables
var_index = c(c_index, d_index)
dose_index = which(names(Data) == dose_var)
outcome_index = which(names(Data) == outcome_var)

#standardize the covariates in the dataset
Data[, var_index] = Data[, var_index] - matrix(rep(colMeans(Data[, var_index]), NN), nrow = NN, byrow = TRUE)
Data[, var_index] = Data[, var_index] / matrix(rep(apply(Data[, var_index], 2, sd), NN), nrow = NN, byrow = TRUE)

#extract variables
Xc = as.matrix(Data[, c_index])
Xd = as.matrix(Data[, d_index])
X = Data[, var_index]
A = Data[, dose_index]
Y = Data[, outcome_index] 

p = length(var_index)  
pc = length(c_index)   #number of continuous variables
pd = length(d_index)   #number of discrete variables

#setting bandwidths
par = c(1.25, 1.75)
constX = par[1:pc]     
constA = par[pc+1]
expo = 4.5
hc = constX*apply(as.matrix(Xc), 2, sd)*NN^(-1/expo)
ha = constA*sd(A)*NN^(-1/expo)

#find the range of the continuous variable
xlim = as.matrix(cbind(apply(Xc, 2, min), apply(Xc, 2, max)))

#inital points
beta0 = rbind(rep(0, p+1), diag(rep(1, p+1)))
scale = c(1, 1/apply(as.matrix(X), 2, sd))
beta0 = (beta0/10)%*%diag(scale) #beta0 = beta0%*%diag(scale/2)

results = finddose2(Y = Y, A = A, Xc = as.matrix(Xc), 
                  Xd = as.matrix(Xd), hc = hc, ha = ha, init = beta0, xlim = xlim, sampling = 4000, parscale = scale,   #parscale = scale*100
                  g = g, g1 = g1, g2 = g2)
beta_hat = as.numeric(results$beta)
#Here we got the solution
#c(0.14455516 ,  0.01001513 ,  0.51562184,  -1.74799745)

dose_hat = g(as.matrix(cbind(1, X)) %*% beta_hat)
#save.image("results/realdata-test/KA-real-Y2.Rdata")
write.csv(dose_hat, "results/realdata-test/KAL-real-dose-Y2.csv", row.names = FALSE)

      

########################################################
###########2. Estimate value function
########################################################

beta0 = t(beta_hat)
###Using one third as testing dataset, two thirds as training dataset
rep = 200
testsize = floor(NN/3)
trainsize = NN-testsize
boot2_beta = matrix(nrow = rep, ncol = (p + 1))
value_hat = rep(NA, rep)
boot2_results = list()
for (i in (1:rep)){
  print(i)
  set.seed(10000+i)
  train_index = sample(NN, size = trainsize)
  #create training dataset and test dataset
  train = Data[train_index, ]
  test = Data[-train_index, ]
  
  trainX = train[, var_index]
  trainY = train[, outcome_index]
  trainA = train[, dose_index]
  trainXc = train[, c_index]
  trainXd = train[, d_index]
  
  testX = test[, var_index]
  testY = test[, outcome_index]
  testA = test[, dose_index]
  testXc = test[, c_index]
  testXd = test[, d_index]
  
  #bandwidths
  hc = constX*apply(as.matrix(trainXc), 2, sd) * trainsize^(-1 / expo)
  ha = constA*sd(trainA)*trainsize^(-1 / expo)
  
  #dose finding
  boot2_results[[i]] = finddose2(
                         Y = trainY, 
                         A = trainA, 
                         Xc = as.matrix(trainXc), 
                         Xd = as.matrix(trainXd), 
                         hc = hc, ha = ha, 
                         sampling = 4000, 
                         xlim = xlim, 
                         init = beta0, 
                         parscale = scale, 
                         g = g, g1 = g1, g2 = g2, 
                         calc_cov = FALSE)
  boot2_beta[i, ] = as.numeric(boot2_results[[i]]$beta)
  
  #bandwidths for estimating value with test dataset
  hc = constX *apply(as.matrix(testXc), 2, sd) * testsize^(-1/expo)
  ha = constA * sd(testA) * testsize^(-1/expo)
  
  #estimate value function
  value_hat[i] = testvalue(betahat = boot2_beta[i, ], 
                           testY = testY, 
                           testA = testA, 
                           testXc = testXc, 
                           testXd = testXd, 
                           testsize = testsize, 
                           hc = hc, ha = ha)
  write.csv(value_hat, "results/realdata-test/KAL-real-value-Y2.csv", row.names = FALSE)
}



#########################################
#########################################

