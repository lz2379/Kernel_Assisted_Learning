source('KAl-functions.r')

solution0=c(1,0.5,0.5,0)

set.seed(1000)
numberofsets=500
for (samplesize in c(400,800)){
  result_05=matrix(rep(NA,4*numberofsets),nrow=numberofsets)
  cov_05=matrix(rep(NA,4*numberofsets),nrow=numberofsets)
  cov2_05=matrix(rep(NA,4*numberofsets),nrow=numberofsets)
  results_05=list()
  d=10
  value=rep(NA,numberofsets)
  for (i in 1:numberofsets){
    print(i)
    set.seed(1000+i)
    C=-25
    N=800
    beta=c(1,0.5,0.5,0)
    mbeta=c(8,4,-2,-2)
    sd0=1
    p=length(beta)-1
    A=runif(N,min=0,max=2)
    X=matrix(runif(N*p,min=-1,max=1),nrow=N)
    Q=C*(A-cbind(1,X)%*% beta)^2
    mu=Q+cbind(1,X)%*% mbeta
    Y=rnorm(rep(1,N),mean=mu,sd=rep(sd0,N))
    Ad=floor(A*d/2)
    X2=X^2
    X=cbind(X,X2)
    colnames(X)=c("X1","X2","X3","X1_2","X2_2","X3_2")
    lm_fit=lm(Y~X+as.factor(Ad):X)
    beta_est=lm_fit$coefficients
    
    set.seed(5000+i)
    C_new=-25
    N_new=1000
    beta_new=c(1,0.5,0.5,0)
    mbeta_new=c(8,4,-2,-2)
    p_new=length(beta_new)-1
    A_new=runif(N_new,min=0,max=2)
    X_new=matrix(runif(N_new*p_new,min=-1,max=1),nrow=N_new)
    X_new2=cbind(X_new,X_new^2)
    
    predY=matrix(rep(NA,d*N_new),nrow=N_new)
    predY[,1]=beta_est[1]+X_new2%*%beta_est[2:7]
    for (j in 2:10){
      predY[,j]=beta_est[1]+X_new2%*%beta_est[2:7]+X_new2%*%beta_est[(6*j-4):(6*j+1)]
    }
    outdose=(apply(predY,1,which.max)-0.5)/5
    Yhat=C*(outdose-cbind(1,X_new)%*% beta_new)^2+cbind(1,X_new)%*% mbeta
    value[i]=mean(Yhat)
    write.csv(value,paste("results/discreteQ-sim05-value", samplesize, ".csv"),row.names=FALSE)
    }
}


