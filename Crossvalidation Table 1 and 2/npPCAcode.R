library("RobStatTM")

#Compute the dimension estimator. It allows choosing between three loss functions (average, ratio, and identity).
# ratio is the default
#Split the sample into ngrupos numbers of folds (default =5)
#The penalty is a power, where the default is 1


npPCA <- function(datos, ngrupos=5, gamma_max=1,potencia_dim=1, mm=c("ratio") , tipo="cl") 
{
  n <- dim(datos)[1]
  p <- dim(datos)[2]
  
  grilla_gamma <- seq(0,gamma_max, by=0.01)
  lgrilla<-length(grilla_gamma)
  dd_gamma <- rep(NA,lgrilla )
  perdida <- rep(0, lgrilla)
  
  for(j in 1:ngrupos)
  {
    nobs<-round(n/ngrupos,0)
    indices<-(j-1)*nobs+seq(1:nobs)
    datos_train<- datos[-indices,]
    datos_test<- datos[indices,]
    Sigma_test_rob <- cov_rob_todas(datos_test, tipo=tipo)
    Sigma_train_rob<- cov_rob_todas(datos_train, tipo=tipo)
    total.svd <- svd(Sigma_train_rob)
    lambda_est <- total.svd$d
    autovectores <-   total.svd$v 
    
    if(mm=="ratio")  {mm_hat <- lambda_est/lambda_est[1]  }
    if(mm=="identity")  { mm_hat <- lambda_est  }
    if(mm=="average")  { mm_hat <- lambda_est/sum(lambda_est)}
    
    for(i in 1:lgrilla) {
      gamma_n <- grilla_gamma[i]
      pen <- mm_hat+gamma_n*seq(1:p)^potencia_dim
      dd_gamma[i] <-  which.min(pen)-1
    }
    
    for(i in 1:lgrilla) {
      perdida[i] <- perdida[i]+sum((Sigma_test_rob-Sigma_est_modelo_rob(datos_train,dd_gamma[i],tipo))^2) }
    
  }
  
  
  gamma_final <- grilla_gamma[which.min(perdida)]
  C_est <- cov_rob_todas(datos, tipo=tipo)
  lambda_est <- eigen(C_est)$values 
  
  if(mm=="ratio")  {mm_hat <- lambda_est/lambda_est[1]  }
  if(mm=="identity")  { mm_hat <- lambda_est  }
  if(mm=="average")  { mm_hat <- lambda_est/sum(lambda_est)}
  
  pen <-   mm_hat+gamma_final*seq(1:p)^potencia_dim 
  ddhat <- which.min(pen)-1
  sigmahat<-Sigma_est_modelos(datos, ddhat, tipo=tipo)
  result <- list(dimention= ddhat, cov= sigmahat)
  return(result)
}
  


# Estimator of the covariance matrix under the model
#Given an estimated d, estimate the covariance matrix
# we can choose classical estimators (default) of several options of robust procesures

Sigma_est_modelo <- function(datos,dd,tipo="cl")
{
  p <- dim(datos)[2]
  Sigma_est <- cov_rob_todas(datos, tipo=tipo)
  total.svd <- svd(Sigma_est)
  autovalores <- total.svd$d
  autovectores <-   total.svd$v
  sigma_cuad_est <- mean(autovalores[(dd+1):p])
  if(dd>0)  {diagonal_est <- diag(c(autovalores[1:dd], rep(sigma_cuad_est, (p-dd))  ))}
  if(dd==0) {diagonal_est <- diag( rep(sigma_cuad_est, p) ) }
  Sigma_est_modelo <- autovectores%*%diagonal_est%*%t(autovectores)
  Sigma_est_modelo 
}

#cl (default) de empirical covariance matrix
#covRob (is an MM estimator for p<10 and the Rocke proprosal for p>=10)
#fastmve
#initPP
cov_rob_todas <- function(datos, tipo="cl")
{
  if(tipo=="fastmve") { covhat <- fastmve(datos)$cov  }
  if(tipo=="initPP")  { covhat <- initPP(datos)$cova  }
  if(tipo=="covRob")  { covhat <- covRob(datos)$cov  }
  if(tipo=="cl")  { covhat <- cov(datos,use="pairwise.complete.obs")  }
  covhat
}



npPCAcross <- function(datos, ngrupos=5,  tipo="cl")
{
  n <- dim(datos)[1]
  p <- dim(datos)[2]
  dd <- 1:p
  perdida <- rep(0, p-1)
  for(j in 1:ngrupos)
  {
    nobs<-round(n/ngrupos,0)
    indices<-(j-1)*nobs+seq(1:nobs)
    datos_train<- datos[-indices,]
    datos_test<- datos[indices,]
    Sigma_test_rob <- cov_rob_todas(datos_test, tipo=tipo)
    Sigma_train_rob<- cov_rob_todas(datos_train, tipo=tipo)
    total.svd <- svd(Sigma_train_rob)
    lambda_est <- total.svd$d
    autovectores <-   total.svd$v
    for(i in 1:p) {
      perdida[i] <- perdida[i]+sum((Sigma_test_rob-Sigma_est_modelo(datos_train,i,tipo))^2)
    }
  }
  ddhat <- dd[which.min(perdida)]
  sigmahat<-Sigma_est_modelo(datos, ddhat, tipo=tipo)
  result <- list(dimention= ddhat, cov= sigmahat)
  return(result)
}
