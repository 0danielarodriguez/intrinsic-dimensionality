Sigma_est_modelo <- function(datos,dd)# dim<=p-1; dim=p-1 es como no tener modelo. 
{
  p <- dim(datos)[2]
  Sigma_est <- cov(datos,use="pairwise.complete.obs")
  total.svd <- svd(Sigma_est)
  autovalores <- total.svd$d#autovalores.   eigen(Sigma_est)$values 
  autovectores <-   total.svd$v #autovectores.
  sigma_cuad_est <- mean(autovalores[(dd+1):p])
  diagonal_est <- diag(c(autovalores[1:dd], rep(sigma_cuad_est, (p-dd))  ))
  Sigma_est_modelo <- autovectores%*%diagonal_est%*%t(autovectores)
  Sigma_est_modelo 
}


dim_est_cross <- function(datos, gamma_max=1,potencia_dim=1, mm=c("cociente") ) #, cociente original promedio
  {
#spliteamos
  n <- dim(datos)[1]
  p <- dim(datos)[2]
  n_train <- round(0.8*n,0)
  datos_train<- datos[1:n_train,]
  datos_test<- datos[(n_train+1):n,]
  ########################################
#2 
  Sigma_test <- cov(datos_test,use="pairwise.complete.obs")
#3
  Sigma_train <- cov(datos_train,use="pairwise.complete.obs")
  total.svd <- svd(Sigma_train)
  lambda_est <- total.svd$d#autovalores.   eigen(Sigma_est)$values 
  autovectores <-   total.svd$v #autovectores.
#esto es para elegir el estimador
  if(mm=="cociente")
  {
  mm_hat <- lambda_est/lambda_est[1]
  }

  if(mm=="original")
  {
    mm_hat <- lambda_est
  }
  
#4
  grilla_gamma <- seq(0,gamma_max, by=0.01)
  lgrilla<-length(grilla_gamma)
  dd_gamma <- rep(NA,lgrilla ) #aca vamos a guardar los d(gamma)
 
for(i in 1:lgrilla)
{
  gamma_n <- grilla_gamma[i]
  #dividimos por el mas grande lambda_est
  pen <- mm_hat+gamma_n*seq(1:p)^potencia_dim #funcion a minimizar divido por l
  dd_gamma[i] <-  which.min(pen)-1
}
  
  perdida <- rep(NA, lgrilla) # loss 
  for(i in 1:lgrilla)
  {
    if(dd_gamma[i]>0)
    {
      perdida[i] <- sum((Sigma_test-Sigma_est_modelo(datos_train,dd_gamma[i]))^2) 
    }
    
  }
  
  gamma_final <- grilla_gamma[which.min(perdida)]
  C_est <- cov(datos)
  lambda_est <- eigen(C_est)$values 

  
  #aca hay que elegir el estimador 
  if(mm=="cociente")
  {
    mm_hat <- lambda_est/lambda_est[1]
  }
  
  if(mm=="original")
  {
    mm_hat <- lambda_est
  }
  
  
  pen <-   mm_hat+gamma_final*seq(1:p)^potencia_dim #funcion a minimizar divido por l
  dd_estimado <- which.min(pen)-1
  dd_estimado

  }



dim_est_cross_grupos <- function(datos,ngrupos=10, gamma_max=1,potencia_dim=1, mm=c("cociente") ) #, cociente original promedio
{
  #spliteamos
  n <- dim(datos)[1]
  p <- dim(datos)[2]
  
  
  grilla_gamma <- seq(0,gamma_max, by=0.01)
  lgrilla<-length(grilla_gamma)
  dd_gamma <- rep(NA,lgrilla ) #aca vamos a guardar los d(gamma)
  perdida <- rep(0, lgrilla) # loss 
  
  for(j in 1:ngrupos)
  {
    nobs<-n/ngrupos
    indices<-(j-1)*nobs+seq(1:nobs)
    datos_train<- datos[-indices,]
    datos_test<- datos[indices,]
    
    
  ########################################
  #2 
  Sigma_test <- cov(datos_test,use="pairwise.complete.obs")
  #3
  Sigma_train <- cov(datos_train,use="pairwise.complete.obs")
  total.svd <- svd(Sigma_train)
  lambda_est <- total.svd$d#autovalores.   eigen(Sigma_est)$values 
  autovectores <-   total.svd$v #autovectores.
  #esto es para elegir el estimador
  if(mm=="cociente")
  {
    mm_hat <- lambda_est/lambda_est[1]
  }
  
  if(mm=="original")
  {
    mm_hat <- lambda_est
  }
  
  #4
  
  for(i in 1:lgrilla)
  {
    gamma_n <- grilla_gamma[i]
    #dividimos por el mas grande lambda_est
    pen <- mm_hat+gamma_n*seq(1:p)^potencia_dim #funcion a minimizar divido por l
    dd_gamma[i] <-  which.min(pen)-1
  }
  
  for(i in 1:lgrilla)
  {
    if(dd_gamma[i]>0)
    {
      perdida[i] <- perdida[i]+sum((Sigma_test-Sigma_est_modelo(datos_train,dd_gamma[i]))^2) 
    }
    
  }
  
  }
  
  ######################################
  
  gamma_final <- grilla_gamma[which.min(perdida)]
  C_est <- cov(datos)
  lambda_est <- eigen(C_est)$values 
  
  
  #aca hay que elegir el estimador 
  if(mm=="cociente")
  {
    mm_hat <- lambda_est/lambda_est[1]
  }
  
  if(mm=="original")
  {
    mm_hat <- lambda_est
  }
  
  
  pen <-   mm_hat+gamma_final*seq(1:p)^potencia_dim #funcion a minimizar divido por l
  dd_estimado <- which.min(pen)-1
  dd_estimado
  
}

