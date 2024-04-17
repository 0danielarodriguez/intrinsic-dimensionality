library("RobStatTM")

#Generacion de datos
vectores_canonicos <- function(p)
{
  canonicos <- vector("list", p)
  for(j in 1:p)
  {
    ee <- rep(0,p)
    ee[j] <- 1
    canonicos[[j]] <- ee
  }
  canonicos
}


genero_datos_SNR <- function(nsamp,p=50,d=20,SNR,sigma=1)
{
  alpha <- SNR*(p-d)*sigma^2/d
  lambda_generacion <-rep(alpha, d)-sigma^2 # descuento el ruido - o pongo directo el input sigma 
  canonicos <- vectores_canonicos(p)
  W <- matrix(NA, p,d) #matriz para mandar a R^p
  for(i in 1:d)
  {
    W[,i] <- sqrt(lambda_generacion[i])*canonicos[[i]]
  }
  
  C <- W%*%t(W)+diag(rep(sigma^2,p),p,p)
  #eigen(C)$values / verificacion interna
  
  #generacion de datos
  vv <-matrix(rnorm(d*nsamp,0,1), nsamp,d)
  
  ee <- matrix(rnorm(p*nsamp,0,sigma), nsamp,p)
  
  xx <- t(W%*%t(vv))+ee
  
  salida <- xx  
  salida
  # C_est <- cov(xx)
  # lambda_est <- eigen(C_est)$values
  # lambda_est 
}


Sigma_est_modelo <- function(datos,dd)# dim<=p-1; dim=p-1 es como no tener modelo. 
{
  p <- dim(datos)[2]
  Sigma_est <- cov(datos)
  total.svd <- svd(Sigma_est)
  autovalores <- total.svd$d#autovalores.   eigen(Sigma_est)$values 
  
  autovectores <-   total.svd$v #autovectores.
  
  sigma_cuad_est <- mean(autovalores[(dd+1):p])
  
  if(dd>0)
  {
    diagonal_est <- diag(c(autovalores[1:dd], rep(sigma_cuad_est, (p-dd))  ))
  }
  
  if(dd==0)
  {
    diagonal_est <- diag( rep(sigma_cuad_est, p) )
  }
  
  Sigma_est_modelo <- autovectores%*%diagonal_est%*%t(autovectores)
  Sigma_est_modelo 
  
}


dim_est_cross <- function(datos, gamma_max=1,paso_grilla=0.01 , potencia_dim=1, mm=c("cociente") ) #, cociente original promedio
{
  #spliteamos
  n <- dim(datos)[1]
  p <- dim(datos)[2]
  n_train <- round(0.8*n,0)
  datos_train<- datos[1:n_train,]
  datos_test<- datos[(n_train+1):n,]
  ########################################
  #2
  Sigma_test <- cov(datos_test)
  
  #3
  Sigma_train <- cov(datos_train)
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
  
  if(mm=="promedio")
  {
    mm_hat <- lambda_est/sum(lambda_est)
  }
  
  
  
  
  #4
  grilla_gamma <- seq(0,gamma_max, by=paso_grilla)
  dd_gamma <- rep(NA, length(grilla_gamma)) #aca vamos a guardar los d(gamma)
  
  for(i in 1:length(grilla_gamma))
  {gamma_n <- grilla_gamma[i]
  #dividimos por el mas grande lambda_est
  pen <- mm_hat+gamma_n*seq(1:p)^potencia_dim #funcion a minimizar divido por l
  dd_gamma[i] <-  which.min(pen)-1
  
  }
  
  
  perdida <- rep(0, length(grilla_gamma)) # loss 
  for(i in 1:length(grilla_gamma))
  {
    
    perdida[i] <- sum((  Sigma_test-Sigma_est_modelo(datos_train,dd_gamma[i]))^2) 
    
    
  }
  
  gamma_final <- grilla_gamma[which.min(perdida)]
  C_est <- cov(datos)
  lambda_est <- eigen(C_est)$values 
  
  
  #esto es para elegir el estimador
  if(mm=="cociente")
  {
    mm_hat <- lambda_est/lambda_est[1]
  }
  
  if(mm=="original")
  {
    mm_hat <- lambda_est
  }
  
  if(mm=="promedio")
  {
    mm_hat <- lambda_est/sum(lambda_est)
  }
  
  
  
  pen <-   mm_hat+gamma_final*seq(1:p)^potencia_dim #funcion a minimizar divido por l
  dd_estimado <- which.min(pen)-1
  salida <- dd_estimado
  salida
}

#este devuelve el d del gamma que minimiza la perdida. 
dim_est_cross_bis <- function(datos, gamma_max=1,paso_grilla=0.01 , potencia_dim=1, mm=c("cociente") ) #, cociente original promedio
{
  #spliteamos
  n <- dim(datos)[1]
  p <- dim(datos)[2]
  n_train <- round(0.8*n,0)
  datos_train<- datos[1:n_train,]
  datos_test<- datos[(n_train+1):n,]
  ########################################
  #2
  Sigma_test <- cov(datos_test)
  
  #3
  Sigma_train <- cov(datos_train)
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
  
  if(mm=="promedio")
  {
    mm_hat <- lambda_est/sum(lambda_est)
  }
  
  
  
  
  #4
  grilla_gamma <- seq(0,gamma_max, by=paso_grilla)
  dd_gamma <- rep(0, length(grilla_gamma)) #aca vamos a guardar los d(gamma)
  
  for(i in 1:length(grilla_gamma))
  {gamma_n <- grilla_gamma[i]
  #dividimos por el mas grande lambda_est
  pen <- mm_hat+gamma_n*seq(1:p)^potencia_dim #funcion a minimizar divido por l
  dd_gamma[i] <-  which.min(pen)-1
  
  }
  
  
  perdida <- rep(0, length(grilla_gamma)) # loss 
  for(i in 1:length(grilla_gamma))
  {
    perdida[i] <- sum((  Sigma_test-Sigma_est_modelo(datos_train,dd_gamma[i]))^2) 
    
    
    
  }
  
  gamma_final <-which.min(perdida)
  salida <-dd_gamma[gamma_final]
  salida
}




#aca minimizamos la suma de las perdidas sobre los k folds. 
dim_est_cross_grupos <- function(datos, gamma_max=1,potencia_dim=1, mm=c("cociente") ) #, cociente original promedio
{
  #spliteamos
  n <- dim(datos)[1]
  p <- dim(datos)[2]
  
  
  grilla_gamma <- seq(0,gamma_max, by=0.01)
  lgrilla<-length(grilla_gamma)
  dd_gamma <- rep(NA,lgrilla ) #aca vamos a guardar los d(gamma)
  perdida <- rep(0, lgrilla) # loss 
  
  ngrupos<-5
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
      perdida[i] <- perdida[i]+sum((Sigma_test-Sigma_est_modelo(datos_train,dd_gamma[i]))^2) 
      
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


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}



# esta calcula un d para kada fold (tomando el d del gamma que minimiza) y elige el d mas votado. 
dim_est_cross_grupos_bis <- function(datos, gamma_max=1,potencia_dim=1, mm=c("cociente") ) #, cociente original promedio
{
  #spliteamos
  n <- dim(datos)[1]
  p <- dim(datos)[2]
  
  
  grilla_gamma <- seq(0,gamma_max, by=0.01)
  lgrilla<-length(grilla_gamma)
  dd_gamma <- rep(NA,lgrilla ) #aca vamos a guardar los d(gamma)
#  perdida <- rep(0, lgrilla) # loss 
  ngrupos<-5
    dd_k_folds <- rep(NA, ngrupos)
  
  ngrupos<-5
  for(j in 1:ngrupos)
  {
    perdida <- rep(0, lgrilla) # loss 
    
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
      
      
      perdida[i] <- sum((Sigma_test-Sigma_est_modelo(datos_train,dd_gamma[i]))^2) 
      
      
    }

        gamma_final <- which.min(perdida)
    dd_k_folds[j] <- dd_gamma[gamma_final]
    
  }
  
  ######################################
  
  
  salida<- getmode(dd_k_folds)
  salida
  
}




##############################################
#Agregamos la parte robusta
library("RobStatTM")

cov_rob_todas <- function(datos, tipo="covRob")#covRobRocke, fastmve, initPP, covRob
 # covRobMM se recomienda para dimension <10, y covRobRocke para >=10. covRob contiene a las dos y llama a la que corresponda según la dimensiòn.
{
  # if(tipo=="covRobMM")
  # {
  #   salida <- covRobMM(datos)$cov
  # } 
  # if(tipo=="covRobRocke")
  # {
  #   salida <- covRobRocke(datos)$cov
  # }
  if(tipo=="fastmve")
  {
    salida <- fastmve(datos)$cov
  }
  if(tipo=="initPP")
  {
    salida <- initPP(datos)$cova
  }
  if(tipo=="covRob")
  {
    salida <- covRob(datos)$cov
  }
  salida
}




Sigma_est_modelo_rob<- function(datos,dd, tipo="covRob")# dim<=p-1; dim=p-1 es como no tener modelo. 
{
  p <- dim(datos)[2]
  Sigma_est_rob <- cov_rob_todas(datos, tipo=tipo) #aca ponemos metodo robusto
  
  
  total.svd <- svd(Sigma_est_rob)
  autovalores <- total.svd$d#autovalores.   eigen(Sigma_est)$values 
  
  autovectores <-   total.svd$v #autovectores.
  
  sigma_cuad_est <- mean(autovalores[(dd+1):p])
  
  if(dd>0)
  {
    diagonal_est <- diag(c(autovalores[1:dd], rep(sigma_cuad_est, (p-dd))  ))
  }
  
  if(dd==0)
  {
    diagonal_est <- diag( rep(sigma_cuad_est, p) )
  }
  
  Sigma_est_modelo <- autovectores%*%diagonal_est%*%t(autovectores)
  Sigma_est_modelo 
  
}



#aca minimizamos la suma de las perdidas sobre los k folds. 
dim_est_cross_grupos_rob <- function(datos, gamma_max=1,potencia_dim=1, mm=c("cociente") , tipo="covRob") #, cociente original promedio
{
  #spliteamos
  n <- dim(datos)[1]
  p <- dim(datos)[2]
  
  
  grilla_gamma <- seq(0,gamma_max, by=0.01)
  lgrilla<-length(grilla_gamma)
  dd_gamma <- rep(NA,lgrilla ) #aca vamos a guardar los d(gamma)
  perdida <- rep(0, lgrilla) # loss 
  
  ngrupos<-5
  for(j in 1:ngrupos)
  {
    #   print(j)
    nobs<-n/ngrupos
    indices<-(j-1)*nobs+seq(1:nobs)
    datos_train<- datos[-indices,]
    datos_test<- datos[indices,]
    
    
    ########################################
    #2 
    Sigma_test_rob <- cov_rob_todas(datos_test, tipo=tipo)
    #3
    Sigma_train_rob<- cov_rob_todas(datos_train, tipo=tipo)
    total.svd <- svd(Sigma_train_rob)
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
 
      #print(i)
      perdida[i] <- perdida[i]+sum((Sigma_test_rob-Sigma_est_modelo_rob(datos_train,dd_gamma[i]))^2) 
      
    }
    
  }
  
  ######################################
  
  gamma_final <- grilla_gamma[which.min(perdida)]
  C_est <- cov_rob_todas(datos, tipo=tipo)
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
