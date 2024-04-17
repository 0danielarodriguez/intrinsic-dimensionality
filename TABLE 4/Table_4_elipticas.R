rm(list = ls())
#poner directorio de trabajo
#setwd("~/Dropbox/mariela-dani/otra-vez-autovalores/simus_articulo/p_10__normal_nrep_1000")
setwd("~/Dropbox/mariela-dani/otra-vez-autovalores/statistical papers/Revision/simus_articulo_lili/TABLE 4")
getwd()
source("funciones_cross_validation_revisadas_dani-2.R")


library("varclust")
library("SPAC2")
library("pesel")
#library("Matrix")


corrida_normalmix <- function(Nrep=10,nsamp=100,p=10,sigma=sqrt(1), signal=c(10,8,6,4,2),df=3)

  
{    
    d <- length(signal)
    canonicos <- vector("list", p)
    for(j in 1:p)
    {
      ee <- rep(0,p)
      ee[j] <- 1
      canonicos[[j]] <- ee
    }
    
    mu <- rep(0, p)
    lambda_generacion <-signal-sigma^2 # descuento el ruido 
    
    
    W <- matrix(NA, p,d) #matriz para mandar a R^p
    for(i in 1:d)
    {
      W[,i] <- sqrt(lambda_generacion[i])*canonicos[[i]]
    }
    

##################
    
ARCHIVORESULTADOS<-paste("Simu_n_",nsamp, "_Nrep_", Nrep,"signal", signal[1],".txt",sep="")  
set.seed(1200)
  for(a in 1:Nrep)
  {
    #genero datos
    #  print(a)
    
    vv <-matrix(rnorm(d*nsamp,0,1), nsamp,d)
    ee_ppca <- matrix(rnorm(p*nsamp,0,sigma), nsamp,p)
    ss <- rt(nsamp,df)^2 #genero student para las elipticas
    xx <- ss*(t(W%*%t(vv))+ee_ppca)+mu


    #Calculo estimadores 
    
    lambda_est <- svd(cov(xx))$d
    sin_penalizar<- which.min(lambda_est)-1
    
   d_ef_nuestra <- dim_est_cross_grupos(xx)
    
    
    d_ef_nuestra_covRob  <- dim_est_cross_grupos_rob(xx, tipo="covRob")
    lambda_Deng <-lambda <- eigen(as.matrix(Matrix::nearPD(stats::cov(t(scale(t(xx)))))$mat))$val

   # pPPCA(lambda = lambda_est)
    d_ef_pPCA_pen_1 <- pPPCA(lambda=lambda_Deng,Tvotes = 1000,verbose = FALSE,penalty = 1,tau = 0.001,beta = NULL )
    d_ef_pPCA_pen_2<- pPPCA(lambda=lambda_Deng, Tvotes = 1000,verbose = FALSE,penalty = 2,tau = 0.001,beta = NULL  )
    d_ef_pPCA_pen_3<- pPPCA(lambda=lambda_Deng,Tvotes = 1000,verbose = FALSE,penalty = 3,tau = 0.001,beta = NULL)
    
    # minka
    d_ef_minka_1 <- minka2001(lambda=lambda_est, M=nsamp, verbose=FALSE, tau = 0.001, BIC = FALSE)
    d_ef_minka_2 <- minka2001(lambda=lambda_est, M=nsamp, verbose=FALSE, tau = 0.001, BIC = TRUE)
    
    
    # Bayesiano 
        d_ef_bayesiano <- mlcc.bic(xx, numb.clusters = 1, numb.runs = 30, stop.criterion = 1,
                               max.iter = 30, max.dim = p, scale = FALSE, numb.cores = NULL,
                               greedy = TRUE, estimate.dimensions = TRUE, verbose = FALSE,
                               flat.prior = FALSE, show.warnings = FALSE)$subspacesDimensions[[1]]
    
    
    resultados <-     c( sin_penalizar,d_ef_nuestra, d_ef_nuestra_covRob,
                         d_ef_pPCA_pen_1,d_ef_pPCA_pen_2,d_ef_pPCA_pen_3,
                         d_ef_minka_1, d_ef_minka_2,d_ef_bayesiano)
                         
    salida<-c(a,resultados)
    n.columnas<-length(salida)  
    write(t(salida),file=ARCHIVORESULTADOS,ncolumns=n.columnas,append=T)    
  }
  
}    


####################
#simulacion: 
  enes <-c(100,300, 500,1000)
  Nrep <- 1000
  signal_1 <- rep(1,5)
  signal_2 <- c(10,8,6,4,2)
  signal_3 <- c(7,6,5)
  ss_tres_escenarios <- list(signal_1,signal_2, signal_3)
  ss_sigmas <- c(0.5,1,1)
  
  for(s in 3:3)  
  {
    signal <- ss_tres_escenarios[[s]] #aca estan los lambdas
    sigma <- ss_sigmas[s] #aca esta el ruido

      for(i in 1:length(enes))
      {
        
        nsamp <- enes[i]
        {
          print(c(s,i))
          
          corrida_normalmix(Nrep,nsamp=nsamp,p=10,sigma=sigma, signal=signal)
            
          
        }
      }
      }
    
###############################
#proceso salidas todo 

  
  enes <- c(100,300, 500,1000)
  Nrep <- 1000
  #errores <- c("normal", "student")
  signal_1 <- rep(1,5)
  signal_2 <- c(10,8,6,4,2)
  signal_3 <- c(7,6,5)
  ss_tres_escenarios <- list(signal_1,signal_2, signal_3)
  ss_sigmas <- c(0.5,1,1)
  d_true_varios <- c(5,5,3)
  salida_grilla <- c()
  
    for(s in 1:3)  
  {
    signal <- ss_tres_escenarios[[s]] #aca estan los lambdas
    sigma <- ss_sigmas[s] #aca esta el ruido
    d_true <- d_true_varios[s]
    for(i in 1:length(enes))
    {
      nsamp <- enes[i]
      {
        ARCHIVORESULTADOS<-paste("Simu_n_",nsamp, "_Nrep_", Nrep,"signal", signal[1],".txt",sep="")  
        en_grilla <- read.table(ARCHIVORESULTADOS)
        salida_grilla <- rbind(salida_grilla, c(s,nsamp,apply(en_grilla[,-1]==d_true, 2,mean)))
      }
    }
    }
   
    
  

#Con esto armamos tabla con frecuencia relativa de acierto para cada estimador. 
pirulo_bishop<-round(salida_grilla,2) # saco columna que indexa la replicacion
pirulo_bishop<-pirulo_bishop[, c(1,2,4, 5,6,10,11)]# saco columna que indexa la replicacion
colnames(pirulo_bishop)<-c( "signal", "nsamp",  "npPCA", "robPCA", "pPCA", "L","Bayes")


library(xtable)
xtable(pirulo_bishop, digits = c(c(0,0,0,0),rep(2,dim(pirulo_bishop)[2]-3)))