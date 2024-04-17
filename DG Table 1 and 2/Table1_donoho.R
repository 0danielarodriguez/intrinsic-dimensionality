rm(list = ls())
#poner directorio de trabajo
#setwd("~/Dropbox/mariela-dani/otra-vez-autovalores/simus_articulo/p_10__normal_nrep_1000")
getwd()
source("C:/Users/Usuario/Dropbox/mariela-dani/otra-vez-autovalores/statistical papers/Revision/simus_articulo_lili/donoho/funciones.R") 


corrida_normal <- function(Nrep=10,nsamp=100,p=10,sigma=sqrt(1), signal=c(10,8,6,4,2))
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
    
ARCHIVORESULTADOS<-paste("normalSimu_n_",nsamp, "_Nrep_", Nrep,"signal", signal[1],".txt",sep="")  
set.seed(1200)
  for(a in 1:Nrep)
  {
    #genero datos
    #  print(a)
    
    vv <-matrix(rnorm(d*nsamp,0,1), nsamp,d)
    ee_ppca <- matrix(rnorm(p*nsamp,0,sigma), nsamp,p)
    xx <- t(W%*%t(vv))+ee_ppca+mu

        
    #Calculo estimadores 
     
    resultados <-    elije_donoho_el_d(xx)                     
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
  
  for(s in 1:3)  
  {
    signal <- ss_tres_escenarios[[s]] #aca estan los lambdas
    sigma <- ss_sigmas[s] #aca esta el ruido

      for(i in 1:length(enes))
      {
        
        nsamp <- enes[i]
        {
          print(c(s,i))
          
          corrida_normal(Nrep,nsamp=nsamp,p=10,sigma=sigma, signal=signal)
            
          
        }
      }
      }
    
###############################
#proceso salidas todo 

  
  enes <- c(100,300, 500,1000)
  Nrep <- 1000
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
        ARCHIVORESULTADOS<-paste("normalSimu_n_",nsamp, "_Nrep_", Nrep,"signal", signal[1],".txt",sep="")  
        en_grilla <- read.table(ARCHIVORESULTADOS)
        salida_grilla <- rbind(salida_grilla, c(s,nsamp,mean(en_grilla[,-1]==d_true)))#apply(en_grilla[,-1]==d_true, 2,mean)))
      #  salida_resumen <- table(en_grilla[,2] )
      #  salida_resumen/1000
      }
    }
    }

#Con esto armamos tabla con frecuencia relativa de acierto para cada estimador. 
pirulo_donoho<-round(salida_grilla,2) # saco columna que indexa la replicacion
library(xtable)
xtable(pirulo_donoho)

