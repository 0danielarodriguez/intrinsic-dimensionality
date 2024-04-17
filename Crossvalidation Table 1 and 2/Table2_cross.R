rm(list = ls())
#poner directorio de trabajo
#setwd("~/Dropbox/mariela-dani/otra-vez-autovalores/simus_articulo/p_10__normal_nrep_1000")
getwd()
library("mvtnorm")
source("npPCAcode.R")


corrida_student <- function(Nrep=10,nsamp=100,p=10,sigma=sqrt(1), 
                            signal=c(10,8,6,4,2),  nu=3)
  # error can be=c( normal,student)  nu degree of freedo of student
{ 
  
  d <- length(signal)
  
  #########
  #to generated multivariate student
  Sigma_student <- diag(c(signal,rep(sigma^2,(p-d))))*(nu-2)/nu #corregido para tener la matriz de covarianza dada por la diagonal
  
  
  
  
  ARCHIVORESULTADOS<-paste("tSimu_n_p_10",nsamp,"Student_multi", "_Nrep_", Nrep,"signal",
                           signal[1],".txt",sep="")  
  set.seed(1200)
  for(a in 1:Nrep)
  {
    #generate data
    
    xx <- rmvt(n=nsamp,Sigma_student, df=nu)
    
    
        
    #Calculo estimadores 
    resultados <-    npPCAcross(xx, ngrupos=5, tipo="cl")$dimention
                         
    salida<-c(a,resultados)
    n.columnas<-length(salida)  
    write(t(salida),file=ARCHIVORESULTADOS,ncolumns=n.columnas,append=T)    
  }
  
}    

####################
#simulation: 
enes <-c(100,300, 500,1000)
Nrep <- 1000
signal_1 <- rep(1,5)
signal_2 <- c(10,8,6,4,2)
signal_3 <- c(7,6,5)
ss_tres_escenarios <- list(signal_1,signal_2, signal_3)
ss_sigmas <- c(0.5,1,1)

for(s in 1:3)  
{
  signal <- ss_tres_escenarios[[s]] #  lambdas
  sigma <- ss_sigmas[s] #noise
  for(i in 1:length(enes))
  {
    print(c(s,i))
    nsamp <- enes[i]
    corrida_student(Nrep=Nrep,nsamp=nsamp,p=10,sigma=sigma, 
                    signal=signal,  nu=3)
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
        ARCHIVORESULTADOS<-paste("tSimu_n_p_10",nsamp,"Student_multi", "_Nrep_", Nrep,"signal",
                                 signal[1],".txt",sep="")  
  #      ARCHIVORESULTADOS<-paste("tSimu_n_",nsamp, "_Nrep_", Nrep,"signal", signal[1],".txt",sep="")  
        en_grilla <- read.table(ARCHIVORESULTADOS)
        salida_grilla <- rbind(salida_grilla, c(s,nsamp,mean(en_grilla[,-1]==d_true)))#apply(en_grilla[,-1]==d_true, 2,mean)))
        }
    }
    }

#Con esto armamos tabla con frecuencia relativa de acierto para cada estimador. 
pirulo_bishop<-round(salida_grilla,2) # saco columna que indexa la replicacion
library(xtable)
xtable(pirulo_bishop)
