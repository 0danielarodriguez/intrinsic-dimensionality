rm(list = ls())
setwd("~/Dropbox/mariela-dani/otra-vez-autovalores/statistical papers/Revision/simus_articulo_lili/TABLE 2")

##needed packages
library("mvtnorm")
library("varclust")
install.packages("devtools") # if you have not installed already

devtools::install_github("WeiAkaneDeng/SPAC2")


  library("SPAC2")
library("pesel")
library("RobStatTM")
library(latex2exp)
library("extrafont")
font_import() 
loadfonts(device = "pdf")
##needed functions
source("functions.R")

 


corrida_student <- function(Nrep=10,nsamp=100,p=10,sigma=sqrt(1), 
                    signal=c(10,8,6,4,2),  nu=3)
  # error can be=c( normal,student)  nu degree of freedo of student
{ 
  
    d <- length(signal)
    
    #########
    #to generated multivariate student
    Sigma_student <- diag(c(signal,rep(sigma^2,(p-d))))*(nu-2)/nu #corregido para tener la matriz de covarianza dada por la diagonal
    
    
    
    
    ARCHIVORESULTADOS<-paste("Simu_n_p_10",nsamp,"Student_multi", "_Nrep_", Nrep,"signal",
                             signal[1],".txt",sep="")  
    set.seed(1200)
  for(a in 1:Nrep)
  {
    #generate data
     
        xx <- rmvt(n=nsamp,Sigma_student, df=nu)
    
        
    #Compute estimators 
    
    lambda_est <- svd(cov(xx))$d
    
    sin_penalizar<- which.min(lambda_est)-1
    d_ef_nuestra <- dim_est_cross_grupos(xx)
    
    d_ef_nuestra_covRob  <- dim_est_cross_grupos_rob(xx, tipo="covRob")
    

    # Deng estimators
    lambda_Deng <-lambda <- eigen(as.matrix(Matrix::nearPD(stats::cov(t(scale(t(xx)))))$mat))$val

    # pPPCA(lambda = lambda_est)
    d_ef_pPCA_pen_1 <- pPPCA(
      lambda=lambda_Deng,
      Tvotes = 1000,
      verbose = FALSE,
      penalty = 1,
      tau = 0.001,
      beta = NULL
    )
    
    
    
    d_ef_pPCA_pen_2<- pPPCA(
      lambda=lambda_Deng,
      Tvotes = 1000,
      verbose = FALSE,
      penalty = 2,
      tau = 0.001,
      beta = NULL
    )
    
    
    d_ef_pPCA_pen_3<- pPPCA(
      lambda=lambda_Deng,
      Tvotes = 1000,
      verbose = FALSE,
      penalty = 3,
      tau = 0.001,
      beta = NULL
    )
    
    # minka
    d_ef_minka_1 <- minka2001(lambda=lambda_est, M=nsamp, verbose=FALSE, tau = 0.001, BIC = FALSE)
    
    d_ef_minka_2 <- minka2001(lambda=lambda_est, M=nsamp, verbose=FALSE, tau = 0.001, BIC = TRUE)
    
    
    # Bayes 
      d_ef_bayesiano <- mlcc.bic(xx, numb.clusters = 1, numb.runs = 30, stop.criterion = 1,
                               max.iter = 30, max.dim = p, scale = FALSE, numb.cores = NULL,
                               greedy = TRUE, estimate.dimensions = TRUE, verbose = FALSE,
                               flat.prior = FALSE, show.warnings = FALSE)$subspacesDimensions[[1]]
    
    
    
    
        
    resultados <-     c( sin_penalizar,d_ef_nuestra,
                         d_ef_nuestra_covRob,
                         d_ef_pPCA_pen_1,d_ef_pPCA_pen_2,d_ef_pPCA_pen_3,
                         d_ef_minka_1,
                         d_ef_minka_2, 
                          d_ef_bayesiano)
                         
    
    
    
    salida<-c(a,resultados)
    n.columnas<-length(salida)  
    write(t(salida),file=ARCHIVORESULTADOS,ncolumns=n.columnas,append=T)    
  }
  
}    


#
#corrida_student(Nrep=5,nsamp=100,p=10,sigma=sqrt(1), 
                           # signal=c(10,8,6,4,2),  nu=3)
  


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
#Proces the runs 

  
  enes <-c(100,300, 500,1000)
  Nrep <- 1000
  signal_1 <- rep(1,5)
  signal_2 <- c(10,8,6,4,2)
  signal_3 <- c(7,6,5)
  ss_tres_escenarios <- list(signal_1,signal_2, signal_3)
  ss_sigmas <- c(0.5,1,1)
  salida_grilla <- c()
    for(s in 1:3)  
  {
    signal <- ss_tres_escenarios[[s]] # lambdas
    sigma <- ss_sigmas[s] #noise
    d_true <- length(signal)
    for(i in 1:length(enes))
    {
      nsamp <- enes[i]
      ARCHIVORESULTADOS<-paste("Simu_n_p_10",nsamp,"Student_multi", "_Nrep_", Nrep,"signal",
                                 signal[1],".txt",sep="")  
        en_grilla <- read.table(ARCHIVORESULTADOS)
        salida_grilla <- rbind(salida_grilla, c(s,nsamp,apply(en_grilla[,-1]==d_true, 2,mean)))
        
      }
    }
    

  #Relative frequence table of how many times the destimated is equal to d
  
pirulo_stud<-round(salida_grilla,2) # without the index column
colnames(pirulo_stud)<-c("signal", "nsamp", "sin_pen", "nuestro",
                           "robusto", 
                         "pPCA1","pPCA2","pPCA3", "minka1",
                         "minka2", 
                         "bayesiano")




pirulo_stud<-pirulo_stud[, c(1,2,4, 5,6,10,11)]# saco columna que indexa la replicacion
colnames(pirulo_stud)<-c( "signal", "nsamp",  "npPCA",
                            "robPCA", 
                            "pPCA",
                            "L","Bayes")

 



library(xtable)
xtable(pirulo_stud)
