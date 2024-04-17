rm(list = ls())
setwd("~/Dropbox/mariela-dani/otra-vez-autovalores/statistical papers/Revision/simus_articulo_lili/Deng_error_norm_student")setwd("~/Dropbox/mariela-dani/otra-vez-autovalores/simus_articulo/Deng_error_norm_student")

##needed packages
library("varclust")
library("SPAC2")
library("pesel")
library("RobStatTM")
library(latex2exp)
library("extrafont")
font_import() 
loadfonts(device = "pdf")
##needed functions
source("functions.R")


## Running

corrida <- function(Nrep=5,nsamp=1000,p=100,sigma_cuad=0.8, 
                    d=5,  error="normal", nu=3)
                    
 
# error can be error=c( normal,student)  nu degree of freedom of the  student
{ 

  lambda_1 <- (p-(p-d)*sigma_cuad)/d
  signal <- rep(lambda_1,d)
  canonicos <- vector("list", p)
  for(j in 1:p)
  {
    ee <- rep(0,p)
    ee[j] <- 1
    canonicos[[j]] <- ee
  }
  
  mu <- rep(0, p)
  lambda_generacion <-signal-sigma_cuad # minus the noise 
  
  
  W <- matrix(NA, p,d) #matrix to go to R^p
  for(i in 1:d)
  {
    W[,i] <- sqrt(lambda_generacion[i])*canonicos[[i]]
  }
   
  
  
  ARCHIVORESULTADOS<-paste("Simu_errores",error, "_d_", d,"_xi_cua",round(sigma_cuad,2),".txt",sep="")  
  set.seed(123)
  for(a in 1:Nrep)
  {
    #generating the data
     
    vv <-matrix(rnorm(d*nsamp,0,1), nsamp,d)
    
    if(error=="normal")
    {
    ee_ppca <- matrix(rnorm(p*nsamp,0,sqrt(sigma_cuad)), nsamp,p)
    }
    
    if(error=="student")
    {
      ee_ppca <- sqrt(sigma_cuad)*matrix(rt(p*nsamp,nu), nsamp,p)/sqrt(nu/(nu-2)) #errores student
    }
    
    xx <- t(W%*%t(vv))+ee_ppca+mu
    
  #Compute the estimators
    lambda_est <- svd(cov(xx))$d
    
    sin_penalizar<- which.min(lambda_est)-1
    d_ef_nuestra <- dim_est_cross_grupos(xx)
    
    
    
    lambda_Deng <-lambda <- eigen(as.matrix(Matrix::nearPD(stats::cov(t(scale(t(xx)))))$mat))$val
    
    
    
    
     
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
    
    # Minka
    d_ef_minka_1 <- minka2001(lambda=lambda_est, M=nsamp, verbose=FALSE, tau = 0.001, BIC = FALSE)
    
    d_ef_minka_2 <- minka2001(lambda=lambda_est, M=nsamp, verbose=FALSE, tau = 0.001, BIC = TRUE)
    
    
    #OTHERS
    
     
    
    
    d_ef_pesel_het <- pesel(xx, npc.min = 1, npc.max = p-1, scale = TRUE, method = "heterogenous")$nPCs # results same either t() or ()
    d_ef_pesel_homo <-  pesel(xx, npc.min = 1, npc.max = p-1, scale = TRUE, method = "homogenous")$nPCs
    
    
    d_ef_Law <- LawleyTest(lambda = lambda_est) #new
    
      
    lambda <- as.numeric(lambda_est)
    N <- sum(lambda > 0, na.rm = T)
    
    adjD <- which.min(c(lambda[-1]/lambda[-N], 1))
    cumD <- which.max(cumsum(lambda)/(1:N)/((sum(lambda)-cumsum(lambda))/(N-1:N)))
    varD <- which.max(sapply(1:N, function(x) stats::var(lambda[1:x])))
    cumlog <- which.min(log(cumsum(lambda))- cumsum(log(lambda)))
    
    # Bayes 
    
    
    d_ef_bayesiano <- mlcc.bic(xx, numb.clusters = 1, numb.runs = 30, stop.criterion = 1,
                               max.iter = 30, max.dim = p, scale = FALSE, numb.cores = NULL,
                               greedy = TRUE, estimate.dimensions = TRUE, verbose = FALSE,
                               flat.prior = FALSE, show.warnings = FALSE)$subspacesDimensions[[1]]
    
    
    
    
    resultados <-     c( sin_penalizar,d_ef_nuestra,d_ef_pPCA_pen_1,d_ef_pPCA_pen_2,d_ef_pPCA_pen_3, d_ef_minka_1,
                         d_ef_minka_2, 
                         d_ef_pesel_het, d_ef_pesel_homo,
                         d_ef_Law,
                         adjD,
                         cumD,varD,cumlog,d_ef_bayesiano)
    
    
    
    
    
    salida<-c(a,resultados)
    n.columnas<-length(salida)  
    write(t(salida),file=ARCHIVORESULTADOS,ncolumns=n.columnas,append=T)    
  }
  
}   


 


####################
#SIMULATION: 
nsamp <- 5000
Nrep <- 100
errores <- c("normal","student")
xi_cuad <-seq(0.8,0.99,by=0.01)
n_xi <- length(xi_cuad)
kes <- c(5,10,20)

for(j in 1:2)
{ 
  error_simu <- errores[j]
  
  for(l in 1:3)
  {
    print(c(j,l))
    d_true <- kes[l]
    for(i in 1:n_xi)
    {
      xi2 <- xi_cuad[i]
  
      corrida(Nrep=Nrep,nsamp=nsamp,p=100,sigma_cuad =xi2 , 
              d=d_true,  error=error_simu, nu=3)
    }
  }
}


###############################
#ONCE we run we process the data 


nsamp <- 5000
Nrep <- 100
errores <- c("normal","student")
xi_cuad <-seq(0.8,0.99,by=0.01)
n_xi <- length(xi_cuad)
kes <- c(5,10,20)
salida_grilla <- c()

for(j in 1:1)
{ 
  error_simu <- errores[j]
  
  for(l in 1:3)
  {
    d_true <- kes[l]
    for(i in 1:n_xi)
    {
      xi2 <-xi_cuad[i]
      
      ARCHIVORESULTADOS<-paste("Simu_errores",error_simu, "_d_", d_true,"_xi_cua",xi2,".txt",sep="")  
      
      en_grilla <- read.table(ARCHIVORESULTADOS)
      salida_grilla <- rbind(salida_grilla, c(round(xi2,2),j,d_true,apply(en_grilla[,-1]==d_true, 2,mean)))
      
          }
  }
}





#Relative frequence table of how many times the destimated is equal to d

pirulo_Deng<-round(salida_grilla,2)  

pirulo_Deng <- pirulo_Deng[, c(1,2,3,5,6,10)]
colnames(pirulo_Deng)<-c("xi_cuad","error","d_tru",   "nuestro",
                         "pPCA1",
                         "minka2")



nombres_estimadores <- c( "npPCA",
                         "pPCA",
                         "laplace")


res<- as.vector(pirulo_Deng[,c(-1,-2,-3)])
M <- pirulo_Deng[,c(1,2,3)]

primeras_col_plot <- matrix(rep(t(M),length(nombres_estimadores)),ncol=ncol(M),byrow=TRUE)

estimadores <- rep(nombres_estimadores, each=60)

salplot <- data.frame(primeras_col_plot,estimadores,res)

salplot<-data.frame(salplot)
names(salplot)


tam=30

##FIGURE 1 

colnames(salplot)<-c("xi_cuad","error", "d","estimadores","res")
library(ggplot2)
pp1 <- ggplot(data = salplot, aes(xi_cuad, res)) + 
  geom_point(aes(colour=as.factor(estimadores)))+ 
  geom_line(aes(colour=as.factor(estimadores)), size=1.1)+  
  facet_wrap(~d)+
  labs(colour = "Estimator",y="Fraction of correct estimations")+xlab(TeX("$\\sigma^2$"))+
  theme_classic()+theme(legend.title=element_text(size=tam),
                       legend.text = element_text(size=tam, family = "Arial"))+
  
 theme(legend.text = element_text(family = "Arial"))+
  theme(axis.text=element_text(size=tam),
        axis.title=element_text(size=tam ))
pp1



 
#############################################
#Repetimos para errores student

nsamp <- 5000
Nrep <- 100
errores <- c("normal","student")
xi_cuad <-seq(0.8,0.99,by=0.01)
n_xi <- length(xi_cuad)
kes <- c(5,10,20)
salida_grilla <- c()

j=2
  error_simu <- errores[j]
  
  for(l in 1:3)
  {
    d_true <- kes[l]
    for(i in 1:n_xi)
    {
      xi2 <-xi_cuad[i]
      
      ARCHIVORESULTADOS<-paste("Simu_errores",error_simu, "_d_", d_true,"_xi_cua",xi2,".txt",sep="")  
      
      en_grilla <- read.table(ARCHIVORESULTADOS)
      salida_grilla <- rbind(salida_grilla, c(round(xi2,2),j,d_true,apply(en_grilla[,-1]==d_true, 2,mean)))
      
    }
  }

 