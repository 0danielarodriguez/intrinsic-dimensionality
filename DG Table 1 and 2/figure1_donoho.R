rm(list = ls())
source("~/Dropbox/mariela-dani/otra-vez-autovalores/statistical papers/Revision/simus_articulo_lili/donoho/funciones.R")
library("RobStatTM")
library(latex2exp)
library("extrafont")
font_import() 
loadfonts(device = "pdf")
##needed functions
 
source("C:/Users/Usuario/Dropbox/mariela-dani/otra-vez-autovalores/statistical papers/Revision/simus_articulo_lili/donoho/funciones.R")

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
  
  
  
  ARCHIVORESULTADOS<-paste("CrossSimu_errores",error, "_d_", d,"_xi_cua",round(sigma_cuad,2),".txt",sep="")  
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
    
    
    resultados <-   elije_donoho_el_d(xx) 
    
    
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

for(j in 1:2)
{ 
  error_simu <- errores[j]
  
  for(l in 1:3)
  {
    d_true <- kes[l]
    for(i in 1:n_xi)
    {
      xi2 <-xi_cuad[i]
      
      ARCHIVORESULTADOS<-paste("crossSimu_errores",error_simu, "_d_", d_true,"_xi_cua",xi2,".txt",sep="")  
      
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

