#rm(list = ls())
setwd("~/Dropbox/mariela-dani/otra-vez-autovalores/statistical papers/Revision/simus_articulo_lili/Deng_error_norm_student")setwd("~/Dropbox/mariela-dani/otra-vez-autovalores/simus_articulo/Deng_error_norm_student")

setwd("C:/Users/0dani/Dropbox/daniela/trabajos/mariela-dani/otra-vez-autovalores/statistical papers/Revision/simus_articulo_lili/subir")

##needed packages
library(latex2exp)
library("extrafont")
font_import() 
loadfonts(device = "pdf")



nsamp <- 5000
Nrep <- 100
errores <- c("normal","student")
xi_cuad <-seq(0.8,0.99,by=0.01)
n_xi <- length(xi_cuad)
kes <- c(5,10,20)
salida_grillatodo <- c()

# Primero normales
j<-1
  error_simu <- errores[j]
  
  for(l in 1:3)
  {
    d_true <- kes[l]
    for(i in 1:n_xi)
    {
      xi2 <-xi_cuad[i]
      
      ARCHIVORESULTADOS<-paste("FIGURE 1/Simu_errores",error_simu, "_d_", d_true,"_xi_cua",xi2,".txt",sep="")  
      ARCHIVORESULTADOSCROSS<-paste("Cross/CrossSimu_errores",error_simu, "_d_", d_true,"_xi_cua",xi2,".txt",sep="")  
      ARCHIVORESULTADOSGD<-paste("Donoho/CrossSimu_errores",error_simu, "_d_", d_true,"_xi_cua",xi2,".txt",sep="")  
      en_grilla <- read.table(ARCHIVORESULTADOS)
      en_grillaCROSS <- read.table(ARCHIVORESULTADOSCROSS)
      en_grillaGD <- read.table(ARCHIVORESULTADOSGD)
      en_grillatodo<-cbind(en_grilla,en_grillaCROSS[,2],en_grillaGD[,2])
      salida_grillatodo <- rbind(salida_grillatodo, c(round(xi2,2),j,d_true,apply(en_grillatodo[,-1]==d_true, 2,mean)))
      
          }
  }



#Relative frequence table of how many times the destimated is equal to d

pirulo_Deng<-round(salida_grillatodo,2)  

pirulo_Deng <- pirulo_Deng[, c(1,2,3,5,6,10,19,20)]
colnames(pirulo_Deng)<-c("xi_cuad","error","d_tru","nuestro","pPCA1", "minka2","cross","GD")



nombres_estimadores <- c( "npPCA", "pPCA","laplace","cross","GD")


res<- as.vector(pirulo_Deng[,c(-1,-2,-3)])
M <- pirulo_Deng[,c(1,2,3)]

primeras_col_plot <- matrix(rep(t(M),length(nombres_estimadores)),ncol=ncol(M),byrow=TRUE)

estimadores <- rep(nombres_estimadores, each=60)

salplot <- data.frame(primeras_col_plot,estimadores,res)

salplot<-data.frame(salplot)
names(salplot)


tam=10

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
salida_grillatodo <- c()

j=2
  error_simu <- errores[j]
  
  for(l in 1:3)
  {
    d_true <- kes[l]
    for(i in 1:n_xi)
    {
      xi2 <-xi_cuad[i]
      
      ARCHIVORESULTADOS<-paste("FIGURE 1/Simu_errores",error_simu, "_d_", d_true,"_xi_cua",xi2,".txt",sep="")  
      ARCHIVORESULTADOSCROSS<-paste("Cross/CrossSimu_errores",error_simu, "_d_", d_true,"_xi_cua",xi2,".txt",sep="")  
      ARCHIVORESULTADOSGD<-paste("Donoho/CrossSimu_errores",error_simu, "_d_", d_true,"_xi_cua",xi2,".txt",sep="")  
      en_grilla <- read.table(ARCHIVORESULTADOS)
      en_grillaCROSS <- read.table(ARCHIVORESULTADOSCROSS)
      en_grillaGD <- read.table(ARCHIVORESULTADOSGD)
      en_grillatodo<-cbind(en_grilla,en_grillaCROSS[,2],en_grillaGD[,2])
      salida_grillatodo <- rbind(salida_grillatodo, c(round(xi2,2),j,d_true,apply(en_grillatodo[,-1]==d_true, 2,mean)))
      
      #ARCHIVORESULTADOS<-paste("Simu_errores",error_simu, "_d_", d_true,"_xi_cua",xi2,".txt",sep="")  
      # en_grilla <- read.table(ARCHIVORESULTADOS)
      
    }
  }

  
  #Relative frequence table of how many times the destimated is equal to d
  
  pirulo_Deng<-round(salida_grillatodo,2)  
  
  pirulo_Deng <- pirulo_Deng[, c(1,2,3,5,6,10,19,20)]
  colnames(pirulo_Deng)<-c("xi_cuad","error","d_tru","nuestro","pPCA1", "minka2","cross","GD")
  
  
  
  nombres_estimadores <- c( "npPCA", "pPCA","laplace","cross","GD")
  
  
  res<- as.vector(pirulo_Deng[,c(-1,-2,-3)])
  M <- pirulo_Deng[,c(1,2,3)]
  
  primeras_col_plot <- matrix(rep(t(M),length(nombres_estimadores)),ncol=ncol(M),byrow=TRUE)
  
  estimadores <- rep(nombres_estimadores, each=60)
  
  salplot <- data.frame(primeras_col_plot,estimadores,res)
  
  salplot<-data.frame(salplot)
  names(salplot)
  
  
  tam=10
  
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
  
  