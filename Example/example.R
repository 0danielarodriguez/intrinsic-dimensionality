setwd("C:/Users/0dani/Dropbox/daniela/trabajos/mariela-dani/otra-vez-autovalores/pruebas_R/datos")
bloodpressure <- read.delim("C:/Users/0dani/Dropbox/daniela/trabajos/mariela-dani/otra-vez-autovalores/pruebas_R/datos/presion.txt")

#primero correr lo que esta al final
# que son las funciones

#EJEMPLO
bloodpressure <- read.delim("presion.txt")

# BP respuesta
# 6 predictoras que las guardo en presion
presion<-bloodpressure[,-c(1,2)]
#esta funcion esta en el archivo adjunto.
dim_est_cross_grupos(presion,ngrupos=10, gamma_max=0.5,potencia_dim=1, mm=c("cociente") ) #, cociente original promedio

#lambdan= es una constante (C) por una potencia de n (potencia_n)
#(potencia_n) que se puede elegir por defaul es 0.1
# C tambien se puede elegir por defaul es 1 (C=0 corresponde al estimador sin penalizar)
# la funcion que penaliza la dimension es una potencia de la dimension (potencia_dim) por default es 2
#mm puede ser #, cociente original promedio, segun como se 
#quiera la funcion objetivo cociente del lambda o solo un lambda o  promedios de lambdas 


#en el paper hacemos 100 repeticiones de  cv en grupos 
#para  ver en que medida depende de como reparte los datos
Nrep <- 100
dim_colin <- rep(NA, Nrep)
set.seed(999)
for (i in 1: Nrep)
{
  npre<-dim(presion)[1]
  indran<-sample(1:npre,npre,replace=FALSE)
  presionrandom<-presion[indran,]
  presionrandom1<-presion1[indran,]
  dim_colin[i] <- dim_est_cross_grupos(datos=presionrandom, gamma_max=0.5,potencia_dim=1, mm=c("cociente") ) #, cociente original promedio
  
}

d_colin <- table(dim_colin)
d_colin




