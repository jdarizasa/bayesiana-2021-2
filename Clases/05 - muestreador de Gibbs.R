############################
### MUESTREADOR DE GIBBS ###
############################

###############################################
# EJEMPLO CON EL MODELO NORMAL SEMI-CONJUGADO #
###############################################

#-------------------------------------------------------------------------------
# Muestreador de Gibbs para el modelo normal

# datos 
y <- c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.90,2.08)

# tama?o de la muestra
n <- length(y)

# estadisticos
mean.y <- mean(y)
var.y  <- var(y)

# hiperparametros
mu0 <- 1.9 
t20 <- 0.5^2 ##el 99% de la masa me queda en valores positivos pero es bien aplanada
s20 <- 0.01 
nu0 <- 1      ##mas o menos peso segun lo informativa

S   <- 100000                           # numero de muestras 
PHI <- matrix(data=NA, nrow=S, ncol=2)  # matriz para almacenar las muestras

ncat <- floor(S/10)  # mostrar anuncios cada 10% de las iteraciones

# ALGORITMO (muestreador de Gibbs)
# 1. inicializar la cadena
#    valor inicial: simular de la previa
#    solo es necesario alguno de los valores
set.seed(1)
phi <- c( rnorm(1, mu0, sqrt(t20)), rgamma(1, nu0/2, nu0*s20/2) ) #rgamma porque es precisión y no varianza
PHI[1,] <- phi
# 2. simular iterativamente de las distribuciones condicionales completas
set.seed(2)
for(s in 2:S) {
        # 2.1 actualizar el valor de theta
        mun <- ( mu0/t20 + n*mean.y*phi[2] ) / ( 1/t20 + n*phi[2] )
        t2n <- 1/( 1/t20 + n*phi[2] ) #estamos haciendo el mismo calculo dos veces y se puede recortar eso
        #t2n <- 1/( 1/t20 + n*phi[2] )
        #mun <- ( mu0/t20 + n*mean.y*phi[2] )*t2n
        phi[1] <- rnorm(1, mun, sqrt(t2n) )  #no olvide sacar raiz 
        # 2.2 actualizar el valor de sigma^2
        nun <- nu0+n  #esto es constante y se podria sacar del for
        s2n <- (nu0*s20 + (n-1)*var.y + n*(mean.y-phi[1])^2 ) /nun  #usamos el theta mas reciente
        phi[2] <- rgamma(1, nun/2, nun*s2n/2)
        # 2.3 almacenar
        PHI[s,] <- phi
        # 2.4 progreso
        if (s%%ncat == 0) cat(100*round(s/S, 1), "% completado ... \n", sep = "")
}
#siempre es lo mismo, actualizar-almacenar-progreso

#-------------------------------------------------------------------------------
# grafico del algoritmo
# este grafico no se acostumbra a hacer en la practica

windows(height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
m1<-5
plot(PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
      lty=1,col="gray",xlab=expression(theta),ylab=expression(tilde(sigma)^2))
text(PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )

m1<-15
plot(PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
      lty=1,col="gray",xlab=expression(theta),ylab=expression(tilde(sigma)^2))
text(PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )

m1<-100
plot(PHI[1:m1,],type="l",xlim=range(PHI[1:100,1]), ylim=range(PHI[1:100,2]),
      lty=1,col="gray",xlab=expression(theta),ylab=expression(tilde(sigma)^2))
text(PHI[1:m1,1], PHI[1:m1,2], c(1:m1) )

#-------------------------------------------------------------------------------
# grafico de la distribucion posterior

windows(height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.70,.70,0))
sseq<-1:10000
# distribucion conjunta
plot(PHI[sseq,1], PHI[sseq,2], pch=".", xlim=range(PHI[,1]),ylim=range(PHI[,2]),
     xlab=expression(theta), ylab=expression(tilde(sigma)^2))
# theta
plot(density(PHI[,1], adj=2), xlab=expression(theta), main="", xlim=c(1.55,2.05),
     ylab=expression( paste(italic("p("), theta,"|",italic(y[1]),"...",italic(y[n]),")",sep="")))
abline(v=quantile(PHI[,1],prob=c(.025,.975)),lwd=2,col="gray")
# precision
plot(density(PHI[,2],adj=2), xlab=expression(tilde(sigma)^2),main="",
     ylab=expression( paste(italic("p("),tilde(sigma)^2,"|",italic(y[1]),"...",italic(y[n]),")",sep=""))) 

#-------------------------------------------------------------------------------
# inferencia

# intervalos de credibilidad
quantile(PHI[,1], c(.025,.5,.975))           # media
quantile(PHI[,2], c(.025,.5, .975))          # precision 
quantile(1/sqrt(PHI[,2]), c(.025,.5, .975))  # desviancion estandar
quantile((1/sqrt(PHI[,2]))/PHI[,1], probs = c(0.025, 0.975))  # coeficiente de variacion

# probabilidad posterior
mean( PHI[,1] > 1.8 )

#-------------------------------------------------------------------------------
# diagnosticos de convergencia

stationarity.plot<-function(x,...) 
{
        # grafica una serie de boxplots a lo largo de la cadena
        # x : valores de la cadena
        S <- length(x)
        scan <- 1:S
        ng <- min(round(S/100), 10)
        group<-S*ceiling(ng*scan/S)/ng
        boxplot(x~group,...)              
}

# graficos
windows(height=7,width=7)
par(mfrow=c(2,2))
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
# traceplots
plot(PHI[,1], type = "l", xlab="iteration", ylab=expression(theta))  #se ve muy estacionaria desde el primer momento
plot(1/PHI[,2], type = "l", xlab="iteration", ylab=expression(sigma^2)) 
# boxplots
stationarity.plot(PHI[,1], xlab="iteration", ylab=expression(theta))  
stationarity.plot(1/PHI[,2], xlab="iteration", ylab=expression(sigma^2))  

# autocorrelacion
windows(width = 10, height = 5)
par(mfrow=c(1,2))
acf(PHI[,1])
acf(1/PHI[,2])

# tama?o efectivo de la muestra
#install.packages("coda")
library(coda)

effectiveSize( PHI )
#los resultados están muy bien, peligroso cuando el size sea realmente bajo (10-300)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------