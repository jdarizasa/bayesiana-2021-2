########################
### Regresion Poison ###
########################
########################

# directorio de trabajo
setwd("C:/Users/Juan Camilo/Dropbox/UN/estadistica_bayesiana/codigo")

# cargar funciones 
source("funciones_auxiliares.R")

################################################################################
### Ejemplo: Algoritmo de Metropolis distribucion normal con varianza conocida #
################################################################################
################################################################################

# MODELO (sigma^2 conocido)
# 1.  y_i | theta ~ Normal(theta, sigma^2)
# 2.  theta ~ Normal(mu, tau^2)
# 3.  sigma^2 = 1

# simular datos
n  <- 5
s2 <- 1 
set.seed(1)
y  <- round(rnorm(n,10,1),2)

# previa
t2 <- 10 
mu <- 5

# posterior
mu.n <- ( mean(y)*n/s2 + mu/t2 )/( n/s2 + 1/t2 ) 
t2.n <- 1/( n/s2 + 1/t2 )

######## Metropolis
theta <- 0      # valor inicial
delta <- 2      # parametro de ajuste
S     <- 10000  # numero de muestras
THETA <- NULL   # almacenamiento

# cadena
set.seed(1)
for (s in 1:S) {
        # 1. propuesta
        theta.star <- rnorm(1,theta,sqrt(delta))
        # 2. tasa de aceptacion
        log.r <- ( sum(dnorm(y,theta.star,sqrt(s2),log=TRUE)) + dnorm(theta.star,mu,sqrt(t2),log=TRUE) ) - 
                 ( sum(dnorm(y,theta,     sqrt(s2),log=TRUE)) + dnorm(theta,     mu,sqrt(t2),log=TRUE) ) 
        # 3. actualizar
        if (runif(1) < exp(log.r)) { 
                theta<-theta.star 
        } 
        # 4. almacenar
        THETA <- c(THETA, theta)
}
######## fin MCMC

# grafico
windows(height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
# cadena
# adelgazamiento de la cadena (reducir autocorrelacion)
# skeep<-seq(10,S,by=10)
skeep <- 1:S
plot(skeep,THETA[skeep],type="l",xlab="Iteraci?n",ylab=expression(theta))
# histograma metropolis y posterior analitica
# se omiten las primeras 50 observaciones (periodo de calentamiento)
hist(THETA[-(1:50)],prob=TRUE,main="",xlab=expression(theta),
     ylab="Densidad",col="gray95",border="gray")
th <- seq(min(THETA),max(THETA),length=100)
lines(th, dnorm(th,mu.n,sqrt(t2.n)), col = 2, lty = 2, lwd=2)


######## Metropolis con diferentes parametros de ajuste (delta)

ACR    <- NULL  # tasa de aceptaciones
ACF    <- NULL  # autocorrelaciones
THETAA <- NULL  # muestras

for(delta2 in 2^c(-5,-1,1,5,7) ) {
        # parametros iniciales
        THETA <- NULL
        S     <- 10000
        theta <- 0
        acs   <- 0  # tasa de aceptacion
        # cadena
        set.seed(1)
        for(s in 1:S) {
                # 1. propuesta
                theta.star<-rnorm(1,theta,sqrt(delta2))
                # 2. tasa de aceptacion
                log.r <- sum( dnorm(y,theta.star,sqrt(s2),log=TRUE) - dnorm(y,theta,sqrt(s2),log=TRUE) )  +  
                         dnorm(theta.star,mu,sqrt(t2),log=TRUE) - dnorm(theta,mu,sqrt(t2),log=TRUE) 
                # 3. actualizar
                if(log(runif(1)) < log.r) { 
                        theta <- theta.star 
                        acs   <- acs + 1 
                }
                # 4. almacenar
                THETA <- c(THETA, theta) 
        }
        # fin MCMC
        # almacenar valores de todos los casos (delta2)
        ACR    <- c(ACR, acs/s) 
        ACF    <- c(ACF,acf(THETA,plot=FALSE)$acf[2])
        THETAA <- cbind(THETAA,THETA)
}
########

# graficos
windows(height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
laby <- c(expression(theta),"","","","")
for(k in c(1,3,5)) {
        plot(THETAA[1:500,k],type="l",xlab="Iteraci?n",ylab=laby[k], ylim=range(THETAA) )
        abline(h=mu.n,lty=2)
}

# tasas de aceptacion
ACR

# autocorrelaciones
ACF 

# tama?os efectivos de muestra
coda::effectiveSize(THETAA)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

########################
### Regresion Poison ###
########################
########################

#-------------------------------------------------------------------------------
# Descripcion: 
# Actividades de reproduccion de gorriones en funcion de la edad (Arcese et al, 1992).
# n = 52 gorriones hembras.
# "age"     : edad.
# "fledged" : numero de crias.
#-------------------------------------------------------------------------------

############
### data ###
############

spfage <- structure(c(3, 1, 1, 2, 0, 0, 6, 3, 4, 2, 1, 6, 2, 3, 3, 4, 7, 2, 2, 1, 
                      1, 3, 5, 5, 0, 2, 1, 2, 6, 6, 2, 2, 0, 2, 4, 1, 2, 5, 1, 2, 
                      1, 0, 0, 2, 4, 2, 2, 2, 2, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 
                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                      1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                      1, 1, 1, 1, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 
                      2, 2, 2, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 
                      4, 4, 5, 5, 5, 5, 3, 3, 3, 3, 3, 3, 3, 6, 1, 1, 9, 9, 1, 1, 
                      1, 1, 1, 1, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 25, 25, 16, 16, 
                      16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 25, 16, 16, 16, 16, 
                      25, 25, 25, 25, 9, 9, 9, 9, 9, 9, 9, 36, 1, 1), 
                    .Dim = c(52L, 4L), 
                    .Dimnames = list(NULL, c("fledged", "intercept", "age", "age2")))

spfage <- as.data.frame(spfage)

spf  <- spfage$fledged  # y  = variable dependiente (respuesta)
age  <- spfage$age      # x1 = variable independiente 1
age2 <- age^2           # x2 = variable independiente 2

############################
### analisis descriptivo ###
############################

######## diagrama de caja
windows(height=3.5, width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(spf~as.factor(age), range=0, xlab="Edad (a?os)", ylab="No. Crias", col="gray", border="lightgray")
########

# GLM (frecuentista)
summary(glm(spf~age+age2,family="poisson"))

###################
### Monte Carlo ###
###################

y <- spf                                # variable respuesta
X <- cbind(rep(1,length(y)),age,age^2)  # matriz de dise?o
n <- length(y)                          # tama?o de la muestrsa
p <- dim(X)[2]                          # numero de predictores

# previa
pmn.beta <- rep(0,  p)  # beta0 = 0
psd.beta <- rep(10, p)  # Sigma0 = 100*I

#-------------------------------------------------------------------------------
# En muchos problemas, la varianza posterior es una elecci?n eficiente para la 
# distribuci?n de propuestas. Aunque no se conoce la varianza posterior antes de 
# ejecutar el algoritmo de Metr?polis, a menudo basta con utilizar una aproximaci?n.
# Si esto da como resultado una tasa de aceptaci?n demasiado alta o demasiado baja, 
# siempre es posible ajustar la variabilidad de la propuesta en consecuencia.
#-------------------------------------------------------------------------------
# log: funcion de enlace 
# y + 1/2 evitar problemas en la frontera con 0
var.prop <- var(log(y + 1)) * solve( t(X)%*%X ) # matriz de varianza propuesta
beta <- rep(0,p) # valor inicial beta



######## algoritmo de metropolis
S    <- 250000
BETA <- matrix(NA, nrow=S, ncol=p)
ac   <- 0
ncat <- floor(S/10)

# cadena
library(mvtnorm)
set.seed(1)
for(s in 1:S) {
        # 1. propuesta
        beta.p<- t(rmvnorm(1, beta, var.prop ))
        # 2. tasa de aceptacion
        lhr <- sum(dpois(y,exp(X%*%beta.p),log=T)) -
               sum(dpois(y,exp(X%*%beta),log=T)) +
               sum(dnorm(beta.p,pmn.beta,psd.beta,log=T)) -
               sum(dnorm(beta,pmn.beta,psd.beta,log=T))
        # 3. actualizar
        if (log(runif(1)) < lhr) { 
                beta <- beta.p 
                ac   <- ac + 1 
        }
        # 4. almacenar
        BETA[s,] <- beta
        # 5. Progreso
        if (s%%ncat == 0) cat(100*round(s/S, 1), "% completed ... \n", sep = "" )
}
######### fin mcmc

# tasa de aceptacion
100*ac/S

# diagnosticos
library(coda)
apply(X = BETA, MARGIN = 2, FUN = effectiveSize)

#### grafico diagnostico
blabs <- c(expression(beta[1]),expression(beta[2]),expression(beta[3]))  # etiquetas
thin  <- c(1,(1:1000)*(S/1000))  # muestreo sistematico 

windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
j<-3
plot(thin,BETA[thin,j],type="l",xlab="Iteraci?n",ylab=blabs[j])
abline(h=mean(BETA[,j]) )

acf(BETA[,j],ci.col="gray",xlab="lag")
acf(BETA[thin,j],xlab="lag/10",ci.col="gray")
####


#### grafico posterior
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))

plot(density(BETA[,2],adj=2), type="l", main="",
     xlab=expression(beta[2]),
     ylab=expression(paste(italic("p("),beta[2],"|",italic("y)"),sep="") ) ,
     lwd=2,lty=1,col="black")

plot(density(BETA[,3],adj=2), type="l", main="",
     xlab=expression(beta[3]),
     ylab=expression(paste(italic("p("),beta[3],"|",italic("y)"),sep="") ),
     lwd=2,col="black",lty=1)

Xs       <- cbind(rep(1,6),1:6,(1:6)^2) 
eXB.post <- exp(t(Xs%*%t(BETA )) )
qE       <- apply( eXB.post,2,quantile,probs=c(.025,.5,.975))

plot(c(1,6),range(c(0,qE)),type="n",xlab="Edad (a?os)", ylab="No. Crias")
lines( qE[1,],col="black",lwd=1)
lines( qE[2,],col="black",lwd=2)
lines( qE[3,],col="black",lwd=1)
####

# inferencia posterior

quantile(BETA[,2], c(.025, .975))
quantile(BETA[,3], c(.025, .975))

mean( BETA[,2] > 0  )
mean( BETA[,3] > 0  )

################################################################################

########
# JAGS #
########

library(R2jags)

model <- function() {
    for (i in 1:n) {
        y[i] ~ dpois(theta[i])
        log(theta[i]) <- inprod(X[i,], beta)
    }
    for (j in 1:p) {
        beta[j] ~ dnorm(beta0, phi0)    
    }
}

# previa
beta0 <- 0
phi0  <- 1/100
  
# input
model_data <- list(y = y, X = X, n = n, p = p, beta0 = beta0, phi0 = phi0)

# parameters
model_parameters <- c("beta")

# initial values
initial_values <- list(list("beta" = rep(beta0, p)), 
                       list("beta" = rep(beta0, p)),
                       list("beta" = rep(beta0, p)))

# mcmc settings
niter  <- 26000
nburn  <- 1000
nthin  <- 25
nchain <- length(initial_values)

# mcmc
set.seed(123)
fit <- jags(data = model_data, inits = initial_values, 
            parameters.to.save = model_parameters, model.file = model, 
            n.chains = nchain, n.iter = niter, n.thin = nthin, n.burnin = nburn)

print(fit)

# transformar a objecto MCMC               
fit_mcmc <- coda::as.mcmc(fit)

# plots
windows()
mcmcplots::traplot(fit_mcmc)

windows()
mcmcplots::denplot(fit_mcmc)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

###################################################
#### Modelo lineal con errores correlacionados ####
###################################################
###################################################

#-------------------------------------------------------------------------------
# descripcion:
# Los an?lisis de n?cleos de hielo de la Ant?rtida han permitido a los 
# cient?ficos deducir condiciones atmosf?ricas hist?ricas de los ?ltimos cientos 
# de miles de a?os (Petit et al, 1999). 
#
# Los datos incluyen 200 valores de la temperatura medida en intervalos de tiempo 
# aproximadamente iguales; tiempo entre mediciones consecutivas de aproximadamente 
# 2,000 a?os. 
#
# La temperatura se registra en t?rminos de su diferencia de la temperatura actual
# en grados Celsius, y la concentraci?n de CO2 (dioxido de carbono) se registra 
# en partes por mill?n.
#
# Modelar la temperatura en funci?n del CO2
#
#-------------------------------------------------------------------------------

################
# cargar datos #
################

# directorio de trabajo
setwd("C:/Users/Juan Camilo/Dropbox/UN/estadistica_bayesiana/codigo")

# data
load("dct.RData")

# cargar funciones
source("funciones_auxiliares.R")

########################
# analisis descriptivo #
########################

# series de tiempo de la temperatura y la concentraci?n de di?xido de carbono 
# en una escala estandarizada (centrada y escalada) para tener una media de cero 
# y una varianza de uno).
#
# La gr?fica indica que la historia temporal de temperatura y CO2 siguen patrones 
# muy similares.
#
# La concentraci?n de CO2 parece ser un predictor de la temperatura. 
#
# Un modelo de regresi?n linal simple es razonable?
# Caracterizaci?n de los residuales:
#  * Normalidad? Ok.
#  * IID? No. ACF_1 aprox. 0.5.

####### grafico
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
layout(matrix( c(1,1,2),nrow=1,ncol=3) )

plot(dct[,1], (dct[,3]-mean(dct[,3]))/sd(dct[,3]), 
     type="l",col="black", xlab="A?o",ylab="Medici?n estandarizada",ylim=c(-2.5,3))
legend(-115000,3.2,legend=c("Temp.",expression(CO[2])),bty="n", lwd=c(2,2),col=c("black","gray"))
lines(dct[,1],  (dct[,2]-mean(dct[,2]))/sd(dct[,2]),type="l",col="gray")

plot(dct[,2], dct[,3],xlab=expression(paste(CO[2],"(ppmv)")),ylab="Diferencia de temp. (deg C)")
########


######## grafico residuales
lmfit<-lm(dct$tmp~dct$co2)

windows(height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
hist(lmfit$res,main="",xlab="Residual",ylab="Frecuencia",freq=F)
curve(expr = dnorm(x,mean(lmfit$res),sd(lmfit$res)),add=T,col=2)
acf(lmfit$res,ci.col="gray",xlab="lag")
########


########
# MCMC #
########

# data
n  <-dim(dct)[1]
y  <-dct[,3]
X  <-cbind(rep(1,n),dct[,2])
DY <-abs(outer( (1:n),(1:n) ,"-")) # para construir la matriz de correlacion

# valores iniciales
library(nlme)
lmfit   <- lm(y ~ -1+X)
fit.gls <- gls(y~X[,2], correlation=corARMA(p=1), method="ML")

beta    <-lmfit$coef
s2      <-summary(lmfit)$sigma^2
phi     <-acf(lmfit$res,plot=FALSE)$acf[2]

# previa
nu0  <-1
s20  <-1
T0   <-diag(1/1000,nrow=2)  # SIGMA_0^{-1}

###
set.seed(1)
S     <-1000    # numero de iteraciones
odens <-S/1000  # informacion
OUT   <-NULL    # almacenamiento
ac    <-0       # tasa de aceptacion

# cadena
for (s in 1:S) {
        # simular beta
        Cor    <- phi^DY  
        iCor   <- solve(Cor)
        V.beta <- solve( t(X)%*%iCor%*%X/s2 + T0)
        E.beta <- V.beta%*%( t(X)%*%iCor%*%y/s2  )
        beta   <- t(rmvnorm(1,E.beta,V.beta)  )
        # simular sigma^2
        s2 <- 1/rgamma(1,(nu0+n)/2,(nu0*s20+t(y-X%*%beta)%*%iCor%*%(y-X%*%beta)) /2 )
        # simular rho (metropolis)
        # 1. propuesta
        phi.p <- abs(runif(1,phi-.1,phi+.1))
        phi.p <- min(phi.p, 2-phi.p)
        # 2. tasa de aceptacion
        lr <- -.5*( determinant(phi.p^DY,log=TRUE)$mod - determinant(phi^DY,log=TRUE)$mod + 
                    tr( (y-X%*%beta)%*%t(y-X%*%beta)%*%(solve(phi.p^DY) - solve(phi^DY)) )/s2 )
        # 3. actualizar valor
        if( log(runif(1)) < lr ) { 
                phi<-phi.p
                ac<-ac+1 
        }
        # progreso & almacenar
        if(s%%odens==0) {
                cat(s,ac/s,beta,s2,phi,"\n") 
                OUT <- rbind(OUT,c(beta,s2,phi))
                par(mfrow=c(2,2))
                plot(OUT[,1], type = "l", ylab  = "beta 1")  ; abline(h=fit.gls$coef[1], col = 2, lty = 2)
                plot(OUT[,2], type = "l", ylab  = "beta 2")  ; abline(h=fit.gls$coef[2], col = 2, lty = 2)
                plot(OUT[,3], type = "l", ylab  = "sigma^2") ; abline(h=fit.gls$sigma^2, col = 2, lty = 2)
                plot(OUT[,4], type = "l", ylab  = "rho")     ; abline(h=.8284, col = 2, lty = 2)
        }
}
####### fin MCMC


################ 
# diagnosticos #
################

OUT.1000<-OUT
library(coda)
apply(OUT,2,effectiveSize)


load("OUT25000.RData")
apply(OUT.25000,2,effectiveSize )


mean(OUT.25000[,2])
quantile(OUT.25000[,2], c(0.025,0.5,0.975))


windows(height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
plot(OUT.1000[,4],xlab="scan",ylab=expression(rho),type="l")
acf(OUT.1000[,4],ci.col="gray",xlab="lag")


windows(height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
plot(OUT.25000[,4],xlab="scan/25",ylab=expression(rho),type="l")
acf(OUT.25000[,4],ci.col="gray",xlab="lag/25")


##############
# inferencia #
##############

windows(height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
plot(density(OUT.25000[,2],adj=2),xlab=expression(beta[2]),
   ylab="Densidad marginal posterior",main="")
plot(y~X[,2],xlab=expression(CO[2]),ylab="Temperatura")
abline(mean(OUT.25000[,1]),mean(OUT.25000[,2]),lwd=2)
abline(lmfit$coef,col="gray",lwd=2)
legend(180,2.5,legend=c("GLS","OLS"),bty="n",
      lwd=c(2,2),col=c("black","gray"))

quantile(OUT.25000[,2],probs=c(.025,.975) )

################################################################################

########
# JAGS #
########

# data
load("dct.RData")
y  <- dct[,3]
X  <- cbind(1,dct[,2])
n  <- dim(X)[1]
p  <- dim(X)[2] 
DY <- abs(outer( (1:n),(1:n) ,"-")) # para construir la matriz de correlacion

library(R2jags)

model <- function() {
    y[1:n] ~ dmnorm(X[1:n,1:p]%*%beta[1:p], phi*inverse(ilogit(rho)^DY[1:n,1:n]))
    beta[1:p] ~ dmnorm(beta0[1:p], Omega0[1:p,1:p])
    phi ~ dgamma(a0, b0)
    rho ~ dnorm(c0, d0)
}

# previa
beta0  <- c(rep(0,p))
Omega0 <- diag(1/1000,nrow=p)
a0     <- 1/2
b0     <- 1/2
c0     <- 0
d0     <- 1/1000

# input
model_data <- list(y = y, X = X, DY = DY, n = n, p = p, beta0 = beta0, Omega0 = Omega0, a0 = a0, b0 = b0, c0 = c0, d0 = d0) 

# parameters
model_parameters <- c("beta","phi","rho")

# initial values
initial_values <- list(list("beta" = c(-11.1174396, 0.0284304), "phi" = 1/5.2785807, "rho" = 1.524685), 
                       list("beta" = c(-11.1174396, 0.0284304), "phi" = 1/5.2785807, "rho" = 1.524685),
                       list("beta" = c(-11.1174396, 0.0284304), "phi" = 1/5.2785807, "rho" = 1.524685))

# mcmc settings
niter  <- 1100
nburn  <- 100
nthin  <- 1
nchain <- length(initial_values)

# mcmc
set.seed(123)
fit <- jags(data = model_data, inits = initial_values, 
            parameters.to.save = model_parameters, model.file = model, 
            n.chains = nchain, n.iter = niter, n.thin = nthin, n.burnin = nburn)

print(fit)

# transformar a objecto MCMC               
fit_mcmc <- coda::as.mcmc(fit)

# plots
windows()
mcmcplots::traplot(fit_mcmc)

windows()
mcmcplots::denplot(fit_mcmc)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------