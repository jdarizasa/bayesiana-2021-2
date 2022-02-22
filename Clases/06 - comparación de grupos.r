###########################
### Modelos jerarquicos ###
###########################

###########################
# PUNTAJES DE MATEM?TICAS #
###########################

# directorio de trabajo
setwd("C:/Users/Juan Camilo/Dropbox/UN/estadistica_bayesiana/codigo/")

# cargar algunas funciones
source("funciones_auxiliares.R")

# cargar data
load("math.RData")

# NOTACION
# ids  : identificador de los colegios (c)
# J    : numero de grupos (colegios)
# m    : numero de grupos (colegios)
# n    : numero de estudiantes en cada colegio (c)
# YM   : puntajes  de los estudiantes (matrix)
# Y    : puntajes  de los estudiantes (list)
# ybar : promedios de los puntajes (c)
# ymed : medianas  de los puntajes (c)
# s2   : varianzas de los puntajes (c)
# sv   : varianzas de los puntajes (c)

##################################
# ANALISIS EXPLORATORIO DE DATOS #
###################################

# representacion de los puntajes brutos
windows(height=3.5,width=7)
par(mfrow=c(1,1))
plot(c(1,J), range(Y), type="n", ylab="Puntaje", xlab="Ranking (promedio muestral)")
for(l in 1:J) {
        j<-order(ybar)[l]
        points( rep(l, n[j]), Y[[j]], pch=16, cex=.5 )
        segments( l, min(Y[[j]]), l, max(Y[[j]]) )
}
abline(h=mean(ybar))

# representacion de los promedios de los grupos
windows(height=3.5,width=7)
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0)) 
hist(ybar, freq = F, main="", xlab="Promedio", ylab="Densidad", 
     border="darkgray", col="gray95")
plot(n,ybar, xlab="Tama?o del grupo", ylab="Promedio")

#-------------------------------------------------------------------------------
# No es raro que los grupos con promedios muestrales muy altos o muy bajos sean 
# tambi?n aquellos grupos con tama?os muestrales bajos, ya que var(ybar) = sigma^2/n
#-------------------------------------------------------------------------------

########
# MCMC #
########

# hiper-parametros para obtener previas no informativos
nu0  <- 1  
s20  <- 100
eta0 <- 1  
t20  <- 100
mu0  <- 50 
g20  <- 25

#-------------------------------------------------------------------------------
# ELICITACION DE LA PREVIA
#
# - sigma_0^2 = 100 :
#
#   El examen de matem?ticas se diseñó para dar una varianza nacional de 100. 
#   Dado que esta varianza incluye la varianza tanto dentro de la escuela como 
#   entre escuelas, la varianza dentro de la escuela debe ser como m?ximo de 100, 
#   por lo que tomamos sigma_0^2 = 100.
#
# - nu_0 = 1 :
#
#   tau_0^2 = 100 es probable que sea una sobreestimaci?n, por lo que solo 
#   concentramos d?bilmente la distribuci?n previa alrededor de este valor 
#   al tomar nu_0 = 1. 
#
# - tau_0^2 = 100 and eta_0 = 1 :
#
#   De manera similar, la varianza entre escuelas no debe ser superior a 100, 
#   por lo que tomamos tau_0^2 = 100 y eta_0 = 1. 
#
# - mu_0 = 50 and gamma_0^2 = 25 :
#
#   La media nacional de todas las escuelas es 50. Aunque la media de las grandes
#   escuelas p?blicas urbanas puede ser diferente del promedio nacional, no deber?a 
#   diferir demasiado. Por eso tomamos mu_0 = 50 y gamma = 5, de forma que la 
#   probabiliad previa de que mu_0 este en el intervalo (40, 60) es aprox. 95%.
#-------------------------------------------------------------------------------

# 0. setup
S     <- 50000                            # numero de iteraciones
ncat  <- floor(S/10)                      # progreso
THETA <- matrix(data=NA, nrow=S, ncol=m)  # almacenar thetas 
MST   <- matrix(data=NA, nrow=S, ncol=3)  # almacenar mu, sigma^2, tau^2

# 1. valores inciales
theta  <- ybar
sigma2 <- mean(sv)
mu     <- mean(theta)
tau2   <- var(theta)

# 2. MCMC
set.seed(1)
for (s in 1:S) {
        # 2.1 actualizar theta
        for(j in 1:m) {
                vtheta   <- 1/(n[j]/sigma2+1/tau2)
                etheta   <- vtheta*(ybar[j]*n[j]/sigma2+mu/tau2)
                theta[j] <- rnorm(1,etheta,sqrt(vtheta))
        }
        # 2.2 actualizar sigma^2
        nun <- nu0+sum(n)
        ss  <- nu0*s20
        for (j in 1:m) {
                ss <- ss+sum((Y[[j]]-theta[j])^2)
        }
        sigma2 <- 1/rgamma(1,nun/2,ss/2)
        # 2.3 actualizar mu
        vmu <- 1/(m/tau2+1/g20)
        emu <- vmu*(m*mean(theta)/tau2 + mu0/g20)
        mu  <- rnorm(1,emu,sqrt(vmu)) 
        # 2.4 actualizar tau2
        etam <- eta0+m
        ss   <- eta0*t20 + sum( (theta-mu)^2 )
        tau2 <- 1/rgamma(1,etam/2,ss/2)
        # tau2 <- sample_tau2(m, eta0, t20, mu, theta)
        # 2.5 almacenar valores
        THETA[s,] <- theta
        MST[s,]   <- c(mu,sigma2,tau2)
        # 2.6 progreso
        if (s%%ncat == 0) cat(100*round(s/S, 1), "% completed ... \n", sep = "" )
} 
# FINAL DEL ALGORITMO

################
# DIAGNOSTICOS #
################

######## diagramas de caja
windows(height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
stationarity.plot(MST[,1],xlab="iteration",ylab=expression(mu))
stationarity.plot(MST[,2],xlab="iteration",ylab=expression(sigma^2))
stationarity.plot(MST[,3],xlab="iteration",ylab=expression(tau^2))
########


######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(MST[,1],xlab="iteration",type ="l",ylab=expression(mu))
plot(MST[,2],xlab="iteration",type ="l",ylab=expression(sigma^2))
plot(MST[,3],xlab="iteration",type ="l",ylab=expression(tau^2))
########


######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(THETA[,1],xlab="iteration",type ="l",ylab=expression(theta[1]))
plot(THETA[,2],xlab="iteration",type ="l",ylab=expression(theta[2]))
plot(THETA[,3],xlab="iteration",type ="l",ylab=expression(theta[3]))
########


######## autocorrelaciones
library(coda)
windows(height=1.75,width=5)
par(mfrow=c(1,3))
acf(MST[,1],main=expression(mu)) -> a1
acf(MST[,2],main=expression(sigma^2)) -> a2
acf(MST[,3],main=expression(tau^2)) -> a3
########


######## tama?os efectivos de muestra
effectiveSize(MST)
summary( effectiveSize(THETA) )
########


######## errores estandar de MC : DE(parametro)/sqrt( n_eff )
MCERR <- apply(MST,2,sd) / sqrt( effectiveSize(MST) )
apply(MST,2,mean)

TMCERR <- apply(THETA,2,sd)/sqrt( effectiveSize(THETA) )
summary( TMCERR )
########


########################
# INFERENCIA POSTERIOR #
########################


######## distribuciones posteriores para mu, sigma^2 and tau^2
windows(height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
# mu
plot(density(MST[,1],adj=2),xlab=expression(mu),main="",lwd=2,
ylab=expression(paste(italic("p("),mu,"|",italic(y[1]),"...",italic(y[m]),")")))
abline( v=quantile(MST[,1],c(.025,.5,.975)),col=c(2,4,2),lty=c(3,2,3) )
# sigma^2
plot(density(MST[,2],adj=2),xlab=expression(sigma^2),main="", lwd=2,
ylab=expression(paste(italic("p("),sigma^2,"|",italic(y[1]),"...",italic(y[m]),")")))
abline( v=quantile(MST[,2],c(.025,.5,.975)),col=c(2,4,2),lty=c(3,2,3) )
# tau^2
plot(density(MST[,3],adj=2),xlab=expression(tau^2),main="",lwd=2,
ylab=expression(paste(italic("p("),tau^2,"|",italic(y[1]),"...",italic(y[m]),")")))
abline( v=quantile(MST[,3],c(.025,.5,.975)),col=c(2,4,2),lty=c(3,2,3) )
#######


####### medias posteriores para mu, sigma, tau
mean(MST[,1])
mean(sqrt(MST[,2]))
mean(sqrt(MST[,3]))
#######

#-------------------------------------------------------------------------------
# aproximadamente el 99% de los puntajes dentro de una escuela est?n dentro de 
# 4 ? 9.21 = 37 puntos entre s?, mientras que el 99% de los puntajes promedio 
# de las escuelas est?n dentro de 4 ? 4.97 = 20 puntos el uno del otro.

# ademas se observa que la media global es aproximadamente 50, lo cual coincide 
# con el dise?o de la prueba
#-------------------------------------------------------------------------------


######## shrinkage
#AQUI SE VE LA DIFERENCIA ENTRE COMPARTIR LA INFO ENTRE COLEGIOS        
windows(height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
#
theta.hat<-apply(THETA,2,mean)
plot(ybar,theta.hat,xlab=expression(bar(italic(y))),ylab=expression(hat(theta)))
abline(0,1)
#
plot(n,ybar-theta.hat,ylab=expression( bar(italic(y))-hat(theta) ),xlab="sample size")
abline(h=0)
########

#-------------------------------------------------------------------------------
# la relaci?n sigue aproximadamente una l?nea con una pendiente menor que uno, 
# lo que indica que los valores altos de \bar{y}_j corresponden a valores 
# ligeramente menos altos \hat{\theta}_j
#
# los grupos con tama?os de muestra bajos son los que m?s se "reducen" o "encogen", 
# (shrunk) mientras que los grupos con tama?os de muestra grandes casi no se 
# "reducen" en absoluto.
#
# Cuanto mayor sea el tama?o de la muestra de un grupo, m?s informaci?n tenemos 
# para ese grupo y menos informaci?n necesitamos "tomar prestada" del resto de 
# la poblaci?n.
#-------------------------------------------------------------------------------


######## ranking
theta.order<-order(theta.hat)
theta.order[1:20]

ybar.order<-order(ybar)
ybar.order[1:20]

ybar[c(46,82)]
theta.hat[c(46,82)]

n[c(46,82)]

mean( THETA[,46] < THETA[,82])  # Pr (theta_46 < theta_82 )
########


######## ranking de colegios 46 y 82
windows(height=3.5,width=7)
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,1))
plot(density(THETA[,46],adj=2),col="black",
     xlim=range(c(Y[[46]],Y[[82]],THETA[,c(46,82)])),lwd=2,
     main="",xlab="math score",ylim=c(-.05,.2),ylab="",yaxt="n")
axis(side=2,at=c(0,0.10,0.20) )
lines(density(THETA[,82],adj=2),col="gray",lwd=2)
abline(h=0)
#
points( Y[[46]],rep(-0.01666,n[46]), col="black",pch=16)
points( ybar[46],-.01666,col="black",pch=16 ,cex=1.5)
abline( h=-.01666,col="black")
#
points( Y[[82]],rep(-0.0333,n[82]), col="gray",pch=16)
points( ybar[82],-.0333,col="gray",pch=16 ,cex=1.5)
abline( h=-.0333,col="gray")
#
segments(mean(MST[,1]), 0,mean(MST[,1]),1,lwd=2,lty=2 )
#
legend(52.5,.15,
       legend=c("Colegio 46","Colegio 82",expression(paste("E[", mu,"|",italic(y[1]),"...",italic(y[m]),"]"))),
       lwd=c(2,2),lty=c(1,1,2),col=c("black","gray"),bty="n")
########

#-------------------------------------------------------------------------------
# Los datos brutos para las dos escuelas se muestran en gr?ficos de puntos debajo 
# de las densidades posteriores; los puntos grandes representan las medias del grupo
#
# la densidad posterior del colegio 46 es m?s elevada que la del colegio 82. 
# esto se debe a que el tama?o de la muestra de la escuela 46 es de 21 estudiantes, 
# mientras que la de la escuela 82 es de solo 5 estudiantes. Por lo tanto, nuestro 
# grado de certeza para theta_46 es mucho mayor que el de theta_82.
#
# rankings basados en medias muestrales y medias posteriores no coinciden necesariamente
#
# hay m?s evidencia de que theta_46 es excepcionalmente bajo que hay evidencia 
# de que theta_82 es excepcionalmente bajo.
#
#------------------------------------------------------------------------------- 


###########################################
# MODELO JERARQUICO PARA MEDIA Y VAIRANZA #
###########################################

######## hiper-parametros para obtener una previa no informativa
eta0 <- 1  ; t20 <- 100
mu0  <- 50 ; g20 <- 25
a0   <- 1  ; b0  <-  1/100 ; wnu0 <- 1
########


######## valores iniciales
theta  <- ybar
sigma2 <- sv 
mu     <- mean(theta)
tau2   <- var(theta)
s20    <- 100
nu0    <- 1
########


######## setup
S      <- 50000
ncat   <- floor(S/10)
THETA  <- matrix(data=NA, nrow=S, ncol=m)
SIGMA2 <- matrix(data=NA, nrow=S, ncol=m)
MTSN   <- matrix(data=NA, nrow=S, ncol=4)  # mu, tau2, sigma02, nu0
# s2.pp  <- NULL
nu0s   <- 1:5000 # rango de valores para muestrear en p(nu_0 | rest)
########


######## MCMC
set.seed(1)
for(s in 1:S) {
        # actualizar thetas
        for(j in 1:m) {
                vtheta   <- 1/(n[j]/sigma2[j]+1/tau2)
                etheta   <- vtheta*(ybar[j]*n[j]/sigma2[j]+mu/tau2)
                theta[j] <- rnorm(1,etheta,sqrt(vtheta))
        }
        # actualizar sigma2s
        for(j in 1:m) { 
                nun       <- nu0+n[j]
                ss        <- nu0*s20+ sum((Y[[j]]-theta[j])^2)
                sigma2[j] <- 1/rgamma(1,nun/2,ss/2)
        }
        # actualizar s20
        s20 <- rgamma(1,a0+m*nu0/2,b0+nu0*sum(1/sigma2)/2)
        # actualizar nu0
        lpnu0 <- .5*nu0s*m*log(s20*nu0s/2)-m*lgamma(nu0s/2)+(nu0s/2-1)*sum(log(1/sigma2)) - nu0s*s20*sum(1/sigma2)/2 - wnu0*nu0s
        nu0   <- sample(x = nu0s, size = 1, prob=exp( lpnu0-max(lpnu0) )) #truco de restar por el maximo
        # actualizar mu
        vmu <- 1/(m/tau2+1/g20)
        emu <- vmu*(m*mean(theta)/tau2 + mu0/g20)
        mu  <- rnorm(1,emu,sqrt(vmu))
        # actualizar tau2
        etam <-eta0+m
        ss   <- eta0*t20 + sum( (theta-mu)^2 )
        tau2 <-1/rgamma(1,etam/2,ss/2)
        # almacenar
        THETA[s,]  <- theta
        SIGMA2[s,] <- sigma2
        MTSN[s,]   <- c(mu,tau2,s20,nu0)
        # s2.pp      <- c(s2.pp,1/rgamma(1,nu0/2,nu0*s20/2))
        ### progreso
        if (s%%ncat == 0) cat(100*round(s/S, 1), "% completed ... \n", sep = "" )
}
######## FINAL DEL ALGORITMO



######## comparar valores de la media posterior de sigma^2_j con las varianzas muestrales
apply(SIGMA2,2,mean) -> sigma2.hat

windows(height=3.5,width=7)
par(mfrow=c(1,2),mar=c(3,3,1,1)+.3,mgp=c(1.75,.75,0))
#
plot(sv,sigma2.hat,xlab=expression(s^2),ylab=expression(hat( sigma^2)) )
abline(0,1)
#
plot(n, sv-sigma2.hat,xlab="sample size",ylab=expression(s^2-hat(sigma^2)))  
abline(h=0)
########


########
windows(height=3.5,width=6)
par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
#
plot(density(MST[,1],adj=2),lwd=2,main="",col="gray",xlab=expression(mu),
ylab=expression(paste(italic("p("),mu,"|",italic(y[1]),"...",italic(y[m]),")")))
lines(density(MTSN[,1],adj=2),lwd=2)
#
plot(density(MST[,3],adj=2),lwd=2,main="",col="gray",xlab=expression(tau^2), 
ylab=expression(paste(italic("p("),tau^2,"|",italic(y[1]),"...",italic(y[m]),")")))
lines(density(MTSN[,2],adj=2),lwd=2)
#
plot(table(MTSN[,4]),xlab=expression(nu[0]),
ylab=expression(paste(italic("p("),nu[0],"|",italic(y[1]),"...",italic(y[m]),")")))
#
plot(density(MTSN[,3],adj=2),lwd=2,main="",xlab=expression(sigma[0]^2),
ylab=expression(paste(italic("p("),sigma[0]^2,"|",italic(y[1]),"...",italic(y[m]),")")))
########



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------