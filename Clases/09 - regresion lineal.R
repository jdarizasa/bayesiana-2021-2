####################
# REGRESION LINEAL #
####################

# Estudio del efectos de dos reg?menes de ejercicio sobre la absorci?n de ox?geno.
#
# Seis de doce hombres fueron asignados aleatoriamente a un programa de 
# carrera en terreno plano de 12 semanas, y los seis restantes fueron asignados 
# a un programa de aer?bicos de 12 semanas.
#
# Se midi? el consumo m?ximo de ox?geno de cada sujeto (en litros por minuto) 
# al correr en una cinta de correr inclinada, tanto antes como despu?s del 
# programa de 12 semanas. 
# 
# El objetivo consiste en evaluar c?mo el cambio en la absorci?n m?xima de ox?geno 
# depende del programa de entrenamiento.
#
# y    : cambio en la absorci?n m?xima de ox?geno.
# trat : programa de entrenamiento (1 si aer?bicos, 0 si correr). 
# edad : edad (en a?os).
# n    : tama?o de la muestra.
#
# Modelo: E(y | X) = beta1*x1 + beta2*x2 + beta3*x3 + beta4*x4
#         x1 = 1 
#         x2 = 1 si aer?bicos, 0 si correr
#         x3 = edad
#         x4 = x2*x3
#
# Correr:    E(y | X) = beta1 + beta3*edad
# Aer?bicos: E(y | X) = (beta1 + beta2) + (beta3 + beta4)*edad

# data
trat <- c(0,0,0,0,0,0,1,1,1,1,1,1)
edad <- c(23,22,22,25,27,20,31,23,27,28,22,24)
y    <- c(-0.87,-10.74,-3.27,-1.97,7.50,-7.25,17.05,4.96,10.40,11.05,0.26,2.51)

# data
y <- as.matrix(y)
X <- cbind(1, trat, edad, trat*edad)

# dimensiones
n <- dim(X)[1]
p <- dim(X)[2]
colnames(X) <- paste0("x", 1:p)

# scatterplot
windows(height=3.5, width=7)
par(mfrow=c(1,1), mar=c(3,3,1,1), mgp=c(1.75,.75,0))
plot(y ~ edad, pch=16, xlab = "Edad", ylab = "Cambio en la absorci?n", 
     col = c("black","gray")[trat+1])
legend("topleft", legend=c("Aer?bicos","Correr"), pch=c(16,16), col=c("gray","black"))

# OLS
beta.ols <- solve(t(X)%*%X)%*%t(X)%*%y
sig2.ols <- sum((y - X%*%beta.ols)^2)/(n-p)

# fit.ols <- lm(y ~ -1 + X)
# summary(fit.ols)$coef
# summary(fit.ols)$sigm

# hiperparametros (previa g)
nu0 <- 1
s20 <- sig2.ols
g   <- n

# calculos para la posterior
Hg    <- (g/(g+1))*(X%*%solve(t(X)%*%X)%*%t(X))
SSRg  <- t(y)%*%(diag(1,n) - Hg)%*%y
Vbeta <- (g/(g+1))*solve(t(X)%*%X)
Ebeta <- Vbeta%*%t(X)%*%y

# Monte Carlo
S <- 5000
s2.post   <- matrix(NA, S, 1) 
beta.post <- matrix(NA, S, p) 
set.seed(1)
for(s in 1:S) {
      s2.post[s] <- 1/rgamma(1, (nu0 + n)/2, (nu0*s20 + SSRg)/2)
      beta.post[s,] <- c(mvtnorm::rmvnorm(1, Ebeta, s2.post[s]*Vbeta))
}

# inferencia sigma^2
quantile(s2.post, probs = c(.025,.5,.975))

# inferencia sigma
quantile(sqrt(s2.post), probs = c(.025,.5,.975))

# inferencia beta
apply(X = beta.post, MARGIN = 2, FUN = quantile, probs = c(.025,.5,.975))
apply(X = beta.post, MARGIN = 2, FUN = mean)
apply(X = beta.post, MARGIN = 2, FUN = sd)
colMeans(beta.post > 0)

# posterior beta2 y beta4
windows(height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
# beta2 posterior
x<-seq(-85,130,length=200)
plot(density(beta.post[,2],adj=2), xlab=expression(beta[2]),main="",ylab="",lwd=2)
abline(v=0,col="gray")
abline(v=quantile(beta.post[,2], c(0.025,0.975)),col="blue",lty=2)
# beta4 posterior
x<-seq(-5,5,length=100)
plot(density(beta.post[,4],adj=2),xlab=expression(beta[4]),main="",ylab="",lwd=2)
abline(v=0,col="gray")
abline(v=quantile(beta.post[,4], c(0.025,0.975)),col="blue",lty=2)
# (beta2, beta4) posterior
plot(beta.post[,c(2,4)],type = "p",cex=0.1,xlab=expression(beta[2]),ylab=expression(beta[4]))
abline(h=0,col="gray") ; abline(v=0,col="gray")

# Efecto del tratamiento aerobico: Delta
#
# Delta = E(y | edad, aerobico) - E(y | edad, aerobico)
#       = ( (beta1 + beta2) + (beta3 + beta4)*edad ) - ( beta1 + beta3*edad )
#       = beta2 + beta4*edad
r.edad <- min(edad):max(edad)
n.edad <- length(r.edad)
TE <- matrix(NA, S, n.edad)
for (j in 1:n.edad){
      TE[,j] <- beta.post[,2] + beta.post[,4]*r.edad[j]
}
that <- colMeans(TE)
ic1  <- apply(X = TE, MARGIN = 2, FUN = function(x) quantile(x, c(0.050,0.950)))
ic2  <- apply(X = TE, MARGIN = 2, FUN = function(x) quantile(x, c(0.025,0.975)))
colo <- c("blue","black")[as.numeric((ic2[1,] < 0) & (0 < ic2[2,]))+1]

windows(height=3.5,width=7)
par(mfrow=c(1,1),mar=c(3,3,1.5,1),mgp=c(1.75,.75,0))
plot(NA, NA, xlab = "Edad", ylab = "Efecto tratamiento", main = "", 
     xlim = range(edad), ylim = range(TE), cex.axis = 0.75, xaxt = "n")
axis(side = 1, at = r.edad, labels = r.edad, cex.axis = 0.75)
abline(h = 0, col = "gray", lwd = 2)
abline(v = r.edad, col = "gray95", lwd = 1, lty = 2) 
for (j in 1:n.edad) {
      segments(x0 = r.edad[j], y0 = ic1[1,j], x1 = r.edad[j], y1 = ic1[2,j], lwd = 3, col = colo[j])
      segments(x0 = r.edad[j], y0 = ic2[1,j], x1 = r.edad[j], y1 = ic2[2,j], lwd = 1, col = colo[j])
      lines(x = r.edad[j], y = that[j], type = "p", pch = 16, cex = 0.8, col = colo[j])
}

################################################################################

# descarcar JAGS
# http://www.sourceforge.net/projects/mcmc-jags/files

# user manuals
# https://people.stat.sc.edu/hansont/stat740/jags_user_manual.pdf
# http://www.jkarreth.net/files/bayes-cph_Tutorial-JAGS.pdf 

library(R2jags)

# data
trat <- c(0,0,0,0,0,0,1,1,1,1,1,1)
edad <- c(23,22,22,25,27,20,31,23,27,28,22,24)
y    <- c(-0.87,-10.74,-3.27,-1.97,7.50,-7.25,17.05,4.96,10.40,11.05,0.26,2.51)

# data
y <- c(as.matrix(y))
X <- cbind(1, trat, edad, trat*edad)

# dimensiones
n <- dim(X)[1]
p <- dim(X)[2]
colnames(X) <- paste0("x", 1:p)

# la Normal est? parametrizada en t?rminos de la precisi?n
# los vectores deben tener el mismo formato
# notaci?n:
#     phi    = 1/sigma^2
#     Omega0 = X^T*X
#     a0     = nu0/2
#     b0     = nu0*sigma0^2/2

model <- function() {
      for (i in 1:n) {
            y[i] ~ dnorm(inprod(X[i,], beta), phi)
      }
      beta[1:p] ~ dmnorm(beta0[1:p], (phi/g)*Omega0[1:p,1:p])
      phi ~ dgamma(a0, b0)
}

# previa
beta0  <- c(rep(0,p))
Omega0 <- t(X)%*%X
g      <- n
nu0    <- 1
s20    <- sig2.ols
a0     <- nu0/2
b0     <- nu0*s20/2

# input
model_data <- list(y = y, X = X, n = n, p = p, g = g, beta0 = beta0, Omega0 = Omega0, a0 = a0, b0 = b0)

# parameters
model_parameters <- c("beta", "phi")

# initial values
initial_values <- list(list("beta" = beta0, "phi" = 1/s20), 
                       list("beta" = beta0, "phi" = 1/s20),
                       list("beta" = beta0, "phi" = 1/s20))

# mcmc settings
niter  <- 6000
nburn  <- 1000
nthin  <- 1
nchain <- length(initial_values)

# mcmc
set.seed(123)
fit <- jags(data = model_data, inits = initial_values, 
            parameters.to.save = model_parameters, model.file = model, 
            n.chains = nchain, n.iter = niter, n.thin = nthin, n.burnin = nburn)

print(fit)

# numero de muestras
S <- fit$BUGSoutput$n.sims

# dic
fit$BUGSoutput$DIC

# transformar a objecto MCMC               
fit_mcmc <- coda::as.mcmc(fit)

# diagnosticos
mcmcplots::mcmcplot(fit_mcmc)
superdiag::superdiag(fit_mcmc)

# plots
windows()
mcmcplots::denplot(fit_mcmc)

windows()
mcmcplots::traplot(fit_mcmc)

windows()
mcmcplots::caterplot(fit_mcmc, parms = c("beta[2]","beta[4]"))

# mas opciones ...
fit_mcmc_gg <- ggmcmc::ggs(fit_mcmc)

windows()
ggmcmc::ggs_density(fit_mcmc_gg)

