############################
#### MODELO 1 - POISSON ####
############################

# Cargar datos 
y <- c(1, 3, 2, 12, 1, 1) # variable dependiente (respuesta)
x <- c(33, 14, 27, 90, 12, 17) # variable independiente 

n <- length(y)                    # tamaño de la muestrsa
p <- dim(X)[2]                    # numero de predictores

#hiperparametros
a1 <- 1
b1 <- 1
a2 <- 10
b2 <- 1

# initial values
theta = rep(1, n)
alpha = 1
beta = 1

######## algoritmo de metropolis
S    <- 251000
THETA <- matrix(NA, nrow=S, ncol=n)
AB <- matrix(NA, nrow = S, ncol = 2)
LV1 <- matrix(NA, nrow = S, ncol = 1)
ac   <- 0
ncat <- floor(S/10)

###### samplers
sample_theta <- function(alpha, beta, y, x){
  for (i in 1:n) {
    a_t <- alpha + y[i]
    b_t <- beta + x[i]
    theta[i] <- rgamma(1, a_t, b_t)
  }
  return(theta)
}

sample_beta <- function(a2, b2, n, alpha, theta){
  a_b <- a2 + n*alpha
  b_b <- b2 + sum(theta)
  beta <- rgamma(1, a_b, b_b)
  return(beta)
}

delta <- 1

set.seed(7)
# cadena transformación
for (s in 1:S) {
  # simular theta
  theta <- sample_theta(alpha = alpha, beta = beta, y = y, x = x)
  # simular beta
  beta <- sample_beta(a2 = a2, b2 = b2, n = n, alpha = alpha, theta = theta)
  # simular alpha (metropolis)
  # 1. propuesta
  kappa <- log(alpha)
  kappa.p <- rnorm(1,kappa,sqrt(delta))
  # 2. tasa de aceptacion
  lr <- sum(dgamma(theta, exp(kappa.p), beta, log=T)) + dgamma(exp(kappa.p), a1, b1, log = T) + kappa.p -
    sum(dgamma(theta, exp(kappa), beta, log=T)) - dgamma(exp(kappa), a1, b1, log = T) - kappa
  # 3. actualizar valor
  if( log(runif(1)) < lr ) { 
    alpha <- exp(kappa.p)
    ac <- ac+1 
  }
  # Almacenar
  THETA[s,] <- theta
  AB[s,] <- c(alpha, beta)
  # Log - Verosimilitud
  LV1[s] <- sum(dpois(x = y, lambda = theta*x, log = T))
  # Progreso
  if (s%%ncat == 0) cat(100*round(s/S, 1), "% completed ... \n", sep = "" )
  
}
####### fin MCMC
# tasa de aceptacion
100*ac/S

# diagnosticos
library(coda)
apply(X = THETA, MARGIN = 2, FUN = effectiveSize)
apply(X = AB, MARGIN = 2, FUN = effectiveSize)
effectiveSize(LV1)

######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,1))
plot(LV1,xlab="iteration",type ="l",ylab="deviance", main = "")

######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,1),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(AB[,1],xlab="iteration",type ="l",ylab=expression(alpha), main = "Convergencia Alpha")

windows(height=1.75,width=5)
par(mfrow=c(1,1),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(AB[,2],xlab="iteration",type ="l",ylab=expression(beta), main = "Convergencia Beta")

######## cadenas
windows()
par(mfrow=c(2,3))
plot(THETA[,1],xlab="iteration",type ="l",ylab=expression(theta[1]), main = "Convergencia ")
plot(THETA[,2],xlab="iteration",type ="l",ylab=expression(theta[2]), main = "Convergencia ")
plot(THETA[,3],xlab="iteration",type ="l",ylab=expression(theta[3]), main = "Convergencia ")
plot(THETA[,4],xlab="iteration",type ="l",ylab=expression(theta[4]), main = "Convergencia ")
plot(THETA[,5],xlab="iteration",type ="l",ylab=expression(theta[5]), main = "Convergencia ")
plot(THETA[,6],xlab="iteration",type ="l",ylab=expression(theta[6]), main = "Convergencia ")
####

THETA <- THETA[-c(1:1000),]
AB <- AB[-c(1:1000),]
LV1 <- LV1[-c(1:1000),]

#autocorrelaciones
library(coda)
windows(height=3,width=15)
par(mfrow=c(1,6))
acf(THETA[,1],main=expression(theta1), ci.col="gray",xlab="lag") 
acf(THETA[,2],main=expression(theta2), ci.col="gray",xlab="lag") 
acf(THETA[,3],main=expression(theta3), ci.col="gray",xlab="lag") 
acf(THETA[,4],main=expression(theta4), ci.col="gray",xlab="lag") 
acf(THETA[,5],main=expression(theta5), ci.col="gray",xlab="lag") 
acf(THETA[,6],main=expression(theta6), ci.col="gray",xlab="lag") 

windows(height=3,width=10)
par(mfrow=c(1,3))
acf(AB[,1],main=expression(alpha), ci.col="gray",xlab="lag") 
acf(AB[,2],main=expression(beta), ci.col="gray",xlab="lag") 

#### Adelgazamiento 
blabs <- c(expression(theta[1]),expression(theta[2]),expression(theta[3]),
           expression(theta[4]),expression(theta[5]),expression(theta[6]))  # etiquetas
thin  <- c(1,(1:20000)*(S/25100))  # muestreo sistematico 

#theta1
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,THETA[thin,1],type="l",xlab="Iteración",ylab=blabs[1])
abline(h=mean(THETA[,1]))
acf(THETA[,1],ci.col="gray",xlab="lag")
acf(THETA[thin,1],xlab="lag/10",ci.col="gray")

#theta2
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,THETA[thin,2],type="l",xlab="Iteración",ylab=blabs[2])
abline(h=mean(THETA[,2]))
acf(THETA[,2],ci.col="gray",xlab="lag")
acf(THETA[thin,2],xlab="lag/10",ci.col="gray")

#theta3
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,THETA[thin,3],type="l",xlab="Iteración",ylab=blabs[3])
abline(h=mean(THETA[,3]))
acf(THETA[,3],ci.col="gray",xlab="lag")
acf(THETA[thin,3],xlab="lag/10",ci.col="gray")

#theta4
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,THETA[thin,4],type="l",xlab="Iteración",ylab=blabs[4])
abline(h=mean(THETA[,4]))
acf(THETA[,4],ci.col="gray",xlab="lag")
acf(THETA[thin,4],xlab="lag/10",ci.col="gray")

#theta5
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,THETA[thin,5],type="l",xlab="Iteración",ylab=blabs[5])
abline(h=mean(THETA[,5]))
acf(THETA[,5],ci.col="gray",xlab="lag")
acf(THETA[thin,5],xlab="lag/10",ci.col="gray")

#theta6
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,THETA[thin,6],type="l",xlab="Iteración",ylab=blabs[6])
abline(h=mean(THETA[,6]))
acf(THETA[,6],ci.col="gray",xlab="lag")
acf(THETA[thin,6],xlab="lag/10",ci.col="gray")

#alpha
blabs <- c(expression(alpha))  # etiquetas
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,AB[thin,1],type="l",xlab="Iteración",ylab=blabs[1])
abline(h=mean(AB[,1]))
acf(AB[,1],ci.col="gray",xlab="lag")
acf(AB[thin,1],xlab="lag/10",ci.col="gray")
A <- AB[thin,1]
effectiveSize(A)

#beta
blabs <- c(expression(beta))  # etiquetas
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,AB[thin,2],type="l",xlab="Iteración",ylab=blabs[1])
abline(h=mean(AB[,2]))
acf(AB[,2],ci.col="gray",xlab="lag")
acf(AB[thin,2],xlab="lag/10",ci.col="gray")


####
#PARA THETAS
######## errores estandar de MC : DE(parametro)/sqrt( n_eff )
MCERR <- apply(THETA,2,sd) / sqrt( effectiveSize(THETA) )
MCERR/apply(THETA,2,mean) #dan muy pequeños entonces no hay problema

#PARA ALPHA Y BETA
esA <- sd(A)/sqrt(effectiveSize(A))
esA/mean(A)

MCERRb <- sd(AB[,2]) / sqrt( effectiveSize(AB[,2]) )
MCERRb/mean(AB) #dan muy pequeños entonces no hay problema

#PARA DEVIANCE
esd <- sd(LV1)/sqrt(effectiveSize(LV1))
esd/abs(mean(LV1))

#### RESUMEN POSTERIOR ####
resumen_posterior <- function(x) 
{
  round(c(mean(x), sd(x), quantile(x = x, probs = c(0.025,0.975))), 3)
}

tab <- rbind(resumen_posterior(A),
             resumen_posterior(AB[,2]),
             resumen_posterior(THETA[,1]),
             resumen_posterior(THETA[,2]),
             resumen_posterior(THETA[,3]),
             resumen_posterior(THETA[,4]),
             resumen_posterior(THETA[,5]),
             resumen_posterior(THETA[,6]))
colnames(tab) <- c("Media", "SD", "Q2.5%", "Q97.5%")
rownames(tab) <- c("alpha", "beta",paste0("theta", 1:n))
print(tab)

#### COMPARACION THETA Y OLS ####
tab2 <- rbind(c(mean(THETA[,1]), y[1]/x[1]),
              c(mean(THETA[,2]), y[2]/x[2]),
              c(mean(THETA[,3]), y[3]/x[3]),
              c(mean(THETA[,4]), y[4]/x[4]),
              c(mean(THETA[,5]), y[5]/x[5]),
              c(mean(THETA[,6]), y[6]/x[6]))
colnames(tab2) <- c("Posterior", "OLS")
rownames(tab2) <- c(paste0("theta", 1:n))

print(tab2)

windows()
par(mfrow=c(2,3))
# grafico theta1
#windows()
hist(x =THETA[,1] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(THETA[,1]))
abline(v = y[1]/x[1], col = "red", lwd = 2, lty = 3)
abline(v=c(quantile(THETA[,1], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(THETA[,1]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))


# grafico theta2
#windows()
hist(x =THETA[,2] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(THETA[,2]))
abline(v = y[2]/x[2], col = "red", lwd = 2, lty = 3)
abline(v=c(quantile(THETA[,2], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(THETA[,2]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

# grafico theta3
#windows()
hist(x =THETA[,3] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(THETA[,3]))
abline(v = y[3]/x[3], col = "red", lwd = 2, lty = 3)
abline(v=c(quantile(THETA[,3], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(THETA[,3]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

# grafico theta4
#windows()
hist(x =THETA[,4] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(THETA[,4]))
abline(v = y[4]/x[4], col = "red", lwd = 2, lty = 3)
abline(v=c(quantile(THETA[,4], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(THETA[,4]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

# grafico theta5
#windows()
hist(x =THETA[,5] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(THETA[,5]))
abline(v = y[5]/x[5], col = "red", lwd = 2, lty = 3)
abline(v=c(quantile(THETA[,5], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(THETA[,5]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

# grafico theta6
#windows()
hist(x =THETA[,6] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(THETA[,6]))
abline(v = y[6]/x[6], col = "red", lwd = 2, lty = 3)
abline(v=c(quantile(THETA[,6], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(THETA[,6]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

#### COMPARACION ALPHA Y BETA ####
mean(AB[,1])/mean(AB[,2])
mean(y/x)

propAB <- AB[,1]/AB[,2]

hist(x = propAB, freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(propAB))
abline(v = mean(y/x), col = "red", lwd = 2, lty = 3)
abline(v=c(quantile(propAB, c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(propAB), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación media de y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

#### COMPARACION ENTRE THETAS ####
windows()
par(mfrow=c(2,3))
# dispersión entre thetas
plot(x =THETA[,2], y = THETA[,1], 
     xlab = expression(theta[2]), ylab = expression(theta[1]), main = "", col = "grey", type = "p",
     xlim = c(0,0.8), ylim = c(0,0.45))
points(x=y[2]/x[2], y = y[1]/x[1], pch=15)
abline(0,1, col="red")
legend("topright", legend = c("Observación de y1/x1 vs y2/x2"), 
       fill=1, border = "black")


plot(x =THETA[,2], y = THETA[,3], 
     xlab = expression(theta[2]), ylab = expression(theta[3]), main = "", col = "grey", type = "p",
     xlim = c(0,0.8), ylim = c(0,0.45))
points(x=y[2]/x[2], y = y[3]/x[3], pch=15)
abline(0,1, col="red")
legend("topright", legend = c("Observación de y3/x3 vs y2/x2"), 
       fill=1, border = "black")

plot(x =THETA[,2], y = THETA[,4], 
     xlab = expression(theta[2]), ylab = expression(theta[4]), main = "", col = "grey", type = "p",
     xlim = c(0,0.8), ylim = c(0,0.45))
points(x=y[2]/x[2], y = y[4]/x[4], pch=15)
abline(0,1, col="red")
legend("topright", legend = c("Observación de y4/x4 vs y2/x2"), 
       fill=1, border = "black")

plot(x =THETA[,2], y = THETA[,5], 
     xlab = expression(theta[2]), ylab = expression(theta[5]), main = "", col = "grey", type = "p",
     xlim = c(0,0.8), ylim = c(0,0.45))
points(x=y[2]/x[2], y = y[5]/x[5], pch=15)
abline(0,1, col="red")
legend("topright", legend = c("Observación de y5/x5 vs y2/x2"), 
       fill=1, border = "black")

plot(x =THETA[,2], y = THETA[,6], 
     xlab = expression(theta[2]), ylab = expression(theta[6]), main = "", col = "grey", type = "p",
     xlim = c(0,0.8), ylim = c(0,0.45))
points(x=y[2]/x[2], y = y[6]/x[6], pch=15)
abline(0,1, col="red")
legend("topright", legend = c("Observación de y6/x6 vs y2/x2"), 
       fill=1, border = "black")

## Probabilidades
#Pr(theta2 > thetaj)

mean(THETA[,2] > THETA[,1])
mean(THETA[,2] > THETA[,3])
mean(THETA[,2] > THETA[,4])
mean(THETA[,2] > THETA[,5])
mean(THETA[,2] > THETA[,6])

#Pr(theta2 = max(theta))

mean(THETA[,2] == apply(THETA, 1, max))

y/x
#### DIC ####
theta_hat <- colMeans(THETA)
lpyth_m1 <- sum(dpois(x = y, lambda = theta_hat*x, log = T))
pDIC_m1  <- 2*(lpyth_m1 - mean(LV1))
dic_m1   <- -2*lpyth_m1 + 2*pDIC_m1 

#### BONDAD DE AJUSTE ####
S <- length(THETA[,1])
ybar_hat <- rep(NA, S)
sd_hat <- rep(NA, S)
set.seed(7)
for (i in 1:S) {
  y_pred <- rpois(n, THETA[i,]*x)
  ybar_hat[i] <- mean(y_pred)
  sd_hat[i] <- sd(y_pred)
  if (i%%ncat == 0) cat(100*round(i/S, 1), "% completed ... \n", sep = "" )
  
}

windows()
par(mfrow=c(1,2))
hist(x = ybar_hat, freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(ybar_hat))
abline(v = mean(y), col = "red", lwd = 2, lty = 3)
legend("topright", legend = expression(bar(y)), 
       bty = "n", lwd = 2, lty = 3 ,col = "red")


hist(x = sd_hat, freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(sd_hat))
abline(v = sd(y), col = "red", lwd = 2, lty = 3)
legend("topright", legend = c("sd(y)"), 
       bty = "n", lwd = 2, lty = c(3) ,col = c("red"))

# ppp
mean( ybar_hat > mean(y))
mean( sd_hat > sd(y))

######################################
#### MODELO 2 - BINOMIAL NEGATIVA ####
######################################

# Cargar datos 
y <- c(1, 3, 2, 12, 1, 1) # variable dependiente (respuesta)
x <- c(33, 14, 27, 90, 12, 17) # variable independiente 

X <- cbind(rep(1,length(y)),x)    # matriz de diseño
n <- length(y)                    # tamaño de la muestrsa
p <- dim(X)[2]                    # numero de predictores

#hiperparametros
a1 <- 1
b1 <- 1
a2 <- 10
b2 <- 1
a3 <- 0
b3 <- 10

# initial values
theta = rep(1, n)
alpha = 1
beta = 1
lambda = 1

######## algoritmo de metropolis
S    <- 251000
THETA <- matrix(NA, nrow=S, ncol=n)
ABL <- matrix(NA, nrow = S, ncol = 3)
LV2 <- matrix(NA, nrow = S, ncol = 1)
ac_t   <- rep(0,n)
ac_a <- 0
ac_l <- 0
ncat <- floor(S/10)

###### samplers
sample_beta <- function(a2, b2, n, alpha, theta){
  a_b <- a2 + n*alpha
  b_b <- b2 + sum(theta)
  beta <- rgamma(1, a_b, b_b)
  return(beta)
}

delta_a <- 1
delta_l <- 6
delta_t <- c(1.5, 1.5, 1.5, 1, 1.5, 1.5)

set.seed(7)
# cadena transformacion
for (s in 1:S) {
  ## simular beta
  beta <- sample_beta(a2 = a2, b2 = b2, n = n, alpha = alpha, theta = theta)
  ## simular theta (metropolis)
  for (i in 1:n) {
    # 1. propuesta
    phi <- log(theta[i])
    phi.p <- rnorm(1, phi, delta_t[i])
    # 2. tasa de aceptacion
    lr_t <- dnbinom(x = y[i], mu = exp(phi.p)*x[i], size =  (exp(phi.p)*x[i])/lambda, log = T) + dgamma(exp(phi.p), alpha, beta, log = T) + phi.p -
      dnbinom(x = y[i], mu = exp(phi)*x[i], size = (exp(phi)*x[i])/lambda, log = T) - dgamma(exp(phi), alpha, beta, log = T) - phi
    # 3. actualizar valor
    if( log(runif(1)) < lr_t ) { 
      theta[i] <- exp(phi.p)
      ac_t[i] <- ac_t[i]+1 
    }
  }
  ## simular alpha (metropolis)
  # 1. propuesta
  kappa <- log(alpha)
  kappa.p <- rnorm(1,kappa,sqrt(delta_a))
  # 2. tasa de aceptacion
  lr <- sum(dgamma(theta, exp(kappa.p), beta, log=T)) + dgamma(exp(kappa.p), a1, b1, log = T) + kappa.p -
    sum(dgamma(theta, exp(kappa), beta, log=T)) - dgamma(exp(kappa), a1, b1, log = T) - kappa
  # 3. actualizar valor
  if( log(runif(1)) < lr ) { 
    alpha <- exp(kappa.p)
    ac_a <- ac_a+1 
  }
  ## simular lambda (metropolis)
  # 1. propuesta
  nu <- lambda/10
  delta <- log(nu) - log(1-nu)
  delta.p <- rnorm(1, delta, sqrt(delta_l))
  # 2. tasa de aceptacion
  lr_l <- sum(dnbinom(x = y, mu = theta*x, size = (theta*x)/(1/(1 + exp(-delta.p))) )) + dunif((1/(1 + exp(-delta.p))), a3, b3, log = T) + (-delta.p - 2*log(1 + exp(-delta.p))) -
    sum(dnbinom(x = y, mu = theta*x, size = (theta*x)/(1/(1 + exp(-delta))) )) - dunif((1/(1 + exp(-delta))), a3, b3, log = T)  - (-delta - 2*log(1 + exp(-delta))) 
  # 3. actualizar valor
  if( log(runif(1)) < lr_l ) { 
    nu <- 1/(1 + exp(-delta.p))
    lambda <- nu*10
    ac_l <- ac_l+1 
  }
  ## Almacenar
  THETA[s,] <- theta
  ABL[s,] <- c(alpha, beta, lambda)
  ## Log-Verosimilitud
  LV2[s] <- sum(dnbinom(x = y, mu = theta*x, size = (theta*x)/lambda, log = T))
  ## Progreso
  if (s%%ncat == 0) cat(100*round(s/S, 1), "% completed ... \n", sep = "" )
  
}
####### fin MCMC

# tasa de aceptacion
100*ac_l/S
100*ac_a/S
100*ac_t/S

THETA <- THETA[-c(1:1000),]
ABL <- ABL[-c(1:1000),]
LV2 <- LV2[-c(1:1000)]

# diagnosticos
library(coda)
apply(X = THETA, MARGIN = 2, FUN = effectiveSize)
apply(X = ABL, MARGIN = 2, FUN = effectiveSize)
effectiveSize(LV2)

plot(LV2,xlab="iteration",type ="l",ylab="deviance", main = "")

######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,1),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(ABL[,1],xlab="iteration",type ="l",ylab=expression(alpha), main = "")

windows(height=1.75,width=5)
par(mfrow=c(1,1),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(ABL[,2],xlab="iteration",type ="l",ylab=expression(beta), main = "")

windows(height=1.75,width=5)
par(mfrow=c(1,1),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(ABL[,3],xlab="iteration",type ="l",ylab=expression(lambda), main = "")

######## cadenas
windows()
par(mfrow=c(2,3))
plot(THETA[,1],xlab="iteration",type ="l",ylab=expression(theta[1]), main = "")
plot(THETA[,2],xlab="iteration",type ="l",ylab=expression(theta[2]), main = "")
plot(THETA[,3],xlab="iteration",type ="l",ylab=expression(theta[3]), main = "")
plot(THETA[,4],xlab="iteration",type ="l",ylab=expression(theta[4]), main = "")
plot(THETA[,5],xlab="iteration",type ="l",ylab=expression(theta[5]), main = "")
plot(THETA[,6],xlab="iteration",type ="l",ylab=expression(theta[6]), main = "")
####

#### Adelgazamiento 
blabs <- c(expression(theta[1]),expression(theta[2]),expression(theta[3]),
           expression(theta[4]),expression(theta[5]),expression(theta[6]))  # etiquetas
thin  <- c(1,(1:20000)*(S/25100))  # muestreo sistematico 

#theta1
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,THETA[thin,1],type="l",xlab="Iteración",ylab=blabs[1])
abline(h=mean(THETA[,1]))
acf(THETA[,1],ci.col="gray",xlab="lag")
acf(THETA[thin,1],xlab="lag/10",ci.col="gray")

#theta2
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,THETA[thin,2],type="l",xlab="Iteración",ylab=blabs[2])
abline(h=mean(THETA[,2]))
acf(THETA[,2],ci.col="gray",xlab="lag")
acf(THETA[thin,2],xlab="lag/10",ci.col="gray")

#theta3
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,THETA[thin,3],type="l",xlab="Iteración",ylab=blabs[3])
abline(h=mean(THETA[,3]))
acf(THETA[,3],ci.col="gray",xlab="lag")
acf(THETA[thin,3],xlab="lag/10",ci.col="gray")

#theta4
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,THETA[thin,4],type="l",xlab="Iteración",ylab=blabs[4])
abline(h=mean(THETA[,4]))
acf(THETA[,4],ci.col="gray",xlab="lag")
acf(THETA[thin,4],xlab="lag/10",ci.col="gray")

#theta5
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,THETA[thin,5],type="l",xlab="Iteración",ylab=blabs[5])
abline(h=mean(THETA[,5]))
acf(THETA[,5],ci.col="gray",xlab="lag")
acf(THETA[thin,5],xlab="lag/10",ci.col="gray")

#theta6
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,THETA[thin,6],type="l",xlab="Iteración",ylab=blabs[6])
abline(h=mean(THETA[,6]))
acf(THETA[,6],ci.col="gray",xlab="lag")
acf(THETA[thin,6],xlab="lag/10",ci.col="gray")

THETA <- THETA[thin,]
apply(X = THETA, MARGIN = 2, FUN = effectiveSize)

#alpha
blabs <- c(expression(alpha))  # etiquetas
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,ABL[thin,1],type="l",xlab="Iteración",ylab=blabs[1])
abline(h=mean(ABL[,1]))
acf(ABL[,1],ci.col="gray",xlab="lag")
acf(ABL[thin,1],xlab="lag/10",ci.col="gray")

#beta
blabs <- c(expression(beta))  # etiquetas
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,ABL[thin,2],type="l",xlab="Iteración",ylab=blabs[1])
abline(h=mean(ABL[,2]))
acf(ABL[,2],ci.col="gray",xlab="lag")
acf(ABL[thin,2],xlab="lag/10",ci.col="gray")


#lambda
blabs <- c(expression(lambda))  # etiquetas
windows(height=1.75,width=5)
par(mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
par(mfrow=c(1,3))
plot(thin,ABL[thin,3],type="l",xlab="Iteración",ylab=blabs[1])
abline(h=mean(ABL[,3]))
acf(ABL[,3],ci.col="gray",xlab="lag")
acf(ABL[thin,3],xlab="lag/10",ci.col="gray")

ABL <- ABL[thin,]
apply(X = ABL, MARGIN = 2, FUN = effectiveSize)

####

#PARA THETAS
######## errores estandar de MC : DE(parametro)/sqrt( n_eff )
MCERR <- apply(THETA,2,sd) / sqrt( effectiveSize(THETA) )
MCERR/apply(THETA,2,mean) #dan muy pequeños entonces no hay problema

#PARA ALPHA Y BETA
MCERRb <- apply(ABL, 2, sd) / sqrt( effectiveSize(ABL) )
MCERRb/apply(ABL, 2, mean) #dan muy pequeños entonces no hay problema

#PARA DEVIANCE
esd <- sd(LV2)/sqrt(effectiveSize(LV2))
esd/abs(mean(LV2))

#### RESUMEN POSTERIOR ####
resumen_posterior <- function(x) 
{
  round(c(mean(x), sd(x), quantile(x = x, probs = c(0.025,0.975))), 3)
}

tab <- rbind(resumen_posterior(ABL[,1]),
             resumen_posterior(ABL[,2]),
             resumen_posterior(ABL[,3]),
             resumen_posterior(THETA[,1]),
             resumen_posterior(THETA[,2]),
             resumen_posterior(THETA[,3]),
             resumen_posterior(THETA[,4]),
             resumen_posterior(THETA[,5]),
             resumen_posterior(THETA[,6]))
colnames(tab) <- c("Media", "SD", "Q2.5%", "Q97.5%")
rownames(tab) <- c("alpha", "beta", "lambda",paste0("theta", 1:n))
print(tab)

#### COMPARACION THETA Y OLS ####
tab2 <- rbind(c(mean(THETA[,1]), y[1]/x[1]),
              c(mean(THETA[,2]), y[2]/x[2]),
              c(mean(THETA[,3]), y[3]/x[3]),
              c(mean(THETA[,4]), y[4]/x[4]),
              c(mean(THETA[,5]), y[5]/x[5]),
              c(mean(THETA[,6]), y[6]/x[6]))
colnames(tab2) <- c("Posterior", "OLS")
rownames(tab2) <- c(paste0("theta", 1:n))

print(tab2)

windows()
par(mfrow=c(2,3))
# grafico theta1
#windows()
hist(x =THETA[,1] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "", xlim = c(0,2), ylim = c(0,8))
lines(density(THETA[,1]))
abline(v = y[1]/x[1], col = "red", lwd = 2, lty = 3)
abline(v=c(quantile(THETA[,1], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(THETA[,1]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

# grafico theta2
#windows()
hist(x =THETA[,2] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "", xlim = c(0,4), ylim=c(0,3.5))
lines(density(THETA[,2]))
abline(v = y[2]/x[2], col = "red", lwd = 2, lty = 3)
abline(v=c(quantile(THETA[,2], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(THETA[,2]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

# grafico theta3
#windows()
hist(x =THETA[,3] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "", ylim = c(0,5))
lines(density(THETA[,3]))
abline(v = y[3]/x[3], col = "red", lwd = 2, lty = 3)
abline(v=c(quantile(THETA[,3], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(THETA[,3]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

# grafico theta4
#windows()
hist(x =THETA[,4] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(THETA[,4]))
abline(v = y[4]/x[4], col = "red", lwd = 2, lty = 3)
abline(v=c(quantile(THETA[,4], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(THETA[,4]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

# grafico theta5
#windows()
hist(x =THETA[,5] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "", xlim = c(0,3))
lines(density(THETA[,5]))
abline(v = y[5]/x[5], col = "red", lwd = 2, lty = 3)
abline(v=c(quantile(THETA[,5], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(THETA[,5]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

# grafico theta6
#windows()
hist(x =THETA[,6] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(THETA[,6]))
abline(v = y[6]/x[6], col = "red", lwd = 2, lty = 3)
abline(v=c(quantile(THETA[,6], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(THETA[,6]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

#### POSTERIOR ABL #####
hist(x =ABL[,1] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")#, xlim = c(0,2), ylim = c(0,8))
lines(density(ABL[,1]))
abline(v=c(quantile(ABL[,1], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(ABL[,1]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

hist(x =ABL[,2] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")#, xlim = c(0,2), ylim = c(0,8))
lines(density(ABL[,2]))
abline(v=c(quantile(ABL[,2], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(ABL[,2]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

hist(x =ABL[,3] , freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")#, xlim = c(0,2), ylim = c(0,8))
lines(density(ABL[,3]))
abline(v=c(quantile(ABL[,3], c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(ABL[,3]), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

#### COMPARACION ALPHA Y BETA ####
mean(ABL[,1])/mean(ABL[,2])
mean(y/x)

propAB <- ABL[,1]/ABL[,2]

hist(x = propAB, freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "", xlim = c(0,2), ylim = c(0,6))
lines(density(propAB))
abline(v = mean(y/x), col = "red", lwd = 2, lty = 3)
abline(v=c(quantile(propAB, c(0.025, 0.975))), col="dodgerblue4", pch=3)
abline(v = mean(propAB), col = "mediumorchid4", lwd = 2, lty = 1)
legend("topright", legend = c("Media posterior", "IC 95%", "Observación media y/x"), 
       bty = "n", lwd = 2, lty = c(1,1,3) ,col = c("mediumorchid4", "dodgerblue4", "red"))

#### COMPARACION ENTRE THETAS ####
windows()
par(mfrow=c(2,3))
# dispersión entre thetas
plot(x =THETA[,2], y = THETA[,1], 
     xlab = expression(theta[2]), ylab = expression(theta[1]), main = "", col = "grey", type = "p",
     xlim = c(0,1.4), ylim = c(0,1.2))
points(x=y[2]/x[2], y = y[1]/x[1], pch=15)
abline(0,1, col="red")
legend("topright", legend = c("Observación de y1/x1 vs y2/x2"), 
       fill=1, border = "black")


plot(x =THETA[,2], y = THETA[,3], 
     xlab = expression(theta[2]), ylab = expression(theta[3]), main = "", col = "grey", type = "p",
     xlim = c(0,1.4), ylim = c(0,0.8))
points(x=y[2]/x[2], y = y[3]/x[3], pch=15)
abline(0,1, col="red")
legend("topright", legend = c("Observación de y3/x3 vs y2/x2"), 
       fill=1, border = "black")

plot(x =THETA[,2], y = THETA[,4], 
     xlab = expression(theta[2]), ylab = expression(theta[4]), main = "", col = "grey", type = "p",
     xlim = c(0,1.4), ylim = c(0,1.2))
points(x=y[2]/x[2], y = y[4]/x[4], pch=15)
abline(0,1, col="red")
legend("topright", legend = c("Observación de y4/x4 vs y2/x2"), 
       fill=1, border = "black")

plot(x =THETA[,2], y = THETA[,5], 
     xlab = expression(theta[2]), ylab = expression(theta[5]), main = "", col = "grey", type = "p",
     xlim = c(0,1.5), ylim = c(0,1))
points(x=y[2]/x[2], y = y[5]/x[5], pch=15)
abline(0,1, col="red")
legend("topright", legend = c("Observación de y5/x5 vs y2/x2"), 
       fill=1, border = "black")

plot(x =THETA[,2], y = THETA[,6], 
     xlab = expression(theta[2]), ylab = expression(theta[6]), main = "", col = "grey", type = "p",
     xlim = c(0,1.5), ylim = c(0,1))
points(x=y[2]/x[2], y = y[6]/x[6], pch=15)
abline(0,1, col="red")
legend("topright", legend = c("Observación de y6/x6 vs y2/x2"), 
       fill=1, border = "black")

## Probabilidades
#Pr(theta2 > thetaj)

mean(THETA[,2] > THETA[,1])
mean(THETA[,2] > THETA[,3])
mean(THETA[,2] > THETA[,4])
mean(THETA[,2] > THETA[,5])
mean(THETA[,2] > THETA[,6])

#Pr(theta2 = max(theta))

mean(THETA[,2] == apply(THETA, 1, max))

y/x

#### DIC ####
theta_hat <- colMeans(THETA)
lambda_hat <- mean(ABL[,3])
lpyth_m2 <- sum(dnbinom(x = y, mu = theta_hat*x, size = (theta_hat*x)/lambda_hat, log = T))
pDIC_m2  <- 2*(lpyth_m2 - mean(LV2))
dic_m2   <- -2*lpyth_m2 + 2*pDIC_m2 

#### BONDAD DE AJUSTE ####
S <- length(THETA[,1])
ybar_hat <- rep(NA, S)
sd_hat <- rep(NA, S)
for (i in 1:S) {
  y_pred <- rnbinom(n, mu = THETA[i,]*x, size = (THETA[i,]*x)/ABL[i,3])
  ybar_hat[i] <- mean(y_pred)
  sd_hat[i] <- sd(y_pred)
}

windows()
par(mfrow=c(1,2))
hist(x = ybar_hat, freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(ybar_hat))
abline(v = mean(y), col = "red", lwd = 2, lty = 3)
legend("topright", legend = expression(bar(y)), 
       bty = "n", lwd = 2, lty = 3 ,col = "red")


hist(x = sd_hat, freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(sd_hat))
abline(v = sd(y), col = "red", lwd = 2, lty = 3)
legend("topright", legend = "sd(y)", 
       bty = "n", lwd = 2, lty = 3,col = "red")

# ppp
mean( ybar_hat > mean(y))
mean( sd_hat > sd(y))
