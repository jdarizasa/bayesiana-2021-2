############################
##### PREPROCESAMIENTO ####
###########################
#1 Cargar datos y escoger variables (esto ya no)
SB2020_I <- data.frame(dept = SB11_20201$COLE_COD_DEPTO_UBICACION, sexo = SB11_20201$ESTU_GENERO,
                       htra = SB11_20201$ESTU_HORASSEMANATRABAJA, mat = SB11_20201$PUNT_MATEMATICAS)

SB2020_II <- data.frame(dept = SB11_20202$COLE_COD_DEPTO_UBICACION, sexo = SB11_20202$ESTU_GENERO,
                        htra = SB11_20202$ESTU_HORASSEMANATRABAJA, mat = SB11_20202$PUNT_MATEMATICAS)

rm(SB11_20201)
rm(SB11_20202)

write.csv(SB2020_I, file = "saber2020I.csv")
write.csv(SB2020_II, file = "saber2020II.csv")

#2 Eliminar datos faltantes
sb_1 <- read.csv("saber2020I.csv", header = T, sep = ',', na.strings = '')
sb_2 <- read.csv("saber2020II.csv", header = T, sep = ',', na.strings = '')

sapply(sb_1, function(x) sum(is.na(x)))
sapply(sb_2, function(x) sum(is.na(x)))

sb_1 <- na.omit(sb_1); sb_1 <- sb_1[,-1]
sb_2 <- na.omit(sb_2); sb_2 <- sb_2[,-1]

#3 Codificar sexo
class(sb_1$sexo)
levels(sb_1$sexo)
sexo.t <- as.numeric(sb_1$sexo) - 1
sb_1$sexo <- sexo.t

class(sb_2$sexo)
levels(sb_2$sexo)
sexo.t <- as.numeric(sb_2$sexo) - 1
sb_2$sexo <- sexo.t

#4 Codificar trabajo
class(sb_1$htra)
levels(sb_1$htra)
tra.t <- as.numeric(sb_1$htra)
tra.t[ tra.t > 1] <- 0
sb_1$htra <- tra.t

class(sb_2$htra)
levels(sb_2$htra)
tra.t <- as.numeric(sb_2$htra)
tra.t[ tra.t > 1] <- 0
sb_2$htra <- tra.t

write.csv(sb_1, file = "SB_I.csv")
write.csv(sb_2, file = "SB_II.csv")

################################
# PUNTO 1 ####
##############################
library(rgdal)
library(sp)
sb_1 <- read.csv("SB_I.csv", header = T, sep = ",")
sb_2 <- read.csv("SB_II.csv", header = T, sep = ",")
cod <- sort(unique(sb_2$dept))

puntmed <-NULL
for (i in 1:length(cod)) {
  puntmed[i] <- mean(sb_1$mat[sb_1$dept==cod[i]])
}

mapcol <- readOGR("depto.shp")
plot(mapcol)
box()
mapcol$AREA <- puntmed
mapcol@data <- mapcol@data[,c(-4,-5)]

class(mapcol)
View(mapcol@data)
spplot(mapcol, "AREA")

puntmed2 <-NULL
for (i in 1:length(cod)) {
  puntmed2[i] <- mean(sb_2$mat[sb_2$dept==cod[i]])
}

mapcol2 <- readOGR(file.choose())
mapcol2$AREA <- puntmed2
mapcol2@data <- mapcol2@data[,c(-4,-5)]

class(mapcol2)
spplot(mapcol2, "AREA")

load(file = "data_parcial2.RData")
##########
#sb_1 : base de datos procesada 2020-1
#sb_2 : base de datos procesada 2020-2
#X : matriz de diseño datos entrenamiento
#y : vector de observaciones datos entrenamiento
#Xk : lista de matrices de diseño por departamento
#YK : lista de vectores de observación por departamento
#cod : codigo de los departamentos
#m : cantidad de departamentos
#n: cantidad de observaciones
#nj: cantidad de observaciones por departamento
#p :cantidad de covariables
##########

#############################
### MODELO 1
#############################
# y    : puntaje de matematicas en la prueba saber 11.
# sexo : sexo del estudiante. 
# htra : trabaja (1 si trabajó 0 horas, 0 si trabajó).
#
# Modelo: E(y | X) = beta1*x1 + beta2*x2 + beta3*x3
#         x1 = 1 
#         x2 = 1 si masculino, 0 si femenino
#         x3 = 1 si trabajó 0 horas, 0 si trabajó

#OLS
beta.ols <- solve(t(X)%*%X)%*%t(X)%*%y
sigma2.ols <- sum((y - X%*%beta.ols)^2)/(n-p)

#hiperparametros previa unitaria
beta_0 <- beta.ols
S_0 <- n*sigma2.ols*solve(t(X)%*%X)
nu_0 <- 1
s2_0 <- sigma2.ols

#calculos y funciones para la posterior
S_0i <- solve(S_0)
H <- t(X)%*%X
SiB_0 <- S_0i%*%beta_0
Xty <- t(X)%*%y

a <- (nu_0 + n)/2
SSRb <- function(y, X, b){
  k <- sum((y - X%*%b)^2)
  return(k)
}

#puntos iniciales
set.seed(7)
beta <- c(mvtnorm::rmvnorm(1, beta_0, S_0))
s2 <- 1/rgamma(1, a, (nu_0*s2_0 + SSRb(y=y, X=X, b=beta_0)))

# setup
S    <- 55000                            # numero de muestras 
B   <- matrix(data=NA, nrow=S, ncol=3)  # matriz para almacenar las muestras
S2  <- matrix(data=NA, nrow=S, ncol=1)
LP1 <- matrix(data=NA, nrow = S, ncol = 1)
ncat <- floor(S/10)                      # mostrar anuncios cada 10%

#MCMC
set.seed(7)
for (i in 1:S) {
  
  #actualizar valor de betas
  S2n <- solve(S_0i + (1/s2)*H)
  mun <- S2n%*%(SiB_0 + (1/s2)*Xty)
  beta <- c(mvtnorm::rmvnorm(1, mun, S2n))
  
  #actualizar el valor de s2
  b <- (nu_0*s2_0 + SSRb(y=y, X=X, b= beta))/2
  s2 <- 1/rgamma(1, a, b)
  
  ####
  #almacenar
  B[i,] <- beta 
  S2[i,] <- s2
  
  # log-verosimilitud
  LP1[i] <- sum(dnorm(x = y, mean = X%*%beta, sd = sqrt(s2), log = T))
  
  # progreso  
  if (i%%ncat == 0) cat("MCMC M1, ", 100*round(i/S, 1), "% completado \n", sep = "")
}

###### remover las primeras 5000 iteraciones de las cadenas
B <- B[-c(1:5000),]
S2  <- S2 [-c(1:5000) ]
LP1  <- LP1 [-c(1:5000) ]

# salvar muestras
save(B, S2, LP1, file = "samples_parcial2_punto2.RData")

# diagnostico
S <- 50000
windows(width = 6, height = 4)
par(mar=c(2.75,2.75,.5,.5), mgp=c(1.7,.7,0))
plot(x = 1:S, y = LP1, type = "l", cex.axis = 0.8, 
     main = "Convergencia Log-Verosimilitud", xlab = "iteración", ylab = "log-verosimilitud")


######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,1),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(S2,xlab="iteration",main ="Convergencia Sigma^2",type ="l",ylab=expression(sigma^2))


######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,3))
plot(B[,1],xlab="iteration",type ="l",ylab=expression(beta[1]), main = "Convergencia B1")
plot(B[,2],xlab="iteration",type ="l",ylab=expression(beta[2]), main = "Convergencia B2")
plot(B[,3],xlab="iteration",type ="l",ylab=expression(beta[3]), main = "Convergencia B3")

######## autocorrelaciones
library(coda)
windows()
acf(S2,main=expression(sigma^2)) -> a2


######## tama?os efectivos de muestra
effectiveSize(S2)
summary( effectiveSize(B) )

# resumen posterior
resumen_posterior <- function(x) 
{
  round(c(mean(x), sd(x), quantile(x = x, probs = c(0.025,0.975))), 3)
}

tab <- rbind(resumen_posterior(B[,1]),
             resumen_posterior(B[,2]),
             resumen_posterior(B[,3]),
             resumen_posterior(S2))
colnames(tab) <- c("Media", "SD", "Q2.5%", "Q97.5%")
rownames(tab) <- c(paste0("beta", 1:p), "sig2")
print(tab)

####### DIC
set.seed(7)
B1_hat   <- mean(B[,1])
B2_hat   <- mean(B[,2])
B3_hat   <- mean(B[,3])
beta_hat <- c(B1_hat,B2_hat,B3_hat)
sig2_hat <- mean(S2)
lpyth_m1 <- sum(dnorm(x = y, mean = X%*%beta_hat, sd = sqrt(sig2_hat), log = T))
pDIC_m1  <- 2*(lpyth_m1 - mean(LP1))
dic_m1   <- -2*lpyth_m1 + 2*pDIC_m1 

####### MSE
sb_2 <- read.csv("SB_II.csv", header = T, sep = ",")
sb_2 <- sb_2[order(sb_2$dept),]
cod1 <- sort(unique(sb_1$dept))
cod2 <- sort(unique(sb_2$dept))

sb_2 <- sb_2[!(sb_2$dept==27 | sb_2$dept==70 | sb_2$dept==88 | sb_2$dept==91 | sb_2$dept==94 | sb_2$dept==95 | sb_2$dept==97 | sb_2$dept==99),]
y2 <- as.matrix(sb_2$mat)
X2 <- cbind(1, sb_2$sexo, sb_2$htra)

y_hat <- X2%*%beta_hat
mean((y2 - y_hat)^2)

##############################
### MODELO 2
##############################
# y    : puntaje de matematicas en la prueba saber 11.
# sexo : sexo del estudiante. 
# htra : trabaja (1 si trabajó 0 horas, 0 si trabajó).
#
# Modelo: E(y | X) = beta1*x1 + beta2*x2 + beta3*x3 + theta_j
#         x1 = 1 
#         x2 = 1 si masculino, 0 si femenino
#         x3 = 1 si trabajó 0 horas, 0 si trabajó
#         theta_j = efecto aleatorio

#OLS
beta.ols <- solve(t(X)%*%X)%*%t(X)%*%y
sigma2.ols <- sum((y - X%*%beta.ols)^2)/(n-p)

#hiperparametros previa unitaria
beta_0 <- beta.ols
S_0 <- n*sigma2.ols*solve(t(X)%*%X)
eta_0 <- 1
t2_0 <- sigma2.ols
nu_0 <- 1
s2_0 <- sigma2.ols

#calculos y funciones para la posterior
EPj <- function(y, X, b){
  k <- sum(y - X%*%b)
  return(k)
} 

S_0i <- solve(S_0)
H <- t(X)%*%X
SiB_0 <- S_0i%*%beta_0

a1 <- (eta_0 + m)/2
a2 <- (nu_0 + n)/2

SSRb2 <- function(y, X, v, b){
  k <- sum((y - (X%*%b + v))^2)
  return(k)
}

#puntos iniciales
set.seed(7)
theta <- c(mvtnorm::rmvnorm(1, c(rep(0,25)), diag(x=1/(1/t2_0 + 1/s2_0), m)))
t2 <- 1/rgamma(1, a1, (eta_0*t2_0 + t(theta)%*%theta))
beta <- c(mvtnorm::rmvnorm(1, beta_0, S_0))
s2 <- 1/rgamma(1, a2, (nu_0*s2_0 + SSRb2(y=y, X=X, b=beta_0, v=c(rep(0,n)))))

# setup
S    <- 55000                            # numero de muestras 
B   <- matrix(data=NA, nrow=S, ncol=3)  # matriz para almacenar las muestras
S2  <- matrix(data=NA, nrow=S, ncol=1)
TH  <- matrix(data=NA, nrow=S, ncol=m)
LP2 <- matrix(data=NA, nrow=S, ncol=1)
ncat <- floor(S/10)                      # mostrar anuncios cada 10%

#MCMC
set.seed(7)
for (i in 1:S) {
  #actualizar valor de theta
  for (j in 1:25) {
    y_k <- c(sb_1$mat[sb_1$dept==cod[j]])
    X_k <- cbind(1, sb_1$sexo[sb_1$dept==cod[j]], sb_1$htra[sb_1$dept==cod[j]])
    sn <- 1/(1/t2 + nj[j]/s2)
    mn <- sn*(1/s2)*EPj(y=y_k, X=X_k, b=beta)
    theta[j] <- rnorm(1, mn, sqrt(sn))
  }
  
  #actualizar valor de betas  
  v <- rep(theta[1],nj[1])
  for (j in 2:m) {
    v <- c(v,rep(theta[j],nj[j]))
  }
  S2n <- solve(S_0i + (1/s2)*H)
  mun <- S2n%*%(SiB_0 + (1/s2)*(t(X)%*%(y-v)))
  beta <- c(mvtnorm::rmvnorm(1, mun, S2n))
  
  #actualizar el valor de s2
  b2 <- (nu_0*s2_0 + SSRb2(y=y, X=X, v=v, b= beta))/2
  s2 <- 1/rgamma(1, a2, b2)
  
  
  #actualizar valor de tau2
  b1 <- (eta_0*t2_0 + t(theta)%*%theta)/2
  t2 <- 1/rgamma(1, a1, b1)
  
  ##
  #almacenar
  B[i,] <- beta 
  S2[i,] <- s2
  TH[i,] <- theta
  
  #Log-verosimilitud 
  LP2[i] <- sum(dnorm(x = y, mean = X%*%beta + v, sd = sqrt(s2), log = T))
  
  # progreso
  if (i%%ncat == 0) cat("MCMC M2, ", 100*round(i/S, 1), "% completado \n", sep = "")
}

# salvar muestras
save(B, S2, TH, LP2, file = "samples_parcial2_punto4.RData")


###### remover las primeras 5000 iteraciones de las cadenas
B <- B[-c(1:5000),]
S2  <- S2 [-c(1:5000) ]
TH <- TH[-c(1:5000),]
LP2  <- LP2 [-c(1:5000) ]


# diagnostico
windows(width = 6, height = 4)
S <- 50000
par(mar=c(2.75,2.75,.5,.5), mgp=c(1.7,.7,0))
plot(x = 1:S, y = LP2, type = "l", cex.axis = 0.8, 
     main = "Convergencia Log-Verosimilitud", xlab = "iteración", ylab = "log-verosimilitud")

######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,1),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(S2,xlab="iteration",type ="l",ylab=expression(sigma^2), main = "Convergencia Sigma^2")


######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,3))
plot(B[,1],xlab="iteration",type ="l",ylab=expression(beta[1]), main = "Convergencia B1")
plot(B[,2],xlab="iteration",type ="l",ylab=expression(beta[2]), main = "Convergencia B2")
plot(B[,3],xlab="iteration",type ="l",ylab=expression(beta[3]), main = "Convergencia B3")

######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(TH[,12],xlab="iteration",type ="l",ylab=expression(beta[1]))
plot(TH[,25],xlab="iteration",type ="l",ylab=expression(beta[2]))
plot(TH[,8],xlab="iteration",type ="l",ylab=expression(beta[3]))

######## autocorrelaciones
library(coda)
windows()
acf(S2,main=expression(sigma^2)) -> a2


######## tama?os efectivos de muestra
effectiveSize(S2)
summary( effectiveSize(B) )
summary( effectiveSize(TH) )

# resumen posterior
resumen_posterior <- function(x) 
{
  round(c(mean(x), sd(x), quantile(x = x, probs = c(0.025,0.975))), 3)
}

tab <- rbind(resumen_posterior(B[,1]),
             resumen_posterior(B[,2]),
             resumen_posterior(B[,3]),
             resumen_posterior(S2))
colnames(tab) <- c("Media", "SD", "Q2.5%", "Q97.5%")
rownames(tab) <- c(paste0("beta", 1:p), "sig2")
print(tab)


####### DIC
set.seed(7)
B1_hat   <- mean(B[,1])
B2_hat   <- mean(B[,2])
B3_hat   <- mean(B[,3])
beta_hat <- c(B1_hat,B2_hat,B3_hat)
sig2_hat <- mean(S2)
v_hat <- rep(mean(TH[,1]),nj[1])
for (j in 2:m) {
  v_hat <- c(v_hat,rep(mean(TH[,j]),nj[j]))
}
lpyth_m2 <- sum(dnorm(x = y, mean = X%*%beta_hat + v_hat, sd = sqrt(sig2_hat), log = T))
pDIC_m2  <- 2*(lpyth_m2 - mean(LP2))
dic_m2   <- -2*lpyth_m2 + 2*pDIC_m2 

####### MSE
sb_2 <- read.csv("SB_II.csv", header = T, sep = ",")
sb_2 <- sb_2[order(sb_2$dept),]
sb_2 <- sb_2[!(sb_2$dept==27 | sb_2$dept==70 | sb_2$dept==88 | sb_2$dept==91 | sb_2$dept==94 | sb_2$dept==95 | sb_2$dept==97 | sb_2$dept==99),]

y2 <- as.matrix(sb_2$mat)
X2 <- cbind(1, sb_2$sexo, sb_2$htra)

contdept <- data.frame(table(sb_2$dept))
nj2 <- contdept$Freq

v_hat <- rep(mean(TH[,1]),nj2[1])
for (j in 2:m) {
  v_hat <- c(v_hat,rep(mean(TH[,j]),nj2[j]))
}
y_hat <- X2%*%beta_hat + v_hat
mean((y2 - y_hat)^2)


###############################
#### MODELO 3
###############################
# y    : puntaje de matematicas en la prueba saber 11.
# sexo : sexo del estudiante. 
# htra : trabaja (1 si trabajó 0 horas, 0 si trabajó).
#
# Modelo: E(y_j | X) = beta1_j*x1 + beta2_j*x2 + beta3_j*x3 + theta_j
#         x1 = 1 
#         x2 = 1 si masculino, 0 si femenino
#         x3 = 1 si trabajó 0 horas, 0 si trabajó
#         theta_j = efecto aleatorio

###########OLS
beta.ols <- solve(t(X)%*%X)%*%t(X)%*%y
sigma2.ols <- sum((y - X%*%beta.ols)^2)/(n-p)

######CONDICIONALES COMPLETAS
sample_theta <- function(nj, s2_j, t2, Yk, Xk, B_j, theta)
{
  for (j in 1:25) {
    sn <- 1/(1/t2 + nj[j]/s2_j[j])
    mn <- sn*(1/s2_j[j])*sum(Yk[[j]] - Xk[[j]]%*%B_j[j,])
    theta[j] <- rnorm(1, mn, sqrt(sn))
  }
  return(theta)
}
#set.seed(7)
#sample_theta(nj=nj, theta = c(0,25), B_j = B_j, s2_j=s2_j, t2=sigma2.ols, Xk=Xk, Yk=Yk)

sample_tau2 <- function(eta_0, m, t2_0, theta, t2)
{
  a1 <- (eta_0 + m)/2    
  b1 <- (eta_0*t2_0 + t(theta)%*%theta)/2
  t2 <- 1/rgamma(1, a1, b1)
  return(t2)
}
#sample_tau2(eta_0 = 1, m = 25, t2_0 = 1, theta = theta, t2=0)

sample_Beta_j <- function(Xk, Yk, theta, nj, Sigma, s2_j, beta, B_j)
{
  for (j in 1:25) {
    sig <- solve(solve(Sigma) + (1/s2_j[j])*(t(Xk[[j]])%*%Xk[[j]]))
    mun <- sig%*%(solve(Sigma)%*%beta + (1/s2_j[j])*(t(Xk[[j]])%*%(Yk[[j]] - theta[j])))
    B_j[j,] <- c(mvtnorm::rmvnorm(1, mun, sig))
  }
  return(B_j)
}
#sample_Beta_j(Xk=Xk, Yk=Yk, theta=theta, nj=nj, Sigma=diag(1,3), s2_j=c(rep(1,25)), beta = beta.ols, B_j=matrix(0, 25, 3))

sample_beta <- function(B_j, D_0, m, Sigma, mu_0, beta)
{
  suB_j <- c(apply(B_j, 2, sum))
  S2n <- solve(solve(D_0) + m*solve(Sigma))
  mun2 <- S2n%*%(solve(D_0)%*%mu_0 + solve(Sigma)%*%suB_j)
  beta <- c(mvtnorm::rmvnorm(1, mun2, S2n))
  return(beta)
}
#sample_beta(B_j = B_j, D_0 = diag(1,3), m = m, Sigma = diag(1,3), mu_0 = beta.ols, beta=c(0,3))

library(MCMCpack)
sample_Sigma <- function(B_j, beta, n_0, m, S_0, Sigma)
{
  Baux <- matrix(0, 3, 3)
  for (j in 1:25) {
    Baux <- ((B_j[j,] - beta)%*%t(B_j[j,] - beta)) + Baux 
  }
  aS <- n_0 + m
  bS <- solve(S_0 + Baux) 
  Sigma <- riwish(aS, solve(bS))
  return(Sigma)
}
#sample_Sigma(B_j = B_j, beta = beta.ols, n_0=5, m = m, S_0 = diag(1,3), Sigma=matrix(0,3))

sample_S2_j <- function(Xk, Yk, theta, nj, nu, s2, B_j, s2_j)
{
  for (j in 1:25) {
    a2 <- (nu + nj[j])/2
    b2 <- (nu*s2 + sum((Yk[[j]] - (Xk[[j]]%*%B_j[j,]+ theta[j]))^2))/2
    s2_j[j] <- 1/rgamma(1, a2, b2)
  }
  return(s2_j)
}
#sample_S2_j(Xk=Xk, Yk=Yk, theta = theta, nj = nj, nu = 1, s2 = sigma2.ols, B_j = B_j, s2_j = c(rep(0,25)))

sample_nu <- function(nu0s, m, s2, s2_j, ka_0, nu)
{
  lpnu0 <- .5*nu0s*m*log(s2*nu0s/2)-m*lgamma(nu0s/2)+(nu0s/2)*sum(log(1/s2_j)) - .5*nu0s*s2*sum(1/s2_j) - ka_0*nu0s
  nu   <- sample(x = nu0s, size = 1, prob=exp( lpnu0-max(lpnu0) ))
  return(nu)
}
#sample_nu(nu0s = 1:100, m = m, s2 = sigma2.ols, s2_j = s2_j, ka_0 =1, nu=1)

sample_s2 <- function(a_0, b_0, m, nu, s2_j, s2)
{
  a3 <- a_0 + (m*nu)/2 
  b3 <- b_0 + (nu*sum(1/s2_j))/2
  s2 <- rgamma(1, shape=a3, rate=b3)
  return(s2)
}
#sample_s2(a_0 = 1, b_0 = sigma2.ols, m = m, nu = 3, s2_j = s2_j, s2 = 1)

######MCMC

#hiperparametros previa unitaria
eta_0 <- 1
t2_0 <- sigma2.ols

mu_0 <- beta.ols
D_0 <- n*sigma2.ols*solve(t(X)%*%X)
n_0 <- 5
S_0 <- n*sigma2.ols*solve(t(X)%*%X)

ka_0 <- 1
a_0 <- 1 
b_0 <- 1/sigma2.ols

#puntos iniciales
library(CholWishart)
set.seed(7)
theta <- c(rep(0,m))
t2 <- 1

B_j <- matrix(data=0, nrow = 25, ncol = 3)
beta <- c(0,0,0)
Sigma <- diag(1,3)

s2_j <- c(rep(1,25))
nu <- 1
s2 <- 1

# setup
S   <- 55000                            # numero de muestras 
B   <- matrix(data=NA, nrow=S, ncol=3)# matriz para almacenar las muestras
S2  <- matrix(data=NA, nrow=S, ncol=1)
B1 <- matrix(data = NA, nrow = S, ncol = m)
B2 <- matrix(data = NA, nrow = S, ncol = m)
B3 <- matrix(data = NA, nrow = S, ncol = m)
S2J <- matrix(data = NA, nrow = S, ncol = m)
Si1 <- matrix(data = NA, nrow = S, ncol = 3)
Si2 <- matrix(data = NA, nrow = S, ncol = 3)
Si3 <- matrix(data = NA, nrow = S, ncol = 3)
TH <- matrix(data=NA, nrow = S, ncol = m)
LP3 <- matrix(data = NA, nrow=S, ncol = 1)
NU <- matrix(data = NA, nrow = S, ncol = 1)
ncat <- floor(S/10)                      # mostrar anuncios cada 10%
nu0s <- 1:100

#cadena
set.seed(7)
for (i in 1:S) {
  #actualizar parametros
  B_j <- sample_Beta_j(Xk=Xk, Yk=Yk, theta = theta, nj = nj, Sigma = Sigma, s2_j = s2_j, beta = beta, B_j = B_j)
  
  s2_j <- sample_S2_j(Xk=Xk, Yk=Yk, theta = theta, nj = nj, nu = nu, s2 = s2, B_j = B_j, s2_j = s2_j)
  
  theta <- sample_theta(nj = nj, theta = theta, B_j = B_j, s2_j = s2_j, t2 = t2, Xk=Xk, Yk=Yk)
  
  t2 <- sample_tau2(eta_0 = eta_0, m = m, t2_0 = t2_0, theta = theta, t2 = t2)
  
  beta <- sample_beta(B_j = B_j, D_0 = D_0, m = m, Sigma = Sigma, mu_0 = mu_0, beta = beta)
  
  Sigma <- sample_Sigma(B_j = B_j, beta = beta, n_0 = n_0, m = m, S_0 = S_0, Sigma = Sigma)
  
  nu <- sample_nu(nu0s = nu0s, m = m, s2 = s2, s2_j = s2_j, ka_0 = ka_0, nu = nu)
  
  s2 <- sample_s2(a_0 = a_0, b_0 = b_0, m = m, nu = nu, s2_j = s2_j, s2 = s2)
  
  #almacenar
  B[i,] <- beta 
  S2[i,] <- s2
  B1[i,] <- B_j[,1]
  B2[i,] <- B_j[,2]
  B3[i,] <- B_j[,3]
  S2J[i,] <- s2_j
  Si1[i,] <- Sigma[,1]
  Si2[i,] <- Sigma[,2]
  Si3[i,] <- Sigma[,3]
  NU[i] <- nu
  TH[i,] <- theta
  
  #Log-verosimilitud 
  Xbt <- Xk[[1]]%*%B_j[1,] + theta[1]
  for (j in 2:m) {
    Xbt <- c(Xbt, Xk[[j]]%*%B_j[j,] + theta[j])
  }
  LP3[i] <- sum(dnorm(x = y, mean = Xbt , sd = rep(sqrt(s2_j),nj), log = T))
  
  # progreso
  if (i%%ncat == 0) cat("MCMC M3, ", 100*round(i/S, 1), "% completado \n", sep = "")
}

colMeans(Si1)
colMeans(Si2)
colMeans(Si3)

rbind(t(colMeans(Si1)), t(colMeans(Si2)), t(colMeans(Si3)))

# salvar muestras
save(B, S2, TH, B1, B2, B3, S2J, Si1, Si2, Si3, NU, LP3, file = "samples_parcial2_punto6.RData")

# cargar datos
load(file = "data_parcial2.RData")

# cargar muestras
load(file = "samples_parcial2_punto6.RData")

B <- B[-c(1:5000),]
S2  <- S2 [-c(1:5000) ]
B1 <- B1[-c(1:5000),]
B2 <- B2[-c(1:5000),]
B3 <- B3[-c(1:5000),]
TH <- TH[-c(1:5000),]
S2J <- S2J[-c(1:5000),]
LP3  <- LP3[-c(1:5000) ]

# diagnostico
windows(width = 6, height = 4)
S <- 50000
par(mar=c(2.75,2.75,.5,.5), mgp=c(1.7,.7,0))
plot(x = 1:S, y = LP3, type = "l", cex.axis = 0.8, 
     main = "Convergencia Log-Verosimilitud", xlab = "iteración", ylab = "log-verosimilitud")

mean(LP3)

######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,1),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(S2,xlab="iteration",type ="l",ylab=expression(sigma^2), main = "Convergencia Sigma^2")

######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,3))
plot(B[,1],xlab="iteration",type ="l",ylab=expression(beta[1]), main = "Convergencia B1")
plot(B[,2],xlab="iteration",type ="l",ylab=expression(beta[2]), main = "Convergencia B2")
plot(B[,3],xlab="iteration",type ="l",ylab=expression(beta[3]), main = "Convergencia B3")

######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(TH[,12],xlab="iteration",type ="l",ylab=expression(theta[1]))
plot(TH[,22],xlab="iteration",type ="l",ylab=expression(theta[2]))
plot(TH[,5],xlab="iteration",type ="l",ylab=expression(theta[3]))

######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(B1[,1],xlab="iteration",type ="l",ylab=expression(beta[1]))
plot(B1[,23],xlab="iteration",type ="l",ylab=expression(beta[2]))
plot(B1[,5],xlab="iteration",type ="l",ylab=expression(beta[3]))

######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(B2[,1],xlab="iteration",type ="l",ylab=expression(beta[1]))
plot(B2[,23],xlab="iteration",type ="l",ylab=expression(beta[2]))
plot(B2[,5],xlab="iteration",type ="l",ylab=expression(beta[3]))

######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,3),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(B3[,1],xlab="iteration",type ="l",ylab=expression(beta[1]))
plot(B3[,23],xlab="iteration",type ="l",ylab=expression(beta[2]))
plot(B3[,5],xlab="iteration",type ="l",ylab=expression(beta[3]))

######## cadenas
windows(height=1.75,width=5)
par(mfrow=c(1,1),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
plot(S2J[,25],xlab="iteration",type ="l",ylab=expression(sigma^2))

hist(NU)
mean(NU)

######## autocorrelaciones
library(coda)
windows()
acf(B[,1],main=expression(sigma^2)) -> a2


######## tama?os efectivos de muestra
effectiveSize(S2)
effectiveSize(NU)
apply(B, 2, effectiveSize) 
apply(B1, 2, effectiveSize)
apply(B2, 2, effectiveSize)
apply(B3, 2, effectiveSize)
apply(TH, 2, effectiveSize)
apply(S2J, 2, effectiveSize)

##### resumen posterior #######
resumen_posterior <- function(x) 
{
  round(c(mean(x), sd(x), quantile(x = x, probs = c(0.025,0.975))), 3)
}

tab <- rbind(resumen_posterior(B[,1]),
             resumen_posterior(B[,2]),
             resumen_posterior(B[,3]),
             resumen_posterior(S2))
colnames(tab) <- c("Media", "SD", "Q2.5%", "Q97.5%")
rownames(tab) <- c(paste0("beta", 1:p), "sig2")
print(tab)


####### DIC ##########
set.seed(7)
B1_hat   <- colMeans(B1)
B2_hat   <- colMeans(B2)
B3_hat   <- colMeans(B3)
B_hat <- matrix(data = c(B1_hat, B2_hat, B3_hat), nrow = 25, ncol = 3)
sig2_hat <- colMeans(S2J)
theta_hat <- colMeans(TH)
Xbt <- Xk[[1]]%*%B_hat[1,] + theta_hat[1]
for (j in 2:m) {
  Xbt <- c(Xbt, Xk[[j]]%*%B_hat[j,] + theta_hat[j])
}
lpyth_m3 <- sum(dnorm(x = y, mean = Xbt, sd = rep(sqrt(sig2_hat), nj), log = T))
pDIC_m3  <- 2*(lpyth_m3 - mean(LP3))
dic_m3   <- -2*lpyth_m3 + 2*pDIC_m3 

####### MSE ########
sb_2 <- read.csv("SB_II.csv", header = T, sep = ",")
sb_2 <- sb_2[order(sb_2$dept),]
sb_2 <- sb_2[!(sb_2$dept==27 | sb_2$dept==70 | sb_2$dept==88 | sb_2$dept==91 | sb_2$dept==94 | sb_2$dept==95 | sb_2$dept==97 | sb_2$dept==99),]

y2 <- as.matrix(sb_2$mat)

Yk2 <- list(c(sb_2$mat[sb_2$dept==cod[1]]), c(sb_2$mat[sb_2$dept==cod[2]]), c(sb_2$mat[sb_2$dept==cod[3]]), c(sb_2$mat[sb_2$dept==cod[4]]), c(sb_2$mat[sb_2$dept==cod[5]]), c(sb_2$mat[sb_2$dept==cod[6]]), c(sb_2$mat[sb_2$dept==cod[7]]), c(sb_2$mat[sb_2$dept==cod[8]]), c(sb_2$mat[sb_2$dept==cod[9]]), c(sb_2$mat[sb_2$dept==cod[10]])
            , c(sb_2$mat[sb_2$dept==cod[11]]), c(sb_2$mat[sb_2$dept==cod[12]]), c(sb_2$mat[sb_2$dept==cod[13]]), c(sb_2$mat[sb_2$dept==cod[14]]), c(sb_2$mat[sb_2$dept==cod[15]]), c(sb_2$mat[sb_2$dept==cod[16]]), c(sb_2$mat[sb_2$dept==cod[17]]), c(sb_2$mat[sb_2$dept==cod[18]]), c(sb_2$mat[sb_2$dept==cod[19]]), c(sb_2$mat[sb_2$dept==cod[20]])
            , c(sb_2$mat[sb_2$dept==cod[21]]), c(sb_2$mat[sb_2$dept==cod[22]]), c(sb_2$mat[sb_2$dept==cod[23]]), c(sb_2$mat[sb_2$dept==cod[24]]), c(sb_2$mat[sb_2$dept==cod[25]]))

Xk2 <- list(cbind(1, sb_2$sexo[sb_2$dept==cod[1]], sb_2$htra[sb_2$dept==cod[1]]), cbind(1, sb_2$sexo[sb_2$dept==cod[2]], sb_2$htra[sb_2$dept==cod[2]]), cbind(1, sb_2$sexo[sb_2$dept==cod[3]], sb_2$htra[sb_2$dept==cod[3]]), cbind(1, sb_2$sexo[sb_2$dept==cod[4]], sb_2$htra[sb_2$dept==cod[4]]), cbind(1, sb_2$sexo[sb_2$dept==cod[5]], sb_2$htra[sb_2$dept==cod[5]])
            , cbind(1, sb_2$sexo[sb_2$dept==cod[6]], sb_2$htra[sb_2$dept==cod[6]]), cbind(1, sb_2$sexo[sb_2$dept==cod[7]], sb_2$htra[sb_2$dept==cod[7]]), cbind(1, sb_2$sexo[sb_2$dept==cod[8]], sb_2$htra[sb_2$dept==cod[8]]), cbind(1, sb_2$sexo[sb_2$dept==cod[9]], sb_2$htra[sb_2$dept==cod[9]]), cbind(1, sb_2$sexo[sb_2$dept==cod[10]], sb_2$htra[sb_2$dept==cod[10]])
            , cbind(1, sb_2$sexo[sb_2$dept==cod[11]], sb_2$htra[sb_2$dept==cod[11]]), cbind(1, sb_2$sexo[sb_2$dept==cod[12]], sb_2$htra[sb_2$dept==cod[12]]), cbind(1, sb_2$sexo[sb_2$dept==cod[13]], sb_2$htra[sb_2$dept==cod[13]]), cbind(1, sb_2$sexo[sb_2$dept==cod[14]], sb_2$htra[sb_2$dept==cod[14]]), cbind(1, sb_2$sexo[sb_2$dept==cod[15]], sb_2$htra[sb_2$dept==cod[15]])
            , cbind(1, sb_2$sexo[sb_2$dept==cod[16]], sb_2$htra[sb_2$dept==cod[16]]), cbind(1, sb_2$sexo[sb_2$dept==cod[17]], sb_2$htra[sb_2$dept==cod[17]]), cbind(1, sb_2$sexo[sb_2$dept==cod[18]], sb_2$htra[sb_2$dept==cod[18]]), cbind(1, sb_2$sexo[sb_2$dept==cod[19]], sb_2$htra[sb_2$dept==cod[19]]), cbind(1, sb_2$sexo[sb_2$dept==cod[20]], sb_2$htra[sb_2$dept==cod[20]])
            , cbind(1, sb_2$sexo[sb_2$dept==cod[21]], sb_2$htra[sb_2$dept==cod[21]]), cbind(1, sb_2$sexo[sb_2$dept==cod[22]], sb_2$htra[sb_2$dept==cod[22]]), cbind(1, sb_2$sexo[sb_2$dept==cod[23]], sb_2$htra[sb_2$dept==cod[23]]), cbind(1, sb_2$sexo[sb_2$dept==cod[24]], sb_2$htra[sb_2$dept==cod[24]]), cbind(1, sb_2$sexo[sb_2$dept==cod[25]], sb_2$htra[sb_2$dept==cod[25]]))


contdept <- data.frame(table(sb_2$dept))
nj2 <- contdept$Freq

X2bt <- Xk2[[1]]%*%B_hat[1,] + theta_hat[1]
for (j in 2:m) {
  X2bt <- c(X2bt, Xk2[[j]]%*%B_hat[j,] + theta_hat[j])
}
y_hat <- X2bt
mean((y2 - y_hat)^2)

load(file = "samples_parcial2_MSE.RData")

save(YP2, YP3, file = "samples_parcial2_MSE.RData")

####### representacion de B1 posterior para cada departamento
library(ggplot2)
Val <- B1_hat

LI95 <- NULL
LS95 <- NULL
LI99 <- NULL
LS99 <- NULL

for(l in 1:m){
  LI95[l] <- quantile(x = B1[,l], probs = 0.025)
  LS95[l] <- quantile(x = B1[,l], probs = 0.975)
  LI99[l] <- quantile(x = B1[,l], probs = 0.005)
  LS99[l] <- quantile(x = B1[,l], probs = 0.995)
}


Dat <- data.frame(Dep=c(unique(training_data$depa),unique(training_data$depa)),
                  Val=c(Val,Val),
                  LI=c(LI95,LI99),
                  LS=c(LS95,LS99),
                  Tipo=c(rep("IC95%",25),rep("IC99%",25)))

Dat$Dep[c(19,44)] <- c("NSANTANDER","NSANTANDER")
#windows(height=3.5,width=7)
ggplot(Dat, aes(x = Dep, y = Val)) +
  geom_point(size=1)+
  geom_linerange(aes(ymax = LS, ymin = LI,color=Tipo),
                 position = position_dodge(width = 0.2))+
  geom_hline(yintercept=50,col="black")+
  labs(x="Departamento",y=expression(beta[1]),color="IC")+
  scale_color_manual(values=c("IC95%"="blue","IC99%"="red"),
                     labels = c("IC al 95%", "IC al 99%"))+
  theme_light()+
  theme(axis.text.x=element_text(angle=35, hjust=1),
        plot.title = element_text(hjust = 0.5))




IC1 <- matrix(0, 25, 2)
for (i in 1:25) {
  IC1[i,] <- quantile(x = B1[,i], probs = c(0.025,0.975))
}





########### representacion de B2 posterior para cada departamento
library(ggplot2)
Val <- B2_hat

LI95 <- NULL
LS95 <- NULL
LI99 <- NULL
LS99 <- NULL

for(l in 1:m){
  LI95[l] <- quantile(x = B2[,l], probs = 0.025)
  LS95[l] <- quantile(x = B2[,l], probs = 0.975)
  LI99[l] <- quantile(x = B2[,l], probs = 0.005)
  LS99[l] <- quantile(x = B2[,l], probs = 0.995)
}


Dat <- data.frame(Dep=c(unique(training_data$depa),unique(training_data$depa)),
                  Val=c(Val,Val),
                  LI=c(LI95,LI99),
                  LS=c(LS95,LS99),
                  Tipo=c(rep("IC95%",25),rep("IC99%",25)))

Dat$Dep[c(19,44)] <- c("NSANTANDER","NSANTANDER")
#windows(height=3.5,width=7)
ggplot(Dat, aes(x = Dep, y = Val)) +
  geom_point(size=1)+
  geom_linerange(aes(ymax = LS, ymin = LI,color=Tipo),
                 position = position_dodge(width = 0.2))+
  geom_hline(yintercept=0,col="black")+
  labs(x="Departamento",y=expression(beta[2]),color="IC")+
  #title="Intervalos de confianza del efecto de la covariable 'Sexo del estudiante' \n en Colombia")+
  scale_color_manual(values=c("IC95%"="blue","IC99%"="red"),
                     labels = c("IC al 95%", "IC al 99%"))+
  theme_light()+
  theme(axis.text.x=element_text(angle=35, hjust=1),
        plot.title = element_text(hjust = 0.5))




IC2 <- matrix(0, 25, 2)
for (i in 1:25) {
  IC2[i,] <- quantile(x = B2[,i], probs = c(0.025,0.975))
}




####### representacion de B3 posterior para cada departamento
library(ggplot2)
Val <- B3_hat

LI95 <- NULL
LS95 <- NULL
LI99 <- NULL
LS99 <- NULL

for(l in 1:m){
  LI95[l] <- quantile(x = B3[,l], probs = 0.025)
  LS95[l] <- quantile(x = B3[,l], probs = 0.975)
  LI99[l] <- quantile(x = B3[,l], probs = 0.005)
  LS99[l] <- quantile(x = B3[,l], probs = 0.995)
}


Dat <- data.frame(Dep=c(unique(training_data$depa),unique(training_data$depa)),
                  Val=c(Val,Val),
                  LI=c(LI95,LI99),
                  LS=c(LS95,LS99),
                  Tipo=c(rep("IC95%",25),rep("IC99%",25)))

Dat$Dep[c(19,44)] <- c("NSANTANDER","NSANTANDER")
#windows(height=3.5,width=7)
ggplot(Dat, aes(x = Dep, y = Val)) +
  geom_point(size=1)+
  geom_linerange(aes(ymax = LS, ymin = LI,color=Tipo),
                 position = position_dodge(width = 0.2))+
  geom_hline(yintercept=0,col="black")+
  labs(x="Departamento",y=expression(beta[3]),color="IC")+
  #title="Intervalos de confianza del efecto de la covariable 'Horas trabajadas por el estudiante' \n en Colombia")+
  scale_color_manual(values=c("IC95%"="blue","IC99%"="red"),
                     labels = c("IC al 95%", "IC al 99%"))+
  theme_light()+
  theme(axis.text.x=element_text(angle=35, hjust=1),
        plot.title = element_text(hjust = 0.5))



IC3 <- matrix(0, 25, 2)
for (i in 1:25) {
  IC3[i,] <- quantile(x = B3[,i], probs = c(0.025,0.975))
}

##### BONDAD DE AJUSTE ###########
## CASANARE
S <- 50000
TS3_cas <- matrix(data = NA, nrow = S, ncol = 1)
ncat <- floor(S/10)                      # mostrar anuncios cada 10%

set.seed(7)
for (s in 1:S) {
  beta <- c(B1[s,24], B2[s,24], B3[s,24])
  theta <- TH[s,24]
  sigma2  <- S2J[s,24]
  yrep    <- rnorm(n = nj[24], mean = Xk[[24]]%*%beta + theta, sd = sqrt(sigma2))
  TS3_cas[s] <- mean(yrep)
  # progreso
  if (s%%ncat == 0) cat("PP M3, ", 100*round(s/S, 1), "% completado \n", sep = "")
}

# grafico
windows(height = 3.5, width = 3.5)
par(mar = c(3,3,1,1), mgp = c(1.75,.75,0))
hist(x = TS3_cas, freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "", ylim = c(0,0.07))
lines(density(TS3_cas))
abline(v = mean(Yk[[24]]), col = "red", lwd = 2, lty = 1)

#ppp
mean(TS3_cas > mean(Yk[[24]]))

## BOGOTA
S <- 50000
TS3_bog <- matrix(data = NA, nrow = S, ncol = 1)
ncat <- floor(S/10)                      # mostrar anuncios cada 10%
set.seed(7)
for (s in 1:S) {
  beta <- c(B1[s,3], B2[s,3], B3[s,3])
  theta <- TH[s,3]
  sigma2  <- S2J[s,3]
  yrep    <- rnorm(n = nj[3], mean = Xk[[3]]%*%beta + theta, sd = sqrt(sigma2))
  TS3_bog[s] <- mean(yrep)
  # progreso
  if (s%%ncat == 0) cat("PP M3, ", 100*round(s/S, 1), "% completado \n", sep = "")
}

# grafico
windows(height = 3.5, width = 3.5)
par(mar = c(3,3,1,1), mgp = c(1.75,.75,0))
hist(x = TS3_bog, freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(TS3_bog))
abline(v = mean(Yk[[3]]), col = "red", lwd = 2, lty = 1)

#ppp
mean(TS3_bog > mean(Yk[[3]]))

## VALLE DEL CAUCA
S <- 50000
TS3_val <- matrix(data = NA, nrow = S, ncol = 1)
ncat <- floor(S/10)                      # mostrar anuncios cada 10%
set.seed(7)
for (s in 1:S) {
  beta <- c(B1[s,22], B2[s,22], B3[s,22])
  theta <- TH[s,22]
  sigma2  <- S2J[s,22]
  yrep    <- rnorm(n = nj[22], mean = Xk[[22]]%*%beta + theta, sd = sqrt(sigma2))
  TS3_val[s] <- mean(yrep)
  # progreso
  if (s%%ncat == 0) cat("PP M3, ", 100*round(s/S, 1), "% completado \n", sep = "")
}

# grafico
windows(height = 3.5, width = 3.5)
par(mar = c(3,3,1,1), mgp = c(1.75,.75,0))
hist(x = TS3_val, freq = F, col = "lightgreen", density = 20, border = "lightgreen", 
     xlab = "", ylab = "", main = "")
lines(density(TS3_val))
abline(v = mean(Yk[[22]]), col = "red", lwd = 2, lty = 1)

#ppp
mean(TS3_val > mean(Yk[[22]]))
