###########################
### Modelos jerarquicos ###
###########################

###########################
# PUNTAJES DE MATEM?TICAS #
###########################

# directorio de trabajo
setwd("C:/Users/Juan Camilo/Dropbox/UN/estadistica_bayesiana/codigo/")

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

###########################################
# MODELO 1 : MODELO NORMAL SEMI-CONJUGADO ####
###########################################

g <- YM[,1] # indice del grupo 
y <- YM[,2] # valor de la variable

# estadisticos
N      <- length(y)
mean.y <- mean(y)
var.y  <- var(y)
sum.y  <- N*mean(y)
ss.y   <- (N-1)*var(y)

# estadisticos de prueba
skewness <- function(x) sum((x - mean(x))^3)/(length(x)*sd(x)^3)
kurtosis <- function(x) sum((x - mean(x))^4)/(length(x)*sd(x)^4) - 3
ratrisd  <- function(x) diff(range(quantile(x, c(0.25,0.75))))/sd(x)

ts_obs <- c(mean(y),median(y),sd(y),ratrisd(y),skewness(y),kurtosis(y))
ts_disp <- c("MEDIA","MEDIANA","SD","RAZON RI-SD","SESGO","CURTOSIS")

# hiperparametros
mu0 <- 50
t20 <- 25
s20 <- 100
nu0 <- 1

# setup
S    <- 10000                            # numero de muestras 
MS   <- matrix(data=NA, nrow=S, ncol=2)  # matriz para almacenar las muestras
LP1  <- matrix(data=NA, nrow=S, ncol=1)
ncat <- floor(S/10)                      # mostrar anuncios cada 10%

# inicializar la cadena
mu   <- mean.y
sig2 <- var.y

# MCMC
set.seed(1)
for(s in 1:S) {
        # actualizar el valor de theta
        t2n <- 1/(1/t20 + N/sig2)
        mun <- t2n*(mu0/t20 + sum.y/sig2)
        mu  <- rnorm(1, mun, sqrt(t2n))
        # actualizar el valor de sigma^2
        an  <- (nu0 + N)/2
        bn  <- (nu0*s20 + ss.y + N*(mean.y - mu)^2)/2
        sig2 <- 1/rgamma(1, an, bn)
        # almacenar
        MS[s,] <- c(mu, sig2)
        # log-verosimilitud
        LP1[s] <- sum(dnorm(x = y, mean = mu, sd = sqrt(sig2), log = T))
        # progreso
        if (s%%ncat == 0) cat("MCMC M1, ", 100*round(s/S, 1), "% completado \n", sep = "")
}
# FINAL DEL ALGORITMO

# DIC
mu_hat   <- mean(MS[,1])
sig2_hat <- mean(MS[,2])
lpyth_m1 <- sum(dnorm(x = y, mean = mu_hat, sd = sqrt(sig2_hat), log = T))
pDIC_m1  <- 2*(lpyth_m1 - mean(LP1))
dic_m1   <- -2*lpyth_m1 + 2*pDIC_m1 

# WAIC
lppd_m1  <- 0
pWAIC_m1 <- 0
for (i in 1:N) {
        # lppd
        tmp1    <- dnorm(x = y[i], mean = MS[,1], sd = sqrt(MS[,2]))
        lppd_m1 <- lppd_m1 + log(mean(tmp1))
        # pWAIC
        tmp2 <- dnorm(x = y[i], mean = MS[,1], sd = sqrt(MS[,2]), log = T)
        pWAIC_m1 <- pWAIC_m1 + 2*( log(mean(tmp1)) - mean(tmp2) )
}
waic_m1 <- -2*lppd_m1 + 2*pWAIC_m1

# BIC
k_m1 <- 2
bic_m1 <- -2*lpyth_m1 + k_m1*log(N)

# test-stats
TS1 <- matrix(NA, S, length(ts_obs))    #ESTADISTICOS
YP1 <- matrix(NA, S, N)                 #MUESTRAS SIMULADAS
for (s in 1:S) {
        mu   <- MS[s,1]
        sig2 <- MS[s,2]
        yrep <- rnorm(n = N, mean = mu, sd = sqrt(sig2))
        YP1[s,] <- yrep
        TS1[s,] <- as.numeric(c(mean(yrep), median(yrep), sd(yrep), ratrisd(yrep), skewness(yrep), kurtosis(yrep)))
        # progreso
        if (s%%ncat == 0) cat("PP M1, ", 100*round(s/S, 1), "% completado \n", sep = "")
}

# MSE
mse_m1 <- mean((y - colMeans(YP1))^2)

###########################################
# MODELO 2 : MODELO JERARQUICO PARA MEDIA ####
###########################################

# hiper-parametros para obtener previas no informativos
nu0  <- 1  
s20  <- 100
eta0 <- 1  
t20  <- 100
mu0  <- 50 
g20  <- 25

# setup
S     <- 10000                            # numero de iteraciones
ncat  <- floor(S/10)                      # progreso
THETA <- matrix(data=NA, nrow=S, ncol=m)  # almacenar thetas 
MST   <- matrix(data=NA, nrow=S, ncol=3)  # almacenar mu, sigma^2, tau^2
LP2   <- matrix(data=NA, nrow=S, ncol=1)

# valores inciales
theta  <- ybar
sigma2 <- mean(sv)
mu     <- mean(theta)
tau2   <- var(theta)

# MCMC
set.seed(1)
for (s in 1:S) {
        # actualizar theta
        for(j in 1:m) {
                vtheta   <- 1/(n[j]/sigma2+1/tau2)
                etheta   <- vtheta*(ybar[j]*n[j]/sigma2+mu/tau2)
                theta[j] <- rnorm(1,etheta,sqrt(vtheta))
        }
        # actualizar sigma^2
        nun <- nu0+sum(n)
        ss  <- nu0*s20
        for (j in 1:m) {
                ss <- ss+sum((Y[[j]]-theta[j])^2)
        }
        sigma2 <- 1/rgamma(1,nun/2,ss/2)
        # actualizar mu
        vmu <- 1/(m/tau2+1/g20)
        emu <- vmu*(m*mean(theta)/tau2 + mu0/g20)
        mu  <- rnorm(1,emu,sqrt(vmu)) 
        # actualizar tau2
        etam <- eta0+m
        ss   <- eta0*t20 + sum( (theta-mu)^2 )
        tau2 <- 1/rgamma(1,etam/2,ss/2)
        # almacenar valores
        THETA[s,] <- theta
        MST[s,]   <- c(mu,sigma2,tau2)
        # log-verosimilitud (opcion 1)
        # logver <- 0 
        # for (j in 1:m)
        #         logver <- logver + sum(dnorm(x = Y[[j]], mean = theta[j], sd = sqrt(sigma2), log = T))
        # LP2[s] <- logver
        # log-verosimilitud (opcion 2)
        LP2[s] <- sum(dnorm(x = y, mean = rep(theta, n), sd = sqrt(sigma2), log = T))
        # progreso
        if (s%%ncat == 0) cat("MCMC M2, ", 100*round(s/S, 1), "% completado \n", sep = "")
} 
# FINAL DEL ALGORITMO

# DIC
theta_hat  <- colMeans(THETA)
sigma2_hat <- colMeans(MST)[2]
lpyth_m2   <- sum(dnorm(x = y, mean = rep(theta_hat, n), sd = sqrt(sigma2_hat), log = T))
pDIC_m2    <- 2*(lpyth_m2 - mean(LP2))
dic_m2     <- -2*lpyth_m2 + 2*pDIC_m2 

# WAIC
lppd_m2  <- 0
pWAIC_m2 <- 0
for (i in 1:N) {
        # lppd
        tmp1    <- dnorm(x = y[i], mean = THETA[,g[i]], sd = sqrt(MST[,2]))
        lppd_m2 <- lppd_m2 + log(mean(tmp1))
        # pWAIC
        tmp2 <- dnorm(x = y[i], mean = THETA[,g[i]], sd = sqrt(MST[,2]), log = T)
        pWAIC_m2 <- pWAIC_m2 + 2*( log(mean(tmp1)) - mean(tmp2) )
}
waic_m2 <- -2*lppd_m2 + 2*pWAIC_m2

# BIC
k_m2 <- m + 3
bic_m2 <- -2*lpyth_m2 + k_m2*log(N)

# test-stats
TS2 <- matrix(NA, S, length(ts_obs))
YP2 <- matrix(NA, S, N)
for (s in 1:S) {
        theta   <- THETA[s,]
        sigma2  <- MST[s,2]
        yrep    <- rnorm(n = N, mean = rep(theta, n), sd = sqrt(sigma2))
        YP2[s,] <- yrep
        TS2[s,] <- as.numeric(c(mean(yrep), median(yrep), sd(yrep), ratrisd(yrep), skewness(yrep), kurtosis(yrep)))
        # progreso
        if (s%%ncat == 0) cat("PP M2, ", 100*round(s/S, 1), "% completado \n", sep = "")
}

# MSE
mse_m2 <- mean((y - colMeans(YP2))^2)

#####################################################
# MODELO 3: MODELO JERARQUICO PARA MEDIA Y VAIRANZA ####
#####################################################

# hiper-parametros
eta0 <- 1  ; t20 <- 100
mu0  <- 50 ; g20 <- 25
a0   <- 1  ; b0  <-  1/100 ; wnu0 <- 1

# valores iniciales
theta  <- ybar
sigma2 <- sv 
mu     <- mean(theta)
tau2   <- var(theta)
s20    <- 100
nu0    <- 1

# setup
S      <- 10000
ncat   <- floor(S/10)
THETA  <- matrix(data=NA, nrow=S, ncol=m)
SIGMA2 <- matrix(data=NA, nrow=S, ncol=m)
MTSN   <- matrix(data=NA, nrow=S, ncol=4)  # mu, tau2, sigma02, nu0
LP3    <- matrix(data=NA, nrow=S, ncol=1)
nu0s   <- 1:100 # rango de valores para muestrear en p(nu_0 | rest)

# MCMC
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
        nu0   <- sample(x = nu0s, size = 1, prob=exp( lpnu0-max(lpnu0) ))
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
        # log-verosimilitud (opcion 2)
        LP3[s] <- sum(dnorm(x = y, mean = rep(theta, n), sd = sqrt(rep(sigma2, n)), log = T))
        ### progreso
        if (s%%ncat == 0) cat("MCMC M3, ", 100*round(s/S, 1), "% completado \n", sep = "")
}
######## FINAL DEL ALGORITMO

# DIC
theta_hat  <- colMeans(THETA)
sigma2_hat <- colMeans(SIGMA2)
lpyth_m3   <- sum(dnorm(x = y, mean = rep(theta_hat, n), sd = sqrt(rep(sigma2_hat, n)), log = T))
pDIC_m3    <- 2*(lpyth_m3 - mean(LP3))
dic_m3     <- -2*lpyth_m3 + 2*pDIC_m3 

# WAIC
lppd_m3  <- 0
pWAIC_m3 <- 0
for (i in 1:N) {
        # lppd
        tmp1    <- dnorm(x = y[i], mean = THETA[,g[i]], sd = sqrt(SIGMA2[,g[i]]))
        lppd_m3 <- lppd_m3 + log(mean(tmp1))
        # pWAIC
        tmp2 <- dnorm(x = y[i], mean = THETA[,g[i]], sd = sqrt(SIGMA2[,g[i]]), log = T)
        pWAIC_m3 <- pWAIC_m3 + 2*( log(mean(tmp1)) - mean(tmp2) )
}
waic_m3 <- -2*lppd_m3 + 2*pWAIC_m3

# BIC
k_m3 <- 2*m + 4
bic_m3 <- -2*lpyth_m3 + k_m3*log(N)

# test-stats
TS3 <- matrix(NA, S, length(ts_obs))
YP3 <- matrix(NA, S, N)
for (s in 1:S) {
        theta   <- THETA[s,]
        sigma2  <- SIGMA2[s,]
        yrep    <- rnorm(n = N, mean = rep(theta, n), sd = sqrt(rep(sigma2, n)))
        YP3[s,] <- yrep
        TS3[s,] <- as.numeric(c(mean(yrep), median(yrep), sd(yrep), ratrisd(yrep), skewness(yrep), kurtosis(yrep)))
        # progreso
        if (s%%ncat == 0) cat("PP M3, ", 100*round(s/S, 1), "% completado \n", sep = "")
}

# MSE
mse_m3 <- mean((y - colMeans(YP3))^2)

##########################
# COMPARACION DE MODELOS ####
##########################

# gr?fico log-verosimilitud
windows(height=4,width=8)
par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(x = NA, y = NA, xlab = "Iteraci?n", ylab = "Log-verosimilitud", 
     cex.axis = 0.7, xlim = c(1, S), ylim = range(LP1, LP2, LP3)) 
lines(x = 1:S, y = LP1, type = "l", col = "lightgray")
lines(x = 1:S, y = LP2, type = "l", col = "mistyrose")
lines(x = 1:S, y = LP3, type = "l", col = "lightblue")
abline(h = mean(LP1), lty = 2, col = "gray")
abline(h = mean(LP2), lty = 2, col = "red")
abline(h = mean(LP3), lty = 2, col = "blue")

# gr?fico log-verosimilitud
windows(height=4,width=8)
par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
plot(x = NA, y = NA, ylab = "Densidad", xlab = "Log-verosimilitud", 
     cex.axis = 0.7, xlim = range(LP2, LP3), ylim = c(0, 0.065)) 
hist(LP2, freq = F, add = T, col = "mistyrose", border = "mistyrose")
hist(LP3, freq = F, add = T, col = "lightblue", border = "lightblue")
lines(density(LP2), col = "red")
lines(density(LP3), col = "blue")
abline(v = mean(LP2), lty = 2, col = "red")
abline(v = mean(LP3), lty = 2, col = "blue")

### ppp (posterior predictive p-values)
for (i in 1:3) {
        TS <- get(paste0("TS",i))
        for (j in 1:length(ts_obs)) {
                assign(x = paste0("ppp",j,"_m",i), value = mean(TS[,j] > ts_obs[j]))
        }
        rm(TS)
}

# grafico de los estadisticos de prueba
windows(height=5,width=10)
par(mfrow=c(2,3),mar=c(3,3,1.5,1),mgp=c(1.75,.75,0))
for (i in 1:length(ts_obs)) {
        ts1 <- TS1[,i]
        ts2 <- TS2[,i]
        ts3 <- TS3[,i]
        den1 <- density(ts1, adjust = 1.5)
        den2 <- density(ts2, adjust = 1.5)
        den3 <- density(ts3, adjust = 1.5)
        plot(NA, NA, xlim = range(ts1,ts2,ts3), ylim = c(0, max(den1$y,den2$y,den3$y)), 
             xlab = ts_disp[i], ylab = "Densidad", main = ts_disp[i])
        lines(den1, col = 1)
        lines(den2, col = 2)
        lines(den3, col = 3)
        abline(v = ts_obs[i], col = "gray", lwd = 2)
        legend("topright", legend = c("M1","M2","M3"), col = 1:3, lwd = 2, bty = "n")
}

####### comparacion
tab <- matrix( c(lpyth_m1, lpyth_m2, lpyth_m3,
                 pDIC_m1,  pDIC_m2,  pDIC_m3,
                 dic_m1,   dic_m2,   dic_m3,
                 lppd_m1,  lppd_m2,  lppd_m3,
                 pWAIC_m1, pWAIC_m2, pWAIC_m3,
                 waic_m1,  waic_m2,  waic_m3,
                 bic_m1,   bic_m2,   bic_m3, 
                 mse_m1,   mse_m2,   mse_m3,
                 ppp1_m1,  ppp1_m2,  ppp1_m3,
                 ppp2_m1,  ppp2_m2,  ppp2_m3,
                 ppp3_m1,  ppp3_m2,  ppp3_m3,
                 ppp4_m1,  ppp4_m2,  ppp4_m3,
                 ppp5_m1,  ppp5_m2,  ppp5_m3,
                 ppp6_m1,  ppp6_m2,  ppp6_m3), 14, 3, T)
colnames(tab) <- c("M1","M2","M3")
rownames(tab) <- c("lp","pDIC","DIC","lppd","pWAIC","WAIC","BIC","MSE","ppp media","ppp mediana","ppp sd","ppp ratio","ppp sesgo","ppp curtosis")
round(tab, 2)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------