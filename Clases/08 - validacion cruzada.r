###########################
### Modelos jerarquicos ###
###########################

###########################
# PUNTAJES DE MATEMÁTICAS #
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
# MODELO 1 : MODELO NORMAL SEMI-CONJUGADO #
###########################################

y  <- YM[,2] # valor de la variable
y0 <- y      # replica de y
N  <- length(y)

# folds
L <- 10
set.seed(1)
u <- runif(N)
folds <- rep(NA, N)
lims  <- seq(from = 0, to = 1, len = L+1) 
for (l in 1:L) {
        indices <- (lims[l] < u) & (u < lims[l+1])
        folds[indices] <- l
}
table(folds)

# hiperparametros
mu0 <- 50
t20 <- 25
s20 <- 100
nu0 <- 1

MSE <- NULL
for (l in 1:L) {
        # folds
        cat("Validación cruzada, fold = ", l, "\n", sep = "")
        y <- y0 
        indices <- folds == l
        ytrue <- y[indices]
        y[indices] <- NA
        nmiss <- length(ytrue)
        # inicializar la cadena
        set.seed(l)
        mu   <- rnorm(1, mu0, t20)
        sig2 <- 1/rgamma(1, nu0/2, nu0*s20/2)
        # MCMC
        S <- 10000 
        YMISS <- matrix(NA, S, nmiss)
        set.seed(l)
        for(s in 1:S) {
                # actualizar y
                y[indices] <- rnorm(n = nmiss, mean = mu, sd = sqrt(sig2))
                sum.y  <- sum(y)
                ss.y   <- (N-1)*var(y)
                mean.y <- sum.y/N
                # actualizar el valor de theta
                t2n <- 1/(1/t20 + N/sig2)
                mun <- t2n*(mu0/t20 + sum.y/sig2)
                mu  <- rnorm(1, mun, sqrt(t2n))
                # actualizar el valor de sigma^2
                an  <- (nu0 + N)/2
                bn  <- (nu0*s20 + ss.y + N*(mean.y - mu)^2)/2
                sig2 <- 1/rgamma(1, an, bn)
                # almacenar
                YMISS[s,] <- y[indices]
        }
        # FINAL DEL ALGORITMO
        yhat <- colMeans(YMISS)
        MSE[l] <- mean((ytrue - yhat)^2)
}

MSE
mean(MSE)