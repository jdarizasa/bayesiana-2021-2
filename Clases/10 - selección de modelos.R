####################
# REGRESION LINEAL #
####################

# Diabetes 
# Efron, Hastie, Johnstone and Tibshirani (2003) "Least Angle Regression" 
# http://statweb.stanford.edu/~tibs/ftp/lars.pdf
# 
# Data:
# * n = 442 pacientes con diabetes.
# * variables de referencia: edad, sexo, índice de masa corporal, presión arterial 
#   promedio y seis mediciones de suero sanguíneo.
# * variable respuesta: medida cuantitativa de progresión de la enfermedad tomada 
#   un año después de las mediciones iniciales.
#
# Objetivo:
# construir un modelo predictivo para la respuesta a partir de las variables de 
# referencia.
#
# Variables regresoras:
# * la relación entre la respuesta y las variables de referencia puede ser no lineal,
#   así que se recomienda incluir en el modelo los términos de segundo orden como 
#   x_j^2 y x_j*x_k para mejorar la predicción.
# * 10 efectos principales  : x_j
# * 10C2 = 45 interacciones : x_j*x_k
# * 9 términos cuadráticos  : x_j^2 (sex es binatia asi que sex = sex^2)
# * p = 64 variables regresoras en total.
#
# Estrategia
# * Para evaluar los modelos, se dividen aleatoriamente los 442 individuos en dos 
#   grupos: entrenamiento (342 individuos) y prueba (100 individuos).
# * Se ajusta el modelo usando los datos de entrenamiento y luego se evaluan las 
#   predicciones con los datos de prueba.
# * Se usa el MSE como medida de ajuste. 

# data
# Para ayudar con la interpretación de los parámetros y poner los regresores en 
# una escala común, todas las variables se han centrado y escalado de modo que y 
# y las columnas de X tengan media 0 y una varianza 1.

########
# DATA #
########

# cargar data
data(diabetes, package = "lars")

# respuesta
yf <- diabetes$y
yf <- (yf-mean(yf))/sd(yf)

# regresores transformados
Xf <- diabetes$x2
Xf <- t((t(Xf)-apply(Xf,2,mean))/apply(Xf,2,sd))

# tamaños
n <- dim(Xf)[1]
p <- dim(Xf)[2] 

# datos de entrenamiento (training) y de prueba (test)
# indices
set.seed(123)
i.te <- sample(x = 1:n, size = 100, replace = F)
i.tr <- (1:n)[-i.te]
# training
y    <- yf[i.tr ] 
X    <- Xf[i.tr,]
y.te <- yf[i.te ]
X.te <- Xf[i.te,]

#######
# OLS #
#######

source('C:/Users/Juan Camilo/Dropbox/UN/estadistica_bayesiana/codigo/backselect.r')

# p = 64
# y^_test = X_test beta^_training
olsfit   <- lm(y ~ -1 + X)
y.te.ols <- X.te%*%(olsfit$coef)

# backwards elimination 
vars <- bselect.tcrit(y, X, tcrit = 1.65)
bslfit <- lm(y ~ -1 + X[,vars$remain])
y.te.bsl <- X.te[,vars$remain]%*%(bslfit$coef)
length(vars$remain)

# MSE
round(mean((y.te - y.te.ols)^2), 4)
round(mean((y.te - y.te.bsl)^2), 4)

# correlacion
round(cor(y.te, y.te.ols), 4)
round(cor(y.te, y.te.bsl), 4)

windows(height=3, width=9)
par(mfrow=c(1,3), mar=c(2.75,3,1.5,.5), mgp=c(1.5,.5,0))
# full test
plot(x = y.te, y = y.te.ols, cex = 1, xlab = expression(italic(y)[test]), 
     ylab = expression(hat(italic(y))[test]), main = "p = 64")
abline(a = 0, b = 1)
# full magnitudes beta
plot(x = olsfit$coef, type = "h", lwd = 2, xlab = "Regresor", 
     ylab = expression(hat(beta)[ols]), main = "p = 64")
# backwards test 
plot(x = y.te, y = y.te.bsl, cex = 1, xlab = expression(italic(y)[test]), 
     ylab = expression(hat(italic(y))[test]), ylim = range(c(y.te.bsl,y.te.ols)),
     main = "Backwards elimination")
abline(a = 0, b = 1)

##################################
# SELECCION DE MODELOS BAYESIANO #
##################################

source('C:/Users/Juan Camilo/Dropbox/UN/estadistica_bayesiana/codigo/regression_gprior.r')

# settings
S    <- 1000
BETA <- matrix(NA, S, p)
Z    <- matrix(NA, S, p)

# valor inicial
z <- rep(1,p)
lpy.c <- lpy.X(y, X[,z == 1,drop = F])

# cadena
set.seed(1)
windows(width = 7, height = 3.5)
par(mfrow=c(1,2),mar=c(2.75,3,.5,.5), mgp=c(1.7,.7,0))
for (s in 1:S) {
      for (j in sample(1:p)) {
            zp <- z 
            zp[j] <- 1 - zp[j]
            lpy.p <- lpy.X(y,X[,zp == 1,drop = F])
            r <- (lpy.p - lpy.c)*(-1)^(zp[j]==0)
            z[j] <- rbinom(n = 1, size = 1, prob = 1/(1+exp(-r)))
            if(z[j] == zp[j]) {
                  lpy.c <- lpy.p
            }
      }
      beta <- z
      if (sum(z) > 0) {
            beta[z == 1] <- lm.gprior(y, X[,z==1,drop=FALSE],S=1)$beta 
      }
      # almacenar
      Z[s,]    <- z
      BETA[s,] <- beta
      if (s > 1) {
            # beta media posterior
            bpm <- apply(X = BETA[1:s,], MARGIN = 2, FUN = mean) 
            # cadena Prop. z
            Zcp <- apply(X = Z[1:s,,drop = F], MARGIN = 2, FUN = cumsum)/(1:s)
            # plots
            plot(bpm, type = "h", lwd = 2, xlab = "Regresor", ylab = "Media post.", cex.axis = 0.7)
            plot(x = c(1,s), y = range(Zcp), type="n", xlab = "Iteración", ylab = "Prop. z", cex.axis = 0.7)
            apply(X = Zcp, MARGIN = 2, FUN = lines)
            cat("s = ", s, ", Media z = ", round(mean(z), 3), ", MSE = ", round(mean((y.te - X.te%*%bpm)^2), 3), "\n", sep = "")
      }
}

save(BETA, Z, file = "C:/Users/Juan Camilo/Dropbox/UN/estadistica_bayesiana/codigo/BZ_bma.RData")

# inferencia
load(file = "C:/Users/Juan Camilo/Dropbox/UN/estadistica_bayesiana/codigo/BZ_bma.RData")
beta.bma <- apply(X = BETA, MARGIN = 2, FUN = mean, na.rm = TRUE)
y.te.bma <- X.te%*%beta.bma
mean((y.te - y.te.bma)^2)

# plot
windows(height=1.75,width=5)
par(mar=c(2.75,3,.5,.5), mgp=c(1.7,.7,0))
layout(matrix(c(1,1,2), nrow=1, ncol=3))
#
plot(apply(X = Z, MARGIN = 2, FUN = mean, na.rm = TRUE), type = "h",lwd = 2, xlab="Regresor", 
     ylab=expression(paste( "Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")))
#
plot(x = y.te, y = y.te.bma, xlab=expression(italic(y)[test]), ylab = expression(hat(italic(y))[test]))
abline(0,1)