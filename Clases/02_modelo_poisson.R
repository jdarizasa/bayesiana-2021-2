rm(list = ls (all.names = TRUE))

setwd("C:/Users/Juan Camilo/Dropbox/UN/estadistica_bayesiana/codigo")

##################
# MODELO POISSON #
##################

#-------------------------------------------------------------------------------
# data ///
#-------------------------------------------------------------------------------

# cargar data
Y <- read.csv("GSS1998.txt")
# y1 : sin pregrado
# y2 : con pregrado o mas
y1 <- Y[(Y$FEMALE == 1) & (Y$YEAR >= 1990) & (Y$AGE == 40) & (Y$DEG <3 ),]$CHILDS
y2 <- Y[(Y$FEMALE == 1) & (Y$YEAR >= 1990) & (Y$AGE == 40) & (Y$DEG >=3),]$CHILDS
y1 <- y1[!is.na(y1)]
y2 <- y2[!is.na(y2)] #eliminamos los NaN


# tama?os de muestra
(n1 <- length(y1))
(n2 <- length(y2))


# estadisticos suficientes (sumatoria de y)
(s1 <- sum(y1))
(s2 <- sum(y2))

#ESCOGER EL MODELO...PERO ¿COMO?
# tendecia
summary(y2)
summary(y1)


# variabilidad
var(y2)
var(y1)

#COMENTARIOS
#note que la media y la varianza no cuadra, esto nos da un picor de que el modelo poisson debemos chequearlo
#podriamos pensar en un modelo poisson con sub dispersion

# chequeo informal del modelo poisson
round(abs(table(y1)/n1 - dpois(x = 0:6, lambda = mean(y1))), 3)
round(abs(table(y2)/n2 - dpois(x = 0:4, lambda = mean(y2))), 3)

#Calculamos las frecuencias observadas de cada punto (numero de hijos) y las frecuencias
#teoricas con la distribución Poisson y su estimación de theta (mean(y))

# distribucion de frecuencias: Analisis descriptivo grafico
windows(height = 3.5, width = 7)
par(mfrow = c(1,2), mar=c(3,3,1,1), mgp = c(1.75,.75,0))
# y1
plot(100*table(y1)/n1, type = "h", xlab = "y", ylab = "F. Relativa", col = gray(.5), 
     lwd = 3, ylim = c(0, 50))
mtext("Menos que pregrado", side = 3)
# y2
plot(100*table(y2)/n2, type = "h", xlab = "y", ylab = "F. Relativa", col = gray(0), 
     lwd = 3, ylim = c(0, 50))
mtext("Pregrado o más", side = 3)

#-------------------------------------------------------------------------------
# distribucion Gamma para varios valores de a y b ///
#-------------------------------------------------------------------------------

p <- function(a, b, ...) {
        x <- seq(.001, 10, length = 1000)
        plot(x, dgamma(x, a, b), type = "l", xlab = expression(theta), lwd = 2,
             ylab = expression(paste("p(",theta,")",sep="")), ...)
        mtext(paste("a=", a, ", b=", b, sep=""), side = 3, line = .12, cex = 0.8)
}


windows(width = 15, height = 10)
par(mfrow = c(2,3), mar = c(3,3,1.5,1), mgp = c(1.75,.75,0))
p(a = 1 , b = 1,  col = 2, ylim = c(0, 1.15))
p(a = 2 , b = 2,  col = 2, ylim = c(0, 1.15))
p(a = 4 , b = 4,  col = 2, ylim = c(0, 1.15))
p(a = 2 , b = 1,  col = 4, ylim = c(0, 1.15))
p(a = 8 , b = 4,  col = 4, ylim = c(0, 1.15))
p(a = 32, b = 16, col = 4, ylim = c(0, 1.15))

#-------------------------------------------------------------------------------
# distribuci?n posterior con previa Gamma(2,1) ///
#-------------------------------------------------------------------------------

# hiperparametros
a <- 2
b <- 1


# parametros de la posterior
(ap1 <- a + s1)
(bp1 <- b + n1)
(ap2 <- a + s2)
(bp2 <- b + n2)


# grilla 
theta <- seq(0, 5, length = 1000)
y <- 0:12


# La distinci?n entre los eventos ( theta_1 > theta_2 ) y ( y*_1 > y*_2 ) es
# extremadamente importante: una fuerte evidencia de una diferencia entre dos 
# poblaciones no significa que la diferencia pr?ctica en s? misma sea grande.


# grafico: distribuciones posterior y predictiva posterior
windows(height = 5, width = 10)
par(mfrow = c(1,2), mar = c(3,3,1,1), mgp = c(1.75,.75,0))
# distribucion posterior
#vea que esto si es continuo
xtheta <- seq(0,5,length=1000)
plot(NA, NA, xlim = c(0,5), ylim = c(0,3), xlab=expression(theta),
     ylab = expression(paste("p","(",theta," | ",y,")",sep="")))
lines(theta, dgamma(theta, shape = ap1, rate = bp1), col = 2)
lines(theta, dgamma(theta, shape = ap2, rate = bp2), col = 4)
lines(theta, dgamma(theta, shape = a  , rate = b  ), col = 1, lty = 2)
abline(h = 0, col = 1)
legend("topright", legend = c("Menos que pregrado", "Pregrado o mas", "Previa"), 
       bty = "n", lwd = 2, col = c(2, 4, 1))

# distribucion predictiva posterior
#vea que y es discreto por eso hago type h
plot(y - .07, dnbinom(y, size = ap1, mu = ap1/bp1), col = 2, type = "h",
     ylab = "p(y* | y )", xlab = "y*", ylim = c(0, .35), lwd = 3)
points(y + .07, dnbinom(y, size = ap2, mu = ap2/bp2), col = 4, type = "h", lwd = 3)


# media posterior e intervalo de credibilidad
tab <- cbind(c(ap1/bp1, qgamma(p = c(.025,.975), shape = ap1, rate = bp1)),
             c(ap2/bp2, qgamma(p = c(.025,.975), shape = ap2, rate = bp2)))
colnames(tab) <- c("Menos que pregrado", "Pregrado o mas")
rownames(tab) <- c("Media", "Q2.5%", "Q97.5%")
round(t(tab), 3)

#comprobar y ser flexible si hay diferencias significativas entre los grupos
# muestras de monte carlo
set.seed(1234)
th1_mc <- rgamma(n = 10000, shape = ap1, rate = bp1)
th2_mc <- rgamma(n = 10000, shape = ap2, rate = bp2)


# media posterior e intervalo de credibilidad
tab_mc <- cbind(c(mean(th1_mc), quantile(th1_mc, probs = c(.025,.975))),
                c(mean(th2_mc), quantile(th2_mc, probs = c(.025,.975))))
colnames(tab_mc) <- c("Menos que pregrado", "Pregrado o mas")
rownames(tab_mc) <- c("Media", "Q2.5%", "Q97.5%")
round(t(tab_mc), 3)


# comparacion
round(t(tab),    3)
round(t(tab_mc), 3)


# Pr( theta_i > 2 | y )
tab <- cbind(c(pgamma(q = 2, shape = ap1, rate = bp1, lower.tail = F), mean(th1_mc > 2)),
             c(pgamma(q = 2, shape = ap2, rate = bp2, lower.tail = F), mean(th2_mc > 2)))
colnames(tab) <- c("Menos que pregrado", "Pregrado o mas")
rownames(tab) <- c("Directo", "Simulacion")
round(t(tab), 3)


# Pr( theta_1 > theta_2 | y )
# ?como se lleva a cabo el calculo directo?
mean( th1_mc > th2_mc )


# Pr( y_1^* > y_2^* | y)
set.seed(1234)
y1_mc <- rpois(n = 10000, lambda = th1_mc)
y2_mc <- rpois(n = 10000, lambda = th2_mc)
mean(y1_mc > y2_mc)


# probabilidades de la distribucion predictiva
tab <- rbind(dnbinom(x = 0:8, size = ap1, mu=ap1/bp1), table(y1_mc)/length(y1_mc))
colnames(tab) <- 0:8
rownames(tab) <- c("Directo", "Simulacion")
round(tab, 3)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
