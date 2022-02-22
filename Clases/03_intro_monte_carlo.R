rm(list = ls (all.names = TRUE))

setwd("C:/Users/Juan Camilo/Dropbox/UN/estadistica_bayesiana/codigo")

#############################
# SIMULACION DE MONTE CARLO #
#############################

#-------------------------------------------------------------------------------
# aproximacion de una densidad gamma con diferentes niveles de precision ///
#-------------------------------------------------------------------------------

a <- 68
b <- 45


# muestras (simulacion iid)
set.seed(1)
theta.sim10   <- rgamma(n = 10,   shape = a, rate = b)
theta.sim100  <- rgamma(n = 100,  shape = a, rate = b)
theta.sim1000 <- rgamma(n = 1000, shape = a, rate = b)


theta.support <- seq(from = 0, to = 3,length = 100)
xrange <- c(.75,2.25)
yrange <- c(0,2.5)


# grafico
windows(height = 3, width = 6)
par(mar = c(3,3,.25,1), mgp = c(1.75,.75,0))
par(mfrow = c(2,3))

hist(x = theta.sim10, prob = T, xlim = xrange, ylim = yrange, xlab = "", ylab = "", main = "", col = "gray95")
lines(x = theta.support, dgamma(x = theta.support, shape = a, rate = b), col = "blue")
text(2.1,2.25, paste("B=10", sep=""))

hist(x = theta.sim100, prob = T, xlim = xrange, ylim = yrange, xlab = "", ylab = "", main = "", col = "gray95")
lines(x = theta.support, dgamma(x = theta.support, shape = a, rate = b), col = "blue")
text(2.1,2.25, paste("B=100", sep=""))

hist(x = theta.sim1000, prob = T, xlim = xrange, ylim = yrange, xlab = "", ylab = "", main = "", col = "gray95")
lines(x = theta.support, dgamma(x = theta.support, shape = a, rate = b), col = "blue")
text(2.1,2.25, paste("B=1000", sep=""))

plot(density(theta.sim10),xlim = xrange, ylim = yrange, xlab = expression(theta), ylab = "", main = "")
lines(x = theta.support, dgamma(x = theta.support, shape = a, rate = b), col = "blue")

plot(density(theta.sim100),xlim = xrange, ylim = yrange, xlab = expression(theta), ylab = "", main = "")
lines(x = theta.support, dgamma(x = theta.support, shape = a, rate = b), col = "blue")

plot(density(theta.sim1000),xlim = xrange, ylim = yrange, xlab = expression(theta), ylab = "", main = "")
lines(x = theta.support, dgamma(x = theta.support, shape = a, rate = b), col = "blue")

#la aproximación se hace mejor cuando se aumenta el número de muestras
#no interprete las muestras generadas como observaciones, nunca

#-------------------------------------------------------------------------------
# Ejemplo 1:modelo Gamma-Poisson ///
#-------------------------------------------------------------------------------

# hiperparametros: a = 2, b = 1
a <- 2
b <- 1


# datos: s = 66, n = 44 (con educacion superior)
s <- 66
n <- 44


# muestras (simulacion iid)
set.seed(1)
theta.sim10   <- rgamma(n = 10,  shape = a+s, rate = b+n)
theta.sim100  <- rgamma(n = 100, shape = a+s, rate = b+n)
theta.sim1000 <- rgamma(n = 1000,shape = a+s, rate = b+n)

# media posterior exacta
(a+s)/(b+n) 


# media posterior aproximada
mean(theta.sim10)
mean(theta.sim100)
mean(theta.sim1000)


# probabilidades exacta
pgamma(q = 1.75, shape = a+s, rate = b+n)


# probabilidades aproximada
mean(theta.sim10 < 1.75)
mean(theta.sim100 < 1.75)
mean(theta.sim1000 < 1.75)


# intervalo de credibilidad al 95% exacto
qgamma(p = c(.025,.975), shape = a+s, rate = b+n)


# intervalo de credibilidad al 95% aproximado
quantile(x = theta.sim10, c(.025,.975))
quantile(x = theta.sim100, c(.025,.975))
quantile(x = theta.sim1000, c(.025,.975))

#-------------------------------------------------------------------------------
# Ejemplo 2: modelo Gamma-Poisson ///
#-------------------------------------------------------------------------------

# hiperparametros
a <- 2
b <- 1


# datos sin educacion superior
s1 <- 217
n1 <- 111
# datos con educacion superior
s2 <- 66
n2 <- 44   


# numero de muestras de MC 
B <- 100000


# muestras (simulacion iid)
set.seed(1)
theta1.mc <- rgamma(n = B, shape = a+s1, rate = b+n1)
theta2.mc <- rgamma(n = B, shape = a+s2, rate = b+n2)


# Pr(theta_1 > theta_2 | y)
mean(theta1.mc > theta2.mc)

#esto es booleano entonces tranquilidad, ya hace el 1/B(sum(I(thetha e A)))

# gamma = theta_1/theta_2
gamma.mc <- theta1.mc/theta2.mc


# Pr(gamma > 1 | y) = Pr(theta_1 > theta_2 | y)
mean(gamma.mc > 1)


# intervalo de credibilidad para gamma
quantile(gamma.mc, probs = c(0.025, 0.975))

#quantile saca los cuantiles empiricos

#estimacion puntual del cociente 

mean(gamma.mc)
# grafico
windows(height = 3.5, width = 3.5)
par(mar = c(3,3,1,1), mgp = c(1.75,.75,0))
plot(density(theta1.mc/theta2.mc), main="", xlim=c(.75,2.25),
     xlab = expression(gamma), ylab = expression(paste("p(",gamma," | y)", sep="")))
abline(v = 1, col = 2)

#como hay tanta masa despues de 1 se puede concluir que hay diferencias


# distribucion predictiva posterior
set.seed(1)
y1.mc <- rpois(n = B, lambda = theta1.mc)
y2.mc <- rpois(n = B, lambda = theta2.mc)


# media predictiva posterior exacta (binomal negativa)
(a+s1)/(b+n1)


# media predictiva posteriori aproximada
mean(y1.mc)

ap1 <- a + s1
bp1 <- b + n1
#probabilidades de la distribución predictiva
tab <- rbind(dnbinom(x=0:10, size = ap1, mu=ap1/bp1), table(y1.mc)/B)
colnames(tab) <- 0:10
rownames(tab) <- c("Directo", "Simulación")
round(tab,3)

# Pr(y*_1 > y*_2 | y)
mean(y1.mc > y2.mc)


# d <- y*_1 - y*_2
diff.mc <- y1.mc - y2.mc


# grafico
windows(height = 3.5, width = 7)
par(mar = c(3,3,1,1), mgp = c(1.75,.75,0))
plot(table(diff.mc)/B, type = "h",lwd = 2, xlab = "d", ylab = "p(d | y )", col = "blue")

#-------------------------------------------------------------------------------
# chequeo del modelo ///
#-------------------------------------------------------------------------------

setwd("C:/Users/Juan Camilo/Dropbox/UN/estadistica_bayesiana/codigo")


# cargar data
Y <- read.csv("GSS1998.txt")
# y1 : sin pregrado
# y2 : con pregrado o mas
y1 <- Y[(Y$FEMALE == 1) & (Y$YEAR >= 1990) & (Y$AGE == 40) & (Y$DEG <3 ),]$CHILDS
y2 <- Y[(Y$FEMALE == 1) & (Y$YEAR >= 1990) & (Y$AGE == 40) & (Y$DEG >=3),]$CHILDS
y1 <- y1[!is.na(y1)]
y2 <- y2[!is.na(y2)]


# grupo 1
# estadistico de prueba: promedio
# estadistico observado
t.obs <- mean(y1)

# distribucion posterior estadistico de prueba
t.mc <- rep(NA, B)
set.seed(1)
for(i in 1:B) {
        y_rep   <- rpois(n = n1, lambda = theta1.mc[i])
        t.mc[i] <- mean(y_rep)
        if (i%%floor(0.1*B) == 0) cat(i/B*100, "% completado ...", "\n", sep = "")
}


# grafico
windows(height = 3.5, width = 3.5)
par(mar = c(3,3,1,1), mgp = c(1.75,.75,0))
hist(x = t.mc, freq = F, col = "lightblue", density = 20, border = "lightblue", 
     xlab = "t", ylab = "p(t | y)", main = "")
lines(density(t.mc))
abline(v = t.obs, col = "red", lwd = 2, lty = 1)


# ppp: (posterior predictive p-value)
mean( t.mc > t.obs )

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------