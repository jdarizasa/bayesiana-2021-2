rm(list = ls (all.names = TRUE))

setwd("C:/Users/Juan Camilo/Dropbox/UN/estadistica_bayesiana/codigo")

################################################################################
#                               MODELO BINOMIAL
################################################################################

#-------------------------------------------------------------------------------
# data ///
#-------------------------------------------------------------------------------

# cargar data
Y <- read.csv("GSS1998.txt")

# extraer variable
y <- Y[(Y$YEAR == 1998) & (Y$AGE >= 65) & (Y$FEMALE == 1),]$HAPUNHAP
y[y >  4] <- NA
y[y <= 2] <- 1 
y[y >  2] <- 0
y <- y[!is.na(y)]

# frecuencias
table(y)

# tama?o de muestra
(n <- length(y))

# estadistico suficiente
(sy <- sum(y))

#-------------------------------------------------------------------------------
# distribuci?n posterior con previa Beta(1,1) ///
#-------------------------------------------------------------------------------

# hiperparametros
a <- 1
b <- 1


# parametros de la posterior
(ap <- a + sy )
(bp <- b + n - sy)


# grilla 
theta <- seq(0, 1, length = 1000)


# grafico: verosimilitud y distribucion posterior entre 0 y 1
windows()
par(mfrow = c(2,1), mar = c(3,3,1,1), mgp = c(1.75,.75,0))
# verosimilitud
plot(theta, 10^17 * theta^sy * (1-theta)^(n-sy), type = "l", xlab = expression(theta),
     ylab = expression(paste(10^27, paste("p","(",y," | ",theta,")", sep="")), sep=""))
# posterior
plot(theta, dbeta(theta, ap, bp), type = "l", xlab = expression(theta),
     ylab = expression(paste("p","(",theta," | ",y,")",sep="")))
lines(theta, rep(1, length(theta)), type = "l", col = "darkgray")


# grafico: verosimilitud y distribucion posterior entre 0.8 y 1
windows()
par(mfrow = c(2,1), mar = c(3,3,1,1), mgp = c(1.75,.75,0))
# verosimilitud
plot(theta, 10^17 * theta^sy * (1-theta)^(n-sy), type = "l", xlab = expression(theta),
     ylab = expression(paste(10^27, paste("p","(",y," | ",theta,")", sep="")), sep=""),
     xlim = c(0.8, 1))
# posterior
plot(theta, dbeta(theta, ap, bp), type = "l", xlab = expression(theta),
     ylab = expression(paste("p","(",theta," | ",y,")",sep="")), xlim = c(0.8, 1))
lines(theta, rep(1, length(theta)), type = "l", col = "darkgray")


# moda posterior
(pmd <- (ap-1)/(ap + bp - 2) )


# media posterior
(pmn <- ap/(ap + bp))


# mediana posterior
(pme <- (ap - 1/3)/(ap + bp - 2/3))


# varianza posterior
(pvr <- (ap*bp)/((ap + bp)^2*(ap + bp + 1)))

#cv
sqrt(pvr)/pmn

# otra manera usando simulacion de la posterior
# simulacion
# siempre ponga semilla
set.seed(1234)
theta_pos <- rbeta(n = 10000, shape1 = ap, shape2 = bp)


# media posterior
mean(theta_pos) 


# varianza posterior
var(theta_pos)


# grafico distribucion posterior
windows()
par(mfrow = c(1,1), mar = c(3,3,1,1), mgp = c(1.75,.75,0))
# histograma
hist(theta_pos, freq = F, xlim = c(0.8, 1), col = "gray99", border = "gray", main = "",
     xlab=expression(theta), ylab=expression(paste("p","(",theta," | ",y,")",sep="")))
# estimacion kernel
lines(density(theta_pos), col = "blue")
# distribucion exacta
lines(theta, dbeta(theta, ap, bp), type = "l", col = "red")
# leyenda
legend("topright", legend = c("Distr. emp?rica", "Distr. exacta", "Estimaci?n kernel"), 
       lwd = 2, col = c("gray","red","blue"), bty="n")

#-------------------------------------------------------------------------------
# ejemplos de previas conjugadas y posteriores usando diferentes valores ///
#-------------------------------------------------------------------------------

p <- function(a, b, n, sy, ...) {
        theta <- seq(0, 1, length = 1000)
        plot(theta, dbeta(theta, a+sy, b+n-sy), type = "l", lwd = 2, 
             ylab = expression(paste("p(",theta,"|y)",sep="")), xlab = expression(theta), ...)
        lines(theta, dbeta(theta, a, b), type = "l", col = "gray", lwd = 2)
        mtext(paste("a = ", a, ", b = ", b, ", n = ", n, ", s = ", sy, sep=""), side = 3, line = .1, cex = 0.8)
        legend("topright", legend = c("Previa", "Posterior"), lwd = c(2,2), col = c("gray","black"), bty="n")
}


# graficos
windows()
par(mfrow = c(2,2), mar = c(3,3,1,1), mgp = c(1.75,.75,0), oma = c(0,0,.5,0))
p(a = 1, b = 1, n = 5,   sy = 1)
p(a = 3, b = 2, n = 5,   sy = 1)
p(a = 1, b = 1, n = 100, sy = 20)
p(a = 3, b = 2, n = 100, sy = 20)

#-------------------------------------------------------------------------------
# inferencia sobre theta con previa Beta(1,1) ///
#-------------------------------------------------------------------------------

a <- 1   ;  b <- 1     # previa
n <- 129 ; sy <- 118   # datos


# intervalo de credibilidad usando distribucion exacta
qbeta( c(.025,.975), shape1 = a+sy, shape2 = b+n-sy )


# intervalo de credibilidad usando simulacion
set.seed(1234)
theta_pos <- rbeta(n = 10000, shape1 = a+sy, shape2 = b+n-sy)
quantile(x = theta_pos, probs = c(.025,.975))


# intervalo de confiabilidad basado en percentiles
windows()
par(mar = c(3,3,1,1), mgp = c(1.75,.75,0))
plot(theta, dbeta(theta, sy+1, n-sy+1), type = "l", xlab = expression(theta),
     ylab = expression(paste("p","(",theta," | ",y,")",sep="")), xlim = c(0.8, 1), lwd = 2)
abline(v = qbeta(c(.025,.975), a+sy,b+n-sy), lty = 2, col = "blue")
abline(v = mean(theta_pos), lty = 2, col = "red")


# probabilidad a priori de que theta > 0.9
pbeta(q = 0.9, shape1 = a, shape2 = b, lower.tail = F)


# probabilidad a posteriori de que theta > 0.9 usando distribucion exacta
pbeta(q = 0.9, shape1 = a+sy, shape2 = b+n-sy, lower.tail = F)


# probabilidad a posteriori de que theta > 0.9 usando simulacion
mean(theta_pos > 0.9)

#-------------------------------------------------------------------------------
# inferencia sobre eta = 1-theta con previa Beta(1,1) ///
#-------------------------------------------------------------------------------

# muestras de la distribucion posterior de eta = 1-theta
eta_pos <- 1 - theta_pos


# distribucion posterior de eta
windows()
par(mar = c(3,3,1,1), mgp = c(1.75,.75,0))
hist(eta_pos, freq = F, xlim = c(0, 0.2), col = "gray99", border = "gray", main = "",
     xlab=expression(eta), ylab=expression(paste("p","(",eta," | ",y,")",sep="")))
lines(density(eta_pos), col = "black", lwd = 2)
abline(v = quantile(x = eta_pos, probs = c(.025,.975)), col = "blue", lty = 2)
abline(v = mean(eta_pos), col = "red", lty = 2)

#-------------------------------------------------------------------------------
# cobertura frecuentista del intervalo de credibilidad para theta ///
#-------------------------------------------------------------------------------

conf <- 0
B <- length(theta_pos)
set.seed(1234)
for (i in 1:B) {
        # simulacion de n observaciones de la distribucion predictiva posterior
        y_rep <- rbinom(n = n, size = 1, prob = theta_pos[i])
        # estadistico suficiente
        s_rep <- sum(y_rep)
        # intervalo de credibilidad
        ic <- qbeta( c(.025,.975), shape1 = a+s_rep, shape2 = b+n-s_rep )
        # theta^(b) pertenece al intervalo?
        conf <- conf + as.numeric((ic[1] < theta_pos[i]) & (theta_pos[i] < ic[2]))/B
        if (i%%floor(0.1*B) == 0) cat(i/B*100, "% completado ...", "\n", sep = "")
}

round(100*conf, 2)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------