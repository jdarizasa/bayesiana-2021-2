#################
# Modelo normal #
#################

# Caso de estudio: Longitud de alas de mosquitos
# 
# En 1981, los bi?logos W. L. Grogan y W. W. Wirth descubrieron en las selvas de 
# Brasil dos nuevas variedades de un diminuto insecto picador llamado mosquito 
# (midge). Llamaron a un tipo de mosquito mosquito Apf y al otro mosquito Af. 
# Los bi?logos descubrieron que el mosquito Apf es portador de una enfermedad 
# debilitante que causa inflamaci?n del cerebro cuando un humano est? mordido 
# por un mosquito infectado. Aunque la enfermedad rara vez es fatal, la 
# discapacidad causada por la hinchaz?n puede ser permanente. La otra forma de 
# mosquito, el Af, es bastante inofensiva y un valioso polinizador. 
# En un esfuerzo por distinguir las dos variedades, los bi?logos tomaron medidas 
# en los mosquitos que capturaron. Este es un conjunto de datos valioso para 
# probar m?todos de clasificaci?n.

# Considere los datos de la longitud del ala en mil?metros de nueve miembros de 
# la especie Af de mosquitos. A partir de estas nueve mediciones, se quiere 
# hacer inferencia sobre la media poblacional theta. Otros estudios sugieren que 
# la longitud de las alas suele ser de alrededor de 1.9 mm. Claramente, se tiene 
# que las longitudes deben ser positivas, lo que implica que theta > 0.


########
# data #
########

# install.packages("Flury")
# cargar data de Grogan y Wirth (1981), disponible en la libreria Flury
# esta libraria no se encuentra disonible para versiones resientes de R  
# library(Flury)
# data(midge)  
# y    <- midge[midge[,1]=="Af", 3]

y    <- c(1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)
n    <- length(y)
ybar <- mean(y)
s2   <- var(y)

# Definici?n de la distribuci?n previa (ELICITACIÓN DE LOS HIPERPARAMETROS)
#
# * Otros estudios sugieren que la longitud promedio suele ser de 1.9 mm con una
#   una desviaci?n est?ndar de 0.1 mm.
# * Se usa kappa_0 = nu_0 = 1 para que las previas sean no informativas dado que 
#   esta poblaci?n en particular puede diferir notoriamente de aquellas estudiadas
#   en otras investigaciones
#
# Hiperpar?metros:
# * mu_0      : media a priori 
# * kappa_0   : n. de "observaciones previas" (cantidad de info) asociadas con mu_0  
# * sigma^2_0 : varianza a priori
# * nu_0      : n. de "observaciones previas" (cantidad de info) asociadas con sigma^2_0  
mu0 <- 1.9  
k0  <- 1  ## como para no hacerlo tan informativo, entre más grande más peso
s20 <- 0.1^2  
nu0 <- 1  ## lo mismo, chiquito para no hacerlo tan informativo. Además el tamaño de muestra es pequeño

# inferencia posterior
kn  <- k0 + n 
nun <- nu0 + n
mun <- (k0*mu0 + n*ybar)/kn  
s2n <- (nu0*s20 +(n-1)*s2 + k0*n*(ybar-mu0)^2/kn)/nun

# DISTRIBUCION POSTERIOR
#
# p( theta, sigma^2 | y ) = p( theta | sigma^2, y ) x p(sigma^2 | y)
#
# 1. p( theta | sigma^2, y ) = Normal( theta | 1.814, sigma^2/10 )
# 2. p( sigma^2 | y ) = Gamma-Inversa( sigma^2 | 10/2, 10*0.0153/2 )

# representanci?n gr?fica de la distribucion posterior de ( theta, sigma^2 )

#es como la densidad de GI
dinvgamma <- function(x, a, b, log = FALSE) {
  # calcula la funcion de densidad de una Gamma-Inversa
  out <- a*log(b) - lgamma(a) - (a+1)*log(x) - b/x 
  if (log == FALSE) out <- exp(out)
  return(out)
}

gs <- 250  # n. de puntos a evaluar en un rango de valores

#grillas
theta <- seq(from = 1.6,  to = 2.0,  length = gs)  # theta     : media
is2   <- seq(from = 15,   to = 160,  length = gs)  # 1/sigma^2 : presici?n
s2g   <- seq(from = .001, to = .045, length = gs)  # sigma^2   : varianza  

# evaluar y almacenar la distribuci?n posterior conjunta (escala log) en el rango de valores
# log p( theta, sigma^2 | y ) = log p( theta | sigma^2, y ) + log p(sigma^2 | y)

#ld means log density between theta and sigma or inverse sigma

ld.th.is2 <- matrix(data = NA, nrow = gs, ncol = gs)  # para ( theta, 1/sigma^2 ) 
ld.th.s2  <- matrix(data = NA, nrow = gs, ncol = gs)  # para ( theta, sigma^2 )
for(i in 1:gs) { 
  for(j in 1:gs) {
    ld.th.is2[i,j] <- dnorm(x = theta[i], mean = mun, sd = 1/sqrt(is2[j]*kn), log = T) + dgamma(x = is2[j], shape = nun/2, rate = nun*s2n/2, log = T)
    ld.th.s2 [i,j] <- dnorm(x = theta[i], mean = mun, sd = sqrt(s2g[j]/kn), log = T) + dinvgamma(x = s2g[j], a = nun/2, nun*s2n/2, log = T)
  }
} 

# CURVAS DE NIVEL ESCALA DE GRISES

#exponenciamos para tener una escala real

grays <- gray((10:0)/10)  # paleta de colores (escala de grises)
windows(height=3.5,width=7)
par(mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(1.75,.75,0))
# posterior (theta, 1/sigma^2)
image(x = theta, y = is2, z = exp(ld.th.is2), col = grays, 
      xlab = expression(theta), ylab = expression(tilde(sigma)^2)) 
# posterior (theta, sigma^2)
image(x = theta, y = s2g, z = exp(ld.th.s2), col = grays, 
      xlab = expression(theta), ylab = expression(sigma^2))

# CURVAS DE NIVEL FACHERAS

windows(height=3.5,width=7)
par(mfrow=c(1,2), mar=c(3,3,1,1), mgp=c(1.75,.75,0))
# posterior (theta, 1/sigma^2)
image(x = theta, y = is2, z = exp(ld.th.is2), col = terrain.colors(100), 
      xlab = expression(theta), ylab = expression(tilde(sigma)^2)) 
# posterior (theta, sigma^2)
image(x = theta, y = s2g, z = exp(ld.th.s2), col = heat.colors(100), 
      xlab = expression(theta), ylab = expression(sigma^2))

# GRAFICO 3-D

gs <- 30  # n. de puntos a evaluar en un rango de valores

theta <- seq(from = 1.72,  to = 1.9,  length = gs)  # theta     : media
s2g   <- seq(from = .005, to = .030, length = gs)  # sigma^2   : varianza  

ld.th.s2 <- matrix(data = NA, nrow = gs, ncol = gs)  # para ( theta, sigma^2 )
for(i in 1:gs) { 
  for(j in 1:gs) {
    ld.th.s2[i,j] <- dnorm(x = theta[i], mean = mun, sd = sqrt(s2g[j]/kn), log = T) + dinvgamma(x = s2g[j], a = nun/2, nun*s2n/2, log = T)
  }
}

windows(height=3.5,width=3.5)
par(mfrow=c(1,1), mar=c(3,3,1,1), mgp=c(1.75,.75,0))
persp(x = theta, y = s2g, z = exp(ld.th.s2), theta = 30, phi = 30, expand = 1, 
      xlab = "Media", ylab = "Varianza", zlab = "Densidad", col = "gray95")

#-------------------------------------------------------------------------------
#inferencia
#-------------------------------------------------------------------------------

# GENERACION DE MUESTRAS DE LA DISTRIBUCION POSTERIOR
# * simulacion de monte carlo
# * calcular cualquier cantidad posterior de interes
#   - tendencia
#   - variabilidad
#   - probabilidades
#   - graficos

# n. de simulaciones
S <- 50000

#primero simula los sigma2_n para simular los thetas_n, uno para cada uno

set.seed(1)
s2.postsample    <- 1/rgamma(n = S, shape = nun/2, rate = nun*s2n/2)
theta.postsample <- rnorm(n = S, mean = mun, sd = sqrt(s2.postsample/kn))

# * alternativamente, se podria usan un ciclo for generar los valores uno por uno
# * esta alternativa es mas ineficiente, pero la estructura es util para luego 
#   dise?ar algoritmos mas avanzados
s2.postsample <- NULL
theta.postsample <- NULL
set.seed(1)
for (i in 1:S) {
  # simulacion
  sigma2 <- 1/rgamma(n = 1,  shape = nun/2, rate = nun*s2n/2)
  theta  <- rnorm(n = 1, mean = mun, sd = sqrt(sigma2/kn))
  # almacenamiento
  s2.postsample[i] <- sigma2
  theta.postsample[i] <- theta
  # progreso
  if (i%%5000== 0) cat("Completado ", i/S*100, "% \n", sep = "") 
}

# INFERENCIA POSTERIOR SOBRE theta

# intervalo de credibilidad al 95%
quantile(x = theta.postsample, c(.025,.975))

# intervalo de confianza al 95% bajo normalidad (distribuci?n t)
n    <- length(y)
ybar <- mean(y)
s2   <- var(y)
ybar + qt(p = c(.025,.975), df = n-1)*sqrt(s2/n)

# intervalo de confianza al 95% bajo normalidad (Bootstrap)
library(boot)
f <- function(data, indices) {
  dt <- data[indices]
  mean(dt)
}
set.seed(12345)
Bootstrap <- boot(data = y, statistic = f, R = 50000)
boot.ci(boot.out = Bootstrap)

# media posterior de theta, E( theta | y )
mean(theta.postsample)

# probabilidad posterior de que theta sea mayor que 1.8, Pr( theta > 1.8 | y)
mean( theta.postsample > 1.8 )

# INFERENCIA POSTERIOR SOBRE sigma^2

# intervalo de credibilidad al 95%
quantile(x = s2.postsample, c(.025, .975) )

# media posterior de sigma^2, E( sigma^2 | y )
mean(s2.postsample)

# INFERENICA POSTERIOR SOBRE FUNCIONES DE ( theta, sigma^2 )

# intervalo de credibilidad para la desviaci?n est?ndar sigma
quantile(x = sqrt(s2.postsample), c(.025, .975))

# estimacion puntual del coef. de variaci?n sigma/theta
mean( sqrt(s2.postsample)/theta.postsample )

# Probabilidad de que el CV sea < 7%
mean(sqrt(s2.postsample)/theta.postsample < 0.07)

# intervalo de credibilidad para el coef. de variaci?n sigma/theta
quantile(x = sqrt(s2.postsample)/theta.postsample, c(.025, .975))

# GRAFICO DISTRIBUCION MARGINAL DE THETA
 
windows(height=3.5,width=3.5)
par(mfrow=c(1,1), mar=c(3,3,1,1), mgp=c(1.75,.75,0))
# histograma
hist(x = theta.postsample, freq = F, 
     main="",xlab=expression(theta),xlim=c(1.60,2.0), 
     ylim=c(0,9.5), col = "gray95", border = "darkgray",
     ylab=expression( paste(italic("p("), theta,"| y)",sep="")))
# estimacion kernel
lines(density(theta.postsample,adjust=3)) 
abline(v=quantile(x = theta.postsample, c(.025,.975)), col="gray", lwd=2)

# GRAFICO DISTRIBUCIONES POSTERIORES

gs <- 250
theta <- seq(from = 1.6,  to = 2.0,  length = gs)  # theta     : media
s2g   <- seq(from = .001, to = .045, length = gs)  # sigma^2   : varianza  

ld.th.s2 <- matrix(data = NA, nrow = gs, ncol = gs)  # para ( theta, sigma^2 )
for(i in 1:gs) { 
  for(j in 1:gs) {
    ld.th.s2[i,j] <- dnorm(x = theta[i], mean = mun, sd = sqrt(s2g[j]/kn), log = T) + dinvgamma(x = s2g[j], a = nun/2, nun*s2n/2, log = T)
  }
}

windows(height=7,width=7)
layout(matrix(c(1,1,2,3),2,2,byrow=T))
par(mar=c(3,3,1,1), mgp=c(1.75,.75,0))
# distribucion posterior conjunta
image(x = theta, y = s2g, z = exp(ld.th.s2), col = 'grays', 
      xlab = expression(theta), ylab = expression(sigma^2), 
      xlim = c(1.60,2.0), ylim=c(.001,.07))
# muestras
points(theta.postsample[1:5000], s2.postsample[1:5000], pch = ".", col = "blue")
# distribucion posterior marginal theta
plot(density(theta.postsample,adjust=3), main="", xlab=expression(theta), xlim=c(1.60,2.0),
     ylab=expression(paste(italic("p("), theta,"| y)",sep=""))) 
abline(v=quantile(x = theta.postsample, c(.025,.975)), col="gray", lwd=2)
# distribucion posterior marginal sigma^2
plot(density(s2.postsample,adjust=3), main="", xlab=expression(sigma^2), xlim=c(0,.075),
     ylab=expression(paste(italic("p("), sigma^2,"| y)",sep=""))) 
abline(v=quantile(x = s2.postsample, c(0.025, 0.975)), col="gray", lwd=2)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
