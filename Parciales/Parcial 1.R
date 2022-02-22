#d)
n1 <- 727
n2 <- 583
n3 <- 137
n <- 1447
theta1 <- n1/n
theta2 <- n2/n
theta3 <- n3/n

(theta1 - theta2) - qnorm(0.9995)*sqrt((theta1+theta2-(theta1-theta2)^2)/n)

(theta1 - theta2) + qnorm(0.9995)*sqrt((theta1+theta2-(theta1-theta2)^2)/n)

#e)
set.seed(8)
rdirichlet <- function(a1,a2,a3,b){
  x1 <- rgamma(1,a1,b)
  x2 <- rgamma(1,a2,b)
  x3 <- rgamma(1,a3,b)
  t1 <- x1/(x1+x2+x3)
  t2 <- x2/(x1+x2+x3)
  t3 <- x3/(x1+x2+x3)
  return(c(t1,t2,t3))
}
rdirichlet(1,2,2,1)


rdirichlet2 <- function(x,b){
  g <- NULL
  for (i in 1:length(x)) {
    g[i] <- rgamma(1,x[i],b)
  }
  t <- g/sum(g)
  return(t)
}
set.seed(8)
rdirichlet2(x=c(1,2,2),b=1)
rdirichlet2(x=c(1,2,2,2,4,1,34,43),b=1)

#e.i)
#analisis de la previa
prior1 <- matrix(NA, nrow = 100000, ncol = 3)
prior2 <- matrix(NA, nrow = 100000, ncol = 3)
prior3 <- matrix(NA, nrow = 100000, ncol = 3)
set.seed(8)
for (j in 1:100000) {
    prior1[j,] <- rdirichlet2(x=c(0.01,0.01,0.01),b=1)
    prior2[j,] <- rdirichlet2(x=c(0.01,0.2,0.1),b=1)
    prior3[j,] <- rdirichlet2(x=c(0.5,0.5,0.5),b=1)
}

library(FactoClass)
par(mfrow=c(1,3))
scatterplot3d (prior1, xlab = 'x1', ylab = 'x2', zlab = 'x3', main ="Prior1",type ="p",color =" blue ",box =FALSE ,las =1, angle = 205)
scatterplot3d (prior2, xlab = 'x1', ylab = 'x2', zlab = 'x3', main ="Prior2",type ="p",color =" blue ",box =FALSE ,las =1, angle = 205)
scatterplot3d (prior3, xlab = 'x1', ylab = 'x2', zlab = 'x3', main ="Prior3",type ="p",color =" blue ",box =FALSE ,las =1, angle = 205)
?scatterplot3d

#posterior
N <- c(727,583,137)
alpha <- c(0.01,0.01,0.01)
B <- 100000
theta <- matrix(NA, nrow = B, ncol = 3)
set.seed(8)
for (i in 1:B) {
  theta[i,] <- rdirichlet2(x= alpha + N, 1)
}

scatterplot3d (theta, xlab = 'x1', ylab = 'x2', zlab = 'x3', main =expression(theta),type ="p",color =" blue ",box =FALSE ,las =1, angle = 205)

head(theta)
gamma <- theta[,1] -theta[,2]
head(gamma)

#e.ii)
#descripciones numericas
summary(gamma)

apply(theta, MARGIN = 2, FUN = summary)

#descripciones graficas

#gamma
par(mfrow=c(1,1))
hist(gamma, freq = F, ylab = 'F')
lines(density(gamma), col = "red", lwd = 2)
abline(v = mean(gamma),col="green",lty=2)
abline(v = quantile(gamma, probs = c(0.0005,0.9995)),col="blue", lty=2)

par(mfrow=c(1,3))
#theta1
hist(theta[,1], freq = F)
lines(density(theta[,1]), col= "red", lwd=2)
abline(v=mean(theta[,1]),col="green",lty=2)
abline(v=quantile(theta[,1], probs = c(0.025,0.975)),col="blue", lty=2)
#theta2
hist(theta[,2], freq = F)
lines(density(theta[,2]), col= "red", lwd=2)
abline(v=mean(theta[,2]),col="green",lty=2)
abline(v=quantile(theta[,2], probs = c(0.025,0.975)),col="blue", lty=2)

#theta3
hist(theta[,3], freq = F)
lines(density(theta[,3]), col= "red", lwd=2)
abline(v=mean(theta[,3]),col="green",lty=2)
abline(v=quantile(theta[,3], probs = c(0.025,0.975)),col="blue", lty=2)

#e.iv)
mean(gamma > 0)
sd(gamma>0)
sd(gamma)/sqrt(B)

?gamma
an <- alpha+N
(gamma(sum(an))/(gamma(an[1])*gamma(an[2])*gamma(an[3])))*(1/an[3]*an[2])*(1/(an[1]+an[2]))*300
