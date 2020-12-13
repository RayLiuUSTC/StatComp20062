## ----echo=FALSE---------------------------------------------------------------
library(ggplot2)
myfun <- function(xvar) {
1/(1 + exp(-xvar + 10))
}
ggplot(data.frame(x=c(0, 20)), aes(x=x)) + stat_function(fun=myfun, geom="line")

## ----echo=FALSE---------------------------------------------------------------
data(rock)
x_html <- knitr::kable(cor(rock), "html")
kableExtra::kable_styling(x_html,bootstrap_options = "striped",
                          full_width = F,font_size = 14)

## -----------------------------------------------------------------------------
a <- set.seed(3.3)
n <- 1000
u <- runif(n)
x <- 2/sqrt(1-u)
hist(x, prob = TRUE,breaks = seq(range(x)[1],range(x)[2],length.out = 31),main = bquote(2/sqrt(1-y)),xlim = c(0,60)) #density histogram of sample
y <- seq(0, 50, .01)
lines(y, 8/y^3) #density curve f(x)

## -----------------------------------------------------------------------------
Epan <- function(n){
  a <- runif(n,-1,1)
  b <- runif(n,-1,1)
  c <- runif(n,-1,1)
  d <- ifelse(abs(c)>=abs(a) & abs(c)>=abs(b),b,c)
  return(d)
}
set.seed(3.9)
x <- Epan(1000)
hist(x, prob = TRUE,breaks = seq(range(x)[1],range(x)[2],length.out = 11)) #density histogram of sample
y <- seq(-1,1,.01)
lines(y, 3*(1-y^2)/4) #density curve f(x)

## -----------------------------------------------------------------------------
a <- set.seed(3.13)
n <- 1000
ga <- rgamma(1000,4,2)
y <- rexp(1000,ga)
hist(y, prob = TRUE,breaks = seq(range(y)[1],range(y)[2],length.out = 31),main = "PARETO=GAMMA+exp(gamma)") #density histogram of sample
x <- seq(0, 10, .01)
lines(x, 4*2^4*(2+x)^-5) #density curve f(x)

## -----------------------------------------------------------------------------
a <- set.seed(5.1)
m <- 1e4; x <- runif(m, min=0, max=pi/3)
theta.hat <- mean(sin(x)) * pi/3
print(c(theta.hat,-(cos(pi/3)-cos(0))))

## -----------------------------------------------------------------------------
a <- set.seed(5.7)
m <- 1e5; x <- runif(m)
T1 <- exp(x)
theta1.hat <- mean(T1)#simple MC
x2 <- x[1:m/2]
T2 <- (exp(x2)+exp(1-x2))/2
theta2.hat <- mean(T2)#antithetic variate approach
show <- c(theta1.hat,theta2.hat,exp(1)-1)
names(show) <- c("SMC","antithetic","theoretical integral value")
print(show)
eprv <- (var(T1)-var(T2))/var(T1)
var.g <- (exp(1)^2-1)/2-(exp(1)-1)^2
var.theta <- var.g/m
var.g1 <- ((exp(1)^2-1)/2-(exp(1)-1)^2)*2+2*(exp(1)-(exp(1)-1)^2)
var.theta1 <- var.g1/(2*m)
tprv <- (var.theta- var.theta1)/var.theta
show2 <- c(eprv,tprv)
names(show2) <- c("empirical prv","theoretical prv")
print(show2)


## -----------------------------------------------------------------------------
x <- seq(0, 5, .01)
w <- 2
beta1 <- 2;alpha <- 4
f1 <- dexp(x)
f2 <- dgamma(x,4,2)
g <- x^2/sqrt(2*pi)*exp(-x^2/2)
#figure (a)
plot(x, g, type = "l", main = "", ylab = "",xlim = c(1,5),
ylim = c(0,2), lwd = w)
lines(x, f1, lty = 2, lwd = w)
lines(x, f2, lty = 3, lwd = w)
legend("topright", legend = c("g", 1:2),
lty = 1:3, lwd = w, inset = 0.02)
#figure (b)
plot(x, g/f1, type = "l", main = "", ylab = "",xlim = c(1,5),
ylim = c(0,20), lwd = w, lty = 1)
lines(x, g/f1, lty = 2, lwd = w)
lines(x, g/f2, lty = 3, lwd = w)
legend("topright", legend = c(1:2),
lty = 2:3, lwd = w, inset = 0.02)

## -----------------------------------------------------------------------------
m <- 100000
theta.hat <- var.hat <- numeric(2)
gfun <- function(x) {
x^2/sqrt(2*pi)*exp(-x^2/2) * (x > 1)
}

x <- rexp(m) #using f1
fg <- gfun(x) / dexp(x)
theta.hat[1] <- mean(fg)
var.hat[1] <- var(fg)

x <- rgamma(m,4,2) #using f2
fg <- gfun(x) / dgamma(x,4,2)
theta.hat[2] <- mean(fg)
var.hat[2] <- var(fg)
rbind(theta.hat, var.hat)

## -----------------------------------------------------------------------------
M <- 10000 #number of replicates
k <- 5 #number of strata
r <- M / k #replicates per stratum
N <- 50 #number of times to repeat the estimation
T2 <- numeric(k)
estimates <- matrix(0, N, 2)
g <- function(x) {
exp(-x - log(1+x^2)) * (x > 0) * (x < 1)
}

q <- - log(1 - seq(0,1,0.2) * (1 - exp(-1)))#find the quantile

for (i in 1:N) {
  estimates[i, 1] <- mean(g(runif(M)))
  for (j in 1:k) {
    u <- runif(M/k)   
    x <- - log(exp(-q[j]) - u * (1 - exp(-1)) / 5)
    T2[j] <- mean(g(x)/(5*(exp(-x) / (1 - exp(-1)))))
  }
  estimates[i, 2] <- sum(T2)
}
mean.est <- apply(estimates, 2, mean)
l <- var.est <- apply(estimates, 2, var)
x_html <- knitr::kable(rbind(mean.est,var.est),"html")
kableExtra::kable_styling(x_html,bootstrap_options = "striped",
                          full_width = F,font_size = 14)

## -----------------------------------------------------------------------------
vrp <- as.data.frame((l[1]-l[2])/l[1])#variance reduction
colnames(vrp) <- "vrp"
x_html <- knitr::kable(vrp,"html")
kableExtra::kable_styling(x_html,bootstrap_options = "striped",
                          full_width = F,font_size = 14)

## ----echo=F-------------------------------------------------------------------
example510 <- matrix(c(0.5241140, 0.5313584, 0.5461507, 0.52506988 ,0.5260492,
0.2436559, 0.4181264, 0.9661300, 0.09658794, 0.1427685),2,5,byrow = T)
example510 <- cbind(example510,round(rbind(mean.est,sqrt(var.est)),digits = 7))
rownames(example510) <- c("theta.hat","se")
colnames(example510) <- c(NULL,"f1","f2","f3","f4","f5","stratified","importance stratum")
x_html <- knitr::kable(example510,"html")
kableExtra::kable_styling(x_html,bootstrap_options = "striped",
                          full_width = F,font_size = 14)

## -----------------------------------------------------------------------------
set.seed(6.4)
n <- 20
alpha <- .05
calcCI <- function(n, alpha) {
x <- exp(rnorm(n,0,2))
return(c(mean(log(x)) - sd(log(x)) * qt(1-alpha/2, df = n-1) / sqrt(n),mean(log(x)) + sd(log(x)) * qt(1-alpha/2, df = n-1) / sqrt(n)))
}
CL <- replicate(1000, expr = calcCI(n,alpha) )
mean(0>CL[1,] & 0<CL[2,])

## -----------------------------------------------------------------------------
set.seed(6.5)
n <- 20
alpha <- .05
calcCI <- function(n, alpha) {
x <- rchisq(n,2)
return(c(mean(x) - sd(x) * qt(1-alpha/2, df = n-1) / sqrt(n),mean(x) + sd(x) * qt(1-alpha/2, df = n-1) / sqrt(n)))
}
CL <- replicate(1000, expr = calcCI(n,alpha) )
mean(2>CL[1,] & 2<CL[2,])

## -----------------------------------------------------------------------------
x <- seq(0,1,0.005)
plot(x,dbeta(x,1,1),type = "l",xlim = c(0.1,0.9),ylim = c(0,3.8),ylab = "density")
for (i in 1:10) {
  lines(x,dbeta(x,i,i),type = "l",col=i,lwd=2)
}
legend(0.7,4, legend = paste('a=',1:10), lty = c(1,1), col = 1:10)

## ----warning=FALSE------------------------------------------------------------
sk <- function(x) {
  #computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}#compute skewness of sample

set.seed(6.71)
c.alpha <- .05
n <- 30
m <- 2000
pa.alpha <- seq(1,30,1)

N <- length(pa.alpha)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-c.alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { 
  sktests <- numeric(m)
  for (i in 1:m) {
    x <- rbeta(n, pa.alpha[j], pa.alpha[j])
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}

plot(pa.alpha, pwr, type = "b",
     xlab = bquote(alpha), ylim = c(0,0.1),
     main = "The power estimate of the skewness test of normality against Beta distributions")
abline(h = .05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(pa.alpha, pwr+se, lty = 3)
lines(pa.alpha, pwr-se, lty = 3)

## ----echo=FALSE---------------------------------------------------------------
set.seed(6.72)
c.alpha <- .05
n <- 30
m <- 2000
pa.alpha <- seq(1,30,1)

N <- length(pa.alpha)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(1-c.alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) {
  sktests <- numeric(m)
  for (i in 1:m) {
    x <- rt(n, pa.alpha[j])
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
}

plot(pa.alpha, pwr, type = "b",
     xlab = bquote(alpha),
     main = "The power estimate of the skewness test of normality against t-distributions")
abline(h = .05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(pa.alpha, pwr+se, lty = 3)
lines(pa.alpha, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
set.seed(6.8)
m <- 10000
n <- c(20,50,100)
power.c5 <- power.f <- numeric(3)
sigma1 <- 1
sigma2 <- 1.5
count5test <- function(x, y) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
return(as.integer(max(c(outx, outy)) > 5))
}

# generate samples under H1 to estimate power
for (i in 1:3) {
  power.c5[i] <- mean(replicate(m, expr={
  x <- rnorm(n[i], 0, sigma1)
  y <- rnorm(n[i], 0, sigma2)
  count5test(x, y)
}))
  power.f[i] <- mean(replicate(m, expr={
  x <- rnorm(n[i], 0, sigma1)
  y <- rnorm(n[i], 0, sigma2)
  return(as.integer(var.test(x,y,conf.level = 0.945)$p.value<0.055))
}))
}
power <- cbind(n,power.c5,power.f)
print(power)

## -----------------------------------------------------------------------------
set.seed(6.81)
d <- 1
n <- c(10, 20, 30, 50, 100, 500) #sample sizes
cv <- qchisq(.95,d*(d+1)*(d+2)/6) #crit. values for each n
sk.mt <- function(x) {
  #computes the sample skewness coeff.
  n <- length(x)
  xbar <- mean(x)
  sigma <- var(x)
  xcenter <- x - xbar
  xsum <- sum((outer(xcenter,xcenter,"*")/sigma)^3)
  return( xsum/n^2 )
}#1 dimension situation
#n is a vector of sample sizes
#we are doing length(n) different simulations
p.reject <- numeric(length(n)) #to store sim. results
m <- 10000 #num. repl. each sim.
for (i in 1:length(n)) {
sktests <- numeric(m) #test decisions
for (j in 1:m) {
x <- rnorm(n[i])
#test decision is 1 (reject) or 0
sktests[j] <- as.integer(sk.mt(x)*n[i]/6 >= cv )
}
p.reject[i] <- mean(sktests) #proportion rejected
}
p.reject

## -----------------------------------------------------------------------------
set.seed(6.101)
alpha <- .1
n <- 30
m <- 2500
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qchisq(.9,d*(d+1)*(d+2)/6)
for (j in 1:N) { #for each epsilon
e <- epsilon[j]
sktests <- numeric(m)
for (i in 1:m) { #for each replicate
sigma <- sample(c(1, 10), replace = TRUE,
size = n, prob = c(1-e, e))
x <- rnorm(n, 0, sigma)
sktests[i] <- as.integer(sk.mt(x)*n/6 >= cv )
}
pwr[j] <- mean(sktests)
}
#plot power vs epsilon
plot(epsilon, pwr, type = "b",
xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
library(MASS)
Mardia<-function(mydata){
  n=nrow(mydata)
  c=ncol(mydata)
  central<-mydata
  central <- apply(central,2,function(x) x-mean(x))
  sigmah<-t(central)%*%central/n
  a<-central%*%solve(sigmah)%*%t(central)
  b<-sum(a^{3})/(n*n)
  test<-n*b/6
  chi<-qchisq(0.95,c*(c+1)*(c+2)/6)
  as.integer(test>chi)
}

set.seed(1234)
mu <- c(0,0,0)
sigma <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
m=1000
n<-c(10, 20, 30, 50, 100, 500)
#m: number of replicates; n: sample size
a=numeric(length(n))
for(i in 1:length(n)){
  a[i]=mean(replicate(m, expr={
    mydata <- mvrnorm(n[i],mu,sigma) 
    Mardia(mydata)
  }))
}

## -----------------------------------------------------------------------------
print(a)

## -----------------------------------------------------------------------------
library(MASS)
set.seed(7912)
set.seed(7912)
mu1 <- mu2 <- c(0,0,0)
sigma1 <- matrix(c(1,0,0,0,1,0,0,0,1),nrow=3,ncol=3)
sigma2 <- matrix(c(100,0,0,0,100,0,0,0,100),nrow=3,ncol=3)
sigma=list(sigma1,sigma2)
m=1000
n=50
#m: number of replicates; n: sample size
epsilon <- c(seq(0, .06, .01), seq(.1, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
for (j in 1:N) { #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) { #for each replicate
    index=sample(c(1, 2), replace = TRUE, size = n, prob = c(1-e, e))
    mydata<-matrix(0,nrow=n,ncol=3)
    for(t in 1:n){
      if(index[t]==1) mydata[t,]=mvrnorm(1,mu1,sigma1) 
      else mydata[t,]=mvrnorm(1,mu2,sigma2)
    }
    sktests[i] <- Mardia(mydata)
  }
  pwr[j] <- mean(sktests)
}
plot(epsilon, pwr, type = "b",
     xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .05, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## ----warning=FALSE------------------------------------------------------------
library(bootstrap) #for the law data
data("law")
n <- nrow(law)
theta.hat <- cor(law$LSAT, law$GPA)
##Estimate the bias
theta.jack <- numeric(n)
for (i in 1:n)
theta.jack[i] <- cor(law$LSAT[-i], law$GPA[-i])
bias <- (n - 1) * (mean(theta.jack) - theta.hat)
##Estimate the standard error
se <- sqrt((n-1) * mean((theta.jack - mean(theta.jack))^2))
print(list("bias"=bias,"se"=se))

## ----warning=FALSE------------------------------------------------------------
library(boot)
data("aircondit")
hour <- aircondit$hours
lambda.boot <- function(dat,ind){
  y <- dat[ind]
  mean(y)
}
boot.obj <- boot(data = hour,statistic = lambda.boot, R=2000)
print(boot.obj)
print(boot.ci(boot.obj,conf = 0.95) )
hist(boot.obj$t)


## ----warning=FALSE------------------------------------------------------------
data("scor")
n <- nrow(scor)

Contribution<-function(mydata){
  n=nrow(mydata)
  c=ncol(mydata)
  central<-mydata
  central <- apply(central,2,function(x) x-mean(x))
  sigmah<-t(central)%*%central/n
  eigenh <- eigen(sigmah)
  theta.hat <- eigenh$values[1]/sum(eigenh$values)
  theta.hat
}

theta.hat <- Contribution(scor)
theta.jack <- numeric(n)
for (i in 1:n)
theta.jack[i] <- Contribution(scor[-i,])
bias <- (n - 1) * (mean(theta.jack) - theta.hat)
##Estimate the standard error
se <- sqrt((n-1) * mean((theta.jack - mean(theta.jack))^2))
print(list("bias"=bias,"se"=se))

## ----echo=FALSE,warning=FALSE-------------------------------------------------
B <- 2000 #larger for estimating bias
theta.b <- numeric(B)
tt <- scor
for (b in 1:B) {
i <- sample(1:n, size = n, replace = TRUE)
LSAT <- tt[i,]
theta.b[b] <- Contribution(LSAT)
}
bias <- mean(theta.b - theta.hat)
bias

## ----warning=FALSE,message=FALSE----------------------------------------------
library(DAAG); attach(ironslag)
a <- seq(10, 40, .1)

n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:(n-1)) {
  for (j in (k+1):n) {
    y <- magnetic[-c(k,j)]
    x <- chemical[-c(k,j)]
    
    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[c(k,j)]
    e1[k] <- e1[k] + sum((magnetic[c(k,j)] - yhat1)^2)/2
    
    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[c(k,j)] +
      J2$coef[3] * chemical[c(k,j)]^2
    e2[k] <- e2[k] + sum((magnetic[c(k,j)] - yhat2)^2)/2
    
    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(k,j)]
    yhat3 <- exp(logyhat3)
    e3[k] <- e3[k] + sum((magnetic[c(k,j)] - yhat3)^2)/2
    
    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(k,j)])
    yhat4 <- exp(logyhat4)
    e4[k] <- e4[k] + sum((magnetic[c(k,j)] - yhat4)^2)/2
  }
}
t <- n*(n-1)/2
c(sum(e1)/t,sum(e2)/t,sum(e3)/t,sum(e4)/t)

## -----------------------------------------------------------------------------
library(DAAG)
attach(ironslag)
n <- length(magnetic)
e1 <- numeric(n*(n-1)/2)
e2 <- numeric(n*(n-1)/2)
e3 <- numeric(n*(n-1)/2)
e4 <- numeric(n*(n-1)/2)
count <- 0
for (i in 1:(n-1))
  for (j in (i+1):n) {
    count <- count+1
    y <- magnetic[-c(i,j)]
    x <- chemical[-c(i,j)]
    
    P1 <- lm(y~x)
    y1_1 <- chemical[i]*P1$coef[2] + P1$coef[1]
    y1_2 <- chemical[j]*P1$coef[2] + P1$coef[1]
    e1[count] <- (magnetic[i]-y1_1)^2+(magnetic[j]-y1_2)^2
    
    P2 <- lm(y~x+I(x^2))
    y2_1 <- P2$coef[1] + P2$coef[2] * chemical[i] + P2$coef[3] * chemical[i]^2
    y2_2 <- P2$coef[1] + P2$coef[2] * chemical[j] + P2$coef[3] * chemical[j]^2
    e2[count] <- (magnetic[i]-y2_1)^2+(magnetic[j]-y2_2)^2
    
    P3 <- lm(log(y)~x)
    y3_1 <- exp(P3$coef[1] + P3$coef[2] * chemical[i])
    y3_2 <- exp(P3$coef[1] + P3$coef[2] * chemical[j])
    e3[count] <- (magnetic[i]-y3_1)^2+(magnetic[j]-y3_2)^2
    
    P4 <- lm(log(y)~log(x))
    y4_1 <- exp(P4$coef[1] + P4$coef[2] * log(chemical[i]))
    y4_2 <- exp(P4$coef[1] + P4$coef[2] * log(chemical[j]))
    e4[count] <- (magnetic[i]-y4_1)^2+(magnetic[j]-y4_2)^2
  }

e = c(mean(e1)/2,mean(e2)/2,mean(e3)/2,mean(e4)/2)
matrix(e, nrow=1,
       dimnames=list("prediction error", c("Linear","Quadratic"," Exponential","Log-Log")))
detach(ironslag)

## -----------------------------------------------------------------------------
library(boot)
set.seed(123)
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}

permC <- function(n1,n2,mu1=0,mu2=0,sigma1=1,sigma2=1,R = 100){
  if(n1>n2)stop("n1 should be less than n2")
  reps <- numeric(R)
  x <- rnorm(n1,mu1,sigma1)
  y <- rnorm(n2,mu2,sigma2)
  for (i in 1:R) {
    #generate indices k for the first sample
    k <- sample(1:n2, size = n1, replace = FALSE)
    x1 <- x
    y1 <- y[k] #complement of x1
    reps[i] <- count5test(x1, y1)
  }
  mean(reps)
}

n1 <- 20
n2 <- 40
m <- 1000
alphahat <- mean(replicate(m, expr=permC(n1,n2)))
print(alphahat)


## ----warning=FALSE,message=FALSE----------------------------------------------
library(MASS)
library(RANN)
library(energy)
library(Ball)

alpha <- 0.05
sim <- matrix(0,41,4)

for (i in 0:40) {
  epsilon <- (i+10)*.1
  mu1 <- c(0,0)
sigma1 <- epsilon*diag(nrow = 2)
mu2 <- c(0,0)
sigma2 <- diag(nrow = 2)
X <- mvrnorm(40,mu1,sigma1)
Y <- mvrnorm(40,mu2,sigma2)

#NN test
Z <- rbind(X, Y)
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) # what's the first column?
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}
set.seed(12345); N <- c(nrow(X), nrow(Y))
boot.obj <- boot(data = Z, statistic = Tn, R = 999,
                 sim = "permutation", sizes = N,k=3)
ts <- c(boot.obj$t0,boot.obj$t)
p.value1 <- mean(ts>=ts[1])

#energy test
boot.obs <- eqdist.etest(Z, sizes=N, R=999)
p.value2 <- boot.obs$p.value

#ball test
p.value3 = bd.test(x = X, y = Y, num.permutations = 999)

p <- c(p.value1,p.value2,p.value3$p.value)
names(p) <- c('NN test','energy test','ball test')
sim[i+1,] <- c(epsilon,p)
}
 
plot(sim[,1], sim[,2], ylim = c(0, 1), type = "l",
xlab = bquote(epsilon), ylab = "power")
lines(sim[,1], sim[,3], lty = 2)
lines(sim[,1], sim[,4], lty = 4)
abline(h = alpha, lty = 3)
legend("topright", 1, c('NN test','energy test','ball test'),
lty = c(1,2,4), inset = .02)


## -----------------------------------------------------------------------------

sim <- matrix(0,41,4)

for (i in 0:40) {
  epsilon <- (i+40)*.025
  mu1 <- c(epsilon-1,2*(epsilon-1))
sigma1 <- epsilon*diag(nrow = 2)
mu2 <- c(0,0)
sigma2 <- diag(nrow = 2)
X <- mvrnorm(40,mu1,sigma1)
Y <- mvrnorm(40,mu2,sigma2)

#NN test
Z <- rbind(X, Y)
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) # what's the first column?
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}
set.seed(12345); N <- c(nrow(X), nrow(Y))
boot.obj <- boot(data = Z, statistic = Tn, R = 999,
                 sim = "permutation", sizes = N,k=3)
ts <- c(boot.obj$t0,boot.obj$t)
p.value1 <- mean(ts>=ts[1])

#energy test
boot.obs <- eqdist.etest(Z, sizes=N, R=999)
p.value2 <- boot.obs$p.value

#ball test
p.value3 = bd.test(x = X, y = Y, num.permutations = 999)

p <- c(p.value1,p.value2,p.value3$p.value)
names(p) <- c('NN test','energy test','ball test')
sim[i+1,] <- c(epsilon,p)
}
 
plot(sim[,1], sim[,2], ylim = c(0, 1), type = "l",
xlab = bquote(epsilon), ylab = "p-value")
lines(sim[,1], sim[,3], lty = 2)
lines(sim[,1], sim[,4], lty = 4)
abline(h = alpha, lty = 3)
legend("topright", 1, c('NN test','energy test','ball test'),
lty = c(1,2,4), inset = .02)

## -----------------------------------------------------------------------------

sim <- matrix(0,41,4)

for (i in 0:40) {
  epsilon <- i*.025
X <- rt(40,1)
x1 <- rnorm(40,0,1)
x2 <- rnorm(40,4,2)
u <- runif(40)
k <- as.integer(u>epsilon)
Y <- k*x1+(1-k)*x2

#NN test
Z <- c(X, Y)
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) # what's the first column?
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}
set.seed(12345); N <- c(length(X), length(Y))
boot.obj <- boot(data = Z, statistic = Tn, R = 999,
                 sim = "permutation", sizes = N,k=3)
ts <- c(boot.obj$t0,boot.obj$t)
p.value1 <- mean(ts>=ts[1])

#energy test
boot.obs <- eqdist.etest(Z, sizes=N, R=999)
p.value2 <- boot.obs$p.value

#ball test
p.value3 = bd.test(x = X, y = Y, num.permutations = 999) 

p <- c(p.value1,p.value2,p.value3$p.value)
names(p) <- c('NN test','energy test','ball test')
sim[i+1,] <- c(epsilon,p)
}
 
plot(sim[,1], sim[,2], ylim = c(0, 1), type = "l",
xlab = bquote(epsilon), ylab = "p-value")
lines(sim[,1], sim[,3], lty = 2)
lines(sim[,1], sim[,4], lty = 4)
abline(h = alpha, lty = 3)
legend("topright", 1, c('NN test','energy test','ball test'),
lty = c(1,2,4), inset = .02)

## -----------------------------------------------------------------------------

sim <- matrix(0,21,4)

for (i in 0:20) {
  epsilon <- i+2
  mu1 <- c(0,0)
sigma1 <- 4*diag(nrow = 2)
mu2 <- c(0,0)
sigma2 <- diag(nrow = 2)
X <- mvrnorm(epsilon,mu1,sigma1)
Y <- mvrnorm(epsilon*10,mu2,sigma2)

#NN test
Z <- rbind(X, Y)
Tn <- function(z, ix, sizes,k) {
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) # what's the first column?
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
}
set.seed(12345); N <- c(nrow(X), nrow(Y))
boot.obj <- boot(data = Z, statistic = Tn, R = 999,
                 sim = "permutation", sizes = N,k=3)
ts <- c(boot.obj$t0,boot.obj$t)
p.value1 <- mean(ts>=ts[1])

#energy test
boot.obs <- eqdist.etest(Z, sizes=N, R=999)
p.value2 <- boot.obs$p.value

#ball test
p.value3 = bd.test(x = X, y = Y, num.permutations = 999)

p <- c(p.value1,p.value2,p.value3$p.value)
names(p) <- c('NN test','energy test','ball test')
sim[i+1,] <- c(epsilon,p)
}
 
plot(sim[,1], sim[,2], ylim = c(0, 1), type = "l",
xlab = bquote(epsilon), ylab = "p-value")
lines(sim[,1], sim[,3], lty = 2)
lines(sim[,1], sim[,4], lty = 4)
abline(h = alpha, lty = 3)
legend("topright", 1, c('NN test','energy test','ball test'),
lty = c(1,2,4), inset = .02)

## -----------------------------------------------------------------------------
library(knitr)
set.seed(12345)
dLaplace <- function(x){
  y <- 0.5*exp(-abs(x))
  return(y)
}

rw.Metropolis.L <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (dLaplace(y)/dLaplace(x[i-1])))
      x[i] <- y  
    else {
      x[i] <- x[i-1]
      k <- k + 1
    }
  }
  return(list(x=x, k=k))
}

n <- 4  #degrees of freedom for target Student t dist.
N <- 2000
sigma <- c(.05, .5, 2,  16)

x0 <- 25
rw1 <- rw.Metropolis.L(sigma[1], x0, N)
rw2 <- rw.Metropolis.L(sigma[2], x0, N)
rw3 <- rw.Metropolis.L(sigma[3], x0, N)
rw4 <- rw.Metropolis.L(sigma[4], x0, N)
reject <- c(rw1$k, rw2$k, rw3$k, rw4$k)
no.reject <- NULL
#number of candidate points rejected
no.reject <- data.frame(sigma=sigma,no.reject=reject,accepted.rate=(N-reject)/N)
knitr::kable(no.reject)

#par(mfrow=c(2,2))  #display 4 graphs together
alpha <- 0.025
refline <- c(log(2*alpha),-log(2)-log(alpha))
rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
for (j in 1:4) {
  plot(rw[,j], type="l",
       xlab=bquote(sigma == .(round(sigma[j],3))),
       ylab="X", ylim=range(rw[,j]))
  abline(h=refline)
}
#par(mfrow=c(1,1)) #reset to default

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  
  psi.means <- rowMeans(psi)     #row means
  B <- n * var(psi.means)        #between variance est.
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within est.
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
  r.hat <- v.hat / W             #G-R statistic
  return(r.hat)
}

sigma <- 2
k <- 4          #number of chains to generate
n <- 15000      #length of chains

#choose overdispersed initial values
x0 <- c(-10, -5, 5, 10)

#generate the chains
set.seed(12345)
X <- matrix(0, nrow=k, ncol=n)
for (i in 1:k) {
  dd <- rw.Metropolis.L(sigma=sigma, N=n, x0=x0[i])
  dd3 <- dd$x
  X[i, ] <- dd3
}

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi)) psi[i,] <- psi[i,] / (1:ncol(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0, n)
for (j in 2:n)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[2:n], type="l", xlab="", ylab="R",ylim = c(1,1.3),xlim = c(0,8000))
abline(h=1.2, lty=2)


## -----------------------------------------------------------------------------
k1 <- c(4:25,100,500,1000)
f <- function(a,k=4){
  x <- sqrt((a^2*(k-1))/(k-a^2))
  y <- sqrt((a^2*k)/(k+1-a^2))
  Sk1 <- 1-pt(x,df=k-1)
  Sk2 <- 1-pt(y,df=k)
  Sk1-Sk2
}   
res <- uniroot(f,c(1,2))
res
aa <- seq(-2,2,0.01)
plot(aa,f(aa,k=4),type = "l")
abline(h=0,lty=2)
points(1.49,0,cex=2,col="red")

## -----------------------------------------------------------------------------
res <- numeric(25)
for (i in 1:length(k1)) {
  f <- function(a,k=k1[i]){
    x <- sqrt((a^2*(k-1))/(k-a^2))
    y <- sqrt((a^2*k)/(k+1-a^2))
    Sk1 <- 1-pt(x,df=k-1)
    Sk2 <- 1-pt(y,df=k)
    Sk1-Sk2
  }  
  sol <- round(uniroot(f,c(1,2))$root,5)
  res[i] <- sol
}
sol <- cbind(k1,res)
print(sol)

## ----echo=FALSE,message=FALSE-------------------------------------------------
Genotype <- c("AA","BB","OO","AO","BO","AB","sum")
Frequency <- c("p^2","q^2","r^2","2pq","2pr","2qr","1")
Count <- c("nAA","nBB","nOO","nAO","nBO","nAB","n")
sheet <- cbind(Genotype,Frequency,Count)
library(knitr)
kable(t(sheet))

## -----------------------------------------------------------------------------
bloodtype<-function(p,n.obs){
  n<-sum(n.obs)
  nA<-n.obs[1]
  nB<-n.obs[2]
  nAB<-n.obs[3]
  nOO<-n.obs[4]
  cat(p,"\n")
  pAt<-pBt<-pOt<-Q <- rep(0,20)
  pAt[1]<-p[1]
  pBt[1]<-p[2]
  pOt[1]<-1-p[1]-p[2]
  Q[1] <- 0
  for(i in 2:20){
  pA.old<-pAt[i-1]
  pB.old<-pBt[i-1]
  pO.old<-pOt[i-1]
  
  nAA <- nA * pA.old^2/(pA.old^2+2*pA.old*pO.old)
  nAO <- nA * 2*pA.old*pO.old/(pA.old^2+2*pA.old*pO.old)
  nBB <- nB * pB.old^2/(pB.old^2+2*pB.old*pO.old)
  nBO <- nB * 2*pB.old*pO.old/(pB.old^2+2*pB.old*pO.old)
  
  pAt[i]<-(2*nAA+nAO+nAB)/(2*n)
  pBt[i]<-(2*nBB+nBO+nAB)/(2*n)
  pOt[i]<-(2*nOO+nAO+nBO)/(2*n)
  Q[i] <- 2*nAA*log(pAt[i])+2*nBB*log(pBt[i])+2*nOO*log(pOt[i])+
    nAO*log(2*pAt[i]*pOt[i])+nBO*log(2*pBt[i]*pOt[i])+nAB*log(2*pAt[i]*pBt[i])
  }
  return(list(pAt=pAt,pBt=pBt,pOt=pOt,Q=Q))
}
n.obs<-c(444,132,63,361) # observed data,n_A,n_B,n_AB,n_OO
p<-c(1/3,1/3)
a<-bloodtype(p,n.obs)
print(a)

## ----eval=FALSE---------------------------------------------------------------
#  formulas <- list(
#  mpg ~ disp,
#  mpg ~ I(1 / disp),
#  mpg ~ disp + wt,
#  mpg ~ I(1 / disp) + wt
#  )

## -----------------------------------------------------------------------------
attach(mtcars)
formulas = list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)

for (i in 1:4) {print(summary(lm(formulas[[i]])))}
lapply(formulas, function(x) return (summary(lm(x))))

## ----eval=FALSE---------------------------------------------------------------
#  trials <- replicate(
#  100,
#  t.test(rpois(10, 10), rpois(7, 10)),
#  simplify = FALSE
#  )

## -----------------------------------------------------------------------------
set.seed(12345)
sapply(1:100, function(x) return (t.test(rpois(10, 10), rpois(7, 10))$p.value))

## -----------------------------------------------------------------------------
Mapvapplay<-function (f,n,type, ...) {  
 #n:the length of output  
 f <- match.fun(f) 
 fM=Map(f, ...) 
 if(type=="numeric") return(vapply(fM,cbind,numeric(n))) 
 else if (type=="character") return(vapply(fM,cbind,character(n))) 
 else if (type=="logical") return(vapply(fM,cbind,logical(n))) 
 } 

## -----------------------------------------------------------------------------
library(Rcpp)
sourceCpp("../src/rw_Metropolis.cpp")

dLaplace <- function(x){
  y <- 0.5*exp(-abs(x))
  return(y)
}
rw.Metropolis.L <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (dLaplace(y)/dLaplace(x[i-1])))
      x[i] <- y  
    else {
      x[i] <- x[i-1]
      k <- k + 1
    }
  }
  return(list(x=x, k=k))
}

n <- 4  #degrees of freedom for target Student t dist.
N <- 2000
x0 <- 25
sigma <- c(.05, .5, 2,  16)

rw2L <- rw.Metropolis.L(sigma[3], x0, N)
rw2C <- rw_Metropolis(sigma[3], x0, N)
qqplot(rw2C$x[-(1:800)],rw2L$x[-(1:800)])

library(microbenchmark)
ts <- microbenchmark(rw2L = rw.Metropolis.L(sigma[2], x0, N), rw2C = rw_Metropolis(sigma[2], x0, N))
summary(ts)[,c(1,3,5,6)]



