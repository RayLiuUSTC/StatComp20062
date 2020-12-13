# StatComp20062

The package is a collection of two functions about my recent work. One is related to the knockoffs method used in the variable selection while the other is related to resampling.

## Overview

__StatComp20062__ is a simple R package developed to compute the  standard error of an estimator with jackknife method, and to compute the knockoff matrix which is useful in knockoff method. Two functions are included, namely, _jack_ and _create.equi.X_ . For each function, examples are given.

## Introduction to _jack_ 

The source R code for _jack_ is as follows:

***Jackknife***
```{r}
jack <- function(data,fun=NULL){
  theta.hat <- fun(data)
  n <- length(data)      #sample size
  M <- numeric(n)
  for (i in 1:n) { #leave one out
    y <- data[-i]
    M[i] <- fun(y)
  }
  Mbar <- mean(M)
  se.jack <- sqrt(((n - 1)/n) * sum((M - Mbar)^2))
  return(se.jack)
}
```

## Introduction to _create.equi.X_ 

The source R code for _create.equi.X_ is as follows:

***knockoff***
```{r}
create_equi_X <- function(X){
  n = nrow(X)
  p = ncol(X)
  #The rows of input matrix must be more than the columns num
  if (n < 2*p) 
    stop("Input X must have dimensions n >= 2*p") 
  #Compute l-2 norm of each column
  Xsqure <- apply(X,2,FUN = function(x)sqrt(sum(x^2))) 
  #Standardize the input matrix by column
  X3 <- sweep(X,2,Xsqure,FUN = "/")
  #Compute the Gram matrix of the input matrix, which is named as 'Sigma' in the 'knockoff' paper
  Gram <- t(X3)%*%X3
  GramD <- eigen(Gram)
  sj <- min(2*min(GramD$values),1)#extract the minimum eigenvalue of 'Sigma', to compute the identical s_j for all j
  diags <- diag(rep(sj,p))#construct diag{s}
  #Compute the matrix C in the formula
  CC <- 2*diags- diags%*%solve(Gram)%*%diags
  a <- svd(CC)
  C <- t(a$u%*%diag(sqrt(a$d)))
  #Compute the matrix Utide in the formula, which is orthogonal to the X3
  Utide <- svd(X3,nu=n,nv=0)$u[,(n-p+1):n]
  #Compute the knockoff matrix Xtide 
  Xtide <- X3%*%(diag(nrow = p)-solve(Gram)%*%diags)+Utide%*%C
  structure(list(X = X3, Xk = Xtide, s_j = sj), class = "knockoff.variables")
}
```





 
 
