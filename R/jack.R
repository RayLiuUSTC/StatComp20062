#' @title jackknife method
#' @description Compute the jackknife estimate of standard error using R
#' @param data the data given as a vector
#' @param fun function to be bootstrapped
#' @return the standard error to be estimated
#' @examples
#' \dontrun{
#' data <- 20 * rnorm(500,2,3)
#' jack(data = data, fun = mean)
#' }
#' @export
jack <- function(data,fun=NULL){
  theta.hat <- fun(data)
  #set up the bootstrap
  #B is the number of replicates
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


