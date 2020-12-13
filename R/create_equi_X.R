#' @title Creating knockoff matrix
#' @description Creating knockoff matrix using Equi-correlated knockoffs
#' @param X the data as a matrix requiring n>=2p
#' @references Barber and Candes, Controlling the false discovery rate via knockoffs. Ann. Statist. 43 (2015), no. 5, 2055–2085
#' @return An object of class "knockoff.variables". This is a list containing at least the following components:\code{n}
#' @return X n-by-p matrix of original variables(standardized).
#' @return Xk	n-by-p matrix of knockoff variables.
#' @return s_j 2*the minimum eigenvalue of 'Sigma'
#' @examples
#' \dontrun{
#' n2 <- 200;p2 <- 100
#' X2 = matrix(rnorm(n2*p2),n2)
#' X2q <- create_equi_X(X2)}
#' @export
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


