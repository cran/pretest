#' Forecasting h-steps ahead using Recursive Least Squares Slow
#'
#' Consider the following LS-fitted Model with intercept:
#' y_(t+h) = beta_0 + x_t * beta + u_(t+h)
#' which is used to generate out-of-sample forecasts of y, h-steps ahead (h=1,2,3,. . . ).
#' It calculates the recursive residuals starting from the first (n * pi0) data points, where n is the total number of data points.
#'
#' recursive_hstep_fast is the fast version that avoids the recursive calculation of inverse of the matrix using Sherman-Morrison formula.
#' recursive_hstep_slow is the slow version that calculates the standard OLS recursively.
#'
#' @param y n x 1 Outcome series, which should be numeric and one dimensional.
#' @param x n x p Predictor matrix (intercept would be added automatically).
#' @param pi0 Fraction of the sample, which should be within 0 and 1.
#' @param h Number of steps ahead to predict, which should be a positive integer.
#' @return Series of residuals estimated
#' @author Rong Peng, \email{r.peng@@soton.ac.uk}
#' @examples
#' x<- rnorm(15);
#' y<- x+rnorm(15);
#' temp2 <- recursive_hstep_slow(y,x,pi0=0.5,h=1);
#' @export


recursive_hstep_slow=function(y,x,pi0,h){

  #Stopping criteria
  if(is.null(y)){
    stop("y must be one dimension")
  }

  y <- as.matrix(y)
  ny=dim(y)[1]
  py=dim(y)[2]

  if(py > 1){
    stop("y must be one dimension")
  }

  if(anyNA(as.numeric(y))){
    stop("y must not contain NA")
  }

  if(is.null(x)){
    stop("x must be one dimension")
  }

  x <- as.matrix(x)
  n=dim(x)[1]
  x =  cbind(rep(1,n),x)
  p=dim(x)[2] #[n,p] = size(X);

  if(anyNA(as.numeric(x))){
    stop("x must not contain NA")
  }

  if(ny != n){
    stop("y and x must have same length of rows")
  }

  if(is.null(pi0)){
    stop("pi0 must be between 0 and 1")
  }else if(!is.numeric(pi0)){
    stop("pi0 must be between 0 and 1")
  }

  if(pi0 >= 1 || pi0 <= 0){
    stop("pi0 must be between 0 and 1")
  }

  if(is.null(h)){
    stop("h must be a positive integer")
  }else if(!is.numeric(h)){
    stop("h must be a positive integer")
  }

  h <- as.integer(h)

  if(h <= 0 || h > (n-1)){
    stop("h must be a positive integer")
  }

  #Initialisation of parameters
  k0=round(n * pi0)
  bhat=matrix(NA,p,n-k0+1)
  zehat = matrix(NA,(n-k0-h+1),1)

  #iteractively update Beta
  for (s in k0:n)
  {
    bhat[,(s-k0+1)]=solve(t(x[1:(s-h),])%*%x[1:(s-h),])%*%t(x[1:(s-h),])%*%y[(h+1):s]
  }

  #update recursive residuals
  for (i in (k0+h):n )
  {
    zehat[i-k0-h+1] = y[i]-x[i-h,]%*%bhat[,(i-k0-h+1)]
  }


  return(zehat)
}

