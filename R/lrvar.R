#' Heteroskedastic Long run variance
#'
#' Long-run covariance estimation using Newey-West (Bartlett) weights
#'
#' Copyright: Kevin Sheppard
#' Kevin.sheppard\@economics.ox.ac.uk
#' Revision: 3    Date: 5/1/2007
#'
#' @param u P by K vector of residual series, for which we recommend to use the recursive residuals from larger model.
#' @param nlag Non-negative integer containing the lag length to use. If empty or not included,
#' nleg = min(floor(1.2*T^(1/3)),T) will be used.
#' @param demean Logical true of false (0 or 1) indicating whether the mean should be subtracted when computing.
#' @return K by K vector of Long run variance using Newey-West (Bartlett) weights.
#' @examples
#' x<- rnorm(15);
#' #Newey-West covariance with automatic BW selection
#' lrcov = lr_var(x)
#' #Newey-West covariance with 10 lags
#' lrcov = lr_var(x, 10)
#' #Newey-West covariance with 10 lags and no demeaning
#' lrcov = lr_var(x, 10, 0)
#' @export


lr_var = function(u,nlag=NULL,demean=TRUE){

  #Stopping criteria
  if(is.null(u)){
    stop("u must be a numeric series")
  }else if(!is.numeric(u)){
    stop("u must be a numeric series")
  }

  if(length(dim(u))>2){
    stop('u must be a P by K matrix of data.')
  }

  P = dim(as.matrix(u))[1]
  K = dim(as.matrix(u))[2]

  if(is.null(nlag)){
    nlag = min(floor(1.2*P^(1/3)),P)
  }

  if(is.null(demean)){
    demean = TRUE
  }

  if(demean == 1){
    demean = TRUE
  }
  else if(demean == 0){
    demean = FALSE
  }else if(is.logical(demean)==FALSE){
      stop('DEMEAN must be either logical true or false')
  }

  if(nlag != floor(nlag) || nlag < 0){
    stop('NLAG must be a non-negative integer')
  }


  u = as.matrix(u)

  #Initialisation of parameters
  if(demean){
    u = u - kronecker(matrix(1,P,K),mean(u))
  }



  #NW weights
  w = rep(NA,nlag+1)
  w=(nlag+1-(0:nlag))/(nlag+1)

  #Calculate covariance
  V=t(u)%*%u/P

  for (i in 1:nlag){
    Gammai=(t(u[(i+1):P,])%*%u[1:(P-i),])/P
    GplusGprime=Gammai+t(Gammai)
    V=V+w[i+1]*GplusGprime
  }

return(V)

}
