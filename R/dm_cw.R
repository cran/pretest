#' Diebold-Mariano Test and Clark & West  Test
#'
#' It calculates the original DM statistics and the CW adjusted version of DM statistics, including the versions based on a Newey-West type estimator of the long run variance.
#'
#' @references Clark, T. E., & West, K. D. (2007). Approximately normal tests for equal predictive accuracy in nested models. Journal of econometrics, 138(1), 291-311.
#' @references Diebold, F. X., & Mariano, R. S. (1995). Com paring predictive accu racy. Journal of Business and Economic Statistics, 13(3), 253-263.
#' @param Ehat1 Residual series from Model 1 (the smaller model). One dimension and numeric.
#' @param Ehat2 Residual series from Model 2 (the larger/nested model). One dimension and numeric.
#' @return A list of  statistics and corresponding P values will be produced.
#' @examples
#' e1<- rnorm(15);
#' e2<- rnorm(15);
#' temp1 <- dm_cw(e1,e2)
#'
#' @export



dm_cw=function(Ehat1,Ehat2){

  #Stopping criteria
  if(is.null(Ehat1)){
    stop("Ehat1 must be a numeric series")
  }else if(!is.numeric(Ehat1)){
    stop("Ehat1 must be a numeric series")
  }

  if(is.null(Ehat2)){
    stop("Ehat2 must be a numeric series")
  }else if(!is.numeric(Ehat2)){
    stop("Ehat2 must be a numeric series")
  }

  Ehat1 = as.matrix(Ehat1)
  Ehat2 = as.matrix(Ehat2)

  if(dim(Ehat1)[1] != dim(Ehat2)[1]){
    stop("Ehat1 and Ehat2 must have same length")
  }else if(dim(Ehat1)[2]>1 | dim(Ehat2)[2]>1){
    stop("Ehat1 and Ehat2 must be one dimension")
  }


  #Initialisation of parameters
  n=length(Ehat1)
  Ehat1sq = Ehat1^2
  Ehat2sq = Ehat2^2
  Ehat2sq_adj = Ehat2^2-(Ehat1-Ehat2)^2

  mse_1 = mean(Ehat1sq)
  mse_2 = mean(Ehat2sq)
  mse_2_adj = mean(Ehat2sq_adj)

  #Calculation of DM and CW statistics
  dmspread = Ehat1sq-Ehat2sq
  cwspread = Ehat1sq-Ehat2sq_adj

  psisqdm = t(dmspread-mean(dmspread))%*%(dmspread-mean(dmspread))/n
  psisqcw = t(cwspread-mean(cwspread))%*%(cwspread-mean(cwspread))/n

  dm = (sqrt(n))*mean(dmspread)/sqrt(psisqdm)
  cw = (sqrt(n))*mean(cwspread)/sqrt(psisqcw)

  #Calculation of DM and CW statistics based on Newy-West type of estimator
  psisqdm_nw = lr_var(dmspread-mean(dmspread))
  psisqcw_nw = lr_var(cwspread-mean(cwspread))

  dm_nw = (sqrt(n))*mean(dmspread)/sqrt(psisqdm_nw)
  cw_nw = (sqrt(n))*mean(cwspread)/sqrt(psisqcw_nw)

  #P-values
  pv_dm = 1-stats::pnorm(dm)
  pv_cw = 1-stats::pnorm(cw)
  pv_dm_nw = 1-stats::pnorm(dm_nw)
  pv_cw_nw = 1-stats::pnorm(cw_nw)

  return(list(dm = dm, cw=cw,dm_nw=dm_nw,cw_nw=cw_nw,pv_dm=pv_dm,pv_cw=pv_cw,pv_dm_nw=pv_dm_nw,pv_cw_nw=pv_cw_nw))


}
