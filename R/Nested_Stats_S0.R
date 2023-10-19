#' Predictive Accuracy Testing for Nested Environment S^0
#'
#' It calculates the S^0 statistics for nested models with null hypothesis being the two models having equal predictive power following Pitarakis (2023).
#' There are in total four versions of S^0, based on the assumptions of variance (homo or hete) and residuals (original or adjusted).
#' All S^0 statistics will be standarised to a standard N(0,1) normal distribution, and corresponding P values would be provided.
#'
#' @param Ehat1 Residual series from Model 1 (the smaller model). One dimension and numeric.
#' @param Ehat2 Residual series from Model 2 (the larger/nested model). One dimension and numeric.
#' @param lam10 Fraction of the sample used for Model 1, which should be within 0 and 1.
#' @param lam20 Fraction of the sample used for Model 2, which should be within 0 and 1. Note that lam10 cannot equal to lam20 (c.f. Pitarakis, 2023).
#' @return A list of S^0 statistics and corresponding P values will be produced. "adj" means a Clark and West's (2007) reformulation of sample MSE has been applied
#' , and "NW" means robust Newey-West type estimator (c.f. Deng and Perron, 2008) for heteroskedastic errors has been used.
#' @seealso \code{\link{Nested_Stats_Sbar}}
#' @author Rong Peng, \email{r.peng@@soton.ac.uk}
#' @references Pitarakis, J. Y. (2023). A novel approach to predictive accuracy testing in nested environments. Econometric Theory, 1-44.
#' @references Deng, A., & Perron, P. (2008). The limit distribution of the CUSUM of squares test under general mixing conditions. Econometric Theory, 24(3), 809-822.
#' @references Clark, T. E., & West, K. D. (2007). Approximately normal tests for equal predictive accuracy in nested models. Journal of econometrics, 138(1), 291-311.
#' @examples
#' e1<- rnorm(15);
#' e2<- rnorm(15);
#' temp1 <- Nested_Stats_S0(e1,e2,lam10=0.5,lam20=0.8)
#' temp1$S_lam10_lam20_adj_NW     #S^0_T(lam10, lam^20) with CW adjustment and NW correction
#' temp1$pv_S_lam10_lam20_adj_NW  #P value of it
#' @export



Nested_Stats_S0=function(Ehat1,Ehat2,lam10,lam20){

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

  if(is.null(lam10)){
    stop("lambda10 must be between 0 and 1")
  }else if(!is.numeric(lam10)){
    stop("lambda1 must be between 0 and 1")
  }

  if(lam10 > 1 || lam10 <= 0){
    stop("lambda10 must be between 0 and 1")
  }

  if(is.null(lam20)){
    stop("lambda20 must be between 0 and 1")
  }else if(!is.numeric(lam20)){
    stop("lambda20 must be between 0 and 1")
  }

  if(lam20 > 1 || lam20 <= 0){
    stop("lambda20 must be between 0 and 1")
  }


  if(lam10 == lam20 ){
    stop("lambda10 cannot equal to lambda20")
  }


  #Initialisation of parameters
  n=dim(Ehat1)[1]
  l10 = round(n*lam10)
  l20 = round(n*lam20)

  Ehat1sq = Ehat1^2
  Ehat2sq = Ehat2^2
  Ehat2sq_adj = Ehat2^2-(Ehat1-Ehat2)^2#Clark&West adjusted version

  mse_1 = mean(Ehat1sq)
  mse_2 = mean(Ehat2sq)
  mse_2_adj = mean(Ehat2sq_adj)


  #Calculation of homoskedastic variance for Z0 statistics
  psisq = t((Ehat2sq-mse_2))%*%(Ehat2sq-mse_2)/n
  vcoef =(abs(lam10-lam20))/(lam10%*%lam20)
  lrvar_chom = vcoef*psisq
  s_lrvar_chom = sqrt(lrvar_chom) #sigma_homo

  #Calculation of heteroskedastic variance for Z0 statistics
  nuhat = Ehat2sq-mse_2
  psisq_nw = lr_var(nuhat)#long run variance based on Newey-West type estimator
  lrvar_chet = vcoef*psisq_nw
  s_lrvar_chet = sqrt(lrvar_chet) #sigma_hetero

  #Z0 statistics and its Clark&West adjusted version
  Z = (n/l10)*((sum(Ehat1sq[1:l10])/sqrt(n))-(l10/l20)*(sum(Ehat2sq[1:l20])/sqrt(n)))
  Z_adj = (n/l10)*((sum(Ehat1sq[1:l10])/sqrt(n))-(l10/l20)*(sum(Ehat2sq_adj[1:l20])/sqrt(n)))


  #Output S0 statistics and P-values
  S_lam10_lam20 = Z/s_lrvar_chom #homo S0: standardise Z0 with homo sd
  S_lam10_lam20_adj = Z_adj/s_lrvar_chom #homo S0_adj: standardise Z0_adj with homo sd

  S_lam10_lam20_NW = Z/s_lrvar_chet #hetero S0: standardise Z0 with with NW type sd
  S_lam10_lam20_adj_NW = Z_adj/s_lrvar_chet #hetero S0_adj: standardise Z0_adj with NW type sd

  pv_S_lam10_lam20 = 1-stats::pnorm(S_lam10_lam20)
  pv_S_lam10_lam20_adj = 1-stats::pnorm(S_lam10_lam20_adj)

  pv_S_lam10_lam20_NW = 1-stats::pnorm(S_lam10_lam20_NW)
  pv_S_lam10_lam20_adj_NW = 1-stats::pnorm(S_lam10_lam20_adj_NW)

return(list(S0=S_lam10_lam20,S0_adj=S_lam10_lam20_adj, S0NW=S_lam10_lam20_NW,S0NW_adj=S_lam10_lam20_adj_NW,P_S0=pv_S_lam10_lam20,P_S0adj=pv_S_lam10_lam20_adj,P_S0NW=pv_S_lam10_lam20_NW,P_S0NW_adj=pv_S_lam10_lam20_adj_NW))

}

