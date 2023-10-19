#' Predictive Accuracy Testing for Nested Environment SBar
#'
#' It calculates the SBar statistics for nested models with null hypothesis being the two models having equal predictive power following Pitarakis (2023).
#' There are in total four versions of SBar, based on the assumptions of variance (homo or hete) and residuals (original or adjusted).
#' All SBar statistics will be standarised to a standard N(0,1) normal distribution, and corresponding P values would be provided.
#'
#' @param Ehat1 Residual series from Model 1 (the smaller model). One dimension and numeric.
#' @param Ehat2 Residual series from Model 2 (the larger/nested model). One dimension and numeric.
#' @param lam20 Fraction of the sample used for Model 2, which should be within 0 and 1.
#' @param tau0 Fraction to determine the user-chosen range of lam10 over which the average is taken.
#' @return A list of SBar statistics and corresponding P values will be produced. "adj" means a Clark and West's (2007) reformulation of sample MSE has been applied
#' , and "NW" means robust Newey-West type estimator (c.f. Deng and Perron, 2008) for heteroskedastic errors has been used.
#' @seealso \code{\link{Nested_Stats_S0}}
#' @author Rong Peng, \email{r.peng@@soton.ac.uk}
#' @references Pitarakis, J. Y. (2023). A novel approach to predictive accuracy testing in nested environments. Econometric Theory, 1-44.
#' @references Deng, A., & Perron, P. (2008). The limit distribution of the CUSUM of squares test under general mixing conditions. Econometric Theory, 24(3), 809-822.
#' @references Clark, T. E., & West, K. D. (2007). Approximately normal tests for equal predictive accuracy in nested models. Journal of econometrics, 138(1), 291-311.
#' @examples
#' e1<- rnorm(15);
#' e2<- rnorm(15);
#' temp1 <- Nested_Stats_S0(e1,e2,lam10=0.5,lam20=0.8)
#' temp1$S_lam10_lam20_adj_NW     #\S^0_T(lam10, lam^20) with CW adjustment and NW correction
#' temp1$pv_S_lam10_lam20_adj_NW  #P value of it
#' @export



Nested_Stats_Sbar=function(Ehat1,Ehat2,lam20,tau0){

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

  if(is.null(lam20)){
    stop("lambda20 must be between 0 and 1")
  }else if(!is.numeric(lam20)){
    stop("lambda20 must be between 0 and 1")
  }

  if(lam20 > 1 || lam20 <= 0){
    stop("lambda20 must be between 0 and 1")
  }

  if(is.null(tau0)){
    stop("tau0 must be between 0 and 1")
  }else if(!is.numeric(tau0)){
    stop("tau0 must be between 0 and 1")
  }

  if(tau0 > 1 || tau0 < 0){
    stop("tau0 must be between 0 and 1")
  }


  #Initialisation of parameters
  n=length(Ehat1)

  l20 = round(n*lam20)
  taua = round(n*tau0)
  taub = round(n*(1-tau0))

  zvec =matrix(0, n,1)
  zvec_adj =matrix(0,n,1)

  Ehat1sq = Ehat1^2;
  Ehat2sq = Ehat2^2;
  Ehat2sq_adj = Ehat2^2-(Ehat1-Ehat2)^2#Clark&West adjusted version

  mse_1 = mean(Ehat1sq)
  mse_2 = mean(Ehat2sq)
  mse_2_adj = mean(Ehat2sq_adj)

  #Calculation of homoskedastic variance for Zbar statistics
  psisq = t((Ehat2sq-mse_2))%*%(Ehat2sq-mse_2)/n

  vartau=0

  if ( tau0 == 0 ){
    vartau=(1+2*(lam20)*log(lam20))/(lam20)
  } else if (tau0>0 & tau0<1 & lam20<=tau0 ) {
    vartau = (((1-tau0)^2)+2*lam20*(1-tau0+log(tau0)))/(lam20*((1-tau0)^2))
  } else if (tau0>0 & tau0<1 & lam20>tau0) {
    vartau = (1-(tau0^2)+2*lam20*((1-tau0)*log(lam20)+tau0*log(tau0)))/(lam20*((1-tau0)^2))
  }

  lrvar_chom = vartau*psisq
  s_lrvar_chom = sqrt(lrvar_chom) #sigma_homo


  #Calculation of heteroskedastic variance for Zbar statistics
  nuhat = Ehat2sq-mse_2
  psisq_nw = lr_var(nuhat) #long run variance based on Newey-West type estimator
  lrvar_chet = vartau*psisq_nw
  s_lrvar_chet = sqrt(lrvar_chet) #sigma_hetero


  #Zbar statistics and its Clark&West adjusted version
  for (l1  in (taua+1):n){
    zvec[l1] = (n/l1)*((sum(Ehat1sq[1:l1])/sqrt(n))-(l1/l20)*(sum(Ehat2sq[1:l20])/sqrt(n)))
    zvec_adj[l1] = (n/l1)*((sum(Ehat1sq[1:l1])/sqrt(n))-(l1/l20)*(sum(Ehat2sq_adj[1:l20])/sqrt(n)))
  }

  Zbar = sum(zvec)/taub
  Zbar_adj = sum(zvec_adj)/taub

  #Output Sbar statistics and P-values
  Sbar_tau0_lam20 = Zbar/s_lrvar_chom #homo Sbar: standardise Zbar with homo sd
  Sbar_tau0_lam20_adj = Zbar_adj/s_lrvar_chom #homo Sbar_adj: standardise Zbar_adj with homo sd
  pv_Sbar_tau0_lam20 = 1-stats::pnorm(Sbar_tau0_lam20)
  pv_Sbar_tau0_lam20_adj = 1-stats::pnorm(Sbar_tau0_lam20_adj)

  Sbar_tau0_lam20_NW = Zbar/s_lrvar_chet #hetero Sbar: standardise Zbar with NW type sd
  Sbar_tau0_lam20_adj_NW = Zbar_adj/s_lrvar_chet #hetero Sbar_adj: standardise Zbar_adj with NW type sd
  pv_Sbar_tau0_lam20_NW = 1-stats::pnorm(Sbar_tau0_lam20_NW)
  pv_Sbar_tau0_lam20_adj_NW = 1-stats::pnorm(Sbar_tau0_lam20_adj_NW)

  return(list(SBar=Sbar_tau0_lam20,SBar_adj=Sbar_tau0_lam20_adj, SBarNW=Sbar_tau0_lam20_NW,SBarNW_adj=Sbar_tau0_lam20_adj_NW,P_SBar=pv_Sbar_tau0_lam20,P_SBaradj=pv_Sbar_tau0_lam20_adj,P_SBarNW=pv_Sbar_tau0_lam20_NW,P_SBarNW_adj=pv_Sbar_tau0_lam20_adj_NW))




}
