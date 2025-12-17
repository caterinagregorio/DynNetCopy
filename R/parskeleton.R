#' Prepare skeleton object for initial values
#'
#' @param K 
#' @param nD 
#' @param nL 
#' @param mapping.to.LP 
#' @param mapping.to.LP2 
#' @param vec_ncol_x0n 
#' @param n_col_x 
#' @param nb_RE 
#' @param stochErr 
#' @param indexparaFixeUser 
#' @param paraFixeUser 
#' @param L 
#' @param paras.ini 
#' @param ncolMod.MatrixY 
#' @param link 
#' @param npara_k 
#' @param Survdata 
#' @param basehaz 
#' @param knots_surv 
#' @param assoc 
#' @param truncation 
#' @param data 
#' @param outcomes 
#' @param df 
#' @param nE 
#' @param np_surv 
#' @param fixed.survival.models 
#' @param interactionY.survival.models 
#' @param nYsurv 
#'
#' @returns parameters object
#'
parskeleton <- function(K,
                       nD,
                       nL,
                       mapping.to.LP,
                       mapping.to.LP2,
                       vec_ncol_x0n,
                       n_col_x,
                       nb_RE,
                       stochErr = FALSE,
                       indexparaFixeUser = NULL,
                       paraFixeUser = NULL,
                       L = 1,
                       paras.ini,
                       ncolMod.MatrixY,
                       link,
                       npara_k,
                       Survdata = NULL,
                       basehaz = NULL,
                       knots_surv = NULL,
                       assoc = NULL,
                       truncation = F,
                       data,
                       outcomes,
                       df,
                       nE = 0,
                       np_surv = 0,
                       fixed.survival.models = NULL,
                       interactionY.survival.models = NULL,
                       nYsurv = 0) {
  
  nb_paraD = nb_RE*(nb_RE+1)/2
  
  #alpha_mu0
  alpha_mu0 <- rep(NA, sum(vec_ncol_x0n))
  
  #alpha_mu
  alpha_mu <- rep(NA, n_col_x)
  
  # para of transition matrix vec_alpha_ij
  vec_alpha_ij <- rep(NA, L * nD * nD)
  
  #paraB
  paraB <- NULL
  if (stochErr == TRUE) {
    paraB <- rep(NA, nD)
  }
  
  #Random effects
  alpha_D <- rep(NA, nb_paraD)
  
  #weights for formative structure
  if(!is.null(nL)){
    para_weights <- rep(NA,nL)
  }
  #paraSig
  paraSig <- rep(NA, K)
  
  
  #Links
  ParaTransformY <- rep(NA,ncolMod.MatrixY)
  #Survival
  para_surv <- NULL
  para_basehaz <- NULL
  if(!is.null(Survdata)){
    np_baz <- ifelse(basehaz=="Weibull",2, 0)# changer 0!!
    
    for (jj in 1:nE){
      para_basehaz <- c(para_basehaz,rep(NA,np_baz)) 
      para_surv <- rep(para_surv,  np_surv[jj]) 
      
    }
    if(basehaz=="Splines") cat('add number of parameters for splines in p and para_surv')
    if(basehaz=="Splines") cat('Define knots_surv para_basehaz')
    
  }
  

  
  
  
  #matrix version
  vec_alpha_ij_m <- matrix(vec_alpha_ij,nrow=nD,byrow = T)
  alpha_D_m <- DparChol(nD,alpha_D)
  
  par_obj <- list(baseline_par=alpha_mu0,
                  mixed_par=alpha_mu,
                  random_par=alpha_D_m,
                  transition_par=vec_alpha_ij_m,
                  markers_par=paraSig,
                  links_par=ParaTransformY)
  if(stochErr)par_obj$stocherr_par <- paraB
  if(!is.null(nL))par_obj$formative_weights <- para_weights
  if(!is.null(Survdata)){
    par_obj$survbase_par <- para_basehaz
    par_obj$survcov_par <- para_surv
  }
  par <- list(paras.ini = par_obj, 
              Fixed.para.index = par_obj, 
              Fixed.para.values = par_obj)
}

revertparskeleton <- function(par_obj){
  names <- names(par_obj)
  par <- c(par_obj$baseline_par,
           par_obj$mixed_par,
           as.numeric(par_obj$random_par[lower.tri(par_obj$random_par,diag=T)]),
           t(as.numeric(par_obj$transition_par)))
  if("stocherr_par"%in%names){
    par <- c(par,par_obj$stocherr_par)
  }
  par <- c(par,par_obj$markers_par,par_obj$links_par,par_obj$survbase_par,
           par_obj$survcov_par,par_obj$formative_weights)
  return(par)
}
