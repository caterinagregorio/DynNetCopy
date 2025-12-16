#' Title
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
#' @param data_F 
#'
#' @returns
#' @export
#'
#' @examples
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
                       nYsurv = 0,
                       data_F = data_F) {
  
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
  #paraSig
  paraSig <- rep(NA, K)
  
  
  #Links
  ParaTransformY <- c()
  for (k in 1:K) {
    if (link[k] == 'thresholds') {
      ParaTransformY <- c (ParaTransformY, c(
        2 * qunif(0.98) * (-median(data[, outcomes[k]]) + min(data[, outcomes[k]]) +
                             1) / (length(unique(data[, outcomes[k]])) - 2),
        rep(sqrt(2 * qunif(0.98) /
                   (
                     length(unique(data[, outcomes[k]])) - 2
                   )), length(unique(data[, outcomes[k]])) - 2)
      ))
    } else{
      ParaTransformY <- c(ParaTransformY, rep(1, ncolMod.MatrixY))
    }
  }
  
  para_surv <- NULL
  para_basehaz <- NULL
  knots_surv <- c(0, 0) # changer !!
  
  #Survival
  if (!is.null(Survdata)) {
    cov_surv <- unique(unlist(strsplit(fixed.survival.models, "[+*]")))
    
    if (nE == 2) {
      tmat <- trans.comprisk(2, names = c("event-free", "event1", "event2"))
      
      Survdata$stat1 <- as.numeric(Survdata$StatusEvent == 1)
      Survdata$stat2 <- as.numeric(Survdata$StatusEvent == 2)
      Survdatalong <- msprep(
        time = c(NA, "Event", "Event"),
        status = c(NA, "stat1", "stat2"),
        data = Survdata,
        keep = cov_surv,
        trans = tmat
      )
      dim0 <- dim(Survdatalong)[2]
      Survdatalong <- expand.covs(Survdatalong, cov_surv)
      cov_survexpanded <- paste(names(Survdatalong)[(dim0 + 1):dim(Survdatalong)[2]], collapse =
                                  '+')
      form <- as.formula(
        paste(
          "Surv(time, status) ~ ",
          cov_survexpanded,
          "+ factor(trans)+ strata(trans)",
          sep = ""
        )
      )
      mod_surv <- survreg(form, data = Survdatalong, dist = "weibull")
      
      #   1/(survreg's scale)  =    rweibull shape a
      #   exp(survreg's intercept) = rweibull scale b
      # (a/b) (x/b)^(a-1) shape a, scale b
      para_basehaz <- c(
        1 / (mod_surv$coefficients[["(Intercept)"]]),
        exp(mod_surv$icoef[2]),
        1 / (mod_surv$coefficients[["(Intercept)"]] + mod_surv$coefficients[["factor(trans)2"]]),
        exp(mod_surv$icoef[3])
      )
      n1 <- 1 + np_surv[1] - 1 + ifelse(assoc %in% c(0, 1, 3, 4), 1, 2) *
        nD
      para_surv <- c(
        mod_surv$coefficients[2:(1 + np_surv[1] - 1)],
        rep(0, ifelse(assoc %in% c(0, 1, 3, 4), 1, 2) * nD),
        mod_surv$coefficients[(n1):(n1 - 1 + np_surv[1] -
                                      1)],
        rep(0, ifelse(assoc %in% c(0, 1, 3, 4), 1, 2) * nD)
      )
      
    } else if (nE == 1) {
      form <- as.formula(paste("Surv(Event, StatusEvent) ~ ", cov_surv, sep =
                                 ""))
      mod_surv <- survreg(form, data = Survdata, dist = "weibull")
      para_basehaz <- c(1 / (mod_surv$coefficients[["(Intercept)"]]), exp(mod_surv$icoef[2]))
      para_surv <- c(mod_surv$coefficients[2:(1 + np_surv - 1)], rep(0, ifelse(assoc %in%
                                                                                 c(0, 1, 3, 4), 1, 2) * nD))
    }
    p <- p + length(para_basehaz) + length(para_surv)
    np_baz <- length(para_basehaz) / nE
   
     if (!is.null(interactionY.survival.models)) {
      browser()
    }
  }
  
  vec_alpha_ij_m <- matrix(vec_alpha_ij,nrow=nD,byrow = T)
  alpha_D_m <- DparChol(nD,alpha_D)
  
  par_obj <- list(baseline_par=alpha_mu0,
                  mixed_par=alpha_mu,
                  random=alpha_D_m,
                  transition_par=vec_alpha_ij_m,
                  )
}

