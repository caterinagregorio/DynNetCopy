#' Helper function to give initial model parameters
#'
#' @param structural.model a list of 7 arguments used to specify the structural model: \describe{
#'    \item{\code{fixed.LP0}}{a one-sided linear formula object for specifying the 
#'    fixed effects in the submodel for the baseline level of processes. Note that 
#'    there is no need to specify a random effect model for the baseline level of 
#'    processes as we systematically set a process-specific random intercept (with 
#'    variance fixed to 1 for identifiability purpose). For identifiability purposes, 
#'    the mean intercepts are fixed to 0 (not estimated).}
#' 
#'    \item{\code{fixed.DeltaLP}}{a two-sided linear formula object for 
#'    specifying the response outcomes (one the left part of ~ symbol) 
#'    and the covariates with fixed-effects (on the right part of ~ symbol) 
#'    in the submodel for change over time of latent processes.}
#' 
#'    \item{\code{random.DeltaLP}}{a one-sided linear formula object 
#'    for specifying the random effects in the submodel for change 
#'    over time of latent processes.}
#' 
#'    \item{\code{trans.matrix}}{a one-sided linear formula object 
#'    for specifying a model for elements of the transition matrix, 
#'    which captures the temporal influences between latent processes.}
#' 
#'    \item{\code{fixed.survival}}{a one-sided linear formula object 
#'    for specifying the covariates in the survival sub-model. In competing 
#'    risks model, the specification for the two events should be separated 
#'    by the "|" symbol.}
#' 
#'    \item{\code{interactionY.survival}}{a one-sided linear formula 
#'    object for specifying the covariates in interaction with the 
#'    dynamics of the latent processes, in the survival sub-model. 
#'    In competing risks model, the specification for the two events 
#'    should be separated by the "|" symbol. Only additional terms 
#'    should be included (No "*" symbol). Covariates in 
#'    \code{interactionY.survival} should also be included in \code{fixed.survival.}}
#' 
#'    \item{\code{delta.time}}{indicates the discretisation step to be used for latent processes}
#'    }
#' 
#' @param measurement.model is a list of arguments detailed below used to specify the measurement model: \describe{
#' 
#' \item{\code{link}}{ indicates the link functions to be used to transform the outcomes. It takes values in "linear" for a linear transformation and "n-type-d" for a I-splines transformation where "n" indicates the number of nodes, "type" (which takes values in "quant", "manual", "equi") indicates where the nodes are placed, and "d" indicates the degree of the I-splines.}
#' 
#' \item{\code{knots}}{ argument indicates if necessary the place of knots (when placed manually with "manual"), default value is NULL}
#' }
#' 
#' @param parameters a list of 3 arguments about parameters of the models (e.g., initial parameters, parameters one would like to fix, etc.) possibly created with \link{helperDynNetpar}: \describe{
#' \item{\code{paras.ini}}{ indicates initial values for parameters, default values is NULL.}
#' \item{\code{Fixed.para.indix}}{ indicates the positions of parameters to be constrained.}
#' 
#' \item{\code{Fixed.para.values}}{ indicates the values associated to the index of parameters to be constrained. }
#' }
#'   \item{\code{assocT}}{ Specifies the type of association between the time-to-event(s) and the latent process(es). Values include "r.intercept" for random intercept, "r.slope" for random slope, "r.intercept/slope" for both random intercept and slope, "c.value" for current value. }
#' @param Time indicates the name of the covariate representing the time 
#' @param subject indicates the name of the covariate representing the grouping structure
#' @param data indicates the data frame containing all the variables for estimating the model.
#' @param cholesky logical indicating if the variance covariance matrix is parameterized using the cholesky (TRUE by default) or the correlation (FALSE)
#' @param Tentry name of the variable of entry time
#' @param Event name of the variable of event time
#' @param StatusEvent name of the variable of event status
#' @param basehaz type of baseline hazard function
#'
#' @returns object of class DynNetinit
#' @export
#'
helperDynNetpar <- function(structural.model, measurement.model, Time, Tentry ="Tentry", Event = "Event", StatusEvent = "StatusEvent" ,basehaz = NULL, assocT=NULL,subject, data){
 

  ### check if all component of the model specification are well filled ####
  if(missing(structural.model))stop("The argument structural.model must be specified")
  if(missing(measurement.model))stop("The argument measurement.model must be specified")
  if(missing(subject))stop("The argument subject must be specified")
  if(missing(Time))stop("The argument time must be specified")
  if(missing(data))stop("The argument data must be specified")
  
  if(is.null(structural.model$fixed.DeltaLP))stop("The argument structural.model$fixed.DeltaLP must be specified")
  #if(is.null(structural.model$random.DeltaLP))stop("The argument structural.model$random.DeltaLP must be specified")
  if(is.null(structural.model$trans.matrix))stop("The argument structural.model$trans.matrix must be specified")
  if(is.null(structural.model$delta.time)){
    structural.model$delta.time <- 1
  }
  
  if(is.null(measurement.model$link.functions) || all(is.null(measurement.model$link.functions$links))){
    links <- NULL
    knots <- NULL
    measurement.model$link.functions =list(links = links, knots = knots)
  }else{
    if(all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="quant", x))) ==0)&& 
       all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="manual", x))) ==0)&& 
       all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="equi", x))) ==0)&& 
       all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="linear", x))) ==0)&& 
       all(sapply(measurement.model$link.functions$links, function(x) length(grep(pattern="thresholds", x))) ==0))
      stop("The only available link functions are 'linear', 'splines' and 'thresholds' functions.")
  }
  
  
  survival= FALSE
  if(!is.null(structural.model$fixed.survival)){
    survival = TRUE
    fixed.survival <- structural.model$fixed.survival
  }
  
  ### identification of model components #####
  ## components of structural model
  fixed_X0 <- structural.model$fixed.LP0
  fixed_DeltaX <- structural.model$fixed.DeltaLP
  randoms_DeltaX <- structural.model$random.DeltaLP
  mod_trans <- structural.model$trans.matrix
  DeltaT <- structural.model$delta.time
  
  
  interactionY.survival <- NULL
  if(!is.null(structural.model$interactionY.survival)){
    interactionY.survival <- structural.model$interactionY.survival
  }
  
  # components of measurement model
  link <- measurement.model$link.functions$links
  knots <- measurement.model$link.functions$knots
  
  
  colnames<-colnames(data)
  # if(missing(DeltaT) || DeltaT < 0 ) stop("The discretization step DeltaT cannot be null or negative")
  if(!(subject%in%colnames))stop("Subject should be in the dataset")
  if(!(Time %in% colnames)) stop("Time variable should be indicated and should be in the dataset")
  if(dim(unique(data))[1] != dim(data)[1]) stop("Some rows are the same in the dataset, perhaps because of a too large discretisation step")
  
  
  data <- data[order(data[,subject], data[,Time]),]
  
  #### fixed effects pre-traitement ####
  
  ### for DeltaLP
  # if(missing(fixed_DeltaX)) stop("The argument fixed_DeltaX must be specified in any model")
  if(!inherits(fixed_DeltaX,"formula")) stop("The argument fixed_DeltaX must be a formula")
  
  ### outcomes and latent processes ####
  outcome <- as.character(attr(terms(fixed_DeltaX),"variables"))[2]
  outcomes_by_LP<-strsplit(outcome,"[|]")[[1]]
  if(any(grepl("()",outcomes_by_LP))){
    formative <- TRUE
  }
  nD <- length(outcomes_by_LP) # nD: number of latent process
  
  outcomes <- NULL
  mapping.to.LP <- NULL
  mapping.to.LP2 <- NULL
  for(n in 1:nD){
    outcomes_n <- strsplit(outcomes_by_LP[n],"[+]")[[1]]
    outcomes_n <-as.character(sapply(outcomes_n,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE))
    outcomes_n <- unique(outcomes_n)
    if(is.null(outcomes_n)) stop("at least one marker must be specified for each latent process" )
    outcomes <- c(outcomes, outcomes_n)
    mapping.to.LP <- c(mapping.to.LP, rep(n,length(outcomes_n)))
  }
  
  if(formative){
    outcomes <- as.character(sapply(outcomes,FUN = function(x)gsub("[()+]","",x),simplify = FALSE))
    for (n in 1:nD){
      outcomes_n <-  regmatches(outcomes_by_LP[n], gregexpr("\\([^()]*\\)|Y\\d+", outcomes_by_LP[n]))[[1]]
      outcomes_n <-as.character(sapply(outcomes_n,FUN = function(x)gsub("[()]","",x),simplify = FALSE))
      outcomes_n <-as.character(sapply(outcomes_n,FUN = function(x)gsub("[[:space:]]","",x),simplify = FALSE))
      outcomes_n <- unique(outcomes_n)
      for (l in 1:length(outcomes_n)){
        outcomes_n2 <- regmatches(outcomes_n[l], gregexpr("Y\\d+", outcomes_n[l]))[[1]]
        if(is.null(outcomes_n2)) stop("at least one marker must be specified for each latent exogeneous process" )
        add <- ifelse(n==1,0,max(mapping.to.LP2))
        mapping.to.LP2 <- c(mapping.to.LP2, rep(l+add,length(outcomes_n2)))
      }
      
    }
  }
  if(length(unique(outcomes))!=length(outcomes)) stop("outcomes cannot be mapped to multiple latent processes")
  if(!all(outcomes%in% colnames)) stop("outcomes must be in the data")
  nL <- ifelse(formative==T,max(mapping.to.LP2),NULL)
  K <- length(outcomes)
  all.Y<-seq(1,K)
  
  fixed_DeltaX.model=strsplit(gsub("[[:space:]]","",as.character(fixed_DeltaX)),"~")[[3]]
  fixed_DeltaX.models<-strsplit(fixed_DeltaX.model,"[|]")[[1]]# chaque model d'effet fixe mais en vu de connaitre tous les pred.fixed du modele multi
  
  if(nD !=length(fixed_DeltaX.models)) stop("The number of models does not correspond to the number of latent processes")
  
  if(formative){
    if(nL > K){
      stop("There are too many latent processes compared to the indicated number of markers")
    }
  }else{
    if(nD > K){
      stop("There are too many latent processes compared to the indicated number of markers")
    }
  }
  
  
  ### pre-traitement of fixed effect on initial levels of processes
  if(is.null(fixed_X0)){
    fixed_X0<- ~1
    fixed_X0.models <- rep("1",nD)
  }
  if(!inherits(fixed_X0,"formula")) stop("The argument fixed_X0 must be a formula")
  
  fixed_X0.models =strsplit(gsub("[[:space:]]","",as.character(fixed_X0)),"~")[[2]]
  fixed_X0.models<- as.vector(strsplit(fixed_X0.models,"[|]")[[1]]) 
  for(nd in 1:nD){
    if(fixed_X0.models[nd]=="~-1") fixed_X0.models[nd] <-"~1" # au moins l'intcpt
  }
  
  ### pre-traitement of fixed effect on survival 
  fixed.survival.models <- NULL
  truncation = FALSE
  assocT <- 0
  if(survival){
    if(is.null(fixed.survival)){
      fixed.survival<- ~1
    }else{
      if(!inherits(fixed.survival,"formula")) stop("The argument fixed.survival must be a formula") 
    }
    
    fixed.survival.models <- strsplit(gsub("[[:space:]]","",as.character(fixed.survival)),"~")[[2]]
    covsurv <- unique(as.vector(strsplit(fixed.survival.models,"[|*+]")[[1]]))
    fixed.survival.models <- as.vector(strsplit(fixed.survival.models,"[|]")[[1]]) 
    if(!all(covsurv[-which(covsurv=="1")]%in%colnames))stop("All covariates in fixed.survival should be in the dataset")
    
    if(!is.null(assocT)){
      if(!assocT%in%c("r.intercept", "r.slope", "r.intercept/slope", "c.value"))
        stop("assocT should be defined as r.intercept, r.slope, r.intercept/slope, c.value.")
      assocT <- switch(assocT, "r.intercept"=0, "r.slope"=1, "r.intercept/slope"=2, "c.value"=3, "c.slope"=4, "c.value/slope"=5)
    }else{
      assocT <- "r.intercept/slope"
      assocT <- 2 # random intercept and slope
    }
    
    if(assocT <= 2){
      message(" add interactions ui * X in survival model")
      
    }
    if(!is.null(option$truncation)){
      truncation <- option$truncation
    }
  }
  
  ### pre-traitement of interactions with Y on survival 
  
  interactionY.survival.models <- NULL
  if(survival){
    
    if(!is.null(interactionY.survival)){
      if(!inherits(interactionY.survival,"formula")) stop("The argument interactionY.survival must be a formula") 
      interactionY.survival.models <- strsplit(gsub("[[:space:]]","",as.character(interactionY.survival)),"~")[[2]]
      intYsurv <- (as.vector(strsplit(interactionY.survival.models,"[|*+]")[[1]]))
      interactionY.survival.models <- as.vector(strsplit(interactionY.survival.models,"[|]")[[1]]) 
      if(!all(intYsurv%in%colnames))stop("All covariates in interactionY.survival should be in the dataset")
      if(any(grepl( "*", interactionY.survival, fixed = TRUE)))stop("Only + terms should be included in interactionY.survival, no *.")
    }
  }
  
  
  ### pre-traitement of random effect on processes  intercept and slope
  #### randoms effet on DeltaLP 
  randoms_X0.models <- rep("1",nD)
  #### randoms effet on DeltaX
  
  if(missing(randoms_DeltaX) || is.null(randoms_DeltaX)){
    randoms_DeltaX<- ~1
    randoms_DeltaX.models <- rep("1",nD)
  }
  if(!inherits(randoms_DeltaX,"formula")) stop("The argument random must be a formula")
  randoms_DeltaX.model=strsplit(gsub("[[:space:]]","",as.character(randoms_DeltaX)),"~")[[2]]
  randoms_DeltaX.models<-strsplit(randoms_DeltaX.model,"[|]")[[1]]    
  
  #### traitement of  mod_trans: transition matrix##### 
  if(missing(mod_trans)){
    mod_trans <- ~ 1 # constant transition matrix
  } 
  if(!inherits(mod_trans,"formula")) stop("The argument mod_trans must be a formula")
  mod_trans.model=strsplit(gsub("[[:space:]]","",as.character(mod_trans)),"~")[[2]]
  
  if(nD!=length(fixed_X0.models)){
    stop("The number of models for initial latent processes does not correspond with the number of latent processes")
  }
  if(nD!=length(fixed_DeltaX.models)){
    stop("The number of models for the change over time of latent processes does not correspond with the number of latent processes")
  }
  
  ### traitement of transformation models ##
  if(is.null(link)){
    link <- rep("linear",K)
  }
  else if(length(link)!=K) stop("The number transformation links must be equal to the number of markers")
  
  if(any(link == "thresholds")){
    j     <- which(link == 'thresholds')
    
    Y0    <- data[,outcomes[j]]
    minY0 <- apply(as.matrix(Y0), 2, min, na.rm=TRUE)# min(Y0,rm.na=TRUE)
    maxY0 <- apply(as.matrix(Y0), 2, max, na.rm=TRUE)# max(Y0,rm.na=TRUE)
    
    if(length(j==1))
      ide0 <- matrix(0, nrow = length(j),  ncol = max(sapply(j, function(x) max(Y0, na.rm =T)-min(Y0, na.rm =T) ))) #change dimensions
    else
      ide0 <- matrix(0, nrow = length(j),  ncol = max(sapply(j, function(x) max(Y0[,x], na.rm =T)-min(Y0[,x], na.rm =T) ))) #change dimensions
    
    nbzitr0 <- 2
    #idlink0 <- 3
    #ntrtot0 <- as.integer(maxY0 - minY0)
    if (!(all(is.integer(minY0) | !all(is.integer(maxY0)))))
      stop("With the thresholds link function, the longitudinal outcome must be discrete")
    
    zitr <- c()
    for(i in 1:length(j)){
      if(length(j)>1){
        Y0tmp <- Y0[,i]
      }else{
        Y0tmp <- Y0
      }
      
      if (!all(Y0tmp[which(!is.na(Y0tmp))] %in% minY0[i]:maxY0[i]))
        stop("With the threshold link function, problem with the outcome data, must be discrete")
      
      IND <- sort(unique(Y0tmp))
      IND <- IND[1:(length(IND) - 1)] - minY0[i] + 1
      #ide0 <- rep(0, as.integer(maxY0[i] - minY0[i])) #change dimensions
      ide0[i, IND] <- 1
      
      #if(i==1)
      #  zitr <- matrix(rep(0, length(j)*nbzitr0), length(j), nbzitr0)
      #zitr[i, nbzitr0] <- maxY0[i]
      zitr <- c(zitr, minY0[i], maxY0[i])
    }
    
    if(!is.null(type_int)){
      if(!type_int %in% c("MC", "sobol", "halton", "torus"))
        stop("With the thresholds link function, type_int should be either antithetic, sobol, halton or torus. antithetic not developed yet, sorry.")
    }
    
  }else{
    zitr <- 0
    ide0  <- 0 
  }
  
  
  sequence  <- 0
  ind_seq_i <- 0
  
  unique_id <- unique(data[,subject])
  
  nmes <- sapply(unique_id, function(x) 
    length(which(!is.na(data[which(data[,subject]==x),outcomes]))))
  
  Survdata <- NULL
  knots_surv <- NULL
  ## If joint model -  Event and StatusEvent data
  if(survival){
    if(!(Event%in%colnames)) stop("Event should be in the dataset")
    if(!(StatusEvent%in%colnames)) stop("StatusEvent should be in the dataset")
    
    basehaz <- ifelse(!is.null(option$basehaz), option$basehaz, "Weibull")
    if(Tentry != "Tentry" & !(Tentry%in%colnames)) stop("Tentry should be in the dataset")
    
    if(!basehaz%in%c("Weibull",'splines'))
      stop('basehaz should be either Weibull or splines.')
    
    if(basehaz == "Splines")
      cat("Define knots_surv here")
    #One survival dataframe with one line per individual
    first_line <- sapply(unique(data[,subject]), function(x) which(data[,subject]==x)[1])
    
    if(!(Tentry %in%names(data))) data$Tentry <- 0
    Survdata <- data[first_line, c(Tentry, Event, StatusEvent)]
    names(Survdata) <- c("Tentry", "Event", "StatusEvent")
    
    if(length(which(covsurv!="1"))>0){
      Survdata <- cbind(Survdata, data[first_line, covsurv])
      names(Survdata)[(dim(Survdata)[2]-length(covsurv)+1):dim(Survdata)[2]]<-covsurv
    }
  }

  formula <- list(fixed_X0.models=fixed_X0.models, fixed_DeltaX.models = fixed_DeltaX.models, 
                      randoms_X0.models=randoms_X0.models, randoms_DeltaX.models=randoms_DeltaX.models, 
                      mod_trans.model = mod_trans.model)
  
  
  data_F <- DataFormat(data=data, subject = subject, fixed_X0.models = fixed_X0.models,
                       randoms_X0.models = randoms_X0.models, fixed_DeltaX.models = fixed_DeltaX.models, 
                       randoms_DeltaX.models = randoms_DeltaX.models, mod_trans.model = mod_trans.model, 
                       outcomes = outcomes, nD = nD, link=link, knots = knots, zitr= zitr, ide =  ide0, 
                       Time = Time, Survdata = Survdata, basehaz = basehaz, fixed.survival.models =fixed.survival.models, 
                       interactionY.survival.models = interactionY.survival.models, DeltaT=DeltaT, assoc = assoc, truncation = truncation)
  
  
  K <- data_F$K #  number of markers
  vec_ncol_x0n <- data_F$vec_ncol_x0n # number of parameters on initial level of processes
  n_col_x <- ncol(data_F$x) # number of parameters on processes slope
  nb_RE <- data_F$nb_RE # number of random effects on the processes slope
  L <- ncol(data_F$modA_mat)
  ncolMod.MatrixY <- ncol(data_F$Mod.MatrixY)
  
  assocT <- NULL
  if(!is.null(assoc)){
    assocT <- ifelse(assoc==0, "r.intercept",ifelse(assoc==1, "r.slope",ifelse(assoc==2, "r.intercept/slope",ifelse(
      assoc==3, "c.value",ifelse(assoc==4, "c.slope","c.value/slope")
    ))))
  }
  
  par_obj <- parskeleton(
    K = K,
    nD = nD,
    mapping.to.LP = mapping.to.LP,
    nL = nL,
    mapping.to.LP2 = mapping.to.LP2,
    vec_ncol_x0n,
    n_col_x,
    nb_RE,
    indexparaFixeUser = indexparaFixeUser,
    paraFixeUser = paraFixeUser,
    L = L,
    ncolMod.MatrixY = ncolMod.MatrixY,
    paras.ini = paras.ini,
    link = link,
    npara_k = npara_k,
    Survdata = Survdata,
    basehaz = basehaz,
    knots_surv = knots_surv,
    assoc = assoc,
    truncation = truncation,
    data = data,
    outcomes = outcomes,
    df = data_F$df,
    nE = data_F$nE,
    np_surv = data_F$np_surv,
    fixed.survival.models = fixed.survival.models,
    interactionY.survival.models = interactionY.survival.models,
    nYsurv = data_F$nYsurv
  )
  class(par_obj) <- "DynNetinit"
return(par_obj)
  }