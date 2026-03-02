## -----------------------------------------------------------------------------

#' @title honest_did
#'
#' @description a function to compute a sensitivity analysis
#'  using the approach of Rambachan and Roth (2021)
#' @param object an event study
#' @export
honest_did <- function(object, ...) {
  UseMethod("honest_did", object)
}


#' @title honest_did.AGGTEobj
#'
#' @description a function to compute a sensitivity analysis
#'  using the approach of Rambachan and Roth (2021) when
#'  the event study is estimating using the `did` package
#'
#' @param e_time event time to compute the sensitivity analysis for.
#'  The default value is `e_time=0` corresponding to the "on impact"
#'  effect of participating in the treatment.
#' @param type Options are "smoothness" (which conducts a
#'  sensitivity analysis allowing for violations of linear trends
#'  in pre-treatment periods) or "relative_magnitude" (which
#'  conducts a sensitivity analysis based on the relative magnitudes
#'  of deviations from parallel trends in pre-treatment periods).
#' @inheritParams HonestDiD::createSensitivityResults
#' @inheritParams HonestDiD::createSensitivityResults_relativeMagnitudes
#' 
#' @export
honest_did.AGGTEobj <- function(object,
                                e_time=0,
                                type=c("smoothness", "relative_magnitude"),
                                method=NULL,
                                bound="deviation from parallel trends",
                                Mvec=NULL,
                                Mbarvec=NULL,
                                monotonicityDirection=NULL,
                                biasDirection=NULL,
                                alpha=0.05,
                                parallel=FALSE,
                                gridPoints=10^3,
                                grid.ub=NA,
                                grid.lb=NA,
                                ...) {


  type <- type[1]
  
  # make sure that user is passing in an event study
  if (object$type != "dynamic") {
    stop("'object' must be an AGGTEobj of type 'dynamic' (an event study). Use `aggte(., type = 'dynamic')` to create one.")
  }

  # check if used universal base period and warn otherwise
  if (object$DIDparams$base_period != "universal") {
    warning("It is recommended to use a universal base period for honest_did(). Set `base_period = 'universal'` in your call to att_gt().")
  }

  # recover influence function for event study estimates
  es_inf_func <- object$inf.function$dynamic.inf.func.e

  # recover variance-covariance matrix
  n <- nrow(es_inf_func)
  V <- t(es_inf_func) %*% es_inf_func / (n*n) 

  #Remove the coefficient normalized to zero
  referencePeriodIndex <- which(object$egt == -1)
  V <- V[-referencePeriodIndex,-referencePeriodIndex]
  beta <- object$att.egt[-referencePeriodIndex]
  
  nperiods <- nrow(V)
  npre <- sum(1*(object$egt < 0))
  npost <- nperiods - npre

  baseVec1 <- basisVector(index=(e_time+1),size=npost)

  orig_ci <- constructOriginalCS(betahat = object$att.egt,
                                 sigma = V, numPrePeriods = npre,
                                 numPostPeriods = npost,
                                 l_vec = baseVec1)

  if (type=="relative_magnitude") {
    if (is.null(method)) method <- "C-LF"
    robust_ci <- createSensitivityResults_relativeMagnitudes(betahat = object$att.egt, sigma = V, 
                                                             numPrePeriods = npre, 
                                                             numPostPeriods = npost,
                                                             bound=bound,
                                                             method=method,
                                                             l_vec = baseVec1,
                                                             Mbarvec = Mbarvec,
                                                             monotonicityDirection=monotonicityDirection,
                                                             biasDirection=biasDirection,
                                                             alpha=alpha,
                                                             gridPoints=gridPoints,
                                                             grid.lb=grid.lb,
                                                             grid.ub=grid.ub,
                                                             parallel=parallel)
    
  } else if (type=="smoothness") {
    robust_ci <- createSensitivityResults(betahat = object$att.egt,
                                          sigma = V, 
                                          numPrePeriods = npre, 
                                          numPostPeriods = npost,
                                          method=method,
                                          l_vec = baseVec1,
                                          monotonicityDirection=monotonicityDirection,
                                          biasDirection=biasDirection,
                                          alpha=alpha,
                                          parallel=parallel)
  }

  list(robust_ci=robust_ci, orig_ci=orig_ci, type=type)
}
