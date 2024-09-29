#' @title DIDparams
#'
#' @description Object to hold DiD parameters that are passed across functions
#'
#' @inheritParams att_gt2
#' @inheritParams pre_process_did2
#' @param did_tensor list of outcome tensors that are used in the estimation
#' @param args list of arguments that are used in the estimation
#' @noRd
#' @export
DIDparams2 <- function(did_tensors, args, call=NULL) {
  # get the arguments from args
  yname <- args$yname
  tname <- args$tname
  idname <- args$idname
  gname <- args$gname
  xformla <- args$xformla # formula of covariates
  panel <- args$panel
  est_method <- args$est_method
  bstrap <- args$bstrap
  biters <- args$biters
  cband <- args$cband
  anticipation <- args$anticipation
  control_group <- args$control_group
  allow_unbalanced_panel <- args$allow_unbalanced_panel
  weightsname <- args$weightsname
  base_period <- args$base_period
  clustervars <- args$clustervars
  cores <- args$cores
  pl <- args$pl
  print_details <- args$print_details
  alp <- args$alp
  true_repeated_cross_sections <- args$true_repeated_cross_sections
  time_periods_count <- args$time_periods_count
  time_periods <- args$time_periods
  treated_groups_count <- args$treated_groups_count
  treated_groups <- args$treated_groups
  id_count <- args$id_count

  # get the arguments from did_tensors
  outcomes_tensor <- did_tensors$outcomes_tensor
  time_invariant_data <- did_tensors$time_invariant_data
  cohort_counts <- did_tensors$cohort_counts
  covariates <- did_tensors$covariates # matrix of covariates
  cluster_vector <- did_tensors$cluster
  weights_vector <- did_tensors$weights


  out <- list(yname=yname,
              tname=tname,
              idname=idname,
              gname=gname,
              xformla=xformla,
              panel=panel,
              est_method=est_method,
              bstrap=bstrap,
              biters=biters,
              cband=cband,
              anticipation=anticipation,
              control_group=control_group,
              allow_unbalanced_panel=allow_unbalanced_panel,
              weightsname=weightsname,
              base_period=base_period,
              clustervars=clustervars,
              cores=cores,
              pl = pl,
              print_details=print_details,
              alp=alp,
              true_repeated_cross_sections=true_repeated_cross_sections,
              time_periods_count=time_periods_count,
              time_periods=time_periods,
              treated_groups_count=treated_groups_count,
              treated_groups=treated_groups,
              id_count=id_count,
              outcomes_tensor=outcomes_tensor,
              time_invariant_data=time_invariant_data,
              cohort_counts=cohort_counts,
              covariates=covariates,
              cluster_vector=cluster_vector,
              weights_vector=weights_vector,
              call=call)
  class(out) <- "DIDparams"
  return(out)
}
