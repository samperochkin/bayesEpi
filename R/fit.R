#' Fit a model
#'
#' Function used for fitting a ccModel (case crossover model).
#' See the vignettes for examples of its usage.
#'
#' @param model A ccModel object.
#' @param data The dataset on which the model is fitted. It should contain all the variables involved in the model.
#' @param silent A boolean indicating whether you wish to silence the `aghq::marginal_laplace_tmb()` function (default is `TRUE`).
#' @return An updated ccModel object that contains the fitted model (e.g., in ccModel$fit).
#' @export
fitModel <- function(model, data, silent = F, dll = NULL){
  UseMethod("fitModel", model)
}

#' @method fitModel ccModel
#' @export
#' @useDynLib cc
#' @import TMB
#' @importFrom magrittr %>%
fitModel.ccModel <- function(model, data, silent = F, dll = NULL){

  data <- as.data.frame(data)
  list2env(checkups(model, data), envir = environment())

  # contruct and import case_day and control_days. modifies data.
  list2env(getCaseControl(data, model), envir = environment())
  list2env(setRefValues(data, model), envir = environment())

  if(any(is.na(data[c(names(model$random),names(model$fixed))]))) stop("Problem in preprocessing?")

  # overdispersion and design matrices
  overdispersion <- !is.null(model$overdispersion)
  X <- as.matrix(data[names(model$fixed)])
  U <- as.matrix(data[names(model$random)])

  # discretizes nonlinear random effects, builds polynomial interpolation around reference values,
  # creates As, Xs_int and gamma_dims
  list2env(discretizedRandomDesigns(model, U), envir = environment())

  # prior parameters
  prior_lookup <- c("pc_prec", "gamma")
  beta_prec = c(purrr::map(model$fixed, ~ .x$prior$params$prec), purrr::map(model$random, ~ .x$beta_prior$params$prec)) %>% unlist
  theta_prior <- c(purrr::map(model$random, ~ .x$theta_prior$type), z = model$overdispersion$theta_prior$type) %>% unlist
  theta_prior_id = match(theta_prior , prior_lookup)
  theta_hypers = c(purrr::map(model$random, ~ .x$theta_prior$params), z = model$overdispersion$theta_prior$params) %>% unlist
  if(is.null(theta_hypers)) theta_hypers <- numeric(0)

  # Model fit ---------------------------------------------------------------
  if(overdispersion){
    z_pos <- (1:nrow(data))[-selectFixedOD(data, model, case_day, control_days)]
  }else{
    z_pos <- integer(0)
  }
  tmb_data <- list(count = data[case_day, model$response],
                   case_day = case_day, control_days = control_days,
                   X = cbind(X,Reduce("cbind", Xs_int)), A = Reduce("cbind", As),
                   Q = constructQ(model$random), gamma_dims = gamma_dims,
                   beta_prec = beta_prec, theta_prior_id = theta_prior_id, theta_hypers = theta_hypers,
                   z_pos = z_pos)

  theta_init <- getPriorInit(model)

  parameters <- list(beta = rep(0,ncol(X)+sum(sapply(Xs_int, ncol))),
                     gamma = rep(0, sum(sapply(As, ncol))),
                     z = rep(0,length(z_pos)*overdispersion),
                     # z = rnorm(length(z_pos)*overdispersion, sd = sqrt(exp(-theta_init[length(theta_init)]))),
                     theta = theta_init)


  if(is.null(dll)) dll <- "bayesEpi"
  obj <- TMB::MakeADFun(tmb_data, parameters, random = c("beta","gamma","z"), DLL=dll, hessian=T)
  # obj <- TMB::MakeADFun(tmb_data, parameters, random = c("beta","gamma","z"), DLL=dll, hessian=T, silent=silent)

  if(length(theta_init) == 0){
    # messy code for cases with no random effects. For simulations purposes only; users should refrain from using it.
    invisible(obj$fn())
    opt_beta <- obj$env$last.par.best
    invisible(capture.output(obj <- TMB::MakeADFun(tmb_data, parameters, random = c("gamma","z","theta"), DLL=dll, hessian=T)))
    return(list(beta = opt_beta, sd = 1/sqrt(diag(obj$he(opt_beta, atomic = T))), model = model))
  }

  if(silent){
    invisible(capture.output(quad <- aghq::marginal_laplace_tmb(obj, model$control_aghq$k, theta_init)))
  }else{
    quad <- aghq::marginal_laplace_tmb(obj, model$control_aghq$k, theta_init)
  }

  list(quad = quad, obj = obj, model = model)
}
