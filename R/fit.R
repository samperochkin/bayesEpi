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
fitModel <- function(model, data, silent = F, params_init = NULL){
  UseMethod("fitModel", model)
}

#' @method fitModel ccModel
#' @export
#' @useDynLib bayesEpi
#' @importFrom magrittr %>%
#' @import OSplines
fitModel.ccModel <- function(model, data, silent = F, params_init = NULL){

  data <- as.data.frame(data)
  list2env(checkups(model, data), envir = environment())

  # contruct and import case_day and control_days. modifies data.
  list2env(getCaseControl(data, model), envir = environment())
  list2env(setRefValues(data, model), envir = environment())

  # overdispersion and design matrices
  overdispersion <- !is.null(model$overdispersion)
  X <- as.matrix(data[names(model$fixed)])
  U <- as.matrix(data[names(model$random)])

  # bluids design matrices for random effects, polynomial interpolation around reference values,
  # creates As, Xs_int and gamma_dims
  list2env(createRandomDesigns(model, U), envir = environment())

  # prior parameters
  prior_lookup <- c("pc_prec", "gamma", "log-gamma")
  beta_prec = c(purrr::map(model$fixed, ~ .x$prior$params$prec), purrr::map(model$random, ~ .x$beta_prior$params$prec)) %>% unlist
  theta_prior <- c(purrr::map(model$random, ~ .x$theta_prior$type), z = model$overdispersion$theta_prior$type) %>% unlist
  theta_prior_id = match(theta_prior , prior_lookup)
  theta_hypers = c(purrr::map(model$random, ~ .x$theta_prior$params), z = model$overdispersion$theta_prior$params) %>% unlist
  if(is.null(theta_hypers)) theta_hypers <- numeric(0)


  #############################
  if(length(unique(sapply(model$random, function(ran) ran$model$type))) > 1)
    stop("Random effects must all be of the same type (either 'random walk' or 'integrated Wiener process'). Mixing to be implemented.")
  #############################


  # Model fit ---------------------------------------------------------------
  tmb_data <- list(count = data[case_day, model$response],
                   case_day = case_day, control_days = control_days,
                   X = cbind(X,Reduce("cbind", Xs_int)), A = Reduce("cbind", As),
                   random_effect_id = match(model$random[[1]]$model$type,
                                            c("random walk", "integrated Wiener process")),
                   Q_rw = constructQ_rw(model$random), Q_iwp = constructQ_iwp(model$random),
                   gamma_dims = gamma_dims, beta_prec = beta_prec,
                   theta_prior_id = theta_prior_id, theta_hypers = theta_hypers)

  theta_init <- getPriorInit(model)
  parameters <- list(beta = rep(0,ncol(X)+sum(sapply(Xs_int, ncol))),
                     gamma = rep(0, sum(sapply(As, ncol))),
                     z = rep(0,nrow(X)*overdispersion),
                     theta = theta_init)

  for(param_name in intersect(names(params_init), names(parameters))){
    if(param_name == "theta") theta_init <- params_init$theta
    parameters[[param_name]] <- params_init[[param_name]]
  }

  # dll <- ifelse(model$random[[1]]$model$type == "random walk", "cc_rw", "cc_iwp")
  dll <- "bayesEpi"
  obj <- TMB::MakeADFun(tmb_data, parameters, random = c("beta","gamma","z"), DLL=dll, hessian=T, silent = silent)
  if(silent){
    capture.output(quad <- aghq::marginal_laplace_tmb(ff = obj,
                                                      k = model$aghq_input$k,
                                                      startingvalue =  theta_init,
                                                      control = model$aghq_input$control))

  }else{
    quad <- aghq::marginal_laplace_tmb(ff = obj, k = model$aghq_input$k,
                                       startingvalue =  theta_init, control = model$aghq_input$control)
  }

  list(quad = quad, obj = obj, model = model, U = U)
}
