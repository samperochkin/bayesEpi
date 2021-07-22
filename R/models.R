#' Constructor for ccModel (S3) objects.
#' @param . TO DO
#' @return . TO DO
#' @examples
#' TO DO
#' @export
ccModel <- function(response, time_index, fixed, random, overdispersion, design, control_aghq){
  model <- mget(names(formals()),sys.frame(sys.nframe()))
  attr(model, "class") <- "ccModel"
  model
}


#' Specifies the case-crossover design to be used by fitModel.
#'
#' @param model_class Either "cc" (case crossover) or "ts" (time series).
#' @param ... See details.
#' @return A list of parameters for use by aghq::aghq.
#' @examples
#' control_aghq(model_class = "cc", k = 3)
#' control_aghq(model_class = "ts")
#' @export
ccDesign <- function(...){

  params = list(...)

  # default cc parameters
  design <- list(scheme = "bidirectional",
                 n_control = 2)

  design[names(params)] <- params
  return(design)
}


#' Constructor for tsModel (S3) objects.
#' @param . TO DO
#' @return . TO DO
#' @examples
#' TO DO
#' @export
tsModel <- function(response, time_index, fixed, random, overdispersion, control_aghq){
  model <- as.list(match.call()[-1])
  attr(model, "class") <- "tsModel"
  model
}


#' Specify a nonlinear random effect.
#'
#' @param variable_name
#' @param prior_dist A
#' @examples
#' @export
# fixedEffect <- function(variable, prior = gaussian_prior(prec = .01)){
#   list(variable = variable, prior = prior)
# }
fixedEffect <- function(prior = gaussian_prior(prec = .01)){
  list(prior = prior)
}


#' Specify a nonlinear random effect.
#'
#' @param variable_name
#' @param prior_dist A
#' @examples
#' @export
# randomEffect <- function(variable, model, prior){
#   list(variable = variable, model = model, prior = prior)
# }
randomEffect <- function(model, prior, reference_value){
  list(model = model, prior = prior)
}


#' Sets aghq parameters for fitting a model of model_type.
#'
#' @param model_class Either "cc" (case crossover) or "ts" (time series).
#' @param ... See details.
#' @return A list of parameters for use by aghq::aghq.
#' @examples
#' control_aghq(model_class = "cc", k = 3)
#' control_aghq(model_class = "ts")
#' @export
controlAGHQ <- function(...){

  params = list(...)

  # default cc parameters
  control <- list(k = 3)

  control[names(params)] <- params
  return(control)
}


#' Fits a model (either a ccModel or a tsModel)
#'
#' @param model A model, either a ccModel or a tsModel object.
#' @param data
#' @return An updated ccModel or tsModel object that contains the fitted model (e.g., in ccModel$fit).
#' @examples
#' @importFrom magrittr %>%
fitModel <- function(model, data){
  UseMethod("fitModel", model)
}

fitModel.ccModel <- function(model, data){
  # name local environment and add the parameters in parameter_list to it
  ccModel_env <- environment()
  list2env(model_parameters, envir = ccModel_env)

  data <- as.data.frame(data)

  # identify fixed and random effects
  fixed_names <- names(model$fixed)
  random_names <- names(model$random)

  # create response vector and design matrices
  y <- unlist(data[model$response], use.names = F)
  X <- as.matrix(data[fixed_names])
  U <- as.matrix(data[random_names])

  # discretize nonlinear random effects
  As <- discretizedRandomDesigns(model$random, U)
  # add polynomial interpolation around reference values (for discretized random effects) as fixed effects
  Xs_int <- interpolationFixedEffects(model$random, U)


  # Model fit ---------------------------------------------------------------



  # set parameters (beta, gamma, z, theta) and hyperparameters (u, alpha, beta_prec)
  source("R/models/case-crossover-AGHQ-dev/scripts/setup-parameters.R", local = T)

  # construct objective functions with derivatives (TMB)
  source("R/models/case-crossover-AGHQ-dev/functions/constructObjective.R", local = T)
  ff <- constructObjective(design, model_parameters, param_init)

  # compute quadrature
  quad <- aghq::marginal_laplace_tmb(ff,3,theta_init)


  list(quad = quad, ff = ff, design = design, model_parameters = model_parameters)
}




getRoundedRefValue <- function(random_model, u){

    if(is.function(random_model$ref_value)){
      if(missing(u)) stop("A function was specified as 'ref_value' but no data 'u' was provided to 'getRefValue' for computing it")
      ref_value <- random_model$ref_value(U[,name])
    }else{
      ref_value <- random_model$ref_value
    }

    round(ref_value/random_model$binwidth)*random_model$binwidth
}

discretizedRandomDesigns <- function(random, U){

  random_names <- names(random)

  sapply(random_names, function(name){

    if("binwidth" %in% names(random[[name]]$model)){
      random_model <- random[[name]]$model

      new_u <- round(U[,name]/random_model$binwidth)
      fac <- factor(new_u, levels = seq(min(new_u), max(new_u),1), labels = paste0(name,"__",seq(min(new_u), max(new_u),1)*random_model$binwidth))
      A <- Matrix::t(Matrix::fac2sparse(fac, drop.unused.levels = F))

      # Set reference value by setting corresponding column of A (and neighbours) to zero.
      rounded_ref_value <- getRoundedRefValue(random_model, U[,name])
      ref_value_col <- which(colnames(A) == paste0(name,"__",rounded_ref_value))
      removed_cols <- ref_value_col
      if(random_model$order >= 2) removed_cols <- c(removed_cols,ref_value_col+1)
      if(random_model$order >= 3) removed_cols <- c(ref_value_col-1,removed_cols)
      if(random_model$order > 3) stop("Random walks of order > 3, which you specified for ", name, " are not yet implemented.")
      A[, removed_cols] <- 0

      attr(A, "rounded_ref_value") <- rounded_ref_value
      attr(A, "removed_cols") <- removed_cols

      return(A)

    }else{
      stop("You specified a model for the random effect ", name, " that's not yet implemented. Only general Gaussian random walk for general (discretized) nonlinear effect are implemented... ")
      return(U[,name])
    }

  }, simplify = FALSE, USE.NAMES = TRUE)
}
#


interpolationFixedEffects <-  function(random, U){

  random_names <- names(random)
  sapply(random_names, function(name){

    if("poly_degree" %in% names(random[[name]]$model)){

      random_model <- random[[name]]$model
      poly_degree <- random_model$poly_degree
      if(poly_degree %in% c(0,1)){
        cat("poly_degree for the random walk ", name, " was set to ", poly_degree,". Such poly_degree means no interpolation.")
        return(NULL)
      }
      rounded_ref_value <- attr(As[[name]], "rounded_ref_value")

      # SHOULD WE NORMALIZE THE DATA HERE?
      X_new <- poly(U[,name] - rounded_ref_value, degree = random_model$poly_degree, raw = TRUE)
      colnames(X_new) <- paste0(name,"__", attr(X_new, "degree"))
      return(X_new)

    }else{
      return(NULL)
    }
  }, simplify = FALSE, USE.NAMES = TRUE)
}
#
