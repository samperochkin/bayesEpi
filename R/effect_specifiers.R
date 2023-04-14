#' Specify a fixed effect.
#'
#' @param prior Prior distribution on the (beta) coefficient. Only the gaussian prior is implemented.
#' @examples
#' fixedEffect()
#' @export
fixedEffect <- function(model){
  if(model$type == "poly"){
    prior <- gaussian_prior(prec = rep(.01, model$params$degree))
  }else if(model$type == "bs"){
    prior <- gaussian_prior(prec = rep(.01, length(model$params$knots)+model$params$degree-2))
  }
  list(model = model, prior = prior)
}


#' Specify a nonlinear random effect.
#'
#' @param model A model for the random effect: `rw_effect()` or `iwp_effect()`.
#' @param theta_prior Prior distribution for the parameters theta of the random effect model.
#' @param beta_prior Prior distribution for the fixed effects created as a by-product. Only the Gaussian prior is implemented.
#' The number of beta_prior parameters must match the number of coefficients created by the chosen model, i.e., `order-1`.
#' For example, a random walk of order 3 creates (by default) 2 fixed effects. In most cases, this argument should be left unspecified.
#' @examples
#' randomEffect(model = rw_effect(), theta_prior = pc_prec())
#' @export
randomEffect <- function(model, theta_prior, beta_prior){
  if(!is.null(model$params$poly_degree) & missing(beta_prior)){
    beta_prior <- gaussian_prior(prec = rep(.01, model$params$poly_degree))
    list(model = model, theta_prior = theta_prior, beta_prior = beta_prior)
  }else if(!is.null(model$params$poly_degree)){
    # DO A CHECK FOR THE LENGTH OF BETA_PREC AND POLY_DEGREE
    list(model = model, theta_prior = theta_prior, beta_prior = beta_prior)
  }else{
    # CHECK if beta_prior if provided a message that it's not going to be used.
    list(model = model, theta_prior = theta_prior)
  }
}

