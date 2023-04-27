#' Specify a polynomial fixed effect.
#'
#' @param poly_degree Degree of the polynomial.
#' @return A list specifying a polynomial fixed effect.
#' @export
poly_effect <- function(degree = 1,
                        ref_value = median){

  list(type="poly", params=mget(names(formals()),sys.frame(sys.nframe())))
}


#' Specify a polynomial fixed effect.
#'
#' @param poly_degree Degree of the polynomial.
#' @return A list specifying a polynomial fixed effect.
#' @export
bs_effect <- function(knots,
                      degree = 3,
                      ref_value = knots[[2]][1]){

  list(type="bs", params=mget(names(formals()),sys.frame(sys.nframe())))
}


#' Specify a Gaussian random walk distribution for random effect.
#'
#' This is used to approximate Gaussian processes random effects.
#' When a Gaussian random walk is specified, the range of observed values is discretized into
#' bins and the random effect is modeled as a random walk on these. For technical reasons, the
#' bin that contain the reference value is removed, as well as some of its neighbors, so that
#' the number of bins removed equals `order`. By default, fixed effects (x - x_ref)^i, i=1,...,`order`-1
#' are then added to the model.
#'
#' @param order Order of the random walk. Restricted to 1, 2 and 3.
#' @param poly_degree Degree of the polynomial interpolation at reference value; see details.
#' @param ref_value Either a value or a function to be used for computing the reference value.
#' @param binwidth Binwidth to be used for discretizing the continuous variable; see details.
#' @return A list specifying a random walk model for a nonlinear random effect.
#' @examples
#' rw_effect()
#' @export
rw_effect <- function(order = 2,
                      poly_degree = order-1,
                      ref_value = median,
                      binwidth = 1){

  list(type="random walk", params=mget(names(formals()),sys.frame(sys.nframe())))
}



#' Specify an integrated Wiener process for random effect.
#'
#' This is used to approximate Gaussian processes random effects.
#' By default, fixed effects (x - x_ref)^i, i=1,...,`order`-1
#' are added to the model.
#'
#' @param order Order of the first non-zero derivative.
#' @param poly_degree Degree of the polynomial for additional fixed effects (by default `order`-1).
#' Note that the intercept is excluded.
#' @param ref_value Either a value or a function to be used for computing the reference value.
#' @param knots Number of (equally spaced) knots to use.
#' @return A list specifying a integrated Wiener process model for a nonlinear random effect.
#' @examples
#' iwp_effect()
#' @export
iwp_effect <- function(order = 2,
                       poly_degree = order-1,
                       ref_value = median,
                       knots){

  list(type="integrated Wiener process", params=mget(names(formals()),sys.frame(sys.nframe())))
}



#' Specify a Gaussian random effect.
#'
#' This is used to include overdispersion in the model. This is the only model implemented so far.
#' It does practically nothing except holding the character "gaussian".
#'
#' @return A list specifying a centered Gaussian distribution.
#' @export
gaussian_effect <- function(){
  list(type="gaussian")
}



#' Specify a Gaussian prior distribution.
#'
#' @param prec TO DO.
#' @return A list specifying a Gaussian prior distribution.
#' @examples
#' gaussian_prior()
#' @export
gaussian_prior <- function(prec = .01){
  list(type="gaussian", params = mget(names(formals()),sys.frame(sys.nframe())))
}



#' Specify an exponential (penalized complexity) prior distribution.
#'
#' PC prior distribution; see \insertCite{Simpson:2017;textual}{bayesEpi}. The parameters are such that \eqn{P(\sigma > u) = \alpha}, where \eqn{\sigma = -2log(\theta)}.
#'
#' @importFrom Rdpack reprompt
#' @references{
#'   \insertRef{Simpson:2017}{bayesEpi}
#' }
#' @param alpha A probability.
#' @param u A value greater than 0.
#' @return A list specifying a pc prep prior distribution.
#' @examples
#' pc_prec_prior(alpha = .5, u = .01)
#' @export
pc_prec_prior <- function(alpha = .5, u = .1){
  list(type="pc_prec", params=mget(names(formals()),sys.frame(sys.nframe())))
}



#' Specify an log-gamma prior distribution.
#'
#' @param alpha shape parameter.
#' @param beta rate parameter.
#' @return A list specifying a log-gamma prior distribution.
#' @examples
#' log_gamma_prior(alpha = .3, beta = 1e-6)
#' @export
log_gamma_prior <- function(shape = .3, rate = 1e-6){
  list(type="log_gamma", params=mget(names(formals()),sys.frame(sys.nframe())))
}
