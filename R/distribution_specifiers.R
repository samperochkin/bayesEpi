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
### SAM : TALK WITH ALEX FOR MIDMAT THING!
#' @export
rw_effect <- function(order = 2,
                      poly_degree = order-1,
                      ref_value = median,
                      binwidth = 1){

  list(type="random walk", params=mget(names(formals()),sys.frame(sys.nframe())))
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
gaussian_prior <- function(prec = .001){
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
#' pc_prep_prior(alpha = .5, u = .01)
#' @export
pc_prep_prior <- function(alpha = .5, u = .001){
  list(type="pc_prep", params=mget(names(formals()),sys.frame(sys.nframe())))
}



#' Specify a gamma prior distribution.
#'
#' @param shape Shape parameter (\eqn{\alpha} or \eqn{k}) of the gamma distribution.
#' @param rate Rate parameter (\eqn{\beta}) of the gamma distribution.
#' @param scale Scale parameter (\eqn{\theta = \beta^{-1}}) of the gamma distribution. Cannot be specified if `rate` is.
#' @return A list specifying a gamma prior distribution.
#' @examples
#' gamma_prior(shape=1, rate=.1)
#' gamma_prior(shape=1, scale=10)
#' @export
gamma_prior <- function(shape=1, rate=.1, scale=NULL){
  if(!is.null(rate) & !is.null(scale)) stop("Only one of shape and scale can be passed to gamma_prior")
  if(!is.null(scale)) rate <- 1/scale
  list(type="gamma", params=list(shape=shape, rate=rate))
}
