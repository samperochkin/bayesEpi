#' Specify a Gaussian random walk distribution for random effect.
#'
#' @param order Order of the random walk.
#' @param poly_degree Degree of the polynomial interpolation at reference value; see details.
#' @param ref_value Either a value or a function to be used for computing the reference value.
#' @param binwidth Binwidth to be used for discretizing the continuous variable; see details.
#' @return A list specifying a model for a nonlinear random effect.
#' @examples
### SAM : TALK WITH ALEX FOR MIDMAT THING!
#' @export
rw_effect <- function(order = 3,
                      poly_degree = order-1,
                      ref_value = median,
                      binwidth = 1){

  c(type="random_walk", mget(names(formals()),sys.frame(sys.nframe())))
}

#' Specify a Gaussian prior distribution.
#'
#' @param prec TO DO.
#' @return A list specifying a Gaussian prior distribution.
#' @examples
#' @export
gaussian_prior <- function(prec = .001){
  c(type="gaussian", mget(names(formals()),sys.frame(sys.nframe())))
}


#' Specify an exponential (Penalized Complexity) prior distribution.
#'
#' @param prec TO DO.
#' @return A list specifying a pc prep prior distribution.
#' @examples
#' @export
pc_prep_prior <- function(alpha = .5, u = .01){
  c(type="pc_prep", mget(names(formals()),sys.frame(sys.nframe())))
}


#' Specify a gamma prior distribution.
#'
#' @param prec TO DO.
#' @return A list specifying a gamma prior distribution.
#' @examples
#' @export
gamma_prior <- function(shape, rate = 1, scale = 1/rate){
  c(type="gamma", mget(names(formals()),sys.frame(sys.nframe())))
}
