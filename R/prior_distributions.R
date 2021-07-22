#' Specify a Gaussian random walk prior distribution.
#'
#' @param order Order of the random walk.
#' @param poly_degree Degree of the polynomial interpolation at reference value; see details.
#' @param ref_value Either a value or a function to be used for computing the reference value.
#' @param binwidth Binwidth to be used for discretizing the continuous variable; see details.
#' @return A list specifying a model for the nonlinear random effect.
#' @examples
### SAM : TALK WITH ALEX FOR MIDMAT THING!
#' @export
rw_prior <- function(order = 3,
                     poly_degree = order-1,
                     ref_value = median,
                     binwidth = 1){
  as.list(match.call()[-1])
}
