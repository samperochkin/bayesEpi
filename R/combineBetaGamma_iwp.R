#' Contruct o-spline from knots and boundary conditions.
#'
#' Function used to evaluate random effects (o-splines, including fixed effects for boundary conditions) on a given grid.
#'
#' @param samps Corresponding entry of the output of `aghq::marginal_laplace_tmb()`.
#' @param model Object of class "ccModel" with iwp random effects (fitted with o-splines).
#' @param grids Named list of grids on which to evaluate the o-splines.
#' @return A modification of the samps object containing the o-splines evaluated at the given grid points.
#'
#' @details If `grids` is NULL, then the grids for all random effects are constructed automatically.
#' If grids is a named list, then only the named random effects are considered;
#' if their corresponding grid is NULL, then it is constructed automatically.
#'
#' @import OSplines
#' @export
combineBetaGamma_iwp <- function(samps, model, grids = NULL, silent = FALSE){

  fixed_names <- names(model$fixed)
  random_names <- names(model$random)
  if(is.null(random_names)){
    cat("No random effects to work with. Returning samps as is.\n")
    return(samps)
  }

  K <- length(grids)
  grids <- grids[[intersect(random_names, names(grids))]]
  if(length(grids) < K & !silent) cat("Some name(s) of grids could not be matched. They were dropped.\n")

  if(is.null(names(grids))){
    if(!silent) cat("No grid was specified. All random effects are considered.\n")
    grids <- replicate(length(random_names), NULL)
    names(grids) <- random_names
  }

  grids_names <- names(grids)

  if(any(sapply(grids, is.null)) & !silent) cat("Some named grid(s) were set to `NULL`.
  Using stepsizes of length 1 and min/max knots to determine the range.\n")

  for(nam in grids_names[sapply(grids, is.null)]){
    grids[[nam]] <- seq(min(model$random[[nam]]$model$params$knots),
                     max(model$random[[nam]]$model$params$knots), 1)
  }


  new_samps <- NULL
  counter_beta <- length(fixed_names)
  counter_gamma <- sum(rownames(samps) == "beta") - counter_beta

  for(nam in grids_names){

    random_params <- model$random[[nam]]$model$params
    knots <- random_params$knots
    ref_value <- random_params$ref_value

    u <- grids[[nam]]
    ref_pos <- which(knots == ref_value)
    if(length(ref_pos) == 0){
      stop("The reference value provided for random effect", nam, "could not be matched in the corresponding set of knots. Abort.\n")
    }
    A <- NULL
    if(ref_pos != 1) A <- cbind(A, OSplines::local_poly(knots = rev(ref_value - knots[1:ref_pos]),
                                                        refined_x = ref_value - u,
                                                        p = random_params$order))
    if(ref_pos != length(knots)) A <- cbind(A, OSplines::local_poly(knots = knots[ref_pos:length(knots)] - ref_value,
                                                                    refined_x = u - ref_value,
                                                                    p = random_params$order))

    id_gamma <- counter_gamma + 1:ncol(A)
    counter_gamma <- counter_gamma + ncol(A)
    y <- A %*% samps[id_gamma,]
    if(random_params$poly_degree > 0){
      X <- poly(u - ref_value,
                degree = random_params$poly_degree,
                raw = TRUE)
      id_beta <- counter_beta + 1:ncol(X)
      counter_beta <- counter_beta + ncol(X)
      y <- y + X %*% samps[id_beta,]
    }

    new_samps <- c(new_samps, list(y))
  }
  names(new_samps) <- grids_names

  return(new_samps)
}
