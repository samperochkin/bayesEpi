#' Constructor for ccModel (S3) objects.
#' @param response Name of the response variable (e.g., count).
#' @param time_index Name of the time variable (e.g., date or time)
#' @param fixed Named list of fixed effects. Each name should correspond to a variable in your data.
#' @param random Named list of random effects. Each name should correspond to a variable in your data.
#' @param overdispersion List specifying the model for the overdispersion terms. This is limited to `gaussian_effect()` for now.
#' @param design List specifying the case control referent frames.
#' @param control_aghq List of parameters to be passed to the aghq::quad function.
#' @return A list of lists of class "ccModel".
#' @examples
#' model <- ccModel(response = "count", time_index = "date",
#'                  fixed = list("hum_mean" = fixedEffect(gaussian_prior())),
#'                  random = list("o3_lag" = randomEffect(rw_effect(), pc_prep_prior()),
#'                                "temp_lag" = randomEffect(rw_effect(), pc_prep_prior())),
#'                  overdispersion = randomEffect(gaussian_effect(), pc_prep_prior(u = 10, alpha = .9)),
#'                  design = ccDesign(), control_aghq = controlAGHQ())
#' @export
ccModel <- function(response, time_index, fixed, random, overdispersion, design, control_aghq){
  model <- mget(names(formals()),sys.frame(sys.nframe()))
  attr(model, "class") <- "ccModel"
  model
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



#' Specify the case-crossover design to be used by fitModel.
#'
#' Function for specifying the case control design.
#'
#' A default design is returned by `ccDesign()`. This latter can be changed, e.g., with `ccDesign(lag = 2)`.
#' Arguments of interest are `scheme` ("unidirectional", "bidirectional", "time stratified"),
#' `n_control` (number of control days, for uni- and bidirectional schemes), `lag` (e.g., 7 for ensuring that
#' case and control fall on the same weekday) and `stratum_rule` (for time stratified schemes, the only rule
#' implemented is "month").
#'
#' @param ... See details.
#' @examples
#' ccDesign()
#' @export
ccDesign <- function(...){

  params = list(...)

  # default cc parameters
  design <- list(scheme = "bidirectional",
                 n_control = 2,
                 lag = 7,
                 stratum_rule = NULL)
  design[names(params)] <- params
  if(design$scheme == "time stratified" & is.null(design$stratum_rule)) design$stratum_rule <- "month"
  return(design)
}



#' Set aghq::quad parameters for model fitting.
#'
#' .
#'
#' So far, controlAGHQ() holds only one parameter, `k`, the number of quadrature points.
#'
#' @param ... See details.
#' @return A list of parameters for use by aghq::aghq.
#' @examples
#' controlAGHQ()
#' controlAGHQ(k=5)
#' @export
controlAGHQ <- function(...){

  params = list(...)

  # default cc parameters
  control <- list(k = 3)

  control[names(params)] <- params
  return(control)
}