#' Creates a new_variable so that new_variables[t,] == variables[t-lag,], where t is the time_index (e.g. a variable "date") of your data frame.
#'
#' @param df A data.frame.
#' @param variables A character vector of variable names.
#' @param lags A integer vector giving the lags to use (must be of length 1 or same as variables).
#' @param remove_original A boolean indicating whether the variables should be removed from the data.frame. Default is FALSE.
#' @param ext A character vector indicating giving the names of the new_variables. Default is paste0("_", lags, "lag").
#' @return An updated data.frame.
#' @examples
# createLaggedVariable <- function(df, variable, lags, remove_original = F, ext = paste0("_", lags, "lag")){
#   stop("The function preprocessing_functions/createLaggedVariable is not yet written.")
# }

