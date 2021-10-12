#' Fit a model
#'
#' Function used for fitting either a ccModel (case crossover model) or a tsModel (time series model, yet to be implemented).
#' See the vignettes for examples of its usage.
#'
#' @param model A model, either a ccModel or a tsModel object.
#' @param data The dataset on which the model is fitted. It should contain all the variables involved in the model.
#' @param silent A boolean indicating whether you wish to silence the `aghq::marginal_laplace_tmb()` function (default is `TRUE`).
#' @return An updated ccModel or tsModel object that contains the fitted model (e.g., in ccModel$fit).
#' @export
fitModel <- function(model, data, silent = F, dll = NULL){
  UseMethod("fitModel", model)
}

#' @method fitModel ccModel
#' @export
#' @useDynLib cc
#' @import TMB
#' @importFrom magrittr %>%
fitModel.ccModel <- function(model, data, silent = F, dll = NULL){

  data <- as.data.frame(data)
  list2env(checkups(model, data), envir = environment())

  # contruct and import case_day and control_days. modifies data.
  list2env(getCaseControl(data, model), envir = environment())
  list2env(setRefValues(data, model), envir = environment())

  if(any(is.na(data[c(names(model$random),names(model$fixed))]))) stop("Problem in preprocessing?")

  # overdispersion and design matrices
  overdispersion <- !is.null(model$overdispersion)
  X <- as.matrix(data[names(model$fixed)])
  U <- as.matrix(data[names(model$random)])

  # discretizes nonlinear random effects, builds polynomial interpolation around reference values,
  # creates As, Xs_int and gamma_dims
  list2env(discretizedRandomDesigns(model, U), envir = environment())

  # prior parameters
  prior_lookup <- c("pc_prec", "gamma")
  beta_prec = c(purrr::map(model$fixed, ~ .x$prior$params$prec), purrr::map(model$random, ~ .x$beta_prior$params$prec)) %>% unlist
  theta_prior <- c(purrr::map(model$random, ~ .x$theta_prior$type), z = model$overdispersion$theta_prior$type) %>% unlist
  theta_prior_id = match(theta_prior , prior_lookup)
  theta_hypers = c(purrr::map(model$random, ~ .x$theta_prior$params), z = model$overdispersion$theta_prior$params) %>% unlist


  # Model fit ---------------------------------------------------------------
  z_pos <- unlist(lapply(1:model$design$lag, function(x) seq(x,nrow(data), (model$design$n_control + 1) * model$design$lag)))
  z_pos <- (1:nrow(data))[-z_pos]
  tmb_data <- list(count = data[case_day, model$response],
                   case_day = case_day, control_days = control_days,
                   X = cbind(X,Reduce("cbind", Xs_int)), A = Reduce("cbind", As),
                   Q = constructQ(model$random), gamma_dims = gamma_dims,
                   beta_prec = beta_prec, theta_prior_id = theta_prior_id, theta_hypers = theta_hypers,
                   z_pos = z_pos)

  theta_init <- getPriorInit(model)
  # theta_init <- rep(10,length(theta_init))

  parameters <- list(beta = rep(0,ncol(X)+sum(sapply(Xs_int, ncol))),
                     gamma = rep(0, sum(sapply(As, ncol))),
                     z = rep(0,length(z_pos)*overdispersion),
                     # z = rep(0,nrow(data)*overdispersion),
                     theta = theta_init)


  if(is.null(dll)) dll <- "bayesEpi"
  obj <- TMB::MakeADFun(tmb_data, parameters, random = c("beta","gamma","z"), DLL=dll, hessian=T)
  # obj <- TMB::MakeADFun(tmb_data, parameters, random = c("beta","gamma","z"), DLL=dll, hessian=T, silent=silent)

  if(silent){
    invisible(capture.output(quad <- aghq::marginal_laplace_tmb(obj, model$control_aghq$k, theta_init)))
  }else{
    quad <- aghq::marginal_laplace_tmb(obj, model$control_aghq$k, theta_init)
  }

  list(quad = quad, obj = obj, model = model)
}



# Helpers -----------------------------------------------------------------

# Checkups
checkups <- function(model, data){
  findVariables(model, data)
  list2env(renameResponseTimeIndex(model, data), envir = environment())
  list2env(removeNA(model, data), envir = environment())
  list(model = model, data = data)
}

# check if all variables are actually in the data.
findVariables <- function(model, data){
  var_names <- c(model$response, model$time_index, names(model$fixed),names(model$random))
  if(!all(var_names %in% names(data))){
    stop("One of the names provided for the response, time_index, fixed effects and/or random effects could not be matched in the data provided.")
  }
}

# For convenience
renameResponseTimeIndex <- function(model, data){

  if(!(model$response == "count")){
    if("count" %in% names(data)){
      message("Renaming the response 'count' and removing the original 'count' variable.")
      data$count <- NULL
    }
    names(data)[names(data) == model$response] <- "count"
    model$response <- "count"
  }

  if(!(model$time_index == "time")){
    if("time" %in% names(data)){
      message("Renaming the time_index 'time' and removing the original 'time' variable.")
      data$time <- NULL
    }
    names(data)[names(data) == model$time_index] <- "time"
    model$time_index <- "time"
  }

  list(model = model, data = data)
}

# check for NA and remove them.
removeNA <- function(model, data){
  var_names <- c(model$response, model$time_index, names(model$fixed),names(model$random))
  NA_rows <- which(apply(is.na(data[var_names]),1,any))

  if(length(NA_rows) > 0){
    data <- data[-NA_rows,]
    message(length(NA_rows), "row(s) with one or more variables of interest with value NA were found. They were removed.")
    return(list(data = data))
  }else{
    return(list())
  }
}

# Set reference values, in case a function was provided as `ref_value`.
setRefValues <- function(data, model){

  for(name in names(model$random)){
    ref_value <- model$random[[name]]$model$params$ref_value
    if(is.function(ref_value)) model$random[[name]]$model$params$ref_value <- ref_value(data[,name])
  }

  list(model = model)
}

# Identify rows/columns that need to be removed in A and Q.
getColsToRemove <- function(ref_value_pos, order){
  removed_cols <- ref_value_pos
  if(order >= 2) removed_cols <- c(removed_cols,ref_value_pos+1)
  if(order >= 3) removed_cols <- c(ref_value_pos-1,removed_cols)
  if(order > 3) stop("Random walks of order > 3, which you specified for ", name, " are not yet implemented.")
  removed_cols
}

# Discretize range(U[,random_effect]) into bins so that a random walk model can be fitted.
discretizedRandomDesigns <- function(model, U){

  random <- model$random
  if(is.null(random)){
    return(list(As = list(as(matrix(nrow=nrow(U), ncol=0), "dgTMatrix")),
                Xs_int = list(matrix(nrow=nrow(U), ncol=0)),
                gamma_dims = integer(0),
                model = model))
  }

  random_names <- names(random)
  As <- list()

  for(name in random_names){

    if(!("binwidth" %in% names(random[[name]]$model$params))) stop("binwidth is not specified for random effect ", name,".")

    random_params <- model$random[[name]]$model$params

    new_u <- round(U[,name]/random_params$binwidth)
    bin_values <- min(new_u):max(new_u) * random_params$binwidth
    fac <- factor(new_u, levels = seq(min(new_u), max(new_u),1), labels = paste0(name,"__",bin_values))
    A <- Matrix::t(Matrix::fac2sparse(fac, drop.unused.levels = F))

    # Set reference value by setting corresponding column of A (and neighbours) to zero.
    ref_value_pos <- which.min(abs(random_params$ref_value - bin_values))
    rounded_ref_value <- bin_values[ref_value_pos]
    removed_cols <- getColsToRemove(ref_value_pos, random_params$order)

    As[[length(As)+1]] <- A[, -removed_cols]

    # for later
    model$random[[name]]$model$extra$bin_values <- bin_values
    model$random[[name]]$model$extra$rounded_ref_value <- rounded_ref_value
    model$random[[name]]$model$extra$ref_value_pos <- ref_value_pos
    model$random[[name]]$model$extra$removed_cols <- removed_cols

    if(random[[name]]$model$params$poly_degree > 0){
      model$random[[name]]$model$extra$bin_values_int <- poly(bin_values - rounded_ref_value, #****************************************
                                                              degree = model$random[[name]]$model$params$poly_degree,
                                                              raw = TRUE)
      # model$random[[name]]$model$extra$bin_values_int <- poly((bin_values - rounded_ref_value)/diff(range(bin_values)),
      #                                                         degree = model$random[[name]]$model$params$poly_degree,
      #                                                         raw = TRUE)
    }
  }

  names(As) <- random_names
  Xs_int <- interpolationFixedEffects(model$random, U)
  gamma_dims <- sapply(As, ncol)
  gamma_dims <- gamma_dims[gamma_dims != 0]

  list(As = As, gamma_dims = gamma_dims, Xs_int = Xs_int, model = model)
}
#

# Create new fixed effects as by-products of random effects.
interpolationFixedEffects <-  function(random, U){

  random_names <- names(random)
  sapply(random_names, function(name) {
    if ("poly_degree" %in% names(random[[name]]$model$params)) {
      random_extra <- random[[name]]$model$extra
      poly_degree <- random[[name]]$model$params$poly_degree
      if (poly_degree == 0) return(matrix(nrow=nrow(U), ncol=0))

      X_new <- poly(U[, name] - random_extra$rounded_ref_value, #***********************************************
                    degree = poly_degree, raw = TRUE)
      # X_new <- poly((U[, name] - random_extra$rounded_ref_value)/diff(range(U[, name])),
      #               degree = poly_degree, raw = TRUE)
      colnames(X_new) <- paste0(name, "__", attr(X_new, "degree"))
      return(X_new)
    }
    else {
      return(matrix(nrow=nrow(U), ncol=0))
    }
  }, simplify = FALSE, USE.NAMES = TRUE)
}
#

# Create the case_day vector and the corresponding control_days matrix.
#' @import purrr
getCaseControl <- function(data, model){

  design <- model$design
  time <- as.integer(data$time)
  case_day <- time[data$count > 0]

  if(design$scheme == "unidirectional"){

    control_days <- purrr::map(-(design$n_control:1)*design$lag, ~ case_day + .x) %>% Reduce(f="cbind")
    if(design$n_control == 1) control_days <- as.matrix(control_days)

  }else if(design$scheme == "bidirectional"){

    if(design$n_control %% 2 == 0){a <- design$n_control/2; a <- design$lag*(-a:a)[-(a+1)]}
    else{a <- (design$n_control+1)/2; a <- (-a:a)[-c(a+1,2*a+1)]}
    control_days <- purrr::map(a, ~ case_day + .x) %>% Reduce(f="cbind")

  }else if(design$scheme == "time stratified"){

    case_day_id <- match(case_day, time)

    if(design$stratum_rule == "sequential"){
      t0 <- min(time)
      id <- paste(floor((time - t0)/(design$lag * (design$n_control+1))),
                  (time - t0) %% design$lag, sep = "-")

    }else if(design$stratum_rule == "month"){
      id <- paste(format(data$time, "%Y-%m"), time %% design$lag, sep=".")

    }else stop("The stratum rule", design$stratum_rule, "is not implemented.")

    stratum <- split(time, id)
    max_len <- max(sapply(stratum, length)) - 1
    control_days <- lapply(case_day_id, function(c_day_id){
      con <- setdiff(stratum[[id[c_day_id]]], time[c_day_id])
      con <- c(con, rep(0, max_len-length(con)))
      con
    }) %>% Reduce(f="rbind")

  }else{stop("The scheme", design$scheme, "is not implemented.")}


  # filter out case day with no control days
  keep <- apply(matrix(control_days %in% time, nrow=nrow(control_days)),1,any)
  case_day <- case_day[keep]
  control_days <- control_days[keep,,drop=F]

  # filter out days that are neither case nor control days
  keep <- time %in% unique(c(case_day,control_days))
  time <- time[keep]
  data <- data[keep,]

  case_day <- (1:nrow(data))[match(case_day, time)]
  control_days <- matrix((1:nrow(data))[match(control_days, time, nomatch = NA)], nrow(control_days))
  control_days[is.na(control_days)] <- 0
  if(any(rowSums(control_days) == 0)) stop("Error in selecting the control days")

  list(data = data, case_day = case_day, control_days = control_days)
}
#

# Construct the precision matrix Q for the random effects (Gaussian random walks).
constructQ <- function(random){

  if(is.null(random)) return(as(matrix(nrow=0,ncol=0), "dgTMatrix"))

  createD <-function(d,p){
    if(p==0) return(Diagonal(d,1))
    D <- Matrix::bandSparse(d,k =c(0,1),diagonals =list(rep(-1,d),rep(1,d-1)))[-d, ]
    if(p==1) return(D)
    else return(createD(d,p-1)[-1,-1] %*% D)
  }

  Qs <- lapply(random, function(ran){
    if(ran$model$type == "random walk"){

      removed_cols <- ran$model$extra$removed_cols
      order <- ran$model$params$order
      Matrix::crossprod(createD(length(ran$model$extra$bin_values), order)[,-removed_cols])

    }else stop("Not sure how to construct the Q matrix associated to a ", ran$model$type, " random effect model")
  })

  Matrix::bdiag(Qs)
}

# Compute initial theta parameter to be passed to aghq::quad.
getPriorInit <- function(model){
  random_priors <- purrr::map(model$random, ~ .x$theta_prior)
  if(!is.null(model$overdispersion)){
    random_priors <- c(random_priors, list(model$overdispersion$theta_prior))
  }

  sapply(random_priors, function(ran_prior){
    if(ran_prior$type == "pc_prec"){
      return(-2*log(-ran_prior$params$u/log(ran_prior$params$alpha)))
    }else if(ran_prior$type == "gamma"){
      return(ran_prior$params$shape/ran_prior$params$rate)
    }else{
      stop("Can't get inital values for some theta parameter... Check getPriorInit.")
    }
  })
}
