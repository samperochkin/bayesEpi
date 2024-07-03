# Helpers -----------------------------------------------------------------

# Checkups
checkups <- function(model, data){
  findVariables(model, data)
  # list2env(renameResponseTimeIndex(model, data), envir = environment())
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
# renameResponseTimeIndex <- function(model, data){
#
#   if(model$response != "count"){
#     if("count" %in% names(data)){
#       message("Renaming the response 'count' and removing the original 'count' variable.")
#       data$count <- NULL
#     }
#     names(data)[names(data) == model$response] <- "count"
#     model$response <- "count"
#   }
#
#   if(!(model$time_index == "time")){
#     if("time" %in% names(data)){
#       message("Renaming the time_index 'time' and removing the original 'time' variable.")
#       data$time <- NULL
#     }
#     names(data)[names(data) == model$time_index] <- "time"
#     model$time_index <- "time"
#   }
#
#   list(model = model, data = data)
# }

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
  if(order >= 4) removed_cols <- c(ref_value_pos-2,removed_cols)
  if(order > 4) stop("Random walks of order > 4, which you specified for ", name, " are not yet implemented.")
  removed_cols
}

# builds design matrices for random effects (and those for the associated fixed effects)
createFixedDesigns <- function(model, X){

  fixed <- model$fixed

  # If no fixed effects
  if(is.null(fixed)) return(list(Xs_exp = list(matrix(nrow=nrow(X), ncol=0))))


  fixed_names <- names(fixed)
  Xs_exp <- list()

  for(name in fixed_names){

    fixed_params <- model$fixed[[name]]$model$params

    if(fixed[[name]]$model$type == "poly"){
      new_cols <- stats::poly(X[,name] - fixed_params$ref_value, degree = fixed_params$degree, raw = T)
      names(new_cols) <- paste0(name, "_", 1:fixed_params$degree)
      model$fixed[[name]]$model$extra$range <- range(X[,name])

    }else if(fixed[[name]]$model$type == "bs"){


      knots <- fixed_params$knots
      degree <- fixed_params$degree
      ref_value <- fixed_params$ref_value

      if(!(ref_value %in% knots[[1]] & ref_value %in% knots[[2]])) stop("ref_value of ", name, "cannot be found in the corresponding knots vector. \n")
      new_cols <- constructBS(x = X[,name], knots = knots, degree = degree, ref_value = ref_value)

    }else{
      stop("Invalid fixed effect model")

    }

    Xs_exp <- c(Xs_exp, list(new_cols))
  }

  names(Xs_exp) <- fixed_names
  X <- do.call("cbind", Xs_exp)

  list(X = X, Xs_exp = Xs_exp, model = model)
}
#

# builds design matrices for random effects (and those for the associated fixed effects)
createRandomDesigns <- function(model, U){

  random <- model$random

  # If no random effects
  if(is.null(random)){
    return(list(As = list(methods::as(matrix(nrow=nrow(U), ncol=0), "dgTMatrix")),
                Xs_int = list(matrix(nrow=nrow(U), ncol=0)),
                gamma_dims = integer(0),
                model = model))
  }


  random_names <- names(random)
  As <- list()

  for(name in random_names){

    random_params <- model$random[[name]]$model$params

    if(random[[name]]$model$type == "random walk"){
      if(!("binwidth" %in% names(random[[name]]$model$params))) stop("binwidth is not specified for random effect ", name,".")

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
        model$random[[name]]$model$extra$bin_values_int <- stats::poly(bin_values - rounded_ref_value,
                                                                degree = model$random[[name]]$model$params$poly_degree,
                                                                raw = TRUE)
      }

    }else if(random[[name]]$model$type == "integrated Wiener process"){

      knots <- random_params$knots
      ref_value <- random_params$ref_value
      ran <- model$random[[name]]$model$extra$range <- range(U[,name])

      if(!(ref_value %in% knots)) stop("ref_value of", name, "cannot be found in the corresponding knots vector. \n")
      if(!(ran[1] >= knots[1] & ran[2] <= knots[length(knots)])) warning("knots for ", name, " do not span its range. Continuing anyway. \n")
      if(length(knots) <= 2) stop("knots for ", name, " is too small")

      ref_pos <- which(knots == ref_value)
      A <- NULL
      if(ref_pos != 1) A <- cbind(A, methods::as(local_poly(knots = rev(ref_value - knots[1:ref_pos]),
                                                   refined_x = ref_value - U[,name],
                                                   p = random_params$order), "sparseMatrix"))
      if(ref_pos != length(knots)) A <- cbind(A, methods::as(local_poly(knots = knots[ref_pos:length(knots)] - ref_value,
                                                               refined_x = U[,name] - ref_value,
                                                               p = random_params$order), "sparseMatrix"))
      As[[length(As)+1]] <- A
    }

  }

  names(As) <- random_names
  Xs_int <- interpolationFixedEffects(model$random, U)
  gamma_dims <- sapply(As, ncol)
  gamma_dims <- gamma_dims[gamma_dims != 0]

  list(As = As, gamma_dims = gamma_dims, Xs_int = Xs_int, model = model)
}
#

# findKnotsPlacement <- function(k, ran, ref_value){
#
#   knots <- seq(ran[1], ran[2], length.out = k)
#   nearest_knot <- which.min(abs(knots - ref_value))
#   if(nearest_knot == 1) nearest_knot <- 2
#   if(nearest_knot == k) nearest_knot <- k-1
#   shift <- ref_value - knots[nearest_knot]
#
#   if(shift > 0) knots <- seq(ran[1], ran[2] + shift*(1+1/(nearest_knot-1)*(k-nearest_knot)), length.out = k)
#   if(shift < 0) knots <- seq(ran[1] + shift*(1+1/(k-nearest_knot)*(nearest_knot-1)), ran[2], length.out = k)
#
#   return(knots)
# }
#

interpolationFixedEffects <-  function(random, U){

  random_names <- names(random)
  sapply(random_names, function(name) {
    if ("poly_degree" %in% names(random[[name]]$model$params)) {

      poly_degree <- random[[name]]$model$params$poly_degree
      if (poly_degree == 0) return(matrix(nrow=nrow(U), ncol=0))

      if(random[[name]]$model$type == "random walk") cen <- random[[name]]$model$extra$rounded_ref_value
      if(random[[name]]$model$type == "integrated Wiener process") cen <- random[[name]]$model$params$ref_value

      X_new <- stats::poly(U[, name] - cen, degree = poly_degree, raw = TRUE)
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
  if(!(is.null(model$design$stratum_var))){
    design <- model$design
    if(design$scheme == "time stratified" & design$stratum_rule == "sequential"){
      time_stratum <- data[, c(model$time_index, design$stratum_var)]
      time_stratum[,model$time_index] <- sapply(data[, model$time_index], as.integer)
      time_stratum1 = apply(time_stratum , 1 , paste , collapse = "," )
      case_day_stratum <- time_stratum[data[, model$response] > 0,]
      case_day_stratum1 = apply(case_day_stratum , 1 , paste , collapse = "," )
      case_day_id <- match(case_day_stratum1, time_stratum1)
      t0 <- min(time_stratum[,model$time_index])

      stratum_variable <- as.matrix(data[, model$design$stratum_var])
      strata_1 =  as.matrix(floor((time_stratum[,model$time_index] - t0)/(design$lag * (design$n_control+1))))#block
      strata_2 = as.matrix((time_stratum[,model$time_index] - t0)%%(design$lag)) #day of week
      id = paste(strata_1, strata_2,stratum_variable,  sep = "-")

      # stata (case and control days togeteher)
      stratum <- split(time_stratum1, id)

      # number of columns of control_days matrix
      max_len <- max(sapply(stratum, length)) - 1

      # for each case day, enumerates control days (0 means empty)
      control_days <- lapply(case_day_id, function(c_day_id){
        con <- setdiff(stratum[[id[c_day_id]]], time_stratum1[c_day_id])
        con <- c(con, rep(0, max_len-length(con)))
        con
      }) %>% Reduce(f="rbind")
      # filter out case day with no control days
      keep <- apply(matrix(control_days %in% time_stratum1, nrow=nrow(control_days)),1,any)
      case_day_stratum1 <- case_day_stratum1[keep]
      control_days <- control_days[keep,,drop=F]

      # filter out days that are neither case nor control days
      keep <- time_stratum1 %in% unique(c(case_day_stratum1 ,control_days))
      time_stratum1 <- time_stratum1[keep]
      data <- data[keep,]

      case_day_stratum1 <- (1:nrow(data))[match(case_day_stratum1, time_stratum1)]
      control_days <- matrix((1:nrow(data))[match(control_days, time_stratum1, nomatch = NA)], nrow(control_days))
      control_days[is.na(control_days)] <- 0
      if(any(rowSums(control_days) == 0)) stop("Error in selecting the control days")

      list(data = data, case_day = case_day_stratum1, control_days = control_days)
    }
    else stop("The stratum rule ", design$stratum_rule, " is not implemented if a stratum variable is supplied.")
  }
  else{
  design <- model$design
  time <- as.integer(data[, model$time_index])
  case_day <- time[data[, model$response] > 0]
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
      # do something with model$design$stratum_var --- data[,model$design$stratum_var]
      # stop("error)
      # id for the stratum (window_id, dow_id)
      id <- paste(floor((time - t0)/(design$lag * (design$n_control+1))),
                  (time - t0) %% design$lag, sep = "-")
      # id <- paste(floor((time - t0)/(design$lag * (design$n_control+1))),
      #             (time - t0) %% design$lag,
      #             stratum_var, sep = "-")

    }else if(design$stratum_rule == "month"){
      id <- paste(format(data[, model$time_index], "%Y-%m"), time %% design$lag, sep=".")

    }else stop("The stratum rule", design$stratum_rule, "is not implemented.")

    # stata (case and control days togeteher)
    stratum <- split(time, id)

    # number of columns of control_days matrix
    max_len <- max(sapply(stratum, length)) - 1

    # for each case day, enumerates control days (0 means empty)
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
}
#

# Identify where to remove (not to put) overdispersion terms
selectFixedOD <- function(data, model, case_day, control_days){

  if(model$design$scheme == "time stratified"){
    z_rem <- unique(apply(cbind(case_day, control_days),1, function(days){
      sort(days[days != 0])[1]
    }))
  }else{
    lag_group <- as.integer(data[, time_index]) %% model$design$lag

    z_rem <- lapply(1:model$design$lag - 1, function(l){
      mem <- which(lag_group == l)

      max_lag <- if(model$design$scheme == "unidirectional"){
        model$design$n_control*model$design$lag
      }else if(model$design$scheme == "bidirectional"){
        # model$design$lag # if only complete ref. frames are allowed..
        model$design$n_control/2*model$design$lag
      }

      mem[c(1,which(diff(data$z[mem]) > max_lag) + 1)]
    }) %>% unlist
  }

  return(z_rem)
}


# Construct the precision matrix Q for the random effects (Gaussian random walks).
# NEED TO BE GENERALIZED TO ALLOW BOTH TYPES (i.e., return sparse diag matrix for iwp too)
# BUT FOR THIS, the TMB templates need to be merged...
# constructQ <- function(random){
#
#   if(is.null(random)) return(as(matrix(nrow=0,ncol=0), "dgTMatrix"))
#
#   createD <-function(d,p){
#     if(p==0) return(Diagonal(d,1))
#     D <- Matrix::bandSparse(d,k =c(0,1),diagonals =list(rep(-1,d),rep(1,d-1)))[-d, ]
#     if(p==1) return(D)
#     else return(createD(d,p-1)[-1,-1] %*% D)
#   }
#
#   Qs <- lapply(random, function(ran){
#     if(ran$model$type == "random walk"){
#       removed_cols <- ran$model$extra$removed_cols
#       order <- ran$model$params$order
#       return(Matrix::crossprod(createD(length(ran$model$extra$bin_values), order)[,-removed_cols]))
#     }else if(ran$model$type == "integrated Wiener process"){
#       return(diag(compute_weights_precision(x = ran$model$params$knots)))
#     }
#   })
#
#   # HERE ****** TO GENERALIZE
#   if(random[[1]]$model$type == "random walk") return(Matrix::bdiag(Qs))
#   if(random[[1]]$model$type == "integrated Wiener process") return(unlist(Qs))
# }

constructQ_rw <- function(random){

  if(is.null(random)) return(methods::as(matrix(nrow=0,ncol=0), "dgTMatrix"))
  ids <- which(sapply(random, function(ran) ran$model$type == "random walk"))
  if(length(ids) == 0) return(methods::as(matrix(nrow=0,ncol=0), "dgTMatrix"))

  createD <-function(d,p){
    if(p==0) return(Matrix::Diagonal(d,1))
    D <- Matrix::bandSparse(d,k =c(0,1),diagonals =list(rep(-1,d),rep(1,d-1)))[-d, ]
    if(p==1) return(D)
    else return(createD(d,p-1)[-1,-1] %*% D)
  }

  Qs <- lapply(random[ids], function(ran){
    removed_cols <- ran$model$extra$removed_cols
    order <- ran$model$params$order
    Matrix::crossprod(createD(length(ran$model$extra$bin_values), order)[,-removed_cols])
  })

  methods::as(methods::as(Matrix::bdiag(Qs), "generalMatrix"), "TsparseMatrix")
}

#' @import OSplines
constructQ_iwp <- function(random){

  if(is.null(random)) return(numeric(0))
  ids <- which(sapply(random, function(ran) ran$model$type == "integrated Wiener process"))
  if(length(ids) == 0) return(numeric(0))

  unlist(lapply(random[ids], function(ran){
      diag(OSplines::compute_weights_precision(x = ran$model$params$knots))
  }))
}



# Compute initial theta parameter to be passed to aghq::quad.
getPriorInit <- function(model, init_od_to_none = F){
  random_priors <- purrr::map(model$random, ~ .x$theta_prior)
  if(!is.null(model$overdispersion) & !init_od_to_none){
    random_priors <- c(random_priors, list(model$overdispersion$theta_prior))
  }

  if(length(random_priors) == 0) return(numeric(0))

  theta_init <- sapply(random_priors, function(ran_prior){
    if(ran_prior$type == "pc_prec"){
      return(-2*log(-ran_prior$params$u/log(ran_prior$params$alpha)))
    }else if(ran_prior$type == "log_gamma"){
      # return(0)
      return(digamma(ran_prior$params$shape) - log(ran_prior$params$rate))
    }else{
      stop("Can't get inital values for some theta parameter... Check getPriorInit.")
    }
  })

  if(!is.null(model$overdispersion) & init_od_to_none){
    # theta_init = 10 means that var = exp(-10) = .0000454 (i.e. almost no overdispersion at init)
    theta_init <- c(theta_init, 100)
  }

  return(theta_init)
}
