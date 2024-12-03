###########################################################################
# Helpers -----------------------------------------------------------------
###########################################################################





# Checkups ----------------------------------------------------------------
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





# General setup -----------------------------------------------------------
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
      }) |> Reduce(f="rbind")
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
      control_days <- purrr::map(-(design$n_control:1)*design$lag, ~ case_day + .x) |> Reduce(f="cbind")
      if(design$n_control == 1) control_days <- as.matrix(control_days)
    }else if(design$scheme == "bidirectional"){
      if(design$n_control %% 2 == 0){a <- design$n_control/2; a <- design$lag*(-a:a)[-(a+1)]}
      else{a <- (design$n_control+1)/2; a <- (-a:a)[-c(a+1,2*a+1)]}
      control_days <- purrr::map(a, ~ case_day + .x) |> Reduce(f="cbind")
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
      }) |> Reduce(f="rbind")

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





# Design matrices ---------------------------------------------------------

# helper to apply transformation if needed (used to apply reference value too)
applyTransformation <- function(frmodel, xORu){

  ref_value <- frmodel$params$ref_value
  lambda <- frmodel$params$lambda
  c <- frmodel$params$c

  # if(is.null(lambda)) return(xORu - ref_value) # REFREF
  # if(lambda == 0) return(log(xORu + c) - log(ref_value + c)) # REFREF
  # if(lambda != 0) return((xORu + c)^lambda - (ref_value + c)^lambda) # REFREF

  if(is.null(lambda)) return(xORu - ref_value)

  if(is.null(c)) c <- 0
  tc <- c + ref_value

  if(lambda == 0) return(tc*(log(xORu + c) - log(tc)))
  if(lambda != 0) return(((xORu + c)^lambda - tc^lambda)/(lambda * tc^(lambda - 1))) # REFREF
}

applyInvTransformation <- function(frmodel, xORu){

  ref_value <- frmodel$params$ref_value
  lambda <- frmodel$params$lambda
  c <- frmodel$params$c

  if(is.null(lambda)) return(xORu + ref_value)

  if(is.null(c)) c <- 0
  tc <- c + ref_value

  if(lambda == 0) return(exp(xORu/tc + log(tc)) - c)
  if(lambda != 0) return((xORu*(lambda * tc^(lambda - 1)) + tc^lambda)^(1/lambda) - c)
}

# helper to distribute appropriately distribute knots around
# reference value for mgp random effects (note: output standardized so that ref_value = 0)
# also used for iwp
# NOTE: IT IS BEST TO SPECIFY ONES OWN KNOTS. THIS IS SHAKY
splitKnots <- function(region, ref_value, range, stepsize, type = "mgp",
                       extra_left=0, extra_right=0){

  if(range[1] < region[1] || range[2] > region[2])
    stop("For random effect ", name, ": Make sure that region covers all of the data.")
  if(ref_value < region[1] || ref_value > region[2])
    stop("For random effect ", name, ": Make sure ref_value is inside the region provided.")

  region[1] <- ref_value - ceiling((ref_value - region[1])/stepsize)*stepsize
  region[2] <- ref_value + ceiling((region[2]-ref_value)/stepsize)*stepsize

  region <- region + c(extra_left, extra_right)*stepsize
  knots_neg <- seq(region[1]-stepsize*extra_left, ref_value, stepsize)
  knots_pos <- seq(ref_value, region[2]+stepsize*extra_right, stepsize)

  if(length(knots_neg) < 3 | length(knots_pos) < 3)
    stop("For random effect ", name, ": region and stepsize leads to too few knots.")

  if(type == "iwp") return(c(knots_neg, knots_pos[-1]))
  if(type == "mgp"){
    knots_neg <- length(knots_neg)
    knots_pos <- length(knots_pos)
    region_neg <- c(0, ref_value-min(region))
    region_pos <- c(0, max(region)-ref_value)
    return(list(knots = list(neg = knots_neg, pos = knots_pos),
                regions = list(neg = region_neg, pos = region_pos)))
  }
}


# builds design matrices for fixed effects
createFixedDesigns <- function(model, X){

  fixed <- model$fixed

  # If no fixed effects
  if(is.null(fixed)) return(list(Xs_exp = list(matrix(nrow=nrow(X), ncol=0))))

  fixed_names <- names(fixed)
  Xs_exp <- list()

  for(name in fixed_names){

    fmodel <- model$fixed[[name]]$model
    fixed_params <- fmodel$params
    ref_value <- fixed_params$ref_value
    model$fixed[[name]]$model$extra$range <- range(X[,name])

    # This is where transformations are handled!
    x <- applyTransformation(fmodel, X[,name])

    if(fixed[[name]]$model$type == "poly"){
      new_cols <- stats::poly(x, degree = fixed_params$degree, raw = T)
      names(new_cols) <- paste0(name, "_", 1:fixed_params$degree)

    }else if(fmodel$type == "bs"){

      knots <- fixed_params$knots
      degree <- fixed_params$degree
      new_knots <- list(applyTransformation(fmodel, knots[[1]]),
                        applyTransformation(fmodel, knots[[2]]))

      # if(!(ref_value %in% knots[[1]] & ref_value %in% new_knots[[2]])) stop("ref_value of ", name, "cannot be found in the corresponding knots vector. \n")
      # new_cols <- constructBS(x = X[,name], knots = knots, degree = degree, ref_value = ref_value)
      if(!(0 %in% new_knots[[1]] & 0 %in% new_knots[[2]])) stop("ref_value of ", name, "cannot be found in the corresponding knots vector. \n")
      new_cols <- constructBS(x = x, knots = new_knots, degree = degree, ref_value = 0)

    }else{
      stop("Invalid fixed effect model")

    }

    Xs_exp <- c(Xs_exp, list(new_cols))
  }

  names(Xs_exp) <- fixed_names
  X <- do.call("cbind", Xs_exp)

  list(X = X, Xs_exp = Xs_exp, model = model)
}


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

    rmodel <- model$random[[name]]$model
    random_params <- rmodel$params
    ref_value <- random_params$ref_value
    ran <- rmodel$extra$range <- range(U[,name])
    random_type <- rmodel$type

    # This is where transformations are handled!
    u <- applyTransformation(rmodel, U[,name])
    tref_value <- applyTransformation(rmodel, random_params$ref_value)

    if(random_type == "random walk"){
      if(!("binwidth" %in% names(random[[name]]$model$params))) stop("binwidth is not specified for random effect ", name,".")

      # construct bins (binwidth understood on the transformed scale)
      new_u <- round(u/random_params$binwidth)
      bin_tvalues <- min(new_u):max(new_u) * random_params$binwidth
      fac <- factor(new_u, levels = seq(min(new_u), max(new_u),1), labels = paste0(name,"__",bin_tvalues))
      A <- Matrix::t(Matrix::fac2sparse(fac, drop.unused.levels = F))

      # Set reference value by setting corresponding column of A (and neighbours) to zero.
      # Note that the reference value (on the transformed scale) is zero...
      # ref_value_pos <- which.min(abs(bin_values - ref_value))
      ref_value_pos <- which.min(abs(bin_tvalues))
      rounded_ref_tvalue <- bin_tvalues[ref_value_pos]
      if(rounded_ref_tvalue != 0) stop("rounded_ref_tvalue != 0 for ", name, ". Check that Sam...")
      removed_cols <- getColsToRemove(ref_value_pos, random_params$order)
      As[[length(As)+1]] <- A[, -removed_cols]

      # for later # REFREF
      model$random[[name]]$model$extra$tref_value <- 0
      model$random[[name]]$model$extra$bin_tvalues <- bin_tvalues
      model$random[[name]]$model$extra$bin_values <- applyInvTransformation(rmodel, bin_tvalues)
      model$random[[name]]$model$extra$rounded_ref_value <- applyInvTransformation(rmodel, 0)
      model$random[[name]]$model$extra$rounded_ref_tvalue <- rounded_ref_tvalue
      model$random[[name]]$model$extra$ref_value_pos <- ref_value_pos
      model$random[[name]]$model$extra$removed_cols <- removed_cols

      if(random[[name]]$model$params$poly_degree > 0){
        model$random[[name]]$model$extra$bin_values_int <- stats::poly(bin_tvalues,
                                                                degree = model$random[[name]]$model$params$poly_degree,
                                                                raw = TRUE)
      }

    }else if(random_type == "integrated Wiener process"){

      list2env(random_params, envir = environment())
      if(is.null(region)) model$random[[name]]$model$params$region <- region <- range(u)
      knots <- splitKnots(region=region, range = ran, ref_value=ref_value, stepsize=stepsize,
                          type = "iwp", extra_left=extra_left, extra_right=extra_right)
      model$random[[name]]$model$params$knots <- knots
      ref_pos <- which(knots == ref_value)

      # should not happen
      if(length(ref_pos) == 0) stop("ref_value of", name, "cannot be found in the corresponding knots vector. \n")
      if(ran[1] < knots[1] & ran[2] > rev(knots)[1]) warning("knots for ", name, " do not span its range. Continuing anyway. \n")

      new_knots <- applyTransformation(rmodel, knots)
      model$random[[name]]$model$extra$tknots <- new_knots
      A <- methods::as(local_poly(knots = new_knots, refined_x = u, p = random_params$order), "sparseMatrix")
      As[[length(As)+1]] <- A

    }else if(random_type == "monotone Gaussian process"){

      list2env(random_params, envir = environment())
      a <- ifelse(lambda == 0, 1, 1/(1-lambda))

      # (this is always zero now)
      # model$random[[name]]$model$extra$tref_value <- applyTransformation(rmodel, ref_value)


      As[[length(As)+1]] <- if(method == "SS"){
        stop("Method State-Space (SS) not implemented (omitted when reference value was included)")
        Diagonal(1, n = 2*length(u))[seq(1, 2*length(u), by = 2),]

      }else if(method == "FEM"){

        # For mgp, we dont need to transform the data
        u <- U[,name]

        # position of knots (split for left or right of ref value)
        if(is.null(region)) model$random[[name]]$model$params$region <- region <- range(u)
        splitK <- splitKnots(region=region, range = ran, ref_value=ref_value, stepsize=stepsize,
                             type = "mgp", extra_left=extra_left, extra_right=extra_right)

        model$random[[name]]$model$extra$region_split <- splitK$regions
        model$random[[name]]$model$extra$knots_split <- splitK$knots
        model$random[[name]]$model$extra$knots <- sum(unlist(splitK$knots))
        # model$random[[name]]$model$extra$oregion <- range(u)

        # forward and backward design matrices
        u_pos <- pmax(u - ref_value, 0); u_neg <- pmax(ref_value - u, 0)
        A_pos <- A_neg <- NULL
        if(any(u_neg != 0)) A_neg <- Compute_Design(x = u_neg, k = splitK$knots$neg, region = splitK$regions$neg) |> as("dgTMatrix")
        if(any(u_pos != 0)) A_pos <- Compute_Design(x = u_pos, k = splitK$knots$pos, region = splitK$regions$pos) |> as("dgTMatrix")
        Matrix::cbind2(A_neg[, ncol(A_neg):1], A_pos)
      }

    }else{
      stop("model type (", rmodel$type, ") for random effect ", name, " is not valid.")
    }
  }

  names(As) <- random_names
  Xs_int <- interpolationFixedEffects(model$random, U)
  gamma_dims <- sapply(As, ncol)
  gamma_dims <- gamma_dims[gamma_dims != 0]

  list(As = As, gamma_dims = gamma_dims, Xs_int = Xs_int, model = model)
}
#


interpolationFixedEffects <-  function(random, U){

  random_names <- names(random)
  sapply(random_names, function(name) {
    if ("poly_degree" %in% names(random[[name]]$model$params)) {

      poly_degree <- random[[name]]$model$params$poly_degree
      if(poly_degree == 0) return(matrix(nrow=nrow(U), ncol=0))

      # no point in using "buffer" c here. But I do to be coherent with Ziang. REFREF
      # if(random[[name]]$model$params$c > 0) random[[name]]$model$params$c <- 0
      u <- applyTransformation(random[[name]]$model, U[,name])

      X_new <- stats::poly(u, degree = poly_degree, raw = TRUE)
      colnames(X_new) <- paste0(name, "__", attr(X_new, "degree"))
      return(X_new)
    }
    else {
      return(matrix(nrow=nrow(U), ncol=0))
    }
  }, simplify = FALSE, USE.NAMES = TRUE)
}
#


# Identify where to remove (not to put) overdispersion terms (sort of design for OD)
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
    }) |> unlist()
  }

  return(z_rem)
}





# Precision matrices for random effects -----------------------------------
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
    diag(compute_weights_precision(knots=ran$model$extra$tknots))
  }))
}


constructQ_mgp <- function(random, U){

  random_types <- sapply(random, \(ran) ran$model$type)
  random_mgp <- random[random_types == "monotone Gaussian process"]
  if(length(random_mgp) == 0) return(list(Q_mgp = methods::as(matrix(nrow=0,ncol=0), "dgTMatrix"),
                                     log_det_Q_mgp = 0))


  Qs <- lapply(names(random_mgp), function(nam){

    ran_extras <- random_mgp[[nam]]$model$extra
    ran_params <- random_mgp[[nam]]$model$params
    list2env(ran_extras, envir = environment())
    list2env(ran_params, envir = environment())
    a <- ifelse(lambda == 0, 1, 1/(1-lambda))

    if(method == "SS"){
      stop("method SS for mGP random effects was dropped when reference values were included.")
      oo <- order(U[,nam])
      ooo <- rep(2*oo, each=2) - rep(1:0, times=nrow(U))
      Q <- Matrix(data = 0, nrow = 2*nrow(U)+1, ncol = 2*nrow(U), sparse=T)[-1,] |> as("dgTMatrix")
      # stop("SS option needs to be fixed")
      # below there is something to manage with the data/reference value
      # It assumed here that tdata = t(u - ref_value), but its t(u) - t(ref_value)
      Q[ooo,ooo] <- mGP_joint_prec(t_vec = U[oo,nam] - ref_value, alpha = a, c = c) # REFREF

    }else if(method == "FEM"){
      # if region was null, it has been set in createRandomDesigns

      # No need to transform data form mgp
      u <- U[,nam]
      u_pos <- pmax(u - ref_value, 0); u_neg <- pmax(ref_value - u, 0)
      tc <- c + ref_value

      # Define B and penalty matrices based on non-zero regions for training data
      Q_list <- list()
      if(any(u_neg != 0)) Q_list[[1]] <- Compute_Prec(k=knots_split$neg, region = region_split$neg, a = a, c = tc, rev = TRUE)[(knots_split$neg-2):1,(knots_split$neg-2):1]
      if(any(u_pos != 0)) Q_list[[length(Q_list)+1]] <- Compute_Prec(k=knots_split$pos, region = region_split$pos, a = a, c = tc, rev = FALSE)
      Q <- Matrix::bdiag(Q_list)

    }else{
      stop("Unknown method (", method, ") for mGP for covariate ", nam)
    }

    return(Q)
  })

  log_det <- sapply(Qs, \(Q) Matrix::determinant(Q)$modulus)
  Q_mgp <- methods::as(methods::as(Matrix::bdiag(Qs), "generalMatrix"), "TsparseMatrix")
  return(list(Q_mgp=Q_mgp, log_det_Q_mgp=log_det))
}





# Initialisation for priors  ----------------------------------------------
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
