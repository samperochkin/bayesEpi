#' Extract results from a fitted ccModel.
#' @param fit A fitted ccModel (obtained via `fitModel(model)`).
#' @param quantiles Quantiles of the posterior (cumulative) distribution of fixed and random effects (and overdispersion terms).
#' @return A data.frame with values of the beta and gamma coefficients, as well as their credible intervals.
#' @export
getResults <- function(fit, probs = .5, M = 1e4, stepsizes = NULL, u_osplines = "knots+stepsize+u"){
  if(is.null(fit$model$random) & is.null(fit$model$overdispersion) & fit$model$fixed[[1]]$model$type == "poly") return(getResults_poly(fit, probs, stepsizes))
  if(is.null(fit$model$random) & is.null(fit$model$overdispersion) & fit$model$fixed[[1]]$model$type == "bs") return(getResults_bs(fit, probs, stepsizes))
  if(is.null(fit$model$random) & !is.null(fit$model$overdispersion)) return(getResults_od(fit, probs, M))
  if(fit$model$random[[1]]$model$type == "random walk") return(getResults_rw(fit, probs, M))
  if(fit$model$random[[1]]$model$type == "integrated Wiener process") return(getResults_iwp(fit, probs, M, stepsizes, u_osplines))
}


getResults_poly <- function(fit, probs = .5, stepsizes){

  if(is.null(stepsizes)) stepsizes <- 1
  lpb <- fit$obj$env$last.par.best
  H <- fit$obj$env$spHess(lpb)
  S <- solve(H)

  k_cum <- 0
  df <- NULL
  for(nam in names(fit$model$fixed)){
    ran <- fit$model$fixed[[nam]]$model$extra$range

    xx <- seq(ran[1], ran[2], stepsizes)
    XX <- as.matrix(xx, ncol=1)
    # XX <- poly(xx - fit$model$fixed[[nam]]$model$params$ref_value,
    #            degree=fit$model$fixed[[nam]]$model$params$degree)
    colnames(XX) <- nam
    XX <- createFixedDesigns(fit$model, XX)$Xs_exp[[nam]]

    yy <- c(XX %*% lpb[k_cum + 1:ncol(XX)])
    SS <- XX %*% S[k_cum + 1:ncol(XX), k_cum + 1:ncol(XX)] %*% t(XX)

    df0 <- data.frame(parameter_type = as.factor("beta*"),
                      variable_name = as.factor(nam),
                      variable_value = xx,
                      mean = yy,
                      median = yy,
                      sd = sqrt(diag(SS)))

    quants <- sapply(probs, qnorm, mean=df0$mean, sd=df0$sd)
    colnames(quants) <- paste0("perc_", probs*100)
    rownames(quants) <- NULL
    df0 <- cbind(df0, quants)


    df <- rbind(df, df0)
    k_cum <- k_cum + ncol(XX)
  }

  return(df)
}


getResults_bs <- function(fit, probs = .5, stepsizes){

  lpb <- fit$obj$env$last.par.best
  H <- fit$obj$env$spHess(lpb)
  S <- solve(H)

  k_cum <- 0
  df <- NULL
  for(nam in names(fit$model$fixed)){

    kn <- sort(unique(unlist(fit$model$fixed[[nam]]$model$params$knots)))
    xx <- seq(kn[1], kn[length(kn)], stepsizes)

    # OBSOLETE
    # if(is.null(stepsizes)){
    #   xx <- fit$model$fixed[[nam]]$model$params$knots
    # }else{
    #   kn <- fit$model$fixed[[nam]]$model$params$knots
    #   xx <- seq(kn[1], kn[length(kn)], stepsizes)
    # }

    XX <- matrix(xx, ncol=1)
    colnames(XX) <- nam
    XX <- createFixedDesigns(fit$model, XX)$Xs_exp[[nam]]

    yy <- c(XX %*% lpb[k_cum + 1:ncol(XX)])
    SS <- XX %*% S[k_cum + 1:ncol(XX), k_cum + 1:ncol(XX)] %*% t(XX)

    df0 <- data.frame(parameter_type = as.factor("beta*"),
                      variable_name = as.factor(nam),
                      variable_value = xx,
                      mean = yy,
                      median = yy,
                      sd = sqrt(diag(SS)))

    quants <- sapply(probs, qnorm, mean=df0$mean, sd=df0$sd)
    colnames(quants) <- paste0("perc_", probs*100)
    rownames(quants) <- NULL
    df0 <- cbind(df0, quants)


    df <- rbind(df, df0)
    k_cum <- k_cum + ncol(XX)
  }

  return(df)
}



getResults_od <- function(fit, probs = .5, M = 1e4, stepsizes){

  # most of what we need
  quad_samples <- aghq::sample_marginal(fit$quad, M)
  model <- fit$model


  k_cum <- 0
  df <- NULL
  for(nam in names(fit$model$fixed)){

    if(fit$model$fixed[[nam]]$model$type == "bs"){
      xx <- unlist(fit$model$fixed[[nam]]$model$params$knots
      bb <- "beta*"
      set_vals <- NA
    }else if(model$fixed[[nam]]$model$type == "poly" && model$fixed[[nam]]$model$params$degree > 1){
      ran <- model$fixed[[nam]]$model$extra$range
      xx <- seq(ran[1], ran[2], length.out = 100) # this should be taken care of
      bb <- "beta*"
      set_vals <- 1
    }else{
      xx <- 1 # this should be taken care of
      bb <- "beta"
      set_vals <- NA
    }
    XX <- as.matrix(xx, ncol=1)
    # XX <- poly(xx - fit$model$fixed[[nam]]$model$params$ref_value,
    #            degree=fit$model$fixed[[nam]]$model$params$degree)
    colnames(XX) <- nam
    sub_mod <- model
    sub_mod$fixed <- model$fixed[nam]
    XX <- createFixedDesigns(sub_mod, XX)$Xs_exp[[nam]]
    yy <- XX %*% quad_samples$samps[k_cum + 1:ncol(XX),]

    df0 <- data.frame(parameter_type = as.factor(bb),
                      variable_name = as.factor(nam),
                      variable_value = xx * set_vals,
                      mean = rowMeans(yy),
                      median = apply(yy,1,median),
                      sd = apply(yy,1,sd))

    quants <- t(matrix(c(apply(yy, 1, quantile, probs = probs)), nrow = length(probs)))
    colnames(quants) <- paste0("perc_", probs*100)
    rownames(quants) <- NULL
    df0 <- cbind(df0, quants)

    df <- rbind(df, df0)
    k_cum <- k_cum + ncol(XX)

  }

  #--- THETAS ---#

  variable_name <- c("z")

  res_frame_theta <- data.frame(parameter_type = as.factor("theta"),
                                variable_name = as.factor(variable_name),
                                variable_value = as.numeric(NA),
                                mean = sapply(quad_samples$thetasamples, mean),
                                median = sapply(quad_samples$thetasamples, median),
                                sd = sapply(quad_samples$thetasamples, sd))

  quants <- do.call("rbind", lapply(quad_samples$thetasamples, function(samp) t(quantile(samp, probs = probs))))
  colnames(quants) <- paste0("perc_", probs*100)

  df <- rbind(df, cbind(res_frame_theta, quants))
  df <- df[order(df$parameter_type, df$variable_name, df$variable_value),]

  return(df)
}




getResults_rw <- function(fit, probs = .5, M = 1e4){

  # most of what we need
  quad_samples <- aghq::sample_marginal(fit$quad, M)
  model <- fit$model

  # ###############
  # #--- BETAS ---#
  # ###############
  # beta_fixed_frame <- data.frame(parameter_type = "beta",
  #                                variable_name = names(fit$model$fixed),
  #                                bin_value = NA,
  #                                power = 1)
  #
  # beta_poly_frame <- data.frame(parameter_type = "beta",
  #                                variable_name = unlist(lapply(names(fit$model$random), function(nam){
  #                                  rep(nam, model$random[[nam]]$model$params$poly_degree)
  #                                })),
  #                                bin_value = NA,
  #                                power = unlist(lapply(names(fit$model$random), function(nam){
  #                                  seq(1, model$random[[nam]]$model$params$poly_degree,1)
  #                                })))
  # beta_poly_ids <- sapply(random, function(ran){
  #   random[[nam]]$model$params$poly_degree)},
  #   simplify=F, )
  #
  # ################
  # #--- GAMMAS ---#
  # ################
  # bin_values <- lapply(model$random, function(ran) ran$model$extra$bin_values)
  # removed_cols <- lapply(model$random, function(ran) ran$model$extra$removed_cols)
  # gamma_x_values <- lapply(seq_along(bin_values), function(k) bin_values[[k]][-removed_cols[[k]]])
  # names(gamma_x_values) <- names(gamma_ghost_x_values) <- names(fit$model$random)
  #
  # gamma_frame <- data.frame(parameter_type = "gamma",
  #                           variable_name = unlist(lapply(names(fit$model$random), function(nam){
  #                             rep(nam, length(gamma_x_values[[nam]]))
  #                           })),
  #                           bin_value = unlist(gamma_x_values),
  #                           power = NA)
  # parameter_frame <- rbind(beta_fixed_frame, beta_poly_frame, gamma_frame)
  #
  # quantiles_block <- t(apply(quad_samples$samps, 1, quantile, probs = probs))
  # colnames(quantiles_block) <- paste0("perc_", probs*100)
  #
  #
  # ##########################################
  # #--- results involving raw parameters ---#
  # ##########################################
  # parameter_frame <- cbind(parameter_frame,
  #                          data.frame(mean = rowMeans(quad_samples$samps),
  #                                     sd = apply(quad_samples$samps, 1, sd)),
  #                          quantiles_block)
  #
  #
  # #############################################################
  # #--- results involving linear combinations of parameters ---#
  # #############################################################
  #
  # for(name in names(model$random)){
  #   X_int <- model$random[[name]]$model$extra$bin_values_int
  #   A <- diag(nrow(X_int))[,-model$random[[name]]$model$extra$removed_cols]
  #   XA <- as(cbind(X_int,A),"sparseMatrix")
  #
  #   XA %*% c(beta_fixed_ids)
  #
  # }
  #
  #
  # cum_len <- cumsum(c(m, sapply(gamma_x_values, length)))
  # gamma_ids <- lapply(seq_along(cum_len[-1]), function(k){
  #   seq(cum_len[k] + 1, cum_len[k+1], 1)
  # })
  #
  #
  #

  #### setup ####
  # ghost gammas are the gamma parameters set to zero by the algorithm
  bin_values <- lapply(model$random, function(ran) ran$model$extra$bin_values)
  removed_cols <- lapply(model$random, function(ran) ran$model$extra$removed_cols)
  gamma_x_values <- lapply(seq_along(bin_values), function(k) bin_values[[k]][-removed_cols[[k]]])
  ghost_gamma_x_values <- lapply(seq_along(bin_values), function(k) bin_values[[k]][removed_cols[[k]]])
  names(gamma_x_values) <- names(ghost_gamma_x_values) <- names(model$random)


  #### standard outputs ####
  poly_degrees <- sapply(model$random, function(ran) ran$model$params$poly_degree)
  variable_value <- c(rep(as.numeric(NA), length(model$fixed) + sum(poly_degrees)), unlist(gamma_x_values))
  variable_name <- c(names(model$fixed),
                     unlist(lapply(names(model$random), function(nam) rep(nam, poly_degrees[[nam]]))),
                     unlist(lapply(names(model$random), function(nam) rep(nam, length(gamma_x_values[[nam]])))))
  parameter_type <- c(rownames(quad_samples$samps))


  #### Add beta and gamma to get final output ####
  # Note that we do this to the ghost gammas as well, which is why we want them in the output...
  ghost_gamma_samps <- replicate(length(model$random),NULL)
  names(ghost_gamma_samps) <- names(model$random)

  for(nam in names(model$random)){


    # find corresponding gamma (ids) and variable_value
    is_gamma_type <- parameter_type[!grepl("z",parameter_type)] == "gamma" # excluding overdispersion
    gamma_id <- which(is_gamma_type & grepl(nam, x = variable_name))


    # ************************************************************************************************************** WATCH OUT
    if(model$random[[nam]]$model$params$poly_degree == 0){
      ghost_gamma_samps[[nam]] <- matrix(0, length(removed_cols[[nam]]), ncol(quad_samples$samps))

    }else{
      x_values <- c(model$random[[nam]]$model$extra$bin_values_int[-removed_cols[[nam]],1]) # variable_value[gamma_id] - model$random[[nam]]$model$params$ref_value

      # find corresponding betas (ids)
      is_beta_type <- parameter_type[!grepl("z",parameter_type)] == "beta"
      beta_id <- which(is_beta_type & grepl(nam, x = variable_name))

      # multiply beta_i with (x - x0)^i, where x0 is the reference value
      quad_samples$samps[gamma_id,] <- quad_samples$samps[gamma_id,] +
        lapply(seq_along(beta_id), function(i){
          tcrossprod(x_values^i, quad_samples$samp[beta_id[i],])
        }) %>% Reduce(f = "+")

      # do the same with "ghost" gammas
      x_values <- c(model$random[[nam]]$model$extra$bin_values_int[removed_cols[[nam]],1]) # x_values <- bin_values[[nam]][removed_cols[[nam]]] - model$random[[nam]]$model$params$ref_value
      ghost_gamma_samps[[nam]] <- lapply(seq_along(beta_id), function(i){
        tcrossprod(x_values^i, quad_samples$samp[beta_id[i],])
      }) %>% Reduce(f = "+")
    }

  }

  ghost_gamma_samps <- Reduce("rbind", ghost_gamma_samps)

  #### overdispersion output ####
  if(!is.null(model$overdispersion)){
    parameter_type[grepl("z", parameter_type)] <- "gamma" # for consistency with PTS-INLA models
    variable_name <- c(variable_name, rep("z", length(parameter_type) - length(variable_name)))
    variable_value <- c(variable_value, seq(1, length(parameter_type) - length(variable_value),1))
  }

  #### ghost gammas output ####
  # These are the gamma parameters set to zero by the algorithm. We still want them in the output
  parameter_type <- c(parameter_type, rep("gamma", length(unlist(ghost_gamma_x_values))))
  variable_name <- c(variable_name, unlist(lapply(names(model$random), function(nam) rep(nam, length(ghost_gamma_x_values[[nam]])))))
  variable_value <- c(variable_value, unlist(ghost_gamma_x_values))


  res_frame <- data.frame(parameter_type = as.factor(parameter_type),
                          variable_name = as.factor(variable_name),
                          variable_value = variable_value,
                          mean = rowMeans(rbind(quad_samples$samps, ghost_gamma_samps)),
                          median = apply(rbind(quad_samples$samps, ghost_gamma_samps), 1, median),
                          sd = apply(rbind(quad_samples$samps, ghost_gamma_samps), 1, sd))

  quants <- t(matrix(c(apply(rbind(quad_samples$samps, ghost_gamma_samps), 1, quantile, probs = probs)), nrow = length(probs)))
  colnames(quants) <- paste0("perc_", probs*100)
  rownames(quants) <- NULL

  res_frame <- cbind(res_frame, quants)



  #--- THETAS ---#

  variable_name <- names(model$random)
  if(!is.null(model$overdispersion)) variable_name <- c(variable_name, "z")

  res_frame_theta <- data.frame(parameter_type = as.factor("theta"),
                                variable_name = as.factor(variable_name),
                                variable_value = as.numeric(NA),
                                mean = sapply(quad_samples$thetasamples, mean),
                                median = sapply(quad_samples$thetasamples, median),
                                sd = sapply(quad_samples$thetasamples, sd))

  quants <- do.call("rbind", lapply(quad_samples$thetasamples, function(samp) t(quantile(samp, probs = probs))))
  colnames(quants) <- paste0("perc_", probs*100)

  res_frame <- rbind(res_frame, cbind(res_frame_theta, quants))
  res_frame <- res_frame[order(res_frame$parameter_type, res_frame$variable_name, res_frame$variable_value),]

  return(res_frame)
}






#' Extract results from a fitted ccModel with IWP random effects.
#' @param fit A fitted ccModel (obtained via `fitModel(model)`).
#' @param probs Quantiles of the posterior (cumulative) distribution of fixed and random effects (and overdispersion terms).
#' @param stepsizes Named vector specifying the (reported) stepsize for each random effect.
#' @return A data.frame with values of the beta and gamma coefficients, as well as their credible intervals.
#' @import OSplines
getResults_iwp <- function(fit, probs, M, stepsizes, u_osplines){

  # most of what we need
  quad_samples <- aghq::sample_marginal(fit$quad, M)
  model <- fit$model

  fixed_names <- names(model$fixed)
  random_names <- names(model$random)

  # points where to evaluate the o-splines (or list of method to use)
  # u_osplines should be a list (of numeric and/or string) or a single string
  if(!is.list(u_osplines)){
    u_osplines <- as.list(rep(u_osplines, length(random_names)))
    names(u_osplines) <- random_names
  }

  if(is.null(fixed_names)){
    fixed_df <- NULL
  }else{
    fixed_df <- data.frame(parameter_type = "beta",
                           variable_name = fixed_names,
                           variable_value = NA,
                           mean = rowMeans(quad_samples$samps[1:length(fixed_names),,drop=F]),
                           median = apply(quad_samples$samps[1:length(fixed_names),,drop=F], 1, median),
                           sd = apply(quad_samples$samps[1:length(fixed_names),,drop=F], 1, sd))

    quants <- t(matrix(c(apply(quad_samples$samps[1:length(fixed_names),,drop=F], 1, quantile, probs = probs)), nrow = length(probs)))
    colnames(quants) <- paste0("perc_", probs*100)
    rownames(quants) <- NULL
    fixed_df <- cbind(fixed_df, quants)
  }


  random_df <- NULL
  if(!is.null(random_names)){

    if(is.null(names(stepsizes))){
      cat("stepsizes is not specified (or ill-specified). Using 1 for all random effects.\n")
      stepsizes <- rep(1, length(random_names))
      names(stepsizes) <- random_names
    }

    counter_beta <- length(fixed_names)
    counter_gamma <- sum(rownames(quad_samples$samps) == "beta") - counter_beta

    for(nam in random_names){

      random_params <- model$random[[nam]]$model$params
      knots <- random_params$knots
      ref_value <- random_params$ref_value
      ran <- model$random[[nam]]$model$extra$range

      sz <- stepsizes[nam]
      if(is.numeric(u_osplines[[nam]])){
        u <- u_osplines[[nam]]
      }else if(u_osplines[[nam]] == "knots"){
        u <- knots
      }else if(u_osplines[[nam]] == "knots+stepsize"){
        u <- seq(knots[1], knots[length(knots)], stepsizes[nam])
      }else if(u_osplines[[nam]] == "u"){
        u <- seq(min(fit$U[,nam]*step), max(fit$U[,nam]), sz)
      }else if(u_osplines[[nam]] == "knots+stepsize+u"){
        u_ran <- c(knots[1], knots[length(knots)])

        d_lo <- (min(fit$U[,nam])-u_ran[1])/sz
        d_hi <- (max(fit$U[,nam])-u_ran[2])/sz
        if(d_lo < 0) u_ran[1] <- u_ran[1] + floor(d_lo)*sz
        if(d_hi > 0) u_ran[2] <- u_ran[2] + ceiling(d_hi)*sz

        u <- seq(u_ran[1], u_ran[2], stepsizes[nam])
      }else{
        stop("Not a valid u_splines parameter for", nam, "\n")
      }


      ref_pos <- which(knots == ref_value)
      A <- NULL
      if(ref_pos != 1) A <- cbind(A, OSplines::local_poly(knots = rev(ref_value - knots[1:ref_pos]),
                                                   refined_x = ref_value - u,
                                                   p = random_params$order))
      if(ref_pos != length(knots)) A <- cbind(A, OSplines::local_poly(knots = knots[ref_pos:length(knots)] - ref_value,
                                                                      refined_x = u - ref_value,
                                                                      p = random_params$order))

      id_gamma <- counter_gamma + 1:ncol(A)
      counter_gamma <- counter_gamma + ncol(A)
      y <- A %*% quad_samples$samps[id_gamma,,drop=F]
      if(random_params$poly_degree > 0){
        X <- poly(u - ref_value,
                  degree = random_params$poly_degree,
                  raw = TRUE)
        id_beta <- counter_beta + 1:ncol(X)
        counter_beta <- counter_beta + ncol(X)
        y <- y + X %*% quad_samples$samps[id_beta,,drop=F]
      }

      quants <- t(matrix(c(apply(y, 1, quantile, probs = probs)), nrow = length(probs)))
      colnames(quants) <- paste0("perc_", probs*100)
      rownames(quants) <- NULL

      random_df <- rbind(random_df,
                         cbind(data.frame(parameter_type = "gamma",
                                          variable_name = nam,
                                          variable_value = u,
                                          mean = rowMeans(y),
                                          median = apply(y, 1, median),
                                          sd = apply(y, 1, sd)),
                               quants))
    }
  }


  #--- THETAS ---#
  if(!is.null(model$overdispersion)) random_names <- c(random_names, "z")

  theta_df <- data.frame(parameter_type = as.factor("theta"),
                         variable_name = as.factor(random_names),
                         variable_value = as.numeric(NA),
                         mean = sapply(quad_samples$thetasamples, mean),
                         median = sapply(quad_samples$thetasamples, median),
                         sd = sapply(quad_samples$thetasamples, sd))

  quants <- do.call("rbind", lapply(quad_samples$thetasamples, function(samp) t(quantile(samp, probs = probs))))
  colnames(quants) <- paste0("perc_", probs*100)

  theta_df <- cbind(theta_df, quants)
  theta_df <- theta_df[order(theta_df$parameter_type, theta_df$variable_name, theta_df$variable_value),]

  return(rbind(fixed_df, random_df, theta_df))
}
