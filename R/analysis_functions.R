#' Extract results from a fitted ccModel.
#' @param fit A fitted ccModel (obtained via `fitModel(model)`).
#' @param quantiles Quantiles of the posterior (cumulative) distribution of fixed and random effects (and overdispersion terms).
#' @return A data.frame with values of the beta and gamma coefficients, as well as their credible intervals.
#' @export
getResults <- function(fit, probs = .5, M = 1e4){

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

