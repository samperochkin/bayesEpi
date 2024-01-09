#' Extract results from a fitted ccModel.
#' @param fit A fitted ccModel (obtained via `fitModel(model)`).
#' @param probs_pw Nominal levels for computing pointwise credible intervals for fixed and random effects (and overdispersion terms).
#' @param probs_g Nominal levels for computing global envelopes for fixed and random effects (and overdispersion terms).
#' @param M Number of samples used for inference.
#' @param values List of vectors containing the values at which to evaluate the splines (fixed effect), if included.
#' @return A data.frame with values of the beta and gamma coefficients, as well as their credible intervals.
#' @export
getResults <- function(fit, probs_pw = c(.8, .9), probs_g = c(.8, .9), M = 1e4, values = NULL){

  # Two main possibilities:
  #  - fixed effects only (no overdispersion)
  #  - with random effects (incl. overdispersion)

  if(is.null(fit$model$random) & is.null(fit$model$overdispersion)){
    # fixed effects only
    return(getResults_fixed_only(fit, probs_pw, probs_g, M, values))
  }else{
    # with random effects
    return(getResults_general(fit, probs_pw, probs_g, M, values))
  }
}




# Helpers
# pointwise cred. int and global env funs
computePCI_fixed_only <- function(mean, sd, probs_pw){
  probs_pw <- sort(1-probs_pw)/2
  probs_pw <- c(probs_pw, 1-rev(probs_pw))

  quants <- sapply(probs_pw, qnorm, mean=mean, sd=sd)
  colnames(quants) <- paste0("perc_", probs_pw*100)
  rownames(quants) <- NULL
  quants
}

computeGE_fixed_only <- function(mean, sigma, probs_g, M){
  probs_g <- sort(probs_g)
  samps <- mvtnorm::rmvnorm(M, mean=mean, sigma = sigma)
  cset <- GET::create_curve_set(list(r = 1:ncol(samps), obs = t(samps)))
  credible_int <- GET::central_region(cset, type = "erl", coverage = probs_g, alternative = c("two.sided"),
                                      quantile.type = 7, central = "median", nstep = 2)

  cols <- c(grep("lo", names(credible_int), value = T), rev(grep("hi", names(credible_int), value = T)))
  quants <- credible_int[cols] |> do.call(what = "cbind")

  colnames(quants) <- c(paste0("ge_lo_", rev(probs_g)*100), paste0("ge_hi_", probs_g*100))
  quants
}

computePCI_general <- function(yy, probs_pw){
  if(!is.matrix(yy)) yy <- matrix(yy, nrow=1)
  probs_pw <- sort(1-probs_pw)/2
  probs_pw <- c(probs_pw, 1-rev(probs_pw))

  quants <- t(matrix(c(apply(yy, 1, quantile, probs = probs_pw)), nrow = length(probs_pw)))
  colnames(quants) <- paste0("perc_", probs_pw*100)
  rownames(quants) <- NULL
  quants
}

computeGE_general <- function(yy, probs_g){
  if(!is.matrix(yy)) yy <- matrix(yy, nrow=1)
  probs_g <- sort(probs_g)
  cset <- GET::create_curve_set(list(r = 1:nrow(yy), obs = yy))
  credible_int <- GET::central_region(cset, type = "erl", coverage = probs_g, alternative = c("two.sided"),
                                      quantile.type = 7, central = "median", nstep = 2)

  cols <- c(grep("lo", names(credible_int), value = T), rev(grep("hi", names(credible_int), value = T)))
  quants <- credible_int[cols] |> do.call(what = "cbind")

  colnames(quants) <- c(paste0("ge_lo_", rev(probs_g)*100), paste0("ge_hi_", probs_g*100))
  quants
}

getResults_fixed_only <- function(fit, probs_pw, probs_g, M, values = NULL){


  envv <- fit$obj$env
  model <- fit$model
  fixed <- model$fixed

  lpb <- envv$last.par.best
  H <- envv$spHess(lpb)
  S <- solve(H)

  counter_fixed <- 0
  df <- NULL
  for(nam in names(fixed)){

    if(fixed[[nam]]$model$params$degree == 1){
      counter_fixed <- counter_fixed + 1
      next
    }

    # set values at which to evaluate the effect
    if(is.null(values[[nam]])) stop("Please provide values for ", nam)
    xx <- values[[nam]]

    # design matrix (call bayesEpi function to create it)
    XX <- matrix(xx, ncol=1)
    colnames(XX) <- nam
    fake_model <- model
    fake_model$fixed <- fake_model$fixed[nam]
    XX <- bayesEpi:::createFixedDesigns(fake_model, XX)$Xs_exp[[nam]]

    # evaluate effect at xx and compute corresponding variance matrix
    yy <- c(XX %*% lpb[counter_fixed + 1:ncol(XX)])
    SS <- XX %*% S[counter_fixed + 1:ncol(XX), counter_fixed + 1:ncol(XX)] %*% t(XX)

    df0 <- data.frame(parameter_type = as.factor("beta*"),
                      variable_name = as.factor(nam),
                      variable_value = xx,
                      mean = yy,
                      median = yy,
                      sd = sqrt(diag(SS)))


    # include pointwise coverage probs and global envelop
    if(!is.null(probs_pw)) df0 <- cbind(df0, computePCI_fixed_only(df0$mean, df0$sd, probs_pw))
    # if(!is.null(probs_g)) df0 <- cbind(df0, computeGE_fixed_only(df0$mean, as.matrix(SS), probs_g, M))
    if(!is.null(probs_g)) df0 <- cbind(df0, computeGE_fixed_only(df0$mean, as.matrix(SS), probs_g, M))

    df <- rbind(df, df0)
    counter_fixed <- counter_fixed + ncol(XX)
  }

  return(df)
}

getResults_general <- function(fit, probs_pw, probs_g, M, values){

  quad_samples <- aghq::sample_marginal(fit$quad, M)
  model <- fit$model
  fixed <- model$fixed
  random <- model$random

  df <- NULL

  #####################
  ### FIXED EFFECTS ###
  #####################
  counter_fixed <- 0
  for(nam in names(fixed)){

    if(fixed[[nam]]$model$params$degree == 1){
      counter_fixed <- counter_fixed + 1
      next
    }

    if(is.null(values[[nam]])) stop("Please provide values for ", nam)
    xx <- values[[nam]]

    # design matrix (call bayesEpi function to create it)
    XX <- matrix(xx, ncol=1, dimnames = list(NULL, nam))
    fake_model <- model
    fake_model$fixed <- fake_model$fixed[nam]
    XX <- bayesEpi:::createFixedDesigns(fake_model, XX)$Xs_exp[[nam]]

    # evaluate effect at xx
    yy <- XX %*% quad_samples$samps[counter_fixed + 1:ncol(XX),]

    df0 <- data.frame(parameter_type = as.factor("beta*"),
                      variable_name = as.factor(nam),
                      variable_value = xx,
                      mean = rowMeans(yy),
                      median = apply(yy,1,median),
                      sd = apply(yy,1,sd))

    # include pointwise coverage probs and global envelop
    if(!is.null(probs_pw)) df0 <- cbind(df0, computePCI_general(yy, probs_pw))
    if(!is.null(probs_g)) df0 <- cbind(df0, computeGE_general(yy, probs_g))

    df <- rbind(df, df0)
    counter_fixed <- counter_fixed + ncol(XX)
  }




  ######################
  ### RANDOM EFFECTS ###
  ######################
  counter_random <- which(names(fit$obj$env$last.par.best) == "gamma")[1] - 1
  if(is.na(counter_random)) counter_random <- 0
  for(nam in names(random)){

    if(random[[nam]]$model$type != "random walk") stop("Invalid random effect type")

    model0 <- model
    model0$random <- model0$random[nam]
    uu <- model0$random[[nam]]$model$extra$bin_values # get values
    UU <- matrix(uu, ncol=1, dimnames = list(NULL, nam))

    rDesign <- bayesEpi:::createRandomDesigns(model0, UU)
    AA <- rDesign$As[[nam]]
    XX <- rDesign$Xs_int[[nam]]

    # evaluate effect at xx
    yy <- as.matrix(AA %*% quad_samples$samps[counter_random + 1:ncol(AA),]) +
      XX %*% quad_samples$samps[counter_fixed + 1:ncol(XX),]

    df0 <- data.frame(parameter_type = as.factor("gamma*"),
                      variable_name = as.factor(nam),
                      variable_value = uu,
                      mean = rowMeans(yy),
                      median = apply(yy,1,median),
                      sd = apply(yy,1,sd))

    # include pointwise coverage probs and global envelop
    if(!is.null(probs_pw)) df0 <- cbind(df0, computePCI_general(yy, probs_pw))
    if(!is.null(probs_g)) df0 <- cbind(df0, computeGE_general(yy, probs_g))

    df <- rbind(df, df0)
    counter_random <- counter_random + ncol(AA)
    counter_fixed <- counter_fixed + ncol(XX)
  }


  ##############
  ### THETAS ###
  ##############
  variable_name <- names(random)
  if(!is.null(model$overdispersion)) variable_name <- c(variable_name, "z")

  if(length(variable_name) > 0){

    df0 <- data.frame(parameter_type = as.factor("theta"),
                      variable_name = as.factor(variable_name),
                      variable_value = as.numeric(NA),
                      mean = sapply(quad_samples$thetasamples, mean),
                      median = sapply(quad_samples$thetasamples, median),
                      sd = sapply(quad_samples$thetasamples, sd))

    # include pointwise coverage probs and global envelop
    if(!is.null(probs_pw)) df0 <- cbind(df0, lapply(quad_samples$thetasamples, computePCI_general, probs_pw=probs_pw) |> do.call(what = "rbind"))
    if(!is.null(probs_g)) df0 <- cbind(df0, lapply(quad_samples$thetasamples, computeGE_general, probs_g=probs_g) |> do.call(what = "rbind"))

    df <- rbind(df, df0)
  }
}





#' Extract (overdispersion) results from a fitted ccModel.
#' @param fit A fitted ccModel (obtained via `fitModel(model)`).
#' @param probs_pw Nominal levels for computing pointwise credible intervals for fixed and random effects (and overdispersion terms).
#' @param probs_g Nominal levels for computing global envelopes for fixed and random effects (and overdispersion terms).
#' @param M Number of samples used for inference.
#' @return A data.frame with values of the beta and gamma coefficients, as well as their credible intervals.
#' @export
getResults_z <- function(fit, probs_pw, probs_g, M, values){

  samps <- aghq::sample_marginal(fit$quad, M)$samps
  samps <- samps[rownames(samps) == "z",]

  df0 <- data.frame(parameter_type = as.factor("z"),
                    variable_name = fit$model$time_index,
                    variable_value = fit$time,
                    mean = rowMeans(samps),
                    median = apply(samps,1,median),
                    sd = apply(samps,1,sd))

    # include pointwise coverage probs and global envelop
    if(!is.null(probs_pw)) df0 <- cbind(df0, computePCI_general(samps, probs_pw))
    if(!is.null(probs_g)) df0 <- cbind(df0, computeGE_general(samps, probs_g))

    df0
}
