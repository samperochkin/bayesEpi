# library(bayesEpi)
load_all()

# Packages ----------------------------------------------------------------
library(data.table)
library(RhpcBLASctl)
blas_set_num_threads(1)


# Functions ---------------------------------------------------------------

generateCounts <- function(gen_pars, delta_mat, theta_z, with_intercept = T){
  all_covariates <- c(unlist(gen_pars$linear_covariates),
                      unlist(gen_pars$nonlinear_covariates))

  fac_cov <- unlist(gen_pars$factor_covariates)
  all_covariates <- c(all_covariates,
                      unlist(lapply(fac_cov, function(x) grep(x, colnames(delta_mat), value=T))))

  if(gen_pars$time_trend[[1]] != "none") all_covariates <- c(all_covariates, gen_pars$time_trend[[1]])
  if(gen_pars$seasonality[[1]] != "none") all_covariates <- c(all_covariates, grep(gen_pars$seasonality[[1]], colnames(delta_mat), value=T))
  if(with_intercept) all_covariates <- c(all_covariates, "(Intercept)")

  n <- nrow(delta_mat)
  exp_delta <- exp( rowSums(as.matrix(delta_mat[,all_covariates])) + (gen_pars$overdispersion[[1]] != "none")*rnorm(n,0,1/sqrt(theta_z)) )
  rpois(n, exp_delta)
}


# Generate data -----------------------------------------------------------

gen_pars <- fread("~/Git/air-pollution-models/R/simulations/simulation-1/generation_grid.csv")[6,]
gen_pars$nonlinear_covariates <- strsplit(gen_pars$nonlinear_covariates,"|",fixed=T)

# load delta matrix
delta_mat <- readRDS("~/Git/air-pollution-models/R/simulations/simulation-1/PTS-INLA-1/delta_mat.rds")
theta_z <- 4
data <- readRDS("~/Git/air-pollution-models/R/simulations/simulation-1/PTS-INLA-1/data.rds")
data$count <- generateCounts(gen_pars, delta_mat, theta_z)

model <- ccModel(response = "count",
                 time_index = "date",
                 fixed = list("hum_mean" = fixedEffect(gaussian_prior())),
                 random = list("o3_lag" = randomEffect(rw_effect(), pc_prep_prior()),
                               "temp_lag" = randomEffect(rw_effect(), gamma_prior(1,1))),
                 overdispersion = randomEffect(gaussian_effect(), pc_prep_prior()),
                 design = ccDesign(),
                 control_aghq = controlAGHQ())


###############
###############
###############
###############
library(TMB)
compile("templates/cc_debug.cpp")
dyn.load(dynlib("templates/cc_debug"))

data <- as.data.frame(data)
# contruct and import case_day and control_days. modifies and import data.
list2env(getCaseControl(data, model), envir = environment())
list2env(setRefValues(data, model), envir = environment())

# overdispersion, names and design matrices
overdispersion <- !is.null(model$overdispersion)
fixed_names <- names(model$fixed); random_names <- names(model$random)
X <- as.matrix(data[fixed_names]); U <- as.matrix(data[random_names])

# discretize nonlinear random effects add polynomial interpolation around reference values
# create As and Xs_int
list2env(discretizedRandomDesigns(model, U), envir = environment())

# prior parameters
prior_lookup <- c("pc_prep", "gamma")
beta_prec = c(purrr::map(model$fixed, ~ .x$prior$params$prec), purrr::map(model$random, ~ .x$beta_prior$params$prec)) %>% unlist
theta_prior <- c(purrr::map(model$random, ~ .x$theta_prior$type), z = model$overdispersion$theta_prior$type) %>% unlist
theta_prior_id = match(theta_prior , prior_lookup)
theta_hypers = c(purrr::map(model$random, ~ .x$theta_prior$params), z = model$overdispersion$theta_prior$params) %>% unlist


# Model fit ---------------------------------------------------------------
Q <- constructQ(model$random)

tmb_data <- list(count = data[case_day, model$response],
                 case_day = case_day, control_days = control_days,
                 X = cbind(X,Reduce("cbind", Xs_int)), A = Reduce("cbind", As),
                 Q = Q, gamma_dims = sapply(As, ncol),
                 beta_prec = beta_prec, theta_prior_id = theta_prior_id, theta_hypers = theta_hypers,
                 block_to_run = c(0,0,1,0,0))

theta_init <- rep(1,length(model$random)+overdispersion)
parameters <- list(beta = rep(0,ncol(X)+sum(sapply(Xs_int, ncol))),
                   gamma = rep(0, sum(sapply(As, ncol))),
                   z = rep(0,nrow(data)*overdispersion),
                   theta = theta_init)


dyn.load(TMB::dynlib("templates/cc_debug"))
obj <- TMB::MakeADFun(tmb_data, parameters, random = c("beta","gamma","z"), DLL="cc_debug", hessian=T)
obj
# quad <- aghq::marginal_laplace_tmb(obj, model$control$aghq_k, theta_init)
#
# quad
