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
  exp_delta <- exp( rowSums(as.matrix(delta_mat[,all_covariates])) + (gen_pars$overdispersion[[1]] != "none")*rnorm(n,0,1/sqrt(exp(theta_z))) )
  rpois(n, exp_delta)
}


# Generate data -----------------------------------------------------------

gen_pars <- fread("~/Git/air-pollution-models/R/simulations/simulation-1/generation_grid.csv")[6,]
gen_pars$nonlinear_covariates <- strsplit(gen_pars$nonlinear_covariates,"|",fixed=T)

# load delta matrix
delta_mat <- readRDS("~/Git/air-pollution-models/R/simulations/simulation-1/PTS-INLA-1/delta_mat.rds")
theta_z <- 1
data <- readRDS("~/Git/air-pollution-models/R/simulations/simulation-1/PTS-INLA-1/data.rds")
data$count <- generateCounts(gen_pars, delta_mat, theta_z)

model <- ccModel(response = "count",
                 time_index = "date",
                 fixed = list("hum_mean" = fixedEffect(gaussian_prior())),
                 random = list("o3_lag" = randomEffect(rw_effect(), pc_prep_prior()),
                               "temp_lag" = randomEffect(rw_effect(), pc_prep_prior())),
                 overdispersion = randomEffect(gaussian_effect(), pc_prep_prior(u = 10, alpha = .9)),
                 design = ccDesign(),
                 control_aghq = controlAGHQ())

fit <- fitModel(model, data)


response = "count"
time_index = "date"
fixed = list("hum_mean" = fixedEffect(gaussian_prior()))
random = list("o3_lag" = randomEffect(rw_effect(), pc_prep_prior()),
              "temp_lag" = randomEffect(rw_effect(), pc_prep_prior()))
overdispersion = T
design = ccDesign()
control_aghq = controlAGHQ()

fitModel(model, data)
