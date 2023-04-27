devtools::load_all("../bayesEpi")


# Generating parameters and synthetic data --------------------------------

# number of observations and time index
n <- 1000
time <- 1:n


# x <- pmax(0,pmin(50, 10*sin(2*pi*time/365.25) + 22.5 + rnorm(n,0,10)))
x <- pmax(0,pmin(50, rnorm(n,25,10)))
sum(x == 0)
table(round(x,1))
ref_values <- list(x = 20)
values <- list(x = seq(0, round(max(x)), .1))
plot(x)


effect_type <- "bs"
theta_z <- Inf
z <- rnorm(n,0,exp(-theta_z/2))

if(effect_type == "bs"){
  true_knots <- list(c(0,0,0,0,10,20), c(20,30,40,50,50,50))
  true_degree <- 3
  # B <- constructBS(x = x, knots = true_knots, degree = true_degree, ref_value = ref_values$x)
  B <- constructBS(x = x, knots = true_knots, degree = true_degree, ref_value = ref_values$x)

  beta <- sqrt(seq(0,50,length.out=ncol(B)-(true_degree-1)))/sqrt(50) - .25
  beta <- c((-1)^(1 + 1:(true_degree-1))*10^(-2*(1:(true_degree-1))), beta)
  plot(B %*% beta)

  # overdipsersion and response
  y <- rpois(length(x), exp(3 + B %*% beta + z))
  plot(y)
  table(y)

  # data and "true" effect
  data <- data.frame(time = time, x = x, y = y)
  range(data$x)
  range(data$y)

  # B <- constructBS(values$x, knots = true_knots, degree = true_degree, ref_value = ref_values$x)
  B <- constructBS(values$x, knots = true_knots, degree = true_degree)
  true_frame <- data.frame(x = values$x, y = B %*% beta)
  plot(true_frame$x, true_frame$y)
  #
  knots <- true_knots
  degree <- true_degree
  fixed = list(x = fixedEffect(bs_effect(knots = knots, degree = degree, ref_value = ref_values$x)))

}else{
  degree <- 2
  X <- poly(x - ref_values$x, degree = 2, raw = T)

  beta <- c(.05, -.001)
  plot(X %*% beta)

  # overdipsersion and response
  y <- rpois(length(x), exp(3 + X %*% beta + z))
  plot(y)
  table(y)

  # data and "true" effect
  data <- data.frame(time = time, x = x, y = y)
  range(data$x)
  range(data$y)

  X <- poly(values$x - ref_values$x, degree = 2, raw = T)
  true_frame <- data.frame(x = values$x, y = X %*% beta)
  plot(true_frame$x, true_frame$y)
  #
  fixed = list(x = fixedEffect(poly_effect(degree=2, ref_value = ref_values$x)))
}








# Model definition, fit and collection of results -------------------------
# knots <- seq(0,max(data$x)+5,2)
model <- bayesEpi::ccModel(response = "y",
                           time_index = "time",
                           fixed = fixed,
                           random = NULL,
                           design = ccDesign(scheme = "time stratified", n_control = 5),
                           # overdispersion = randomEffect(gaussian_effect(), pc_prec_prior(u = .05)),
                           overdispersion = NULL,
                           aghq_input = aghqInput())


# debug(fitModel.ccModel)
fit <- fitModel(model, data)

# res <- getResults(fit, M = 1e4, probs = c(.025, .5, .975))
# debug(getResults_bs)
res <- getResults(fit, M = 1e4, probs = c(.1, .5, .9), stepsizes = 1)
res[res$parameter_type == "theta",]
res <- res[!(res$parameter_type == "theta"),]

plot(res$variable_value, res$mean, ylim = range(c(res$mean, res$perc_10, res$perc_90, true_frame$y)), type="o", cex=.1)
lines(res$variable_value, res$perc_10, type="o", col=2, cex=.1)
lines(res$variable_value, res$perc_90, type="o", col=2, cex=.1)
lines(true_frame$x, true_frame$y, type="l", col=3, cex=.1)
#
