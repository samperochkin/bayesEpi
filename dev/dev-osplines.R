# devtools::load_all("../bayesEpi")
library(bayesEpi)
# library(OSplines)
# library(aghq)
# library(data.table)
library(ggplot2)


# Generating parameters and synthetic data --------------------------------

# number of observations and time index
n <- (365*4)
time <- 1:n

# reference value and associated random effect "x"
rv <- 10
f <- function(x){
  (x - rv)^2/1000 + sin(2*pi*(x - rv)/30)
}
plot(f(1:75), type="o", cex=.5, pch=19)

set.seed(1)
x <- pmax(10*sin(2*pi*time/365.25) + abs(cumsum(rnorm(n,0,1))) + rnorm(n,0,1), 0)
sum(x == 0)
table(round(x))
ref_values <- list(x = rv)
values <- list(x = seq(0, round(max(x))+1, 1))
plot(x)

# overdipsersion and response
theta_z <- 2
z <- rnorm(n,0,exp(-theta_z/2))
y <- rpois(length(x), exp(1+f(x)+z))
plot(y)

# data and "true" effect
data <- data.frame(time = time, x = x, y = y)
true_frame <- data.frame(x = values$x, y = f(values$x))




# Model definition, fit and collection of results -------------------------
knots <- seq(0,max(data$x)+5,5)
model <- bayesEpi::ccModel(response = "y",
                           time_index = "time",
                           fixed = NULL,
                           random = list("x" = randomEffect(iwp_effect(ref_value = rv, knots = knots),
                                                            pc_prec_prior(u=.05))),
                           design = ccDesign(scheme = "time stratified", n_control = 3),
                           overdispersion = randomEffect(gaussian_effect(), pc_prec_prior(u = .01)),
                           # overdispersion = NULL,
                           aghq_input = aghqInput())



fit <- fitModel(model, data)

fit$obj$env$data$case_day[c(1,7,11)]
(fit$obj$env$data$control_days)[c(1,7,11),]

res <- getResults(fit, probs = c(.025, .5, .975))

samps <- aghq::sample_marginal(fit$quad, 5)$samps
model <- fit$model

new_samps <- combineBetaGamma_iwp(samps, model)

plot(new_samps$x[,1], type = "l", ylim=range(new_samps$x))
for(k in 2:ncol(new_samps$x)) lines(new_samps$x[,k], col=k)


# Plots -------------------------------------------------------------------
gg <- ggplot(res[res$variable_name == "x" & res$parameter_type == "gamma",],
             aes(x=variable_value)) +
  theme_light() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "gray50"),
        strip.text.x = element_text(size = 11.5),
        strip.text.y = element_text(size = 11.5),
        axis.title=element_text(size=12)) +
  xlab("(synthetic) PM") +
  ylab("(log) effect") +
  geom_line(data=true_frame, aes(x = x, y = y), linetype=2) +
  geom_line(aes(y = median), col="gray") +
  geom_line(aes(y = perc_2.5), col = "coral") +
  geom_line(aes(y = perc_97.5), col = "lightblue")
gg

