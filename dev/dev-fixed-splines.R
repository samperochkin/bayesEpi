devtools::load_all("../bayesEpi")


# Generating parameters and synthetic data --------------------------------

# number of observations and time index
n <- (365*1)
time <- 1:n

# reference value and associated random effect "x"
rv <- 10
f <- function(x){
  (x - rv)^2/2000 + sin(2*pi*(x - rv)/30)/log(2+x)
}
plot(seq(0,10,.01), f(seq(0,10,.01)), type="o", cex=.5, pch=19)
plot(f(1:75), type="o", cex=.5, pch=19)
plot(f(1:55), type="o", cex=.5, pch=19)

set.seed(1)
# x <- pmax(10*sin(2*pi*time/365.25) + abs(cumsum(rnorm(n,0,.1))) + rnorm(n,0,.1), 0)
x <- pmin(50, 25*sin(2*pi*time/365.25) + 25 + abs(rnorm(n,.1)))
sum(x == 0)
table(round(x,1))
ref_values <- list(x = rv)
values <- list(x = seq(0, round(max(x))+1, .01))
plot(x)

# overdipsersion and response
theta_z <- 0
z <- rnorm(n,0,exp(-theta_z/2))
y <- rpois(length(x), exp(3+f(x)+z))
plot(y)
table(y)

# data and "true" effect
data <- data.frame(time = time, x = x, y = y)
range(data$x)
range(data$y)
true_frame <- data.frame(x = values$x, y = f(values$x))








# Model definition, fit and collection of results -------------------------
# knots <- seq(0,max(data$x)+5,2)
knots <- seq(0,50,5)
model <- bayesEpi::ccModel(response = "y",
                           time_index = "time",
                           fixed = list(x = fixedEffect(bs_effect(knots = knots, ref_value = ref_values$x))),
                           random = NULL,
                           design = ccDesign(scheme = "time stratified", n_control = 11),
                           # overdispersion = randomEffect(gaussian_effect(), pc_prec_prior(u = .01)),
                           overdispersion = NULL,
                           aghq_input = aghqInput())


# debug(fitModel.ccModel)
fit <- fitModel(model, data)

# res <- getResults(fit, M = 1e4, probs = c(.025, .5, .975))
# debug(getResults_bs)
res <- getResults(fit, M = 1e4, probs = c(.1, .5, .9))

plot(res$variable_value, res$mean, ylim = range(c(res$mean, res$perc_10, res$perc_90, true_frame$y)), type="o")
lines(res$variable_value, res$perc_10, type="o", col=2)
lines(res$variable_value, res$perc_90, type="o", col=2)
lines(true_frame$x, true_frame$y, type="l", col=3)

M <- 5
samps <- aghq::sample_marginal(fit$quad, M)$samps
model <- fit$model
# grids <- list("x" = seq(0,10,2))
grids <- NULL

new_samps <- combineBetaGamma_iwp(samps, model, grids)
head(new_samps)


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

new_samps_long <- reshape(new_samps, varying = paste0("samp_", 1:M), v.names = "sample",
                          idvar = "variable_value", timevar = "sample_id", direction = "long")
gg <- ggplot(new_samps_long, aes(x=variable_value, y=sample, col=as.factor(sample_id))) +
  theme_light() +
  geom_line() +
  facet_wrap(~variable_name)
gg

