# library(bayesEpi)
devtools::load_all("~/Git/packages/bayesEpi")
n <- 500
# toyData <- data.frame(time = seq(n))
toyData <- data.frame(time = as.Date(seq(n), origin = "1990-07-24"))

# explanatory variables and hazard
toyData$hum <- pmin(pmax((sin(as.integer(toyData$time)/365 * 2*pi) + 1)/2 + rnorm(n,0,.75),0),1)
toyData$o3 <- pmax((sin(as.integer(toyData$time)/365 * 2*pi) + 1)*25 + rnorm(n,0,10),0)
loghazard <- abs(1.25-toyData$hum)*1 + (toyData$o3/25)^2

# response variable
toyData$count <- rpois(n,exp(loghazard))

head(toyData)
model <- bayesEpi::ccModel(response = "count",
                           time_index = "time",
                           fixed = list("hum" = fixedEffect(gaussian_prior())),
                           random = list("o3" = randomEffect(iwp_effect(ref_value=10,
                                                                        knots=seq(0,ceiling(max(toyData$o3))-40,2)), pc_prec_prior())),
                           overdispersion = randomEffect(gaussian_effect(), pc_prec_prior()),
                           design = ccDesign(scheme = "time stratified"),
                           aghq_input = aghqInput())

fit <- fitModel(model, toyData)
res <- getResults(fit, probs = c(.05,.5,.95), M = 1000)

samps <- aghq::sample_marginal(fit$quad, 1e3)
z_means <- apply(samps[[1]][rownames(samps[[1]]) == "z",],1,mean)

case_day_ids <- fit$obj$env$data$case_day
control_day_ids <- unique(c(fit$obj$env$data$control_days))
if(control_day_ids[1] == 0) control_day_ids <- control_day_ids[-1]
day_ids <- sort(unique(c(case_day_ids, control_day_ids)))
day_ids <- day_ids[day_ids %in% 1:nrow(toyData)]

length(z_means)
length(toyData$o3)

plot(toyData$time[day_ids], z_means)
plot(toyData$o3[day_ids], z_means)

# res0 <- res[res$parameter_type == "gamma" & res$variable_name == "o3",]
res00 <- res[res$parameter_type == "gamma" & res$variable_name == "o3",]
plot(res0$variable_value, res0$median, type = "o", cex=.25)
lines(res00$variable_value, res00$median, type = "o", cex=.25)
