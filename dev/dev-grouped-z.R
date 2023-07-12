# install.packages("aghq")
# library(bayesEpi)
devtools::load_all()

# Packages ----------------------------------------------------------------
library(data.table)
library(RhpcBLASctl)
blas_set_num_threads(1)



# Dummy setup -------------------------------------------------------------

n <- 2000
x <- rnorm(n,.5)
xx <- zoo::rollmean(rnorm(n+2), k = 3, na.pad = F)
plot(xx, cex=.25, type="l")
dev.off()
z <- rnorm(n,0,.1)

f <- function(xx) sin(2*pi*xx)*xx^2
y <- rpois(n, 5*exp(.1 * x + f(xx) + z))
data <- data.frame(date = rep(as.Date(1:(n/5)), each=5), x = x, xx = xx, y = y)


model <- ccModel(response = "y",
                 time_index = "date",
                 # fixed = NULL,
                 fixed = list("x" = fixedEffect(poly_effect())),
                 random = list("xx" = randomEffect(rw_effect(), pc_prec_prior())),
                 overdispersion = randomEffect(gaussian_effect(), pc_prec_prior()),
                 od_stratum_vars = "date",
                 design = ccDesign(scheme = "time stratified"),
                 aghq_input = aghqInput())

# undebug(fitModel.ccModel)
# debug(createODDesign)
fit <- fitModel(model, data)
# undebug(fitModel.ccModel)

res <- aghq::sample_marginal(fit$quad, M=1000)
sum(rownames(res$samps) == "z")


