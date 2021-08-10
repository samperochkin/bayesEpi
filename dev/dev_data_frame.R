devtools::load_all()

# Data --------------------------------------------------------------------
attach("../../simulations-air-pollution-1/data/AllDatabase1.RData")
data <- as.data.frame(AllDatabase)
detach(pos = grep("data/AllDatabase", search()))

city0 <- "Toronto"
group0 <- "morb"
cause0 <- "all"

data <- data[data$city == city0 & data$group == group0 & data$cause == cause0,]

# Model -------------------------------------------------------------------
model <- ccModel(response = "count",
                 time_index = "time",
                 fixed = list("hum_mean" = fixedEffect(gaussian_prior())),
                 random = list("o3lag" = randomEffect(rw_effect(), pc_prep_prior()),
                               "temp4lag" = randomEffect(rw_effect(), pc_prep_prior())),
                 overdispersion = randomEffect(gaussian_effect(), pc_prep_prior()),
                 design = ccDesign(),
                 control_aghq = controlAGHQ())

model$fixed[["hum_mean"]]
model$random[["o3lag"]]
model$random[["temp4lag"]]

fit <- bayesEpi:::fitModel.ccModel(model, data)
