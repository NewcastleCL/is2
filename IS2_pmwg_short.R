# PMWGS test script
library(is2)

load(here::here("scratch", "forstmann_long.RData"))
load(here::here("scratch", "IS2_pmwg_oldscript.RData"), envir = oldrun <- new.env())
message("Setup")
n_samples <- 100 # number of importance samples
n_particles <- 10 # number of particles

isamples <- is2(sampled, n_samples, n_particles, precomputed_mix = oldrun$mix)

message("Get Maximum Likelihood and Bootstrap for Standard error")

summary_like <- summary(isamples)
print(summary_like)

save.image(here::here("scratch", "IS2_pmwg_testgetlogp.RData"))

plot(sort(oldrun$tmp))
points(sort(isamples$samples) + 1, col = "red")
