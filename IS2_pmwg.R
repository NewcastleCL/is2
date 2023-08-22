# PMWGS test script
# Install IS2
#  devtools::install_github("university-of-newcastle-research/is2", ref="develop")
#  devtools::install_github("NewcastleCL/is2", ref="develop")

library(is2)
#Load the sampler object from memory
load(here::here("scratch", "forstmann_long.RData"))
# IS2 parameters - not sure about good values for these
message("Setup")
importance_samples <- 100 # number of importance samples
n_particles <- 10 # number of particles

# Run importance sampling
importance_samples <- is2(sampled, importance_samples, n_particles)

message("Get Maximum Likelihood and Bootstrap for Standard error")

#Get importance sampling results
summary_like <- summarise(importance_samples)
print(summary_like)

#Save for other works
save.image(here::here("scratch", "IS2_pmwg3.RData"))
