load(here::here("scratch", "IS2_pmwg_oldscript.RData"), envir = oldrun <- new.env())

load(here::here("scratch", "IS2_pmwg_getlogp.RData"), envir = newrun <- new.env())

plot(sort(oldrun$tmp))
points(sort(newrun$importance_samples$samples) + 1, col = "red")
