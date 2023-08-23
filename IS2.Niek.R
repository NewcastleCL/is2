#load packages

library(mvtnorm)
library(MCMCpack)
library(rtdists)
library(invgamma)
library(mixtools)
library(condMVNorm)
library(parallel)
library(corpcor)

# Standard IS2 adjusted by Niek Stevenson
# Rewritten to be easily adapted to use IS2 for other group level distributions
# And updated in terms of speed and vectorized calculations
IS2 <- function(samples, filter = "sample", subfilter = 0, IS_samples = 1000, stepsize_particles = 500, max_particles = 5000, n_cores = 1, df = 5){
  ###### set up variables #####
  info <- add_info_base(samples)
  idx <- which(samples$samples$stage == filter)
  if (length(subfilter)==1) {
    if(subfilter > 0) idx <- idx[-c(1:subfilter)]
  } else idx <- idx[subfilter]
  all_pars <- get_all_pars(samples, idx, info)
  muX<-apply(all_pars$X,2,mean)
  varX<-stats::cov(all_pars$X)
  
  prop_theta=mvtnorm::rmvt(IS_samples,sigma = varX, df=df, delta=muX)
  logw_num <- parallel::mclapply(X=1:IS_samples,
                                 FUN = compute_lw_num,
                                 prop_theta = prop_theta,
                                 stepsize_particles = stepsize_particles,
                                 max_particles = max_particles,
                                 mu_tilde=all_pars$mu_tilde,
                                 var_tilde = all_pars$var_tilde,
                                 info = all_pars$info,
                                 mc.cores = n_cores)

  logw_den <- mvtnorm::dmvt(prop_theta, delta=muX, sigma=varX,df=df, log = TRUE)
  
  finished <- unlist(logw_num) - logw_den
  # max.lw <- max(finished)
  # mean.centred.lw <- mean(exp(finished-max.lw)) #takes off the max and gets mean (avoids infs)
  # lw <- log(mean.centred.lw)+max.lw #puts max back on to get the lw
  return(finished)
}

get_sub_weights <- function(stepsize_particles, condMean, condVar, prop_theta, info, sub){
  wmix = 0.95
  n1=stats::rbinom(n=1,size=stepsize_particles,prob=wmix)
  if (n1<2) n1=2
  if (n1>(stepsize_particles-2)) n1=stepsize_particles-2 ## These just avoid degenerate arrays.
  n2=stepsize_particles-n1
  data <- info$data
  particles1 <- mvtnorm::rmvnorm(n1, condMean,condVar)
  # Group level
  particles2 <- group_dist_IS2(n_samples=n2, parameters = prop_theta,
                               sample=TRUE, info = info)
  particles <- rbind(particles1,particles2)
  # names for ll function to work
  colnames(particles) <- info$par_names
  # do lba log likelihood with given parameters for each subject, gets density of particle from ll func
  lw_first <- apply(particles, 1, info$ll_func, data = data[data$subject == unique(data$subject)[[sub]],])
  # below gets second part of equation 5 numerator ie density under prop_theta
  lw_second <- apply(particles, 1, group_dist_IS2, prop_theta, FALSE, NULL, info)
  # below is the denominator - ie mix of density under conditional and density under pro_theta
  lw_third <- log(wmix*pmax(1e-25 * info$n_randeffect, mvtnorm::dmvnorm(particles, condMean, condVar)) + (1-wmix) * exp(lw_second))
  # does equation 5
  lw <- lw_first+lw_second-lw_third
  return(lw)
}

get_logp=function(prop_theta,stepsize_particles, max_particles, mu_tilde,var_tilde, info){
  # Unload for legibility
  n_subjects <- info$n_subjects
  var_opt_sub <- 1/n_subjects
  
  # for each subject, get N IS samples (particles) and find log weight of each
  lw_subs <- numeric(n_subjects)
  for (j in 1:n_subjects){
    var_z_sub <- 999 # just a place holder > var_opt_sub
    n_total <- 0
    lw <- numeric()
    # generate the particles from the conditional MVnorm AND mix of group level proposals
    conditional = condMVN(mean=mu_tilde[j,],sigma=var_tilde[j,,],
                          dependent.ind=c(1:info$n_randeffect),
                          given.ind=info$given.ind,
                          X.given=prop_theta[info$X.given_ind],
                          check.sigma = F)
    while((var_z_sub > var_opt_sub) & n_total < max_particles){
      n_total <- n_total + stepsize_particles
      lw_tmp <- get_sub_weights(stepsize_particles = stepsize_particles, condMean = conditional$condMean,
                                condVar = conditional$condVar, prop_theta = prop_theta,
                                info = info, sub = j)
      
      # lw <- -weights(psis(-lw, log = F)) # default args are log=TRUE, normalize=TRUE
      lw <- c(lw, lw_tmp)
      max_lw <- max(lw)
      var_z_sub = sum(exp(2*(lw-max_lw)))/(sum(exp(lw-max_lw)))^2-1/n_total
    }
    weight <- exp(lw-max_lw)
    lw_subs[j] <- max_lw+log(mean(weight))
  }
  # sum the logp and return
  return(sum(lw_subs))
}

compute_lw_num=function(i, prop_theta,stepsize_particles, max_particles, mu_tilde,var_tilde,info){
  logp.out <- get_logp(prop_theta[i,], stepsize_particles, max_particles, mu_tilde, var_tilde, info)
  logw_num <- logp.out+prior_dist_IS2(parameters = prop_theta[i,], info)
  return(logw_num)
}

add_info_base <- function(samples){
  info <- list(
    n_randeffect = samples$n_pars,
    n_subjects = samples$n_subjects,
    par_names = samples$par_names,
    data = samples$data,
    ll_func = samples$ll_func,
    prior = samples$prior,
    hyper = attributes(samples)
  )
  return(info)
}


get_all_pars <- function(samples, idx, info){
  n_subjects <- samples$n_subjects
  n_iter = length(samples$samples$stage[idx])
  # Exctract relevant objects
  alpha <- samples$samples$alpha[,,idx]
  theta_mu <- samples$samples$theta_mu[,idx]
  theta_var <- samples$samples$theta_var[,,idx]
  a_half <- log(samples$samples$a_half[,idx])
  theta_var.unwound = apply(theta_var,3,unwind_IS2)
  # Set up
  n_params<- samples$n_pars+nrow(theta_var.unwound)+samples$n_pars
  all_samples=array(dim=c(n_subjects,n_params,n_iter))
  mu_tilde=array(dim = c(n_subjects,n_params))
  var_tilde=array(dim = c(n_subjects,n_params,n_params))
  
  for (j in 1:n_subjects){
    all_samples[j,,] = rbind(alpha[,j,],theta_mu[,],theta_var.unwound[,])
    # calculate the mean for re, mu and sigma
    mu_tilde[j,] =apply(all_samples[j,,],1,mean)
    # calculate the covariance matrix for random effects, mu and sigma
    var_tilde[j,,] = cov(t(all_samples[j,,]))
  }
  
  for(i in 1:n_subjects){ #RJI_change: this bit makes sure that the sigma tilde is pos def
    if(!corpcor::is.positive.definite(var_tilde[i,,], tol=1e-8)){
      var_tilde[i,,]<-corpcor::make.positive.definite(var_tilde[i,,], tol=1e-6)
    }
  }
  X <- cbind(t(theta_mu),t(theta_var.unwound),t(a_half))
  info$n_params <- n_params
  info$given.ind <- (info$n_randeffect+1):n_params
  info$X.given_ind <- 1:(n_params-info$n_randeffect)
  return(list(X = X, mu_tilde = mu_tilde, var_tilde = var_tilde, info = info))
}

unwind_IS2 <- function(x,reverse=FALSE) {
  
  if (reverse) {
    n=sqrt(2*length(x)+0.25)-0.5 ## Dim of matrix.
    out=array(0,dim=c(n,n))
    out[lower.tri(out,diag=TRUE)]=x
    diag(out)=exp(diag(out))
    out=out%*%t(out)
  } else {
    y=t(base::chol(x))
    diag(y)=log(diag(y))
    out=y[lower.tri(y,diag=TRUE)]
  }
  return(out)
}


group_dist_IS2 <- function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, info){
  n_randeffect <- info$n_randeffect
  param.theta.mu <- parameters[1:n_randeffect]
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)]
  param.theta.sig2 <- unwind_IS2(param.theta.sig.unwound, reverse = TRUE)
  if (sample){
    return(rmvnorm(n_samples, param.theta.mu,param.theta.sig2))
  }else{
    logw_second<-max(-1000*info$n_randeffect,  mvtnorm::dmvnorm(random_effect, param.theta.mu,param.theta.sig2,log=TRUE))
    return(logw_second)
  }
}

prior_dist_IS2 <- function(parameters, info){
  n_randeffect <- info$n_randeffect
  prior <- info$prior
  hyper <- info$hyper
  param.theta.mu <- parameters[1:n_randeffect]
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)]
  param.theta.sig2 <- unwind_IS2(param.theta.sig.unwound, reverse = TRUE)
  param.a <- exp(parameters[((length(parameters)-n_randeffect)+1):(length(parameters))])
  log_prior_mu= mvtnorm::dmvnorm(param.theta.mu, mean = prior$theta_mu_mean, sigma = prior$theta_mu_var, log =TRUE)
  log_prior_sigma = log(robust_diwish(param.theta.sig2, v=hyper$v_half+ n_randeffect-1, S = 2*hyper$v_half*diag(1/param.a)))
  log_prior_a = sum(logdinvGamma(param.a,shape = 1/2,rate=1/(hyper$A_half^2)))
  # These are Jacobian corrections for the transformations on these
  logw_den2 <- -sum(log(param.a))
  logw_den3 <- -(log(2^n_randeffect)+sum((n_randeffect:1+1)*log(diag(param.theta.sig2))))
  return(log_prior_mu + log_prior_sigma + log_prior_a - logw_den3 - logw_den2)
}


# Some helper density functions -------------------------------------------

# This one was added to have shape, rate, parametrization
logdinvGamma <- function(x, shape, rate){
  alpha <- shape
  beta <- 1/rate
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha +
                                                        1) * log(x) - (beta/x)
  return(pmax(log.density, -100)) #Roughly equal to 1e-22 on real scale
}

# This one to use the chol2inv rather than the slow inverse the standard condMVN uses.
condMVN <- function (mean, sigma, dependent.ind, given.ind, X.given, check.sigma = TRUE)
{
  if (missing(dependent.ind))
    return("You must specify the indices of dependent random variables in `dependent.ind'")
  if (missing(given.ind) & missing(X.given))
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,
                                                                          dependent.ind])))
  if (length(given.ind) == 0)
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,
                                                                          dependent.ind])))
  if (length(X.given) != length(given.ind))
    stop("lengths of `X.given' and `given.ind' must be same")
  if (check.sigma) {
    if (!isSymmetric(sigma))
      stop("sigma is not a symmetric matrix")
    eigenvalues <- eigen(sigma, only.values = TRUE)$values
    if (any(eigenvalues < 1e-08))
      stop("sigma is not positive-definite")
  }
  B <- sigma[dependent.ind, dependent.ind]
  C <- sigma[dependent.ind, given.ind, drop = FALSE]
  D <- sigma[given.ind, given.ind]
  CDinv <- C %*% chol2inv(chol(D))
  cMu <- c(mean[dependent.ind] + CDinv %*% (X.given - mean[given.ind]))
  cVar <- B - CDinv %*% t(C)
  list(condMean = cMu, condVar = cVar)
}

#  This one is to protect against weird proposals in the diwish function, where sometimes matrices weren't pos def
robust_diwish <- function (W, v, S) { 
  if (!is.matrix(S)) S <- matrix(S)
  if (!is.matrix(W)) W <- matrix(W)
  p <- nrow(S)
  gammapart <- sum(lgamma((v + 1 - 1:p)/2))
  ldenom <- gammapart + 0.5 * v * p * log(2) + 0.25 * p * (p - 1) * log(pi)
  if (corpcor::is.positive.definite(W, tol=1e-8)){
    cholW<-base::chol(W)
  }else{
    return(1e-10)
  }
  if (corpcor::is.positive.definite(S, tol=1e-8)){
    cholS <- base::chol(S)
  }else{
    return(1e-10)
  }
  halflogdetS <- sum(log(diag(cholS)))
  halflogdetW <- sum(log(diag(cholW)))
  invW <- chol2inv(cholW)
  exptrace <- sum(S * invW)
  lnum <- v * halflogdetS - (v + p + 1) * halflogdetW - 0.5 * exptrace
  lpdf <- lnum - ldenom
  out <- exp(lpdf)
  if(!is.finite(out)) return(1e-100)
  if(out < 1e-10) return(1e-100)
  return(exp(lpdf))
}

