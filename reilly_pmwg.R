#######       IS2 Psuedo Code       #########

## set up environment and packages 
rm(list=ls())
library(mvtnorm)
library(MCMCpack)
library(rtdists)
library(invgamma)
library(mixtools)
library(condMVNorm)
library(parallel)
load("forstmann_is2_3b3t.Rdata")

cpus = 1
###### set up variables #####
# number of particles, samples, subjects, random effects etc
n_randeffect=sampled$n_pars
n_subjects = sampled$n_subjects
n_iter = length(sampled$samples$stage[sampled$samples$stage=="sample"])
length_draws = sampled$samples$idx #length of the full transformed random effect vector and/or parameter vector
IS_samples = 100 #number of importance samples
n_particles = 30#number of particles
v_alpha = 2  #?
pars = sampled$pars


# grab the sampled stage of PMwG
# store the random effects
alpha <- sampled$samples$alpha[,,sampled$samples$stage=="sample"]
# store the mu
theta <- sampled$samples$theta_mu[,sampled$samples$stage=="sample"]
# store the cholesky transformed sigma
sig <- sampled$samples$theta_sig[,,sampled$samples$stage=="sample"]
# the a-hlaf is used in  calculating the Huang and Wand (2013) prior. 
# The a is a random sample from inv gamma which weights the inv wishart. The mix of inverse wisharts is the prior on the correlation matrix
a_half <- log(sampled$samples$a_half[,sampled$samples$stage=="sample"])



unwind=function(x,reverse=FALSE) {

  if (reverse) {
    ##        if ((n*n+n)!=2*length(x)) stop("Wrong sizes in unwind.")
    n=sqrt(2*length(x)+0.25)-0.5 ## Dim of matrix.
    out=array(0,dim=c(n,n))
    out[lower.tri(out,diag=TRUE)]=x
    diag(out)=exp(diag(out))
    out=out%*%t(out)
  } else {
    y=t(chol(x))
    diag(y)=log(diag(y))
    out=y[lower.tri(y,diag=TRUE)]
  }
  return(out)
}

#unwound sigma
pts2.unwound = apply(sig,3,unwind)

n.params<- nrow(pts2.unwound)+n_randeffect+n_randeffect
all_samples=array(dim=c(n_subjects,n.params,n_iter))
mu_tilde=array(dim = c(n_subjects,n.params))
sigma_tilde=array(dim = c(n_subjects,n.params,n.params))



for (j in 1:n_subjects){
  all_samples[j,,] = rbind(alpha[,j,],theta[,],pts2.unwound[,])
  # calculate the mean for re, mu and sigma
  mu_tilde[j,] =apply(all_samples[j,,],1,mean)
  # calculate the covariance matrix for random effects, mu and sigma
  sigma_tilde[j,,] = cov(t(all_samples[j,,]))
}



X <- cbind(t(theta),t(pts2.unwound),t(a_half)) 

# do k=2, for a mixture of 2 gaussians (Davids suggestions, could be between 1-5 really)
k = 2 #number of dists

#mvnormalmixEM is a weak point - function can fail. needs a note or output to show if it doesn't work. Should restart if it fails
mix = NULL
while(is.null(mix)) {
  tryCatch(mix<-mvnormalmixEM(X,k=k, maxit = 5000),error=function(e){
  },finally={})
}


mix_weight <- mix$lambda
mix_mu <- mix$mu
mix_sigma <- mix$sigma

###### generate the proposal parameters from the mix of importance samples  #####

# use the weight to get samples for n1. n2 = samples-n1 (i.e 9000 and 1000)
n1=rbinom(n=1, size=IS_samples,prob=max(mix_weight))
n1=pmax(n1,2)
n1=pmin(n1,IS_samples-2)
n2=IS_samples-n1

# generates the 10,000 IS proposals given the mix
proposals1=rmvnorm(n1,mix_mu[[1]],mix_sigma[[1]])
proposals2=rmvnorm(n2,mix_mu[[2]],mix_sigma[[2]])
prop_theta=rbind(proposals1,proposals2)




### main functions

# get_params = function(parameters, n_randeffect){
#   param.theta.mu <- parameters[1:n_randeffect]
#   param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)] ##scott would like it to ask for n(unwind)
#   param.theta.sig2 <- unwind(param.theta.sig.unwound, reverse = TRUE)
#   param.a <- exp(parameters[((length(parameters)-n_randeffect)+1):(length(parameters))])
# }

group_dist = function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, n_randeffect){
  param.theta.mu <- parameters[1:n_randeffect]
  ##scott would like it to ask for n(unwind) rather than doing the calculation for how many it actually needs, you should just input the length of the unwound object
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)] 
  param.theta.sig2 <- unwind(param.theta.sig.unwound, reverse = TRUE)
  if (sample){
    return(mvtnorm::rmvnorm(n=n_samples, mean=param.theta.mu,sigma=param.theta.sig2))
  }else{
    logw_second<-mvtnorm::dmvnorm(random_effect, param.theta.mu,param.theta.sig2,log=TRUE)
    return(logw_second)
  }
}



prior_dist = function(parameters, prior_parameters = sampled$prior, n_randeffect){ ###mod notes: the sampled$prior needs to be fixed/passed in some other time
  param.theta.mu <- parameters[1:n_randeffect]
  param.theta.sig.unwound <- parameters[(n_randeffect+1):(length(parameters)-n_randeffect)] ##scott would like it to ask for n(unwind)
  param.theta.sig2 <- unwind(param.theta.sig.unwound, reverse = TRUE)
  param.a <- exp(parameters[((length(parameters)-n_randeffect)+1):(length(parameters))])
  v_alpha=2
  
  log_prior_mu=mvtnorm::dmvnorm(param.theta.mu, mean = prior_parameters$theta_mu_mean, sigma = prior_parameters$theta_mu_var, log =TRUE)
  log_prior_sigma = log(robust_diwish(param.theta.sig2, v=v_alpha+ n_randeffect-1, S = 2*v_alpha*diag(1/param.a)))  #exp of a-half -> positive only
  log_prior_a = sum(invgamma::dinvgamma(param.a,scale = 0.5,shape=1,log=TRUE))
  
  logw_den2 <- sum(log(1/param.a)) # jacobian determinant of transformation of log of the a-half
  logw_den3 <- log(2^n_randeffect)+sum((n_randeffect:1+1)*log(diag(param.theta.sig2))) #jacobian determinant of cholesky factors of cov matrix
  
  return(log_prior_mu + log_prior_sigma + log_prior_a + logw_den3 - logw_den2)
}





get_logp=function(prop_theta,data,n_subjects,n_particles,n_randeffect,mu_tilde,sigma_tilde,i, group_dist=group_dist){
  #make an array for the density
  logp=array(dim=c(n_particles,n_subjects))
  # for each subject, get 1000 IS samples (particles) and find log weight of each
  for (j in 1:n_subjects){
    #generate the particles from the conditional MVnorm AND mix of group level proposals
    wmix <- 0.95
    n1=rbinom(n=1,size=n_particles,prob=wmix)
    if (n1<2) n1=2
    if (n1>(n_particles-2)) n1=n_particles-2 ## These just avoid degenerate arrays.
    n2=n_particles-n1
    # do conditional MVnorm based on the proposal distribution
    conditional = condMVNorm::condMVN(mean=mu_tilde[j,],sigma=sigma_tilde[j,,],dependent.ind=1:n_randeffect,
                          given.ind=(n_randeffect+1):n.params,X.given=prop_theta[i,1:(n.params-n_randeffect)])
    particles1 <- mvtnorm::rmvnorm(n1, conditional$condMean,conditional$condVar)
    # mix of proposal params and conditional
    particles2 <- group_dist(n_samples=n2, parameters = prop_theta[i,],sample=TRUE, n_randeffect=n_randeffect)
    particles <- rbind(particles1,particles2)
    
    for (k in 1:n_particles){
      x <-particles[k,]
      #names for ll function to work
      #mod notes: this is the bit the prior effects
      names(x)<-sampled$par_names
      #   do lba log likelihood with given parameters for each subject, gets density of particle from ll func
      logw_first=sampled$ll_func(x,data = data[as.numeric(factor(data$subject))==j,]) #mod notes: do we pass this in or the whole sampled object????
      # below gets second part of equation 5 numerator ie density under prop_theta
      # particle k and big vector of things
      logw_second<-group_dist(random_effect = particles[k,], parameters = prop_theta[i,], sample= FALSE, n_randeffect = n_randeffect) #mod notes: group dist
      # below is the denominator - ie mix of density under conditional and density under pro_theta
      logw_third <- log(wmix*dmvnorm(particles[k,], conditional$condMean,conditional$condVar)+(1-wmix)*exp(logw_second)) #mod notes: fine?
      #does equation 5
      logw=(logw_first+logw_second)-logw_third
      #assign to correct row/column
      logp[k,j]=logw 
    }
  }
  #we use this part to centre the logw before addign back on at the end. This avoids inf and -inf values
  sub_max = apply(logp,2,max)
  logw = t(t(logp) - sub_max)
  w = exp(logw)
  subj_logp = log(apply(w,2,mean))+sub_max #means
  
  # this part gets the variance the same as David
  # has been commented out here, but could be an option to use this
  #var_numerator = apply(w^2, 2, sum)
  #var_denominator = apply(w, 2, sum)^2
  #variance
  #logp_variance = (var_numerator/var_denominator) - (1/n_particles) #for each subject
  
  # sum the logp and return 
  return(sum(subj_logp))
}

compute_lw=function(prop_theta,data,n_subjects,n_particles,n_randeffect,mu_tilde,sigma_tilde,i, prior_dist=prior_dist, sampled=sampled){
  
  logp.out <- get_logp(prop_theta,data,n_subjects,n_particles,n_randeffect,mu_tilde,sigma_tilde,i, group_dist=group_dist)
  ##do equation 10
  logw_num <- logp.out[1]+prior_dist(parameters = prop_theta[i,], prior_parameters = sampled$prior, n_randeffect)
  logw_den <- log(mix_weight[1]* mvtnorm::dmvnorm(prop_theta[i,], mix_mu[[1]], mix_sigma[[1]])+ mix_weight[2]* mvtnorm::dmvnorm(prop_theta[i,], mix_mu[[2]], mix_sigma[[2]])) #density of proposed params under the means
  logw <- logw_num-logw_den #this is the equation 10
  return(c(logw))
  #NOTE: we should leave a note if variance is shit - variance is given by the logp function (currently commented out)
}



##### make it work

#makes an array to store the IS samples
tmp<-array(dim=c(IS_samples))

#do the sampling
if (cpus>1){
  tmp <- mclapply(X=1:IS_samples,mc.cores = cpus, FUN = compute_lw, prop_theta = prop_theta,data = data,n_subjects= n_subjects,n_particles = n_particles,
                  n_randeffect = n_randeffect,mu_tilde=mu_tilde,sigma_tilde = sigma_tilde, prior_dist=prior_dist, sampled=sampled)
} else{
  for (i in 1:IS_samples){
    cat(i)
    tmp[i]<-compute_lw(prop_theta,data,n_subjects,n_particles, n_randeffect,mu_tilde,sigma_tilde,i,prior_dist=prior_dist, sampled=sampled)
  }
}


# get the ML value
finished <- tmp
tmp<-unlist(tmp)
max.lw <- max(tmp)
mean.centred.lw <- mean(exp(tmp-max.lw)) #takes off the max and gets mean (avoids infs)
lw <-log(mean.centred.lw)+max.lw #puts max back on to get the lw



##### bootstrapping for SE ######
bootstrap = 10000
log_marglik_boot= array(dim = bootstrap)
for (i in 1:bootstrap){
  log_weight_boot = sample(tmp, IS_samples, replace = TRUE) #resample with replacement from the lw
  max.boot <- max(log_weight_boot)
  centred.boot <- mean(exp(log_weight_boot-max.boot)) #takes off the max and gets mean (avoids infs)
  log_marglik_boot[i] <-log(centred.boot)+max.boot #puts max back on 
}
var(log_marglik_boot)


save.image("IS2_v2.Rdata")

