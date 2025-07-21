
gibbs <- function(seg,      # data structure for segments, see data.simulate.R
                  sources,  # data structure for sources, see data.simulate.R
                  bkgd,     # data structure for bkgd, see data.simulate.R
                  epsilon = .Machine$double.eps, 
                            # source rates (lambda) < epsilon are set to zero
                  epsilon.max.iter = 2000,
                            # no zeroing out source rates after this iteration
                  mcmc.samplesize = 5e4,
                            # number of mcmc samples to draw.
                  start.src.rates.mean = 15,
                  start.src.rates.var  = 50, 
                            # gamma mean and var for starting value of src rates
                  start.zero.prob = 0.5
                            # starting value for zero inflation probability
                  )
{
  
  Time   <- seg$Time     # shorthand for total exposure time
  accept <- 0            # number of random walk metropolis proposals accepted
  reject <- 0.           # number of random walk metropolis proposals rejected
  
  # Model Parameters -- holds MCMC sample of each
  #
  params <- list()
  #  "bkgd.rate",         background rate parameters
  #  "src.rates",         vector of source rates
  #  "src.zero",          vector of indicator variables for zero rate
  #  "pop.shape",         shape parameter for gamma population distribution
  #  "pop.rate",          rate parameter for gamma population distribution
  #  "pop.zero.prob"      zero inflation rate for gamma population distribution
  
  
  # Starting values will be held in first element / row,
  # So number of rows = mcmc.samplesize + 1
  #
  params$bkgd.rate     <- rep(0, mcmc.samplesize +1)
  params$src.rates     <- matrix(0,ncol=sources$n,nrow=mcmc.samplesize+1)
  params$zero          <- matrix(0,ncol=sources$n,nrow=mcmc.samplesize+1)
  params$pop.shape     <- rep(0, mcmc.samplesize +1)
  params$pop.rate      <- rep(0, mcmc.samplesize +1)
  params$pop.zero.prob <- rep(0, mcmc.samplesize +1)
  params$pop.mean      <- rep(0, mcmc.samplesize +1) 
  params$pop.sd        <- rep(0, mcmc.samplesize +1) 
  
  # Set starting values
  #
  params$bkgd.rate[1]     <- bkgd$cnt / bkgd$area / Time
  params$src.rates[1,]    <- rgamma(sources$n, 
                                shape=start.src.rates.mean^2/start.src.rates.var, 
                                rate=start.src.rates.mean/start.src.rates.var)
# params$pop.shape[1]     <- 1   # these are set in Step 5 of first iteration 
# params$pop.rate[1]      <- 1   # using sampled src.rate  
  params$pop.zero.prob[1] <- start.zero.prob
# params$zero[1,]         <-     # sampled in Step 3
  

  for(iter in 1:mcmc.samplesize){
    if(iter%%1000 == 0){cat("Iteration number:", iter, "\n")}
    
    ######################################################
    #                    STEP 1                          #
    # Separate Segment Counts in Sources and Background. #
    ######################################################
    #
    for(s in 1:seg$n){
      # probabilities corresponding to sources
      rates.of.overlap.src <- rep(0, seg$maxoverlap)
      rates.of.overlap.src[1:seg$n.overlap[s]] <- params$src.rates[iter,seg$src.index[s,]] 
      prob <- seg$r[s,] * seg$eff.area[s] * rates.of.overlap.src
      
      # concatenate the probability from background
      prob <- c(prob, seg$area[s] * params$bkgd.rate[iter])
      
      # resample segments counts into corresponding sources and background
      seg$cnt.mis[s,] <- c(rmultinom(n=1, size=seg$cnt.obs[s], prob=prob))
    } # s in 1:seg$n
    # store source counts
    for(i in 1:sources$n){
      # Sum all counts associated with src i
      sources$cnt[i] <- sum((seg$src.index == i) * seg$cnt.mis[,1:seg$maxoverlap])
    }
    
    ######################################################
    #                    STEP 2                          #
    #            Update Background Rate.                 #
    ######################################################  
    # 
    post.shape <- bkgd$prior.shape + bkgd$cnt + sum(seg$cnt.mis[,seg$maxoverlap+1])
    post.rate <- bkgd$prior.rate + Time * (bkgd$area + sum(seg$area)) 
    params$bkgd.rate[iter+1] <- rgamma(1, shape=post.shape, rate=post.rate)
    
    
    #######################################################
    #                    STEP 3                           #
    # Sample source rate parameters (and zero indicators) #
    #######################################################  
    # 
    #
    for(i in 1:sources$n){
      # Sum all counts associated with src i *** NOW IN STEP 1 ***
      # count <- sum((seg$src.index == i) * seg$cnt.mis[,1:seg$maxoverlap])
      
      # Compute parameters of posterior gamma dist'n for src i
      post.shape <- params$pop.shape[iter] + sources$cnt[i]
      post.rate <- params$pop.rate[iter] + sources$exposure[i]
      
      
      # Compute posterior probability that src i has rate zero given count =0
      prob.zero.rate <- params$pop.zero.prob[iter] /
        (params$pop.zero.prob[iter] + 
           (1-params$pop.zero.prob[iter])*
           (params$pop.rate[iter]/post.rate)^params$pop.shape[iter])
       #print(c(sources$cnt[i], post.rate, prob.zero.rate, params$pop.rate[iter], sources$exposure[i]))
      
      # if count = 0, set rate = 0 with computed prob, else sample gamma dist'n
      if((sources$cnt[i]==0) * (runif(1,0,1) < prob.zero.rate)){
        params$zero[iter+1,i] = 1
        params$src.rates[iter+1,i] = 0
      }else{
        params$zero[iter+1,i] = 0
        params$src.rates[iter+1,i] <- rgamma(1, shape=post.shape, rate=post.rate)
        if((params$src.rates[iter+1,i] < epsilon) * (iter<epsilon.max.iter)){ 
          # zero out small values of src.rate in first few iterations
          cat("params$src.rates[", i,"] =", 
              round(params$src.rates[iter+1,i], digits=5),  
              "was zeroed out in iteration", iter,
              ", post.shape: ", post.shape, ", post.rate:", post.rate,
              "\n")
          params$zero[iter+1,i] = 1
          params$src.rates[iter+1,i] = 0
        }
      }
    } # i in 1:sources$n
    
    #####################################################
    #                    STEP 4                         #
    #         Update params$src.zero.prob               #
    #####################################################  
    # 
    #
    dark <- sum(params$zero[iter+1,])
    params$pop.zero.prob[iter+1] <- rbeta(1, dark + 1, sources$n - dark + 1)
    
    ######################################################
    #                    STEP 5                          #
    # Update population dist'n of source rate parameters #
    ######################################################  
    # 
    # 
    not.zero.rates <- params$src.rates[iter+1, params$zero[iter+1,]==0]
    
    # set current to log(mean) and log(variance) of gamma population
    # set starting values in first iteration
    if(iter == 1){
      current <-c(log(mean(not.zero.rates)), log(var(not.zero.rates)))
      }else{
        current <- c(log(params$pop.mean[iter]), 2*log(params$pop.sd[iter]))
      }

    # Use optim to compute hessian for use in random walk Metropolis
    if(iter%%50 == 1){       # recompute only every 50th iteration
      #cat("zero rate length: ", length(not.zero.rates), "\n")
      #cat("params zero: ", params$zero[iter+1,], "\n")
      #cat("src rates: ", params$src.rates[iter+1, ], "\n")
      #cat("current value: ", current, "\n")
      #cat("cond.post.log.value: ", cond.post.log(current, data = not.zero.rates), "\n")
      fit <- optim(current, cond.post.log, data=not.zero.rates, hessian=TRUE)
      variance <- solve(fit$hessian)
    } # iter%%50 == 1
      
    for(inner in 1:10){
      # Random-walk Metropolis
      proposal <- rmvnorm( 1, current, sigma=variance )
      accept.prob <- exp( -cond.post.log(proposal, data=not.zero.rates) 
                          + cond.post.log(current, data=not.zero.rates)  )
      # print(c(current, proposal, accept.prob))
      if(runif(1) < accept.prob){
        current[1] <- proposal[1]
        current[2] <- proposal[2]
        accept <- accept + 1
      } else{reject <- reject + 1}
    }

    # Transform log(mean) and log(var) parameters to gamma shape and rate.
    params$pop.mean[iter+1]  <- exp(current[1])
    params$pop.sd[iter+1]    <- exp(current[2]/2) 
    params$pop.shape[iter+1] <- exp(2*current[1] - current[2])
    params$pop.rate[iter+1]  <- exp(current[1] - current[2])
  } # iter 
  params$random.walk.accept.rate <- accept/(accept+reject)
  params
}
