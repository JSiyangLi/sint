### log.plus(log(a),log(b))=log(a+b) ###
log.plus <- function(x,y){
  if(x>y) x + log(1+exp(y-x))
  else    y + log(1+exp(x-y))
}

### log-likelihood ###
llike = function(n, x=5){
  ifelse(n >= x, 
         dbinom(x, size = n, prob = 0.5, log = TRUE),
         -Inf)
}

### log-prior ###
lprior = function(n){
  dgeom(x = n, prob = 0.5, log = TRUE)
}

### Sampling function ###
samp = function(steps, d, v, worst_llike, worst_obj){
  
  # This function generates samples from the prior with a 
  # likelihood restriction
  
  # steps: number of metropolis iterations to generate the new point
  # d    : dimension of the parameter space
  # v    : variance in the likelihood 
  # worst_llike: log-likelihood restriction
  # worst_obj  : starting point in the metropolis algorithm
  
  u    = matrix(runif(steps*d), nrow=steps, ncol=d);
  umh  = matrix(log(runif(steps*d)), nrow=steps, ncol=d);
  z    = matrix(rnorm(steps*d), nrow=steps, ncol=d); 
  
  di   = matrix(t(replicate(steps, sample(1:d))), ncol = d);
  
  proposal       = worst_obj;
  current        = worst_obj;
  
  current_llike  = worst_llike;
  current_lprior = lprior(worst_obj);
  
  for(i in 1:steps){
    for(j in di[i,]){
      
      proposal[j]     = current[j] + 4^(1.5 - 3*u[i,j]) * z[i,j];
      proposal_llike  = llike(proposal);
      proposal_lprior = lprior(proposal);
      
      if((umh[i,j] < (proposal_lprior - current_lprior)) & 
         proposal_llike > worst_llike){
        
        current[j]     = proposal[j];
        current_llike  = proposal_llike;
        current_lprior = proposal_lprior;
        
      }else{
        
        proposal[j] = current[j]; # Resetting value
        
      }
    }
  }
  return(list(proposal=current, llike=current_llike));
}

### Nested Sampling algorithm ###
ns <- function(n=10, steps=50, d=1, v=0.1, max.iter=2000, tol=0){
  
  # n       : number of active/live points in NS
  # steps   : number of metropolis steps to generate the new points in NS
  # d       : dimension of the parameter space
  # v       : variance in the likelihood 
  # max.iter: maximum number of iterations
  # tol     : tolerance
  
  Obj      = matrix(rnorm(n*d), nrow=n, ncol=d);
  iter     = 0;    
  width1   = log(1-exp(-1/n));
  logZb    = -1.79769e+308; 
  time     = c();
  t0       = proc.time();
  logL     = apply(Obj, 1, function(x)llike(x));    
  H        = 0;
  
  repeat
  {
    iter        = iter + 1;
    Order.logL  = order(logL);
    worst       = Order.logL[1];
    logwidth    = width1 - (iter-1)/n;        
    logZ        = log.plus(logZb, logwidth + logL[worst]);        
    H           = exp(logwidth + logL[worst] - logZ) * logL[worst] - 
      logZ + exp(logZb - logZ) * (H + logZb);
    logZb       = logZ;
    values      = samp(steps, d, v, worst_llike=logL[worst], 
                       worst_obj = Obj[worst,]);
    Obj[worst,] = values[[1]];
    logL[worst] = values[[2]];
    
    if(iter >= max.iter) break    
    #if(iter >= end * n * H) break
    if(max(logL)-iter/n < log(tol) + logZ ) break
    
  }
  
  t1    = proc.time();
  time  = (t1-t0)[1];
  
  logZc = log.plus(logZ, log(mean(exp(logL))) - iter/n); # correction
  
  return(list(logZ = logZc, H = H, iterations=iter, time = time));
}

###########
### RUN ###
###########



# Number of active points for a given number of samples
n_active_d1 =  1800 # 1e6/250 / (3 *  0.74);
n_active_d20 = 90   # 1e6/250 / (3 * 14.89);
n_active_d100 = 20  # 1e6/250 / (3 * 74.44);
OUTDIR = 'results_ns'
# make OUTDIR if it does not exist
if (!dir.exists(OUTDIR)) {
  dir.create(OUTDIR)
}

# save results to file

main  <- function(d, seed)
{
  if(d==1){
    n = n_active_d1;
  } else if(d==20){
    n = n_active_d20;
  } else if(d==100){
    n = n_active_d100;
  } else{
    stop("d must be 1, 20 or 100");
  }
  set.seed(seed);
  filename = paste(OUTDIR, "/results_ns_d", d, "_seed", seed, ".txt", sep="");
  results = ns(n, steps=250, d=d, v=0.1, max.iter = 1e6/250);
  write.table(results, file = filename, row.names = FALSE, col.names = FALSE);
}

main(d=1, seed=42)
