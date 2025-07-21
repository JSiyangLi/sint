
simulate.data <- function(
    Time          = 100000,
    pop.mean      = 10,     # mean and variance of gamma dist'n of non zero  
    pop.sd        = 25,     # source rates
    pop.zero.prob = 0.3,    # proportion of zero source rates
    bkgd.rate     = pop.shape/pop.rate/20,
    n.segments    = 150     # number of segments is rounded to multiple of 15
    )
{
  
  pop.shape     <- (pop.mean / pop.sd)^2
  pop.rate      <- pop.mean / pop.sd^2
  
  ######################################################################
  # Data is generated after details of segments, (effective) areas, etc.
  #
  # The basic data pattern involves 10 sources. Of these, sources 1-5 
  # do not overlap any sources; Sources 6-8 overlap as a triple; and 
  # Source 9-10 overlap as a double. This results in 15 segments. 
  # The basic pattern is repeated reps times. The values of seg.r, seg.area,
  # seg.eff.area, etc also repeate reps time. The counts are drawn 
  # indepednetly for all reps * 15 segments.
  #
  # Thus, the input number of segments is rounded to nearest multiple of 15.
  #
  ######################################################################
  
  # The Segments
  #
  seg <- list()
  #  "Time",           total exposure time
  #  "n",              the number of segments
  #  "max.overlap",    the maximum number of overlapping src regions
  #  "n.overlap",      vector of number of overlapping regions is each segment
  #  "area",           vector of the segment areas
  #  "eff.area",       vector of the segment effective areas
  #  "cnt.obs",        vector of the observed segment counts
  #  "src.index",      matrix (n x max.overlap) of src indices assoc with segment
  #             ***** indices must be 1, 2, ..., sources$n, with no gaps *****
  #  "r",              matrix (n x max.overlap) prob cnt from src lands in segment
  #  "cnt.mis"         matrix (n x max.overlap) of missing segment x src/bkgd cnts
  
  
  seg$Time       <- Time
  seg$n          <- 15*round(n.segments/15, digits=0)
  reps           <- seg$n/15
  seg$maxoverlap <- 3
  seg$n.overlap  <- rep(c(1,1,1,1,1,1,1,1,2,2,2,3,1,1,2), reps)
  seg$area       <- rep(c(8,8,8,8,8,3,3,3,2,2,2,1,5,5,3), reps)
  seg$eff.area   <- rep(c(8,8,7,8,6,6,7,6,4,4,4,2,10,5,3),reps)
  seg$src.index  <- matrix(rep(c(1,0,0, 2,0,0, 3,0,0, 4,0,0, 5,0,0, 6,0,0, 7,0,0, 
                             8,0,0, 6,7,0, 6,8,0, 7,8,0, 6,7,8, 9,0,0, 10,0,0, 
                             9,10,0), reps), byrow=T, ncol=3)
  temp           <- (seg$src.index  > 0)
  seg$src.index  <- seg$src.index + rep(10*(1:reps -1), each=15) * temp
  
  seg$r          <- matrix(rep(c(0.9,0,0, 0.9,0,0, 0.9,0,0, 0.9,0,0, 0.9,0,0, 
                             0.4,0,0, 0.4,0,0, 0.5,0,0, 
                             0.3,0.1,0, 0.1,0.2,0, 0.2,0.1,0, 0.1,0.2,0.1,
                             0.6,0,0, 0.6,0,0, 0.3,0.3,0), reps), byrow=T, ncol=3)
  seg$cnt.mis    <- matrix(0,ncol= seg$maxoverlap + 1, nrow=seg$n)
  
  
  # The Background
  #
  bkgd <- list()
  #  "area",               background area
  #  "cnt",                observed bkgd count
  #  "prior.shape",        prior shape parameter (for gamma prior) 
  #  "prior.rate"          prior rate parameter (for gamma prior) 
  
  bkgd$area          <- 1000
  bkgd$prior.shape   <- 1/10^6.      # prior values form Lazhi's paper
  bkgd$prior.rate    <- Time * bkgd$area / 10^12
  
  
  # The sources
  #
  sources <- list()
  #  "n"                  number of sources
  #                       ***** must be equal to max(seg$src.index) *****
  #. "exposure"           effective exposure: Time * sum(seg$eff.area * seg$r)
  #  "cnt"                the (missing number of cnts from each source)
  
  sources$n          <- reps * 10
  sources$cnt        <- rep(1,sources$n)
  for(i in 1:sources$n){
    sources$exposure[i] <- Time * sum((seg$src.index == i) * seg$r * seg$eff.area)
    if(sum((seg$src.index == i) * seg$r) >1.0) 
      print("Error: element of sum( sources$r ) > 1 for a single source.")
  }
  
  
  
  # GENERATE DATA
  #
  # Generate the source rates
  src.rates          <- rgamma(sources$n, shape=pop.shape, rate=pop.rate)
  src.rates.zero     <- rbinom(sources$n, size=1, prob=pop.zero.prob)
  src.rates          <- src.rates * (1-src.rates.zero)
  
  # Generate the segment counts
  for(s in 1:seg$n){
    # Generate the counts for from each src that contributes to segment s
    seg$cnt.mis[s,1:seg$n.overlap[s]] <- 
      rpois(n=seg$n.overlap[s], 
            lambda=Time * seg$r[s,1:seg$n.overlap[s]] * seg$eff.area[s] * src.rates[seg$src.index[s,]]) 
    # Generate the background count in segment s (stored in last column)
    seg$cnt.mis[s,seg$maxoverlap+1] <-rpois(1,Time * seg$area[s] * bkgd.rate)
  } # s in 1:seg$n 
  # sum over various sources to obtain observed data
  seg$cnt.obs <- apply(seg$cnt.mis, 1, sum)
  
  for(i in 1:sources$n){
    sources$cnt[i] <- 
      sum((seg$src.index == i) * seg$cnt.mis[,-(seg$maxoverlap+1)])
  }
  
  # generate the background counts
  bkgd$cnt    <- rpois(1, Time * bkgd$area * bkgd.rate)
  
  list("Segments" = seg, "Sources" = sources, "Background"= bkgd)
  
} # simulate data
