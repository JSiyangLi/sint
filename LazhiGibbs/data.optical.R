
library(reticulate)
optical.data <- function(
    Time          = 5e4, # from Lazhi, with a temporary scaling of 1e-8
    count_file, # the .txt file of photon counts
    area_file, # the .txt file of segment and background areas
    exposure_file, # the .txt file of segment exposures
    rmat_npz, # the .npz file of r matrix
    bkgd_file # the .txt file of background
)
{
  extract_count_and_area <- function(file_path) {
    # Read the file lines
    lines <- readLines(file_path)
    
    # Initialize variables
    count <- NULL
    area <- NULL
    
    # Loop through lines to find the values
    for (line in lines) {
      # Extract Count (even if line starts with #)
      if (grepl("Count:", line, ignore.case = TRUE)) {
        count <- as.numeric(gsub("[^0-9.]", "", line))
      }
      
      # Extract Area (even if line starts with #)
      if (grepl("Area:", line, ignore.case = TRUE)) {
        area <- as.numeric(gsub("[^0-9.]", "", line))
      }
    }
    
    # Return as a named list
    list(count = count, area = area)
  }
  
  bkgd_info <- extract_count_and_area(bkgd_file)
  # read in the data
  df <- read.table(file = count_file)$V2 # temporary scaling
  segment.count <- df[-length(df)]
  bkgd.count <- bkgd_info$count
  
  Areas <- read.table(file = area_file)$V2
  Aroi <- Areas[-1]
  bkgd.area <- bkgd_info$area
  
  Exposures <- read.table(file = exposure_file)$V2[-1]
  
  np <- import("numpy")
  npzf <- np$load(rmat_npz)
  
  n.segments    <- length(Exposures)     # number of segments
  
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
  sources <- list()
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
  
  
  seg$Time       <- Time
  seg$n          <- n.segments # number of segments
  sources$n      <- npzf$f["shape"][2] # number of sources
  indptr <- npzf$f['indptr'][-1] + 1
  seg$n.overlap  <- diff(indptr)
  seg$area       <- Aroi
  seg$eff.area   <- Exposures
  seg$maxoverlap <- max(diff(indptr))
  seg$cnt.obs <- segment.count
  
  indcol <- npzf$f['indices'] + 1
  src.r <- npzf$f['data']
  seg$r <- seg$src.index <- matrix(rep(0, (seg$n)*(seg$maxoverlap)), ncol = seg$maxoverlap)
  for (i in 1:n.segments) {
    src.pointer <- indptr[i]:(indptr[i+1]-1)
    src.len <- length(src.pointer)
    seg$src.index[i, 1:src.len] <- indcol[src.pointer]
    seg$r[i, 1:src.len] <- src.r[src.pointer]
  }
  
  # The sources
  #
  #  "n"                  number of sources
  #                       ***** must be equal to max(seg$src.index) *****
  #. "exposure"           effective exposure: Time * sum(seg$eff.area * seg$r)
  for(i in 1:sources$n){
    sources$exposure[i] <- Time * sum((seg$src.index == i) * seg$r * seg$eff.area)
    if(sum((seg$src.index == i) * seg$r) >1.0) 
      print("Error: element of sum( sources$r ) > 1 for a single source.")
  }
  
  # The Background
  #
  bkgd <- list()
  #  "area",               background area
  #  "cnt",                observed bkgd count
  #  "prior.shape",        prior shape parameter (for gamma prior) 
  #  "prior.rate"          prior rate parameter (for gamma prior) 
  
  # generate the background counts
  bkgd$cnt    <- bkgd.count
  bkgd$area   <- bkgd.area
  bkgd$prior.shape   <- 1/10^6.      # prior values form Lazhi's paper
  bkgd$prior.rate    <- Time * bkgd$area / 10^12
  
  # latent variables
  seg$cnt.mis    <- matrix(0,ncol= seg$maxoverlap + 1, nrow=seg$n)
  sources$cnt        <- rep(1,sources$n)
  
  output <- list("Segments" = seg, "Sources" = sources, "Background"= bkgd)
  saveRDS(output, "oData.RDS")
  output
} # optical data
