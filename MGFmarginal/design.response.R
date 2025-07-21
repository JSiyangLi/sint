library(reticulate)
library(Matrix)
xray.design.response <- function(
    design_file, # the file containing the design matrix
    response_file # the file containing the response vector
) {
  Time <- 5e4   # from Lazhi
  
  design <- as.matrix(read.table(design_file))
  response <- read.table(response_file)[, 2]
  list("Design" = design, "Response" = response, "Time" = Time)
}

optical.design.response <- function(
    design_file, # the file containing the design matrix
    response_file # the file containing the response vector
) {
  Time <- 5e4   # from Lazhi
  
  scipy_sparse = import("scipy.sparse")
  csr_matrix = scipy_sparse$load_npz("/Users/jl1023/NGC2516/27/cochran/optDesign.txt.npz")
  sparse_design <- Matrix(as.matrix(csr_matrix), sparse = TRUE)
  
  response <- read.table(response_file)[, 2]
  list("Design" = sparse_design, "Response" = response, "Time" = Time)
}


###################
# matrix splitting
##################
matrixSplit <- function(Data) { # splitting Data according to overlapping information
  n.segments <- rowSums(Data$Segments$src.index != 0) # number of sources each segment has
  multipleIndex <- which(n.segments != 1) # which segments have multiple sources
  overlapSource <- Data$Segments$src.index[multipleIndex, ] # the sources that are overlapping
  overlapArray <- sort(unique(c(overlapSource))) # vector of overlapping sources
  isolateSource <- which(!1:Data$Sources$n %in% overlapArray) # the sources that are isolated from all other sources
  isolateSegment <- which(Data$Segments$src.index[, 1] %in% isolateSource)
  
  # verification - the sources that are contained/totally overlapping with other sources should not appear in isolateSources
  singleIndex <- which(n.segments == 1) # which segments have single sources
  singleSegment <- Data$Segments$src.index[singleIndex]
  containedSources <- which(!1:Data$Sources$n %in% singleSegment)
  containing <- !any(containedSources %in% isolateSource) # should be TRUE
  lengthing <- length(isolateSegment) == length(isolateSource) # the length of isolated sources should = the length of isolated segments
  if (!containing) {
    stop("Containation error - Fully contained source considered as isolated")
  }
  if (!lengthing) {
    stop("Length error - the number of isolated sources and segments do not match")
  }
  
  isoData <- list()
  isoData$src.index <- Data$Segments$src.index[isolateSegment, ]
  isoData$area <- Data$Segments$area[isolateSegment]
  isoData$eff.rea <- Data$Segments$eff.area[isolateSegment]
  isoData$cnt.obs <- Data$Segments$cnt.obs[isolateSegment]
  isoData$r <- Data$Segments$r[isolateSegment, ]
  
  multiData <- list()
  multiData$src.index <- Data$Segments$src.index[-isolateSegment, ]
  multiData$area <- Data$Segments$area[-isolateSegment]
  multiData$eff.rea <- Data$Segments$eff.area[-isolateSegment]
  multiData$cnt.obs <- Data$Segments$cnt.obs[-isolateSegment]
  multiData$r <- Data$Segments$r[-isolateSegment, ]
  
  output <- list("Isolate" = isoData, "Multiple" = multiData)
  output
}
oSData <- matrixSplit(oData)
plot(density(oSData$Isolate$r[, 1])) # the density should peak around 0.9

xSData <- matrixSplit(xData)
all(xSData$Isolate$r[, 1] == 0.9) # the density should peak around 0.9
