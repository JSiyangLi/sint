library(scales)
################
# Xray detection
################
# radii
# Simpler version if all lines strictly follow "circle(x, y, r)" format
xfile_path <- "/Users/jl1023/NGC2516/27/repro/b1_xray.xyreg"
xlines <- readLines(xfile_path)

# Extract all numbers from each line
xradii <- sapply(strsplit(xlines, ","), function(x) {
  as.numeric(sub("\\)", "", x[3]))  # Remove closing parenthesis and convert
})

xradius_vector <- unname(xradii)

##################
# area
calculate_source_areas <- function(seg_src_matrix, seg_areas) {
  # Create an empty list to store source indices and areas
  source_data <- list()
  
  # Process each segment (row in the matrix)
  for (seg_idx in 1:nrow(seg_src_matrix)) {
    # Get sources for this segment (ignore zeros)
    sources <- seg_src_matrix[seg_idx, ]
    sources <- sources[sources != 0]
    
    # Get area of current segment
    area <- seg_areas[seg_idx]
    
    # For each source in this segment, add the segment area
    for (src in sources) {
      if (is.null(source_data[[as.character(src)]])) {
        source_data[[as.character(src)]] <- 0
      }
      source_data[[as.character(src)]] <- source_data[[as.character(src)]] + area
    }
  }
  
  # Convert to data frame
  sources <- as.integer(names(source_data))
  areas <- unlist(source_data)
  result <- data.frame(source = sources, area = areas)
  
  # Sort by source index
  result <- result[order(result$source), ]
  
  return(result)
}
xarea_vec <- calculate_source_areas(xData$Segments$src.index, xData$Segments$area)[, 2]
xradii_areac <- sqrt(xarea_vec / pi)

plot(density(xradius_vector), col = alpha("blue", alpha = 0.8), main = "radii in Xray data computation",
     sub = "blue=original, red=segment reconstruction")
lines(density(xradii_areac), col = alpha("red", alpha = 0.8), lty = "dashed")

###################
# optical catalogue
###################
# radii
# Simpler version if all lines strictly follow "circle(x, y, r)" format
ofile_path <- "/Users/jl1023/NGC2516/27/repro/cbind_sub_cleared.xyreg"
olines <- readLines(ofile_path)

# Extract all numbers from each line
oradii <- sapply(strsplit(olines, ","), function(x) {
  as.numeric(sub("\\)", "", x[3]))  # Remove closing parenthesis and convert
})

oradius_vector <- unname(oradii)

##################
# area
oarea_vec <- calculate_source_areas(oData$Segments$src.index, oData$Segments$area)[, 2]
oradii_areac <- sqrt(oarea_vec / pi)

plot(density(oradius_vector), col = alpha("blue", alpha = 0.8), main = "radii in optical data computation",
     sub = "blue=original, red=segment reconstruction")
lines(density(oradii_areac), col = alpha("red", alpha = 0.8), lty = "dashed")
