# Load necessary libraries
library(optparse)
library(spatstat.explore)
library(dplyr)

# Define command-line options
option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Path to the CSV file containing centroids", metavar="character"),
  make_option(c("-p", "--polygon"), type="character", default=NULL, 
              help="Path to the CSV file containing polygon boundaries", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Path to the output file", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="Sample name", metavar="character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are provided
if (is.null(opt$file) || is.null(opt$polygon)) {
  cat("Error: Please provide the paths to the centroid and polygon CSV files using -f and -p.\n")
  print_help(opt_parser)
  quit(status = 1)
}

# Read the centroid data CSV
data <- tryCatch(
  {
    read.csv(opt$file)
  },
  error = function(e) {
    cat("Error: Unable to read the centroid CSV file. Please check the file path.\n")
    quit(status = 1)
  }
)

# Read the polygon CSV file
polygon_data <- tryCatch(
  {
    read.csv(opt$polygon)
  },
  error = function(e) {
    cat("Error: Unable to read the polygon CSV file. Please check the file path.\n")
    quit(status = 1)
  }
)

# Ensure polygon_data has correct structure
if (!all(c("id", "x", "y") %in% colnames(polygon_data))) {
  cat("Error: Polygon CSV must contain 'id', 'x', and 'y' columns.\n")
  quit(status = 1)
}

# Function to check and correct polygon orientation
correct_polygon_orientation <- function(poly) {
  # Shoelace formula to compute signed area
  area <- sum((poly$x[-1] + poly$x[-length(poly$x)]) * (poly$y[-1] - poly$y[-length(poly$y)])) / 2
  
  cat("Original Area:", area, "\n")  # Debugging print
  
  # If the area is negative, the points are clockwise and need to be reversed
  if (area < 0) {
    poly$x <- rev(poly$x)
    poly$y <- rev(poly$y)
  }
  
  # Ensure polygon is closed (append first point at the end if missing)
  if (!(poly$x[1] == poly$x[length(poly$x)] && poly$y[1] == poly$y[length(poly$y)])) {
    poly$x <- c(poly$x, poly$x[1])
    poly$y <- c(poly$y, poly$y[1])
  }
  
  # Recalculate area after correction
  area_after <- sum((poly$x[-1] + poly$x[-length(poly$x)]) * (poly$y[-1] - poly$y[-length(poly$y)])) / 2
  cat("Corrected Area:", area_after, "\n")  # Debugging print
  
  return(poly)
}

# Convert polygon data to spatstat format, correcting orientation if needed
polygon_list <- lapply(unique(polygon_data$id), function(pid) {
  subset_poly <- polygon_data[polygon_data$id == pid, ]
  
  # Ensure the polygon is properly oriented
  correct_polygon_orientation(list(x = subset_poly$x, y = subset_poly$y))
})

# Check if any valid polygons exist
if (length(polygon_list) == 0) {
  cat("Error: No valid polygons found in the polygon CSV.\n")
  quit(status = 1)
}

# Print success message
cat("Successfully read and corrected the polygons!\n")
print(head(data))
print(head(polygon_data))

# Create the observation window using corrected polygon data
ow <- owin(poly = polygon_list)

# Define the point pattern with polygon-based edge correction
P <- ppp(data$centroid.1, data$centroid.0, window = ow, marks=factor(data$labels_2))

# Generate all label combinations for cross-type calculations
combinations <- expand.grid(unique(data$labels_2), unique(data$labels_2))

# Compute cross-type statistics
cross <- tryCatch(
  {
    lapply(1:nrow(combinations), function(k) {
      a <- combinations[k, 1]
      b <- combinations[k, 2]
      
      G1 <- Kcross(P, a, b)
      D1 <- as.data.frame(G1)
      colnames(D1) <- paste0("Kcross.", colnames(D1))
      
      G2 <- Gcross(P, a, b)
      D2 <- as.data.frame(G2)
      colnames(D2) <- paste0("Gcross.", colnames(D2))  
      
      G3 <- Lcross(P, a, b)
      D3 <- as.data.frame(G3)
      colnames(D3) <- paste0("Lcross.", colnames(D3))  
      
      G4 <- Jcross(P, a, b)
      D4 <- as.data.frame(G4)
      colnames(D4) <- paste0("Jcross.", colnames(D4))

      D <- cbind(D1, D2, D3, D4, CA=a, CB=b)
    })
  },
  error = function(e) {
    cat("Error: Unable to calculate Kcross, Gcross, Lcross, Jcross.\n")

    columns <- c("Kcross.r","Kcross.theo","Kcross.border","Gcross.r","Gcross.theo",
                 "Gcross.han","Gcross.rs","Gcross.km","Gcross.hazard","Gcross.theohaz",
                 "Lcross.r","Lcross.theo","Lcross.border","Jcross.r","Jcross.theo",
                 "Jcross.rs","Jcross.han","Jcross.km","Jcross.hazard","CA","CB",
                 "Kcross.trans","Kcross.iso","Lcross.trans","Lcross.iso","sample")
    df <- data.frame(matrix(ncol = length(columns), nrow = 0))
    colnames(df) <- columns
    write.csv(df, opt$output, row.names = FALSE)
    quit(status = 0)    
  }
)

# Combine results and save output
crossDf <- bind_rows(cross)
crossDf <- cbind(crossDf, sample=opt$sample)
write.csv(crossDf, opt$output, row.names = FALSE)
