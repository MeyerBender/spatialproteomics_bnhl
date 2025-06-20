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

# Compute statistics for each cell type
statsList <- tryCatch(
  {
    lapply(unique(data$labels_2), function(ct) {
      subset_data <- data[data$labels_2 == ct,]

      # Create point pattern within the polygon window
      P <- ppp(subset_data$centroid.1, subset_data$centroid.0, window = ow)

      # Compute spatial statistics
      G <- allstats(P)
      D <- as.data.frame(G)

      print("Running cbind with D and ct as follows")
      print(D)
      print(ct)

      cbind(D, CA=ct)
    })
  },
  error = function(e) {
    cat("Error: Creating empty dataframe due to an issue.\n")

    columns <- c("F.function.r","F.function.theo","F.function.cs","F.function.rs","F.function.km",
                 "F.function.hazard","F.function.theohaz","G.function.r","G.function.theo",
                 "G.function.han","G.function.rs","G.function.km","G.function.hazard",
                 "G.function.theohaz","J.function.r","J.function.theo","J.function.rs",
                 "J.function.han","J.function.km","J.function.hazard","K.function.r",
                 "K.function.theo","K.function.border","CA","K.function.trans",
                 "K.function.iso","sample")
    df <- data.frame(matrix(ncol = length(columns), nrow = 0))
    colnames(df) <- columns
    write.csv(df, opt$output, row.names = FALSE)
    quit(status = 0)
  }
)

# Combine results and save output
statsDf <- bind_rows(statsList)
statsDf <- cbind(statsDf, sample=opt$sample)
write.csv(statsDf, opt$output, row.names = FALSE)
