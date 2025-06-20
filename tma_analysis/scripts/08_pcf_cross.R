# Load necessary library
# if(!require(optparse)) install.packages("optparse", dependencies=TRUE, repos = "http://cran.us.r-project.org")
library(optparse)
library(spatstat.explore)
library(dplyr)
# Define command-line options
option_list <- list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="Path to the CSV file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Path to the output file", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="Sample name", metavar="character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if the file argument is provided
if (is.null(opt$file)) {
  print("Error: Please provide the path to the CSV file using the -f or --file argument.")
  print_help(opt_parser)
  quit(status = 1)
}

# Read the CSV file
data <- tryCatch(
  {
    read.csv(opt$file)
  },
  error = function(e) {
    cat("Error: Unable to read the CSV file. Please check the file path.\n")
    quit(status = 1)
  }
)

# Print the first few rows of the CSV data
cat("Successfully read the CSV file!\n")
print(head(data))

cvxhull <- convexhull.xy(cbind(data$centroid.1, y=data$centroid.0))
P <- ppp(data$centroid.1, data$centroid.0, poly = cvxhull$bdry[[1]], marks=factor(data$labels_3))


combinations <- list()
k = 1
for (i in unique(data$labels_3)) {
  for (j in unique(data$labels_3)){
    #if (i != j) {
      combinations[[k]] <- c(i, j)
      k <- k + 1 
    #}
  }
}

cross <- lapply(combinations, function(x) {
  a <- x[[1]]
  b <- x[[2]]
  
  G <- pcfcross(P, a, b)
  D <- as.data.frame(G)
  D <- cbind(D, CA=a, CB=b)
})

crossDf <- bind_rows(cross)
crossDf <- cbind(crossDf, sample=opt$sample)
write.csv(crossDf, opt$output, row.names = FALSE)
