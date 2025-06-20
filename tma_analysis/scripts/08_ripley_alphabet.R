# Load necessary library
#if(!require(optparse)) install.packages("optparse", dependencies=TRUE, repos = "http://cran.us.r-project.org")
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
cat("RIPLEY ALPHABET Successfully read the CSV file!\n")

print(head(data))

cvxhull <- convexhull.xy(cbind(data$centroid.1, y=data$centroid.0))
statsList <- tryCatch(
  {
    lapply(unique(data$labels_3), function(ct) {
      subset_data <- data[data$labels_3 == ct,]
      P <- ppp(subset_data$centroid.1, subset_data$centroid.0, poly = cvxhull$bdry[[1]])
      G <- allstats(P)
      D <- as.data.frame(G)
      print("Running cbind with D and ct as follows")
      print(D)
      print(ct)
      cbind(D, CA=ct)
      })
  },
    # catch errors
    error = function(e) {
      cat("Creating empty df!\n")

      columns <- c("F.function.r","F.function.theo","F.function.cs","F.function.rs","F.function.km","F.function.hazard","F.function.theohaz","G.function.r","G.function.theo","G.function.han","G.function.rs","G.function.km","G.function.hazard","G.function.theohaz","J.function.r","J.function.theo","J.function.rs","J.function.han","J.function.km","J.function.hazard","K.function.r","K.function.theo","K.function.border","CA","K.function.trans","K.function.iso","sample")
      cat("Cols created!\n")
      df <- data.frame(matrix(ncol = length(columns), nrow = 0))
      cat("Renaming cols!\n")
      colnames(df) <- columns
      cat("wrtingcols!\n")
      write.csv(df, opt$output)
      quit(status = 0)
    }
)


statsDf <- bind_rows(statsList)
statsDf <- cbind(statsDf, sample=opt$sample)
write.csv(statsDf, opt$output, row.names = FALSE)