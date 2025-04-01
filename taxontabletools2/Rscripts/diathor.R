#!/usr/bin/env Rscript

# Function to check and install missing packages
install_if_missing <- function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, repos = "https://cloud.r-project.org/")
    }
    library(pkg, character.only = TRUE)
}

# List of required packages
packages <- c("diathor", "openxlsx", "readxl", "optparse")

# Install missing packages
sapply(packages, install_if_missing)

# Parse command-line arguments
library(optparse)
option_list <- list(
    make_option(c("-i", "--input"), type = "character", help = "Input CSV file", metavar = "FILE"),
    make_option(c("-o", "--output"), type = "character", help = "Output XLSX file", metavar = "FILE")
)
opt <- parse_args(OptionParser(option_list = option_list))

# Validate input
if (is.null(opt$input) || is.null(opt$output)) {
    stop("Error: Both -i (input) and -o (output) arguments are required.")
}

# Load input CSV
df1 <- read.csv(opt$input)

# Run DiaThor analysis
allResults <- diaThorAll(df1, resultsPath = dirname(opt$output), vandamReports = FALSE)

# Write output to Excel
write.xlsx(allResults, file = opt$output, sheetName = "Sheet1", rowNames = TRUE)

# Print success message
cat("Analysis complete. Output saved to:", opt$output, "\n")