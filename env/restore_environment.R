# restore_environment.R

# Read the CSV file into a data frame
package_info <- read.csv("R_env.csv", stringsAsFactors = FALSE)

# Install packages from the data frame
for (row in 1:nrow(package_info)) {
  install.packages(
    package_info[row, "Package"], 
    version = package_info[row, "Version"],
    dependencies = TRUE,
    repos = "http://cran.us.r-project.org"
  )
}

BiocManager::install("BiocVersion")
BiocManager::install("BiocStyle")
BiocManager::install("NuPoP")

# Load installed packages
lapply(package_info$Package, library, character.only = TRUE)