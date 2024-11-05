# Set the directory where the files are located
file_directory <- "C:/Users/feiyi/Desktop/github_MEMs/MEMs/Simulation/Results/Equal/MCMC/1_1_1_1/250"  # Replace with the actual path

# Generate the expected file names with ".csv" extension
expected_files <- paste0("mcmc.result.", 1:500, ".csv")

# List the actual files in the directory
actual_files <- list.files(file_directory, pattern = "mcmc.result.*\\.csv")

# Find missing files
missing_files <- setdiff(expected_files, actual_files)

# Print the missing file numbers
if (length(missing_files) > 0) {
  missing_numbers <- as.numeric(gsub("mcmc.result.|\\.csv", "", missing_files))
  print(sort(missing_numbers))
} else {
  print("No files are missing.")
}