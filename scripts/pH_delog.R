# Load required packages
library(tidyverse)
library(readr)

# Read the combined data file - updated for MinMax normalized data
data_path <- "F:/Uni/Bachelorarbeit/R working directory/enmap-soil-correlation/data/combined_data_minmax.csv"
combined_data <- read_csv(data_path)

# Check if pH column exists
if ("pH" %in% colnames(combined_data)) {
  # Convert pH to hydrogen ion concentration [H+] = 10^(-pH)
  combined_data <- combined_data %>%
    mutate(H_concentration = 10^(-pH))
  
  # Save the transformed data - updated filename to indicate MinMax normalization
  output_path <- "F:/Uni/Bachelorarbeit/R working directory/enmap-soil-correlation/data/combined_data_minmax_with_H_concentration.csv"
  write_csv(combined_data, output_path)
  
  cat("Hydrogen ion concentrations have been calculated from pH values.\n")
  cat("New column 'H_concentration' added to the dataset.\n")
  cat("Transformed data saved to:", output_path, "\n")
  
  # Print the first few rows of the transformation to verify
  print(head(combined_data %>% select(pH, H_concentration)))
} else {
  cat("Error: pH column not found in the dataset.\n")
}
