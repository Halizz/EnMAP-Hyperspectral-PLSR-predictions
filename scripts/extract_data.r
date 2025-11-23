# Load necessary libraries
library(tools)

# Function to extract .zip files (nested within .tar.gz)
extract_zip_files <- function(base_dir) {
  # Find all .zip files recursively
  zip_files <- list.files(base_dir, pattern = "\\.zip$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  
  if (length(zip_files) == 0) {
    cat("No .zip files found to extract.\n")
    return(invisible(NULL))
  }
  
  cat("Found", length(zip_files), ".zip files to extract\n")
  
  # Extract each .zip file
  for (zip_file in zip_files) {
    # Create extraction directory (same location as zip file, without .zip extension)
    zip_base_name <- tools::file_path_sans_ext(basename(zip_file))
    zip_extract_path <- file.path(dirname(zip_file), zip_base_name)
    
    # Skip if already extracted (check if directory exists and has files)
    if (dir.exists(zip_extract_path) && length(list.files(zip_extract_path)) > 0) {
      cat("  Skipping already extracted .zip:", basename(zip_file), "\n")
      next
    }
    
    cat("  Extracting .zip file:", basename(zip_file), "\n")
    
    tryCatch({
      # Create extraction directory
      dir.create(zip_extract_path, showWarnings = FALSE, recursive = TRUE)
      
      # Extract zip file
      unzip(zip_file, exdir = zip_extract_path, overwrite = FALSE)
      
      cat("  ✓ Successfully extracted:", basename(zip_file), "\n")
      
      # Optional: Remove .zip file after successful extraction to save space
      # file.remove(zip_file)
      
    }, error = function(e) {
      cat("  ✗ Error extracting .zip file:", basename(zip_file), "\n")
      cat("    Error:", e$message, "\n")
    })
  }
  
  cat("Completed extraction of all .zip files\n\n")
}

# Function to extract .tar.gz files
extract_tar_gz <- function(tar_dir, extract_dir) {
  # Check if the tar_dir exists
  if (!dir.exists(tar_dir)) {
    stop(paste("The directory", tar_dir, "does not exist. Please ensure that the .tar.gz files are located in the specified folder."))
  }

  # List all .tar.gz files
  tar_files <- list.files(tar_dir, pattern = "\\.tar\\.gz$", full.names = TRUE)

  if (length(tar_files) == 0) {
    stop("No .tar.gz files found in the specified directory. Please ensure that the .tar.gz files are located in the specified folder.")
  }

  # Create extraction directory if it doesn't exist
  if (!dir.exists(extract_dir)) {
    dir.create(extract_dir, recursive = TRUE)
  }

  # Extract each .tar.gz file
  for (tar_file in tar_files) {
    tar_base_name <- tools::file_path_sans_ext(basename(tar_file))  # Remove .tar.gz extension
    tar_extract_path <- file.path(extract_dir, tar_base_name)  # Expected extraction folder

    # Skip extraction if already processed
    if (dir.exists(tar_extract_path)) {
      cat("Skipping already extracted .tar.gz file:", tar_file, "\n")
      next
    }

    cat("Extracting .tar.gz file:", tar_file, "\n")
    tryCatch({
      untar(tar_file, exdir = tar_extract_path)
      cat("Successfully extracted:", tar_file, "\n")
    }, error = function(e) {
      cat("Error extracting .tar.gz file:", tar_file, "\n")
      print(e)
    })

    # Pause for a few seconds to avoid reading problems
    cat("Pausing for 2 seconds...\n")
    Sys.sleep(2)
  }
}

# Main extraction function
run_extraction <- function(tar_dir, extract_dir, pause_duration = 5) {
  
  # Step 1: Extract .tar.gz files
  cat("════════════════════════════════════════════════════════════════\n")
  cat("STEP 1: Extracting .tar.gz archives\n")
  cat("════════════════════════════════════════════════════════════════\n")
  extract_tar_gz(tar_dir, extract_dir)
  cat("✓ .tar.gz extraction completed\n\n")
  
  # Step 2: Extract nested .zip files
  cat("════════════════════════════════════════════════════════════════\n")
  cat("STEP 2: Extracting nested .zip files from .tar.gz archives\n")
  cat("════════════════════════════════════════════════════════════════\n")
  extract_zip_files(extract_dir)
  cat("✓ .zip extraction completed\n\n")
  
  cat("════════════════════════════════════════════════════════════════\n")
  cat("✓ EXTRACTION PROCESS COMPLETED\n")
  cat("════════════════════════════════════════════════════════════════\n")
  cat("Data is available in:", extract_dir, "\n\n")
}

# HINWEIS: Auto-Execution auskommentiert - Script wird nur über master_pipeline.R gesteuert
# # Example usage
# run_extraction(
#   tar_dir = "data/ENMAP", 
#   extract_dir = "data/extracted_files", 
#   pause_duration = 5
# )

