# This file acts as a data preparation script for hyperspectral data analysis.
# It includes functions to extract and process .tar.gz files,
# load hyperspectral data, apply sun angle correction, remove cloud points,
# and group/filter data based on specified conditions.
# The final prepared data is saved as an RDS file and a CSV file.

# Load necessary libraries
library(raster)
library(readxl)
library(dplyr)
library(sf)
library(tidyr)
library(terra)
library(xml2)
library(data.table) 


# Define the path to the extracted files
extract_dir <- "data/extracted_files"

# Check if NORMALIZATION_METHOD is defined, otherwise default to "minmax"
if (!exists("NORMALIZATION_METHOD")) {
  NORMALIZATION_METHOD <- "minmax"
  cat("⚠️  NORMALIZATION_METHOD not defined, defaulting to 'minmax'\n")
}
cat("📊 Using normalization method:", toupper(NORMALIZATION_METHOD), "\n")

# Debugging function to list all files in a directory recursively
list_all_files <- function(directory) {
  cat("Listing all files in the directory:", directory, "\n")
  # Check if the directory exists
  if (!dir.exists(directory)) {
    stop("The specified directory does not exist:", directory)
  }

  # List all files recursively
  all_files <- list.files(directory, recursive = TRUE, full.names = TRUE)

  print(all_files)  # Debugging: Print all found files

  # Check if any files were found
  if (length(all_files) == 0) {
    cat("No files found in the directory.\n")
  } else {
    cat("Files found:\n", paste(all_files, collapse = "\n"), "\n")
  }
}

# Function to extract coordinates from the shapefile
extract_coordinates <- function(shapefile_path) {
  cat("Extracting coordinates from the shapefile:", shapefile_path, "\n")
  shp <- st_read(shapefile_path)
  coordinates <- st_coordinates(shp) %>% as.data.frame()
  coordinates$Name <- shp$Name  # Assumption: "Name" is a column in the shapefile
  return(coordinates)
}


# Function to load raster data from processed spectral files
load_raster_data <- function(extract_dir, coordinates) {
  cat("Loading hyperspectral raster data from directory (processed files):", extract_dir, "\n")

  # Build pattern based on normalization method
  if (NORMALIZATION_METHOD == "minmax") {
    pattern <- "_smoothed\\+adaptive_destriped_masked_snv_minmax\\.tif$"
    method_name <- "MinMax"
  } else if (NORMALIZATION_METHOD == "snv") {
    pattern <- "_smoothed\\+adaptive_destriped_masked_snv\\.tif$"
    method_name <- "SNV"
  } else {
    stop("❌ Invalid NORMALIZATION_METHOD: ", NORMALIZATION_METHOD, ". Must be 'snv' or 'minmax'.\n")
  }
  
  cat("Looking for", method_name, "normalized images...\n")
  spectral_files <- list.files(
    extract_dir,
    pattern = pattern,
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  if (length(spectral_files) == 0) {
    stop("❌ No ", method_name, " normalized images found in directory: ", extract_dir, "\n")
  } else {
    cat("✓ Found", length(spectral_files), method_name, "normalized image files\n")
  }

  # Debugging: Print found spectral files
  cat("Found spectral files:\n")
  print(spectral_files)

  # Initialize a data frame for combined valid coordinates
  combined_valid_coordinates <- data.frame()

  for (file in spectral_files) {
    cat("Processing raster file:", file, "\n")

    # Load the raster data
    raster_data <- tryCatch({
      rast(file)
    }, error = function(e) {
      cat("Error loading raster file:", file, "\n", e$message, "\n")
      NULL
    })

    if (is.null(raster_data)) {
      next
    }

    # Compute the extent of the raster
    raster_extent <- tryCatch({
      ext(raster_data)
    }, error = function(e) {
      cat("Error computing extent for raster file:", file, "\n", e$message, "\n")
      NULL
    })

    if (!is.null(raster_extent)) {
      # Filter valid coordinates within the raster extent
      valid_coords <- coordinates %>%
        filter(
          X >= raster_extent$xmin & X <= raster_extent$xmax &
          Y >= raster_extent$ymin & Y <= raster_extent$ymax
        )

      # Debugging: Check the number of valid coordinates
      cat("Number of valid coordinates for file:", file, ":", nrow(valid_coords), "\n")

      if (nrow(valid_coords) > 0) {
        # Extract the base name without suffixes for matching
        base_name <- basename(file)
        
        # Remove suffix based on normalization method
        if (NORMALIZATION_METHOD == "minmax") {
          base_name <- gsub("_smoothed\\+adaptive_destriped_masked_snv_minmax\\.tif$", "", base_name, ignore.case = TRUE)
        } else if (NORMALIZATION_METHOD == "snv") {
          base_name <- gsub("_smoothed\\+adaptive_destriped_masked_snv\\.tif$", "", base_name, ignore.case = TRUE)
        }
        
        # Add valid coordinates along with the cleaned raster file name
        valid_coords <- valid_coords %>%
          mutate(RasterFile = base_name)
        combined_valid_coordinates <- rbind(combined_valid_coordinates, valid_coords)
        cat("Valid coordinates added for file:", file, "\n")
      } else {
        cat("No valid coordinates found for file:", file, "\n")
      }
    }
  }

  # Save combined valid coordinates to a single CSV file
  output_file <- "data/combined_valid_coords.csv"
  write.csv(combined_valid_coordinates, output_file, row.names = FALSE)
  cat("Combined valid coordinates saved to:", output_file, "\n")

  # Return the combined valid coordinates
  return(list(valid_coordinates = combined_valid_coordinates))
}

# Function to load hyperspectral data - Modified to use processed images
load_hyperspectral_data <- function(extract_dir, valid_coordinates) {
  method_name <- ifelse(NORMALIZATION_METHOD == "minmax", "MinMax", "SNV")
  cat("Starting hyperspectral data extraction from", method_name, "normalized images...\n")

  if (is.null(extract_dir) || extract_dir == "") {
    stop("Error: extract_dir is NULL or empty.")
  }

  # Find unique raster files in valid_coordinates
  unique_raster_files <- unique(valid_coordinates$RasterFile)
  
  cat("Number of unique raster files found:", length(unique_raster_files), "\n")
  cat("Unique raster files:", paste(unique_raster_files, collapse = ", "), "\n")

  # Initialize spectral_data with the desired column names
  spectral_data <- data.frame(matrix(ncol = 230, nrow = 0))  # 228 columns
  colnames(spectral_data) <- c(
    "X", "Y", "Name", "ID", paste0("Band", 1:224), "RasterFile", "Redundancy"
  )

  # Ensure Name column in valid_coordinates is character type
  valid_coordinates$Name <- as.character(valid_coordinates$Name)
  
  # Debug: Check the data type of Name column
  cat("Data type of Name column in valid_coordinates:", class(valid_coordinates$Name), "\n")

  for (raster_file in unique_raster_files) {
    cat("Processing raster file:", raster_file, "\n")

    # Remove any quotes from the raster file name
    raster_file_clean <- gsub('"', '', raster_file)

    # Build search pattern based on normalization method
    if (NORMALIZATION_METHOD == "minmax") {
      file_pattern <- paste0(raster_file_clean, ".*_smoothed\\+adaptive_destriped_masked_snv_minmax\\.tif$")
    } else {
      file_pattern <- paste0(raster_file_clean, ".*_smoothed\\+adaptive_destriped_masked_snv\\.tif$")
    }
    
    raster_path <- list.files(
      extract_dir,
      pattern = file_pattern,
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )

    cat("Found", length(raster_path), "matching files for", raster_file_clean, "\n")
    if (length(raster_path) > 0) {
      cat("Using file:", raster_path[1], "\n")
    } else {
      cat("ERROR: No matching file found for:", raster_file_clean, "\n")
      cat("Available files in directory (first 10):\n")
      all_files <- list.files(extract_dir, pattern = "\\.tif$", recursive = TRUE)
      cat(paste(head(all_files, 10), collapse = "\n"), "\n")
    }

    # Load the raster data
    if (length(raster_path) > 0) {
      raster_data <- tryCatch({
        rast(raster_path[1])  # Use first match
      }, error = function(e) {
        cat("Error loading raster file:", raster_path[1], "\n", e$message, "\n")
        NULL
      })
    } else {
      cat("No matching file found for:", raster_file_clean, "\n")
      raster_data <- NULL
    }

    cat("Loaded raster data for file:", raster_file, "\n")
    print(raster_data)  # Debugging: Print raster data object

    if (is.null(raster_data)) {
      next
    }
    # Filter coordinates for the current raster file
    coords_for_raster <- valid_coordinates %>%
      filter(RasterFile == raster_file) %>%
      select(X, Y, Name)  # Keep the 'Name' column

    # Debugging: Check the structure of coords_for_raster
    cat("Coordinates for raster file:\n")
    print(head(coords_for_raster))

    # Extract spectral bands for the coordinates (only pass X and Y)
    extracted_spectra <- tryCatch({
      terra::extract(raster_data, coords_for_raster[, c("X", "Y")], df = TRUE)  # Only pass X and Y
    }, error = function(e) {
      cat("Error extracting spectral data from raster file:", raster_file, "\n", e$message, "\n")
      NULL
    })

    if (is.null(extracted_spectra) || nrow(extracted_spectra) == 0) {
      cat("No spectral data extracted for raster file:", raster_file, "\n")
      next
    }

    # Add the 'Name' column from coords_for_raster
    extracted_spectra <- cbind(coords_for_raster["Name"], extracted_spectra)

    # Debugging: Check the structure of extracted_spectra
    cat("Extracted spectra with Name column:\n")
    print(head(extracted_spectra))

    # Ensure column names are valid and consistent
    colnames(extracted_spectra) <- make.names(colnames(extracted_spectra), unique = TRUE)

    # Rename band columns to a consistent format (Band1, Band2, ..., BandN)
    band_columns <- grep("SPECTRAL_IMAGE|^Band", colnames(extracted_spectra), value = TRUE)
    cat("band_header:\n")
    print(head(band_columns))  # Debugging: Print band columns
    new_band_names <- paste0("Band", seq_along(band_columns))
    colnames(extracted_spectra)[match(band_columns, colnames(extracted_spectra))] <- new_band_names
   # Debugging: Print updated column names

    cat("Updated column names of extracted_spectra:\n")
    print(colnames(extracted_spectra))

    # Combine extracted spectra with coordinates and file name
    combined_data <- cbind(
      coords_for_raster,
      extracted_spectra[, -1],  # Exclude the cell index column
      RasterFile = raster_file
    )

    # Add redundancy column
    combined_data <- combined_data %>%
      group_by(X, Y) %>%
      mutate(Redundancy = n() > 1) %>%
      ungroup()

    print(colnames(combined_data))  # Debugging: Print combined data


    # Ensure column names of combined_data match spectral_data
    if (!all(colnames(combined_data) %in% colnames(spectral_data))) {
      cat("Error: Column names in combined_data do not match spectral_data.\n")
      cat("Columns in combined_data but not in spectral_data:\n")
      print(setdiff(colnames(combined_data), colnames(spectral_data)))
      cat("Columns in spectral_data but not in combined_data:\n")
      print(setdiff(colnames(spectral_data), colnames(combined_data)))
      stop("Column name mismatch detected. Aborting.")
    }
    cat("Column names match successfully.\n")
    # Reorder columns in combined_data to match spectral_data
    combined_data <- combined_data[, colnames(spectral_data), drop = FALSE]

    # Append to the spectral data
    spectral_data <- rbind(spectral_data, combined_data)

    cat("Current spectral data structure:\n")
    print(colnames(spectral_data))
  }

  # Check if spectral_data is empty
  if (nrow(spectral_data) == 0) {
    stop("Error: No spectral data could be loaded or extracted.")
  }

  # Remove rows without spectral data
  cat("Removing rows without spectral data...\n")

  # Check if all spectra (Band1 to Band224) are NA, and remove such rows
  spectral_data <- spectral_data %>%
    filter(rowSums(!is.na(across(starts_with("Band")))) > 0)


  # After the loop over raster files and before writing spectral_data
  cat("Checking redundancy values for RasterFile...\n")

  # Set redundancy value based on multiple occurrences of RasterFile names
  spectral_data <- spectral_data %>%
    group_by(X, Y) %>%
    mutate(Redundancy = n()) %>%
    ungroup()

  cat("Updated redundancy values:\n")
  print(table(spectral_data$Redundancy))

  # Remove specific suffixes from the RasterFile column - based on normalization method
  if (NORMALIZATION_METHOD == "minmax") {
    spectral_data <- spectral_data %>%
      mutate(RasterFile = gsub("_smoothed\\+adaptive_destriped_masked_snv_minmax\\.tif$", "", RasterFile))
  } else if (NORMALIZATION_METHOD == "snv") {
    spectral_data <- spectral_data %>%
      mutate(RasterFile = gsub("_smoothed\\+adaptive_destriped_masked_snv\\.tif$", "", RasterFile))
  }
  
  # Remove generic suffixes that apply to all methods
  spectral_data <- spectral_data %>%
    mutate(RasterFile = gsub("-SPECTRAL_IMAGE.*\\.tif$", "", RasterFile))

  # Debugging: Check the cleaned RasterFile column
  cat("Cleaned RasterFile column:\n")
  print(head(spectral_data$RasterFile))

  # Save the updated spectral data with indicator of normalization method in filename
  output_file <- paste0("data/spectral_data_", NORMALIZATION_METHOD, ".csv")
  write.csv(spectral_data, output_file, row.names = FALSE)
  cat("Spectral data saved to:", output_file, "\n")

  # Return the combined spectral data
  return(spectral_data)
}

# Function to load metadata from the extracted directory - Updated to handle multiple files
load_metadata <- function(extract_dir, spectral_data) {
  cat("Loading metadata from directory:", extract_dir, "\n")
  
  # Find all metadata XML files
  metadata_files <- list.files(
    extract_dir,
    pattern = "-METADATA\\.XML$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  if (length(metadata_files) == 0) {
    cat("No metadata XML files found in directory. Creating empty metadata mapping.\n")
    return(data.frame(
      SpectralFile = character(0),
      sun_elevation = numeric(0),
      sun_azimuth = numeric(0),
      scene_sza = numeric(0),
      scene_aot = numeric(0),
      scene_wv = numeric(0),
      cloud_cover = numeric(0),
      haze_cover = numeric(0),
      dead_pixels_vnir = numeric(0),
      dead_pixels_swir = numeric(0)
    ))
  }
  
  cat("Found", length(metadata_files), "metadata files\n")
  
  # Initialize list to store metadata entries
  metadata_list <- list()
  
  # Process each metadata file
  for (i in 1:length(metadata_files)) {
    metadata_file <- metadata_files[i]
    cat("Processing metadata file [", i, "/", length(metadata_files), "]: ", basename(metadata_file), "\n")
    
    # Extract base name for matching with spectral files
    base_name <- gsub("-METADATA\\.XML$", "", basename(metadata_file), ignore.case = TRUE)
    
    # Load XML
    metadata <- tryCatch({
      read_xml(metadata_file)
    }, error = function(e) {
      cat("Error reading metadata file:", metadata_file, "\n", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(metadata)) {
      next
    }
    
    # Extract metadata parameters
    tryCatch({
      metadata_params <- list(
        SpectralFile = base_name,
        sun_elevation = as.numeric(xml_text(xml_find_first(metadata, ".//sunElevationAngle/center"))),
        sun_azimuth = as.numeric(xml_text(xml_find_first(metadata, ".//sunAzimuthAngle/center"))),
        scene_sza = as.numeric(xml_text(xml_find_first(metadata, ".//sceneSZA"))),
        scene_aot = as.numeric(xml_text(xml_find_first(metadata, ".//sceneAOT"))),
        scene_wv = as.numeric(xml_text(xml_find_first(metadata, ".//sceneWV"))),
        cloud_cover = as.numeric(xml_text(xml_find_first(metadata, ".//cloudCover"))),
        haze_cover = as.numeric(xml_text(xml_find_first(metadata, ".//hazeCover"))),
        dead_pixels_vnir = as.numeric(xml_text(xml_find_first(metadata, ".//deadPixelsVNIR"))),
        dead_pixels_swir = as.numeric(xml_text(xml_find_first(metadata, ".//deadPixelsSWIR")))
      )
      
      metadata_list[[i]] <- metadata_params
    }, error = function(e) {
      cat("Error extracting metadata from", metadata_file, ":", e$message, "\n")
    })
  }
  
  # Combine all metadata entries into a data frame
  if (length(metadata_list) == 0) {
    cat("No metadata could be successfully extracted. Creating empty metadata mapping.\n")
    metadata_mapping <- data.frame(
      SpectralFile = character(0),
      sun_elevation = numeric(0),
      sun_azimuth = numeric(0),
      scene_sza = numeric(0),
      scene_aot = numeric(0),
      scene_wv = numeric(0),
      cloud_cover = numeric(0),
      haze_cover = numeric(0),
      dead_pixels_vnir = numeric(0),
      dead_pixels_swir = numeric(0)
    )
  } else {
    metadata_mapping <- do.call(rbind, lapply(metadata_list, as.data.frame))
  }
  
  # Save metadata mapping to file for future use
  write.csv(metadata_mapping, "data/metadata_mapping.csv", row.names = FALSE)
  cat("Metadata mapping saved to: data/metadata_mapping.csv\n")
  
  return(metadata_mapping)
}


# Function to apply sun angle correction
apply_sun_correction <- function(spectral_data, sun_elevation) {
  sun_elevation_rad <- sun_elevation * pi / 180  # Convert to radians
  correction_factor <- 1 / sin(sun_elevation_rad)  # Calculate correction factor
  corrected_spectral_data <- spectral_data * correction_factor  # Apply correction
  return(corrected_spectral_data)
}

# Function to combine cloud data with spectral data - Works with both SNV and MinMax
load_cloud_data <- function(extract_dir, spectral_data) {
  method_name <- ifelse(NORMALIZATION_METHOD == "minmax", "MinMax", "SNV")
  cat("Loading cloud data for", method_name, "normalized images.\n")


  # Validate that spectral_data contains the required columns
  if (!all(c("X", "Y", "RasterFile") %in% colnames(spectral_data))) {
    stop("Error: spectral_data must contain 'X', 'Y', and 'RasterFile' columns.")
  }

  # Ensure coordinates in spectral_data are rounded to consistent precision
  spectral_data <- spectral_data %>%
    mutate(across(c(X, Y), ~ round(., digits = 6)))  # Adjust precision as needed

  # Add an empty cloud_value column to spectral_data
  if (!"cloud_value" %in% colnames(spectral_data)) {
    spectral_data <- spectral_data %>%
      mutate(cloud_value = NA_real_)  # Initialize with NA (numeric type)
  }

  cat("Spectral data structure before cloud data processing:\n")
  print(colnames(spectral_data))  # Debugging: Print spectral_data structure

  # List all QUALITY_CLOUD.TIF files in the directory
  cloud_files <- list.files(
    extract_dir,
    pattern = "_QUALITY_CLOUD\\.TIF$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  if (length(cloud_files) == 0) {
    cat("Warning: No QUALITY_CLOUD.TIF files found in the directory. Skipping cloud data processing.\n")
    spectral_data$cloud_value <- NA  # Add an empty cloud_score column
    return(spectral_data)
  }

  cat("Found QUALITY_CLOUD.TIF files:\n", paste(cloud_files, collapse = "\n"), "\n")

  # Initialize a column for cloud scores
  spectral_data$cloud_value <- NA

  # Process each QUALITY_CLOUD.TIF file
  for (file_path in cloud_files) {
    cat("Processing QUALITY_CLOUD.TIF file:", file_path, "\n")


    # Bereinige den Namen des Rasterfiles (ohne Endung -QUALITY_CLOUD.TIF)
    raster_file_name <- gsub("-QL_QUALITY_CLOUD\\.TIF$", "", basename(file_path))

    # Filter spectral_data für das aktuelle RasterFile
    coords_in_extent <- spectral_data %>%
      filter(RasterFile == raster_file_name)

    if (nrow(coords_in_extent) == 0) {
      cat("No matching coordinates found for:", raster_file_name, "\n")
      next
    }

    cat("Matching coordinates found for file:", raster_file_name, "\n")
    print(coords_in_extent)

    # Load the raster data
    quality_cloud_raster <- tryCatch({
      rast(file_path)
    }, error = function(e) {
      cat("Error loading QUALITY_CLOUD.TIF file:", file_path, "\n", e$message, "\n")
      NULL
    })

    cat("Loaded raster data for file:", file_path, "\n")
    print(quality_cloud_raster)  # Debugging: Print the raster data object

    if (is.null(quality_cloud_raster)) {
      next
    }

    # Extract cloud values for the valid coordinates
    extracted_cloud_data <- tryCatch({
      extract(quality_cloud_raster, coords_in_extent[, c("X", "Y")], df = TRUE)
    }, error = function(e) {
      cat("Error extracting data from QUALITY_CLOUD.TIF file:", file_path, "\n", e$message, "\n")
      NULL
    })

    if (!is.null(extracted_cloud_data) && nrow(extracted_cloud_data) > 0) {
      # Bind the coordinates (X, Y) and RasterFile to extracted_cloud_data
      extracted_cloud_data <- cbind(
        coords_in_extent[, c("X", "Y", "RasterFile")],  # Add RasterFile column
        extracted_cloud_data[, -1]  # Remove the ID column
      )

      # Rename the column for cloud values
      colnames(extracted_cloud_data)[4] <- "cloud_value"  # The fourth column contains cloud values

      # Debugging: Check the structure after binding coordinates and RasterFile
      cat("Structure of extracted_cloud_data after binding coordinates and RasterFile:\n")
      print(head(extracted_cloud_data))

      # Debugging: Check the structure after removing duplicates
      cat("Structure of extracted_cloud_data after removing duplicates:\n")
      print(head(extracted_cloud_data))

      # Update the cloud_value column in spectral_data
      spectral_data <- spectral_data %>%
        rows_update(extracted_cloud_data, by = c("X", "Y", "RasterFile"))
    }
  }

  # Fill NA values in cloud_score for coordinates not covered by any raster
  spectral_data$cloud_value[is.na(spectral_data$cloud_value)] <- NA

  cat("Spectral_data with cloud values:\n")
  print(colnames(spectral_data))

  cat("Cloud data processing completed. Cloud values added to spectral_data.\n")
  return(spectral_data)
}

# Compute combined raster extent from all raster files
compute_combined_extent <- function(raster_data_list) {
  combined_extent <- NULL
  for (raster_data in raster_data_list) {
    raster_extent <- tryCatch({
      ext(raster_data)
    }, error = function(e) {
      cat("Error computing extent for a raster file: ", e$message, "\n")
      NULL
    })

    if (!is.null(raster_extent)) {
      if (is.null(combined_extent)) {
        combined_extent <- raster_extent
      } else {
        combined_extent <- union(combined_extent, raster_extent)
      }
    }
  }
  return(combined_extent)
}

# Main data preparation function - Update metadata handling
data_preparation <- function(tar_dir, extract_dir, shapefile_path, soil_data_dir) {
  method_name <- ifelse(NORMALIZATION_METHOD == "minmax", "MinMax", "SNV")
  cat("Starting data preparation process using", method_name, "normalized spectral images.\n")

  # Load shapefile and extract coordinates
  coordinates <- extract_coordinates(shapefile_path)  # Extract coordinates from shapefile
  cat("Total coordinates extracted from shapefile:", nrow(coordinates), "\n")
  cat("Coordinates extracted from shapefile:\n")
  print(head(coordinates))

  # Load all soil measurement files
  soil_files <- list.files(soil_data_dir, pattern = "^soil_measurements.*\\.xlsx$", full.names = TRUE)
  soil_data <- do.call(rbind, lapply(soil_files, read_excel))
  cat("Total soil data points loaded:", nrow(soil_data), "\n")

  # Debugging: Check the 'Name' column after replacement
  cat("Updated 'Name' column in valid_coordinates:\n")
  print(head(soil_data$Name))

  # Debugging: Check structure of soil_data before merging
  cat("Structure of soil_data before merging:\n")
  print(head(soil_data))

  # Load raster data using the load_raster_data function
  raster_data_result <- load_raster_data(extract_dir, coordinates)
  valid_coordinates <- raster_data_result$valid_coordinates
  
  cat("Valid coordinates after raster filtering:", nrow(valid_coordinates), "\n")
  cat("Unique raster files in valid coordinates:", length(unique(valid_coordinates$RasterFile)), "\n")

  # Ensure Name column is character type in both dataframes before conversion
  soil_data$Name <- as.character(soil_data$Name)
  valid_coordinates$Name <- as.character(valid_coordinates$Name)
  
  # Apply consistent transformations to Name column
  soil_data <- soil_data %>%
    mutate(Name = gsub(" ", "_", Name) %>% tolower())

  valid_coordinates <- valid_coordinates %>%
    mutate(Name = gsub(" ", "_", Name) %>% tolower())

  # Debug: Check the data types after transformation
  cat("Data type of Name in soil_data after conversion:", class(soil_data$Name), "\n")
  cat("Data type of Name in valid_coordinates after conversion:", class(valid_coordinates$Name), "\n")

  # Check if valid_coordinates is empty
  if (nrow(valid_coordinates) == 0) {
    stop("Error: No valid coordinates found within the raster extents.")
  }

  # Ensure valid_coordinates contains required columns
  if (!all(c("X", "Y", "Name") %in% colnames(valid_coordinates))) {
    stop("Error: valid_coordinates must contain 'X', 'Y', and 'Name' columns.")
  }

  # Debugging: Check structure of valid_coordinates before merging
  cat("Structure of valid_coordinates before merging:\n")
  print(head(valid_coordinates))



  # Check if spectral_data.csv exists and is not empty
  spectral_data_path <- paste0("data/spectral_data_", NORMALIZATION_METHOD, ".csv")
  spectral_data <- NULL
  force_regenerate <- FALSE
  
  if (file.exists(spectral_data_path)) {
    method_name <- ifelse(NORMALIZATION_METHOD == "minmax", "MinMax", "SNV")
    cat(method_name, "normalized spectral data file exists. Loading from:", spectral_data_path, "\n")
    spectral_data <- fread(spectral_data_path)  # Use fread for efficient reading
    cat("Loaded spectral data points:", nrow(spectral_data), "\n")
    
    # Check if the file is empty or has insufficient data
    if (nrow(spectral_data) == 0) {
      cat("⚠️ WARNING: Existing spectral data file is empty! Will regenerate data.\n")
      # Remove the empty file
      file.remove(spectral_data_path)
      cat("✓ Removed empty spectral data file.\n")
      force_regenerate <- TRUE
    } else {
      # Ensure Name is character type when loading from file
      spectral_data$Name <- as.character(spectral_data$Name)
    }
  }
  
  # If spectral data is empty, NULL, or force_regenerate is TRUE, then generate new spectral data
  if (is.null(spectral_data) || nrow(spectral_data) == 0 || force_regenerate) {
    cat("Generating new spectral data from MinMax normalized images...\n")
    spectral_data <- load_hyperspectral_data(extract_dir, valid_coordinates)
    if (is.null(spectral_data) || nrow(spectral_data) == 0) {
      stop("Error: No spectral data could be loaded or extracted from MinMax normalized images.")
    }
    cat("Extracted spectral data points:", nrow(spectral_data), "\n")
    fwrite(spectral_data, spectral_data_path)  # Save for future use
    cat("✓ New spectral data saved to:", spectral_data_path, "\n")
  }

  # Check if metadata_mapping.csv exists
  metadata_mapping_path <- "data/metadata_mapping.csv"
  if (file.exists(metadata_mapping_path)) {
    cat("Metadata mapping file exists. Loading from:", metadata_mapping_path, "\n")
    tryCatch({
      metadata_mapping <- fread(metadata_mapping_path)  # Use fread for efficient reading
      cat("Loaded metadata mapping with", nrow(metadata_mapping), "entries\n")
    }, error = function(e) {
      cat("Error loading metadata mapping file, will regenerate:", e$message, "\n")
      metadata_mapping <- load_metadata(extract_dir, spectral_data)
    })
  } else {
    cat("Metadata mapping file does not exist. Generating metadata mapping from XML files.\n")
    metadata_mapping <- load_metadata(extract_dir, spectral_data)
  }
  
  # Ensure metadata_mapping is not NULL or empty
  if (is.null(metadata_mapping) || nrow(metadata_mapping) == 0) {
    cat("Warning: No metadata available, creating empty metadata mapping.\n")
    metadata_mapping <- data.frame(
      SpectralFile = spectral_data$RasterFile,  # Use RasterFile values as fallback
      SpectralID = tolower(spectral_data$RasterFile)  # Pre-create SpectralID
    )
  }

  # Clean the SpectralFile column in metadata_mapping
  metadata_mapping <- metadata_mapping %>%
    mutate(
      SpectralID = gsub("-spectral_image(\\.(tif|bsq))?$", "", tolower(SpectralFile))
    ) %>%
    distinct(SpectralID, .keep_all = TRUE)  # <-- Eindeutigkeit sicherstellen

  # Bereinige RasterFile für den Join (Suffix entfernen und klein schreiben)
  spectral_data <- spectral_data %>%
    mutate(
      RasterID = gsub("-spectral_image$", "", tolower(RasterFile))
    )

  # Debugging: Check the cleaned values
  cat("Cleaned RasterID values in spectral_data:\n")
  print(unique(spectral_data$RasterID))
  cat("Cleaned SpectralID values in metadata_mapping:\n")
  print(unique(metadata_mapping$SpectralID))



  # Add the columns from soil_data to spectral_data
  # Ensure both Name columns are character type before joining
  spectral_data$Name <- as.character(spectral_data$Name)
  soil_data$Name <- as.character(soil_data$Name)
  
  # Debug: Check data types right before join
  cat("Final check - Data type of Name in spectral_data:", class(spectral_data$Name), "\n")
  cat("Final check - Data type of Name in soil_data:", class(soil_data$Name), "\n")
  cat("Number of rows in spectral_data before join:", nrow(spectral_data), "\n")
  cat("Number of rows in soil_data before join:", nrow(soil_data), "\n")
  
  combined_data <- spectral_data %>%
    left_join(soil_data, by = "Name", relationship = "many-to-many")

  # Debugging: Check the structure of combined_data
  cat("Structure of combined_data after merging soil_data:\n")
  print(colnames(combined_data))

  cat("Number of entries after removing duplicates:\n")
  print(nrow(combined_data))

  # Ensure no list columns remain in combined_data
  if (any(sapply(combined_data, is.list))) {
    stop("Error: combined_data contains unsupported list columns.")
  }

  # Return the prepared data
  cat("Data preparation process completed successfully.\n")

  # Save as CSV for easier access - filename indicates normalization method
  output_csv <- paste0("data/combined_data_", NORMALIZATION_METHOD, ".csv")
  write.csv(combined_data, output_csv, row.names = FALSE)
  cat("Final combined data saved to:", output_csv, "\n")

  return(combined_data)
}

# Run the data preparation function with better error handling
tryCatch({
  # Delete any existing empty spectral data files before starting
  spectral_check_path <- paste0("data/spectral_data_", NORMALIZATION_METHOD, ".csv")
  if (file.exists(spectral_check_path)) {
    file_info <- file.info(spectral_check_path)
    if (file_info$size == 0) {
      cat("Removing existing empty spectral data file...\n")
      file.remove(spectral_check_path)
    }
  }
  
  prepared_data <- data_preparation(
    tar_dir = "data/extracted_files",  # Directory containing .tar.gz files
    extract_dir = "data/extracted_files",  # Directory with extracted and smoothed files
    shapefile_path = "E:/Uni/Bachelorarbeit/Daten/Amazonien/Peru_Bolivien_reprojiziert.shp",  # Directory containing the shapefile
    soil_data_dir = "data/soil_measurements"  # Directory containing soil measurement files
  )

  # Use fwrite for saving combined data with updated filename
  if (!is.null(prepared_data) && nrow(prepared_data) > 0) {
    output_csv <- paste0("data/combined_data_", NORMALIZATION_METHOD, ".csv")
    fwrite(prepared_data, output_csv)
    cat("✓ Data successfully prepared and saved to", output_csv, "\n")
    cat("  Number of rows in final dataset:", nrow(prepared_data), "\n")
  } else {
    cat("❌ ERROR: Prepared data is empty or NULL. Combined data file was not updated.\n")
  }
}, error = function(e) {
  cat("❌ ERROR in data preparation:", e$message, "\n")
})


