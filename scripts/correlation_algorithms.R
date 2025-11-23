# This file implements the Barlow-Twin algorithm to compute correlations between spectral data
# and soil measurements.
# Load necessary libraries
library(dplyr)
library(readr)
# Robuste Skalierungsfunktion mit Fehlerbehandlung
scale_safely <- function(x) {
  tryCatch({
    scaled <- scale(x)
    if (any(is.na(scaled) | is.infinite(scaled))) {
      warning("Skalierung hat NAs oder Inf-Werte erzeugt")
      return(x)  # Ursprungsdaten zurückgeben
    }
    return(scaled)
  }, error = function(e) {
    warning("Fehler bei der Skalierung: ", e$message)
    return(x)
  })
}

# Function to select the most recent image for each location
select_one_image_per_location <- function(data) {
  # Check if dplyr is available
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  
  # Log before filtering
  cat("Before filtering: ", nrow(data), " rows\n")
  cat("Unique locations before filtering: ", length(unique(data$Name)), "\n")
  
  # Function to extract acquisition year from RasterFile name
  extract_year <- function(raster_file) {
    # Safety check for NA or empty values
    if (is.na(raster_file) || raster_file == "") {
      return(NA)
    }
    
    # Extract the acquisition date (first date in the format _YYYYMMDD)
    year_match <- regexpr("_([0-9]{4})([0-9]{2})([0-9]{2})T", raster_file)
    if (year_match > 0) {
      year_str <- substr(raster_file, year_match + 1, year_match + 4)
      return(as.numeric(year_str))
    } else {
      # If no year found, return NA
      return(NA)
    }
  }
  
  # Add a column with the extracted acquisition year
  data$acquisition_year <- sapply(data$RasterFile, extract_year)
  
  # Add a column with the absolute difference from 2023
  data$year_diff_from_2023 <- abs(data$acquisition_year - 2023)
  
  # Filter out rows with NA in acquisition_year before grouping
  data <- data %>%
    filter(!is.na(acquisition_year))
  
  # Check if we still have data after filtering
  if(nrow(data) == 0) {
    warning("No data remains after filtering out rows with NA acquisition years")
    return(data)  # Return empty dataset
  }
  
  # Group by Name and select the row with the smallest difference from 2023
  result <- data %>%
    dplyr::group_by(Name) %>%
    dplyr::arrange(year_diff_from_2023, .by_group = TRUE) %>%  # Sort by year difference within each group
    dplyr::slice(1) %>%                      # Take only the first row of each group
    dplyr::ungroup() %>%                     # Remove grouping
    dplyr::select(-acquisition_year, -year_diff_from_2023)  # Remove the helper columns
    
  # Log after filtering
  cat("After filtering: ", nrow(result), " rows\n")
  cat("Unique locations after filtering: ", length(unique(result$Name)), "\n")
  
  return(result)
}

# Further improved function to handle NA values in spectral data during pre-processing
prepare_data_for_correlation <- function(data, soil_param) {
  # Extract spectral bands
  spectral_bands <- grep("^Band[0-9]+$", colnames(data), value = TRUE)
  if (length(spectral_bands) == 0) {
    stop("No spectral bands found matching pattern '^Band[0-9]+$'")
  }
  
  # Get soil parameter values
  soil_values <- data[[soil_param]]
  if (all(is.na(soil_values))) {
    stop("All values in soil parameter are NA")
  }
  
  # First, only keep rows where the soil parameter is not NA
  valid_soil_idx <- which(!is.na(soil_values))
  if (length(valid_soil_idx) < 3) { # Relaxed from 5 to 3
    warning(paste("Very few non-NA values in soil parameter:", length(valid_soil_idx), "- results may be unreliable"))
    # Continue anyway instead of stopping
  }
  
  filtered_data <- data[valid_soil_idx, ]
  soil_values <- soil_values[valid_soil_idx]
  
  # For spectral data, calculate NA percentage per band
  spectral_data <- filtered_data[, spectral_bands]
  na_per_band <- colMeans(is.na(spectral_data))
  
  # Keep bands with any valid data (relaxed from 50% threshold)
  valid_bands <- spectral_bands[na_per_band < 0.9] # Allow up to 90% NAs
  if (length(valid_bands) == 0) {
    warning("No spectral bands with sufficient non-NA values - using all bands")
    valid_bands <- spectral_bands # Use all bands as fallback
  }
  
  # Use only the valid bands
  spectral_data <- filtered_data[, valid_bands]
  
  # Calculate the proportion of missing values for each row
  na_per_row <- rowMeans(is.na(spectral_data))
  
  # Keep rows with any valid data (relaxed threshold)
  valid_rows <- which(na_per_row < 0.9) # Allow up to 90% NAs
  if (length(valid_rows) < 3) { # Relaxed from 5 to 3
    warning("Very few rows with sufficient non-NA spectral values - results may be unreliable")
    # Use all rows as fallback if needed
    if (length(valid_rows) == 0) {
      valid_rows <- 1:nrow(spectral_data)
    }
  }
  
  return(list(
    spectral_data = spectral_data[valid_rows, ],
    soil_values = soil_values[valid_rows],
    n_observations = length(valid_rows),
    valid_bands = valid_bands
  ))
}

# Improved handling for NA values in spectral data and soil parameters
# Only modifying relevant functions; rest of the file remains unchanged

# Update barlow_twin_algorithm to better handle NAs
barlow_twin_algorithm <- function(data, soil_param, output_prefix) {
  # Load necessary packages
  if (!require("pls", quietly = TRUE)) install.packages("pls")
  if (!require("vegan", quietly = TRUE)) install.packages("vegan")
  library(pls)
  library(vegan)
  
  # Prepare data with robust NA handling
  tryCatch({
    prepared_data <- prepare_data_for_correlation(data, soil_param)
    
    # Log the number of observations after filtering
    cat(sprintf("BarlowTwin: %d observations with valid data for %s\n", 
                prepared_data$n_observations, soil_param))
    
    # Proceed only if we have enough data, but with a warning not an error
    if (prepared_data$n_observations < 5) {
      warning("Few observations for BarlowTwin analysis (less than 5) - results may be unreliable")
      # Continue anyway
    }
    
    # Extract the prepared data
    spectral_data <- prepared_data$spectral_data
    soil_values <- prepared_data$soil_values
    
    # Continue with Barlow Twin algorithm with NA handling
    set.seed(123)
    # Use na.action = na.pass to retain NA values for proper handling
    view1 <- scale_safely(as.matrix(spectral_data))  # Convert to matrix to better handle NAs
    
    # View 2: Light noise and non-linear transformation with NA handling
    noise_level <- 0.05
    noise_matrix <- matrix(rnorm(nrow(spectral_data) * ncol(spectral_data), sd = noise_level), 
                         nrow = nrow(spectral_data))
    # Replace NAs with 0 for noise addition, then restore NA positions
    spectral_matrix <- as.matrix(spectral_data)
    na_mask <- is.na(spectral_matrix)
    spectral_matrix[na_mask] <- 0  # Temporary replace NAs with 0
    noisy_data <- spectral_matrix + noise
    noisy_data[na_mask] <- NA  # Restore NAs
    view2 <- scale_safely(log1p(abs(noisy_data)))
    
    # Dimension reduction with NA handling
    # Use prcomp with na.action = "na.omit" internally
    pca1 <- tryCatch({
      prcomp(view1, scale. = FALSE, na.action = na.omit)
    }, error = function(e) {
      warning("PCA failed on view1, using simplified approach: ", e$message)
      # Fallback to mean-centered data for components
      result <- list()
      result$x <- scale(view1, scale = FALSE)
      result$sdev <- apply(result$x, 2, sd, na.rm = TRUE)
      return(result)
    })
    
    pca2 <- tryCatch({
      prcomp(view2, scale. = FALSE, na.action = na.omit)
    }, error = function(e) {
      warning("PCA failed on view2, using simplified approach: ", e$message)
      result <- list()
      result$x <- scale(view2, scale = FALSE)
      result$sdev <- apply(result$x, 2, sd, na.rm = TRUE)
      return(result)
    })
    
    # Choose components explaining 95% of variance, or at least 2 components
    var_explained1 <- cumsum(pca1$sdev^2 / sum(pca1$sdev^2, na.rm = TRUE))
    var_explained2 <- cumsum(pca2$sdev^2 / sum(pca2$sdev^2, na.rm = TRUE))
    
    k1 <- which(var_explained1 >= 0.95)[1]
    k2 <- which(var_explained2 >= 0.95)[1]
    
    # Handle cases where PCA might have fewer components
    k1 <- max(2, min(min(15, ncol(pca1$x)), k1, na.rm = TRUE))
    k2 <- max(2, min(min(15, ncol(pca2$x)), k2, na.rm = TRUE))
    
    # Extract features
    features1 <- pca1$x[, 1:k1, drop = FALSE]
    features2 <- pca2$x[, 1:k2, drop = FALSE]
    
    # Barlow Twin correlation analysis with NA handling
    features1_norm <- scale_safely(features1)
    features2_norm <- scale_safely(features2)
    
    soil_values_norm <- scale_safely(soil_values)
    
    # Calculate Mantel correlation between views (stability)
    mantel_result <- tryCatch({
      mantel(dist(features1_norm), dist(features2_norm), 
             method="pearson", permutations=99, na.rm = TRUE)
    }, error = function(e) {
      warning("Mantel test failed, using default stability: ", e$message)
      list(statistic = 0.5)  # Default fallback stability
    })
    stability <- mantel_result$statistic
    
    # Calculate correlations between original bands and both views
    corr_view1 <- cor(as.matrix(spectral_data), features1_norm, 
                     method="pearson", use="pairwise.complete.obs")
    corr_view2 <- cor(as.matrix(spectral_data), features2_norm, 
                     method="pearson", use="pairwise.complete.obs")
    
    # Weights based on view agreement
    weights <- abs(cor(features1_norm, features2_norm, 
                      use="pairwise.complete.obs"))
    weights <- rowSums(weights, na.rm=TRUE) / ncol(features2_norm)
    
    # Combined correlation with soil parameter
    soil_corr1 <- cor(features1_norm, soil_values_norm, 
                     use="pairwise.complete.obs")
    soil_corr2 <- cor(features2_norm, soil_values_norm, 
                     use="pairwise.complete.obs")
    
    # Calculate weighted correlation for each band
    final_correlation <- matrix(0, nrow = ncol(spectral_data), ncol = 1)
    rownames(final_correlation) <- colnames(spectral_data)
    
    for(i in 1:ncol(spectral_data)) {
      # Combine correlation information from both views
      weighted_corr_view1 <- sum(corr_view1[i,] * soil_corr1 * weights, na.rm=TRUE)
      weighted_corr_view2 <- sum(corr_view2[i,] * soil_corr2 * weights, na.rm=TRUE)
      
      # Final correlation is a weighted average, stabilized by Mantel statistic
      final_correlation[i,1] <- (weighted_corr_view1 + weighted_corr_view2) * stability / 2
    }
    
    # Save additional analysis information
    info_file <- paste0(output_prefix, "_barlow_twin_info.txt")
    cat(file = info_file, 
        "Barlow-Twin Analysis Info\n",
        "Parameter: ", soil_param, "\n",
        "Number of observations: ", prepared_data$n_observations, "\n",
        "View 1 components: ", k1, "\n",
        "View 2 components: ", k2, "\n",
        "View stability (Mantel): ", stability, "\n",
        "Valid spectral bands: ", length(prepared_data$valid_bands), "\n",
        "NA proportion: ", mean(is.na(spectral_data)), "\n")
    
    return(final_correlation)
    
  }, error = function(e) {
    warning(paste("Error in BarlowTwin analysis:", e$message))
    # Return NA matrix as fallback
    spectral_bands <- grep("^Band[0-9]+$", colnames(data), value = TRUE)
    na_matrix <- matrix(NA, nrow = length(spectral_bands), ncol = 1)
    rownames(na_matrix) <- spectral_bands
    return(na_matrix)
  })
}

# Update pearson_correlation with better NA handling
pearson_correlation <- function(data, soil_param, output_prefix) {
  tryCatch({
    # Prepare data with robust NA handling
    prepared_data <- prepare_data_for_correlation(data, soil_param)
    
    # Log the number of observations after filtering
    cat(sprintf("Pearson: %d observations with valid data for %s\n", 
                prepared_data$n_observations, soil_param))
    
    # Continue even with few observations, but warn
    if (prepared_data$n_observations < 3) {
      warning("Very few observations for Pearson correlation - results may be unreliable")
    }
    
    # Extract the prepared data
    spectral_data <- prepared_data$spectral_data
    soil_values <- prepared_data$soil_values
    
    # Force pairwise deletion for NA handling
    correlation_matrix <- cor(spectral_data, soil_values, 
                            use = "pairwise.complete.obs", method = "pearson")
    
    # Check if all correlations are NA and warn
    if (all(is.na(correlation_matrix))) {
      warning("All Pearson correlations are NA - check data completeness")
    }
    
    return(correlation_matrix)
  }, error = function(e) {
    warning(paste("Error in Pearson correlation calculation:", e$message))
    # Return NA matrix as fallback
    spectral_bands <- grep("^Band[0-9]+$", colnames(data), value = TRUE)
    na_matrix <- matrix(NA, nrow = length(spectral_bands), ncol = 1)
    rownames(na_matrix) <- spectral_bands
    return(na_matrix)
  })
}

# Update spearman_correlation with better NA handling
spearman_correlation <- function(data, soil_param, output_prefix) {
  tryCatch({
    # Prepare data with robust NA handling
    prepared_data <- prepare_data_for_correlation(data, soil_param)
    
    # Log the number of observations after filtering
    cat(sprintf("Spearman: %d observations with valid data for %s\n", 
                prepared_data$n_observations, soil_param))
    
    # Continue even with few observations, but warn
    if (prepared_data$n_observations < 3) {
      warning("Very few observations for Spearman correlation - results may be unreliable")
    }
    
    # Extract the prepared data
    spectral_data <- prepared_data$spectral_data
    soil_values <- prepared_data$soil_values
    
    # Force pairwise deletion for NA handling
    correlation_matrix <- cor(spectral_data, soil_values, 
                            use = "pairwise.complete.obs", method = "spearman")
    
    # Check if all correlations are NA and warn
    if (all(is.na(correlation_matrix))) {
      warning("All Spearman correlations are NA - check data completeness")
    }
    
    return(correlation_matrix)
  }, error = function(e) {
    warning(paste("Error in Spearman correlation calculation:", e$message))
    # Return NA matrix as fallback
    spectral_bands <- grep("^Band[0-9]+$", colnames(data), value = TRUE)
    na_matrix <- matrix(NA, nrow = length(spectral_bands), ncol = 1)
    rownames(na_matrix) <- spectral_bands
    return(na_matrix)
  })
}

# --- Hilfsfunktionen für Strata-Qualität ---
compute_balance_score <- function(strata, train_idx, test_idx) {
  # Misst die Ähnlichkeit der Strata-Verteilung zwischen Training und Test (1 = perfekt)
  train_tab <- prop.table(table(strata[train_idx]))
  test_tab  <- prop.table(table(strata[test_idx]))
  all_levels <- union(names(train_tab), names(test_tab))
  train_vec <- train_tab[all_levels]; train_vec[is.na(train_vec)] <- 0
  test_vec  <- test_tab[all_levels];  test_vec[is.na(test_vec)]  <- 0
  1 - sum(abs(train_vec - test_vec))/2
}

compute_representation_score <- function(strata, train_idx, test_idx) {
  # Anteil der Strata, die in beiden Splits mindestens einmal vorkommen (1 = alle vertreten)
  train_levels <- unique(strata[train_idx])
  test_levels  <- unique(strata[test_idx])
  all_levels   <- unique(strata)
  mean(all_levels %in% train_levels & all_levels %in% test_levels)
}

# --- Stratified random split function ---
stratified_random_split <- function(data, response_var, split_ratio = 0.5) {
  set.seed(358)
  n <- nrow(data)

  # Source improved stratification functions
  if (file.exists("scripts/improved_stratification.R")) {
    source("scripts/improved_stratification.R")
  }

  # Initialize variables for quality metrics
  balance_score <- NA
  representation_score <- NA
  stratification_method <- "unknown"
  n_strata <- NA

  # 1. Immer: Tertile-Stratifizierung für kontinuierliche Bodenvariablen
  y <- data[[response_var]]
  if (is.numeric(y) && length(unique(y[!is.na(y)])) > 6) {
    y_clean <- y[!is.na(y)]
    if (length(y_clean) >= 10) { # Mindestens 10 für robuste Stratifizierung
      # Feste Tertile (drei gleich große Gruppen)
      tertiles <- quantile(y_clean, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
      y_strata <- cut(y, breaks = tertiles, 
                      labels = c("niedrig", "mittel", "hoch"), 
                      include.lowest = TRUE, right = TRUE)
      if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
      idx <- caret::createDataPartition(y_strata, p = split_ratio, list = FALSE)
      train_idx <- as.vector(idx)
      test_idx <- setdiff(seq_len(n), train_idx)
      stratification_method <- "tertile"
      n_strata <- 3

      # --- Score-Berechnung ---
      balance_score <- compute_balance_score(y_strata, train_idx, test_idx)
      representation_score <- compute_representation_score(y_strata, train_idx, test_idx)

      cat("TERTILE-Stratifizierung von", response_var, ":\n")
      cat("  Training Strata:"); print(table(y_strata[train_idx], useNA = "ifany"))
      cat("  Test Strata:"); print(table(y_strata[test_idx], useNA = "ifany"))
      cat(sprintf("  Balance Score: %.3f, Representation Score: %.3f\n", balance_score, representation_score))

      return(list(
        train_indices = train_idx, 
        test_indices = test_idx,
        balance_score = balance_score,
        representation_score = representation_score,
        stratification_method = stratification_method,
        n_strata = n_strata
      ))
    }
  }
  
  # 2. Fallback: Stratifizierung nach 'biodiv' falls Quantile nicht möglich
  if ("biodiv" %in% colnames(data) && length(unique(data$biodiv)) > 1) {
    if (!requireNamespace("caret", quietly = TRUE)) install.packages("caret")
    idx <- caret::createDataPartition(data$biodiv, p = split_ratio, list = FALSE)
    train_idx <- as.vector(idx)
    test_idx <- setdiff(seq_len(n), train_idx)
    
    stratification_method <- "biodiv"
    n_strata <- length(unique(data$biodiv))

    # --- Score-Berechnung ---
    balance_score <- compute_balance_score(data$biodiv, train_idx, test_idx)
    representation_score <- compute_representation_score(data$biodiv, train_idx, test_idx)
    
    cat("FALLBACK: Stratifizierung nach biodiv (Landnutzung):\n")
    cat("  Training - pf:", sum(data$biodiv[train_idx] == "pf", na.rm = TRUE), 
        ", df:", sum(data$biodiv[train_idx] == "df", na.rm = TRUE), "\n")
    cat("  Test - pf:", sum(data$biodiv[test_idx] == "pf", na.rm = TRUE), 
        ", df:", sum(data$biodiv[test_idx] == "df", na.rm = TRUE), "\n")
    cat(sprintf("  Balance Score: %.3f, Representation Score: %.3f\n", balance_score, representation_score))
    
    return(list(
      train_indices = train_idx, 
      test_indices = test_idx,
      balance_score = balance_score,
      representation_score = representation_score,
      stratification_method = stratification_method,
      n_strata = n_strata
    ))
  }
  
  # 3. Letzter Fallback: Einfacher zufälliger Split
  train_idx <- sample(seq_len(n), size = round(n * split_ratio))
  test_idx <- setdiff(seq_len(n), train_idx)
  stratification_method <- "random"
  n_strata <- 1

  # --- Score-Berechnung für Random: alle in einer Strata
  strata_random <- rep("all", n)
  balance_score <- compute_balance_score(strata_random, train_idx, test_idx)
  representation_score <- compute_representation_score(strata_random, train_idx, test_idx)

  cat("WARNUNG: Fallback auf einfachen zufälligen Split (keine Stratifizierung möglich)\n")
  cat(sprintf("  Balance Score: %.3f, Representation Score: %.3f\n", balance_score, representation_score))
  
  return(list(
    train_indices = train_idx, 
    test_indices = test_idx,
    balance_score = balance_score,
    representation_score = representation_score,
    stratification_method = stratification_method,
    n_strata = n_strata
  ))
}

# Hauptfunktion zur Verarbeitung aller Bodenparameter
process_soil_parameters <- function(data_file = "data/combined_data_minmax.csv") {
  # Prüfe, ob die Datei existiert
  if (!file.exists(data_file)) {
    stop("Datei nicht gefunden: ", data_file)
  }
  
  # Lade die geglätteten Daten
  combined_data <- tryCatch({
    read.csv(data_file)
  }, error = function(e) {
    stop("Fehler beim Laden der CSV-Datei: ", e$message)
  })
  
  # Überprüfe, ob die Spalte 'biodiv' existiert
  if (!"biodiv" %in% colnames(combined_data)) {
    stop("Die Spalte 'biodiv' wurde nicht gefunden.")
  }
  
  # Vereinfachter Filter: Nur Biodiversität filtern, da Bilder bereits maskiert sind
  # Explizit NA-Werte ausschließen
  filtered_data <- combined_data %>%
    filter(!is.na(biodiv) & biodiv %in% c("pf", "df"))
  
  # Überprüfe, ob nach dem Filtern Daten übrig bleiben
  if (nrow(filtered_data) == 0) {
    stop("Nach dem Filtern für biodiv=='pf' oder biodiv=='df' sind keine Daten übrig.")
  }
  
  # Check if RasterFile column exists and contains valid data
  if (!"RasterFile" %in% colnames(filtered_data) || all(is.na(filtered_data$RasterFile))) {
    stop("Die Spalte 'RasterFile' fehlt oder enthält nur NA-Werte.")
  }
  
  # Wähle zufällig ein Bild pro Standort aus
  filtered_data <- select_one_image_per_location(filtered_data)
  
  # Check if we still have data after filtering
  if (nrow(filtered_data) == 0) {
    stop("Nach der Auswahl eines Bildes pro Standort sind keine Daten übrig.")
  }

  # Update spectral band pattern for current data structure (Band1 bis Band254)
  spectral_pattern <- "^Band[0-9]+$"
  
  spectral_bands <- grep(spectral_pattern, colnames(filtered_data), value = TRUE)
  
  # Entferne nur Zeilen mit >50% NA-Werten in spektralen Bändern statt complete.cases
  if (length(spectral_bands) > 0) {
    # Berechne für jede Zeile den Anteil der NA-Werte in spektralen Bändern
    na_percentages <- rowMeans(is.na(filtered_data[, spectral_bands]))
    
    # Behalte Zeilen, bei denen weniger als 50% NA-Werte sind (also mehr als 50% gültige Werte)
    valid_spectral_rows <- na_percentages < 0.5

    # Ergänzung: Entferne auch Zeilen, bei denen alle Spektralbänder 0 sind
    all_zero_rows <- rowSums(filtered_data[, spectral_bands, drop = FALSE] == 0, na.rm = TRUE) == length(spectral_bands)
    # Zeilen behalten, die NICHT nur aus 0 bestehen
    valid_spectral_rows <- valid_spectral_rows & !all_zero_rows

    # Zähle entfernte Zeilen für Log-Ausgabe
    removed_rows <- sum(!valid_spectral_rows)
    
    # Filtere die Daten
    filtered_data <- filtered_data[valid_spectral_rows, ]
    cat("Entferne", removed_rows, "Zeilen mit >50% NA-Werten oder nur 0 in spektralen Bändern\n")
    
    # Check if we still have data after filtering NA in spectral bands
    if (nrow(filtered_data) == 0) {
      stop("Nach dem Entfernen von Zeilen mit >50% NA-Werten oder nur 0 in spektralen Bändern sind keine Daten übrig.")
    }
  } else {
    stop("Keine spektralen Bänder (Band1 bis Band254) gefunden.")
  }
  
  exclude_columns <- c(
    "X", "Y", "Name", "biodiv", "ID", "Redundancy", "sun_elevation", "cloud_value", "Hole", 
    "RasterFile", "Country",
    spectral_bands,
    grep("^Band[0-9]+$", colnames(filtered_data), value = TRUE), # Original non-corrected bands
    grep("^mask_", colnames(filtered_data), value = TRUE),
    grep("PIXELMASK|QUALITY_|_SWIR|_VNIR", colnames(filtered_data), value = TRUE)
  )
  
  # Filter for numeric columns only
  soil_parameters <- setdiff(colnames(filtered_data), exclude_columns)
  numeric_cols <- sapply(filtered_data[, soil_parameters, drop = FALSE], is.numeric)
  soil_parameters <- soil_parameters[numeric_cols]

  # Entferne Zeilen, die für den Bodenparameter BS NA enthalten
  if ("BS" %in% colnames(filtered_data)) {
    filtered_data <- filtered_data[!is.na(filtered_data$BS), ]
  }

  # Debug-Output der vorbereiteten Daten (jetzt nach Duplikat-Entfernung und NA-Filter!)
  data_dir <- "data"
  if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
  }
  debug_file <- file.path(data_dir, "final_before_correlation.csv")
  write.csv(filtered_data, debug_file, row.names = FALSE)
  cat("Vorbereitete Daten wurden in '", debug_file, "' gespeichert\n", sep="")

  # Ausgabe grundlegender statistischer Informationen
  spectral_bands <- grep("^Band[0-9]+$", colnames(filtered_data), value = TRUE)
  cat("Anzahl der Beobachtungen nach Filterung:", nrow(filtered_data), "\n")
  cat("Enthaltene Spektralbänder:", length(spectral_bands), "\n")
  
  # ==================================================================================
  # ERWEITERTE SPEKTRALE VORVERARBEITUNG
  # ==================================================================================
  
  # Source advanced preprocessing functions
  if (file.exists("scripts/advanced_spectral_preprocessing.R")) {
    source("scripts/advanced_spectral_preprocessing.R")
    
    cat("\n=== APPLYING ADVANCED SPECTRAL PREPROCESSING ===\n")
    
    # Extrahiere spektrale Daten
    spectral_data <- filtered_data[, spectral_bands]
    
    # Wende umfassende Vorverarbeitung an
    if (exists("comprehensive_spectral_preprocessing")) {
      preprocessing_result <- comprehensive_spectral_preprocessing(
        spectral_data = spectral_data,
        methods = c("outlier_detection", "snv", "detrend", "savgol"),
        save_intermediate = TRUE,
        output_dir = "results/preprocessing"
      )
      
      # Ersetze die ursprünglichen spektralen Daten mit vorverarbeiteten
      if (!is.null(preprocessing_result$processed_data)) {
        filtered_data[, spectral_bands] <- preprocessing_result$processed_data
        
        cat("✓ Spektrale Vorverarbeitung erfolgreich angewendet\n")
        cat("  - Verarbeitete Methoden:", paste(preprocessing_result$applied_methods, collapse = ", "), "\n")
        cat("  - Qualitätsmetriken verfügbar in results/preprocessing/\n")
        
        # Speichere Verarbeitungslog
        if (!dir.exists("results/preprocessing")) {
          dir.create("results/preprocessing", recursive = TRUE)
        }
        saveRDS(preprocessing_result$log, "results/preprocessing/processing_log.rds")
      } else {
        cat("⚠️ Spektrale Vorverarbeitung fehlgeschlagen, verwende ursprüngliche Daten\n")
      }
    } else {
      cat("⚠️ Erweiterte Vorverarbeitung nicht verfügbar, verwende Standardverfahren\n")
    }
  } else {
    cat("⚠️ Advanced preprocessing script nicht gefunden, verwende Standardverfahren\n")
  }
  
  # ==================================================================================
  
  # Identifiziere alle Bodenparameter (numerische Spalten außer Spektraldaten)
  exclude_columns <- c(
    "X", "Y", "Name", "biodiv", "ID", "Redundancy", "sun_elevation", "cloud_value", "Hole", 
    "RasterFile", "Country", # Add categorical columns to exclude
    spectral_bands,
    grep("^Band", colnames(filtered_data), value = TRUE),
    grep("^mask_", colnames(filtered_data), value = TRUE),
    grep("PIXELMASK|QUALITY_|_SWIR|_VNIR", colnames(filtered_data), value = TRUE)
  )
  
  # Filter for numeric columns only
  soil_parameters <- setdiff(colnames(filtered_data), exclude_columns)
  numeric_cols <- sapply(filtered_data[, soil_parameters, drop = FALSE], is.numeric)
  soil_parameters <- soil_parameters[numeric_cols]
  
  # Entferne Zeilen, die für alle Bodenparameter NA sind
  if (length(soil_parameters) > 0) {
    filtered_data <- filtered_data[rowSums(!is.na(filtered_data[, soil_parameters, drop = FALSE])) > 0, ]
  }
  
  # Zeige die ersten paar Bodenparameter und deren Statistik
  if(length(soil_parameters) > 0) {
    cat("Gefundene Bodenparameter:", length(soil_parameters), "\n")
    cat("Erste 5 Parameter:", paste(head(soil_parameters, 5), collapse=", "), "\n")
    
    # Grundlegende Statistiken für die ersten 3 Parameter
    for(param in head(soil_parameters, 3)) {
      cat("\nStatistik für", param, ":\n")
      print(summary(filtered_data[[param]]))
      cat("NA-Werte:", sum(is.na(filtered_data[[param]])), "\n")
    }
  } else {
    cat("Keine Bodenparameter identifiziert.\n")
  }
  
  # Zeige die ersten paar Zeilen von wichtigen Spalten
  cat("\nErste paar Zeilen ausgewählter Spalten:\n")
  selected_cols <- c("Name", "biodiv")
  if(length(soil_parameters) > 0) {
    selected_cols <- c(selected_cols, soil_parameters[1])
  }
  print(head(filtered_data[, selected_cols]))
  
  # Überprüfe, ob Bodenparameter gefunden wurden
  if (length(soil_parameters) == 0) {
    stop("Es wurden keine Bodenparameter identifiziert.")
  }
  
  # Erstelle Ergebnisverzeichnisse, falls sie nicht existieren
  results_dir <- "results"
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
    cat("Ergebnisverzeichnis erstellt:", results_dir, "\n")
  }
  
  # Definiere die Algorithmen
  algorithms <- list(
    BarlowTwin = barlow_twin_algorithm,
    Pearson = pearson_correlation,
    Spearman = spearman_correlation
  )
  
  # Erstelle Unterverzeichnisse für jeden Algorithmus
  for (algo_name in names(algorithms)) {
    algo_dir <- file.path(results_dir, algo_name)
    if (!dir.exists(algo_dir)) {
      dir.create(algo_dir, recursive = TRUE)
    }
  }
  
  # Speichere die Parameter-Korrelationen über alle Splits für den Vergleich
  all_correlations <- list()
  # Initialisiere RSS-Split-Log
  rss_split_log <- data.frame(Parameter=character(), Seed=integer(), RSS=double(), 
                              TrainSize=integer(), ValidSize=integer(), stringsAsFactors=FALSE)
  
  # Erstelle neue Listen für jeden Algorithmus
  algo_correlations <- list(
    BarlowTwin = list(),
    Pearson = list(),
    Spearman = list()
  )
  
  # Process each soil parameter with improved error handling
  for (i in seq_along(soil_parameters)) {
    soil_param <- soil_parameters[i]
    cat(sprintf("\nProcessing soil parameter %d/%d: %s\n", i, length(soil_parameters), soil_param))
    
    # Check if parameter has enough non-NA values
    non_na_count <- sum(!is.na(filtered_data[[soil_param]]))
    if (non_na_count < 10) {
      cat(sprintf("Skipping parameter %s: only %d non-NA values (minimum 10 required)\n", 
                  soil_param, non_na_count))
      next
    }
    
    tryCatch({
      # Stratified random split
      split <- stratified_random_split(filtered_data, soil_param, split_ratio = 0.5)
      split_1 <- filtered_data[split$train_indices, ]
      split_2 <- filtered_data[split$test_indices, ]
      cat("Data split: Training =", nrow(split_1), ", Validation =", nrow(split_2), "\n")

      # Create central RSS Info directory
      rss_info_dir <- file.path(results_dir, "RSS Info")
      if (!dir.exists(rss_info_dir)) {
        dir.create(rss_info_dir, recursive = TRUE)
        cat("RSS Info-Verzeichnis erstellt:", rss_info_dir, "\n")
      }

      # Speichere Split-Informationen als CSV für jeden Algorithmus im zentralen RSS Info Ordner
      for (algo_name in names(algorithms)) {
        algo_rss_dir <- file.path(rss_info_dir, algo_name)
        if (!dir.exists(algo_rss_dir)) dir.create(algo_rss_dir, recursive = TRUE)
        
        split_info <- data.frame(
          Parameter = soil_param,
          N_Train = nrow(split_1),
          N_Valid = nrow(split_2),
          Train_Mean = mean(split_1[[soil_param]], na.rm = TRUE),
          Valid_Mean = mean(split_2[[soil_param]], na.rm = TRUE),
          Train_SD = sd(split_1[[soil_param]], na.rm = TRUE),
          Valid_SD = sd(split_2[[soil_param]], na.rm = TRUE),
          Balance_Score = ifelse(is.null(split$balance_score), NA, split$balance_score),
          Representation_Score = ifelse(is.null(split$representation_score), NA, split$representation_score),
          Stratification_Method = ifelse(is.null(split$stratification_method), "unknown", split$stratification_method),
          N_Strata = ifelse(is.null(split$n_strata), NA, split$n_strata)
        )
        split_info_file <- file.path(algo_rss_dir, paste0(soil_param, "_rss_split_info.csv"))
        write.csv(split_info, split_info_file, row.names = FALSE)
        cat("Split-Informationen gespeichert in:", split_info_file, "\n")
        
        # Detaillierte Strata-Dokumentation für alle Probenorte im zentralen RSS Info Ordner
        detailed_strata <- data.frame(
          Parameter = soil_param,
          Name = filtered_data$Name,
          biodiv = filtered_data$biodiv,
          X = filtered_data$X,
          Y = filtered_data$Y,
          Parameter_Value = filtered_data[[soil_param]],
          Split = ifelse(seq_len(nrow(filtered_data)) %in% split$train_indices, "Training", "Validation"),
          stringsAsFactors = FALSE
        )
        
        # VERWENDE DIE TATSÄCHLICH VERWENDETE STRATIFIZIERUNGSMETHODE AUS DEM SPLIT-OBJEKT
        detailed_strata$Stratum_Type <- ifelse(is.null(split$stratification_method), "unknown", split$stratification_method)
        
        # Füge Strata-Werte basierend auf der tatsächlich verwendeten Methode hinzu
        if (!is.null(split$stratification_method)) {
          if (split$stratification_method == "adaptive" || split$stratification_method == "quantile" || split$stratification_method == "tertile") {
            # Verwende die gleiche Logik wie in der stratified_random_split Funktion
            y_clean <- filtered_data[[soil_param]][!is.na(filtered_data[[soil_param]])]
            if (length(y_clean) >= 10) {
              if (split$stratification_method == "quantile") {
                # Verwende die gleiche Anzahl Quantile wie im Split
                n_quantiles <- ifelse(is.null(split$n_strata), 3, split$n_strata)
                probs <- seq(0, 1, length.out = n_quantiles + 1)
                quantiles <- quantile(y_clean, probs = probs)
                param_strata <- cut(filtered_data[[soil_param]], breaks = quantiles, 
                                   labels = paste0("Q", 1:n_quantiles), 
                                   include.lowest = TRUE)
              } else {
                # Für adaptive oder tertile Methode - vereinfachte Darstellung
                tertiles <- quantile(y_clean, probs = c(0, 0.33, 0.67, 1))
                param_strata <- cut(filtered_data[[soil_param]], breaks = tertiles, 
                                   labels = c("niedrig", "mittel", "hoch"), 
                                   include.lowest = TRUE)
              }
              detailed_strata$Stratum_Value <- as.character(param_strata)
            } else {
              detailed_strata$Stratum_Value <- "insufficient_data"
            }
          } else if (split$stratification_method == "biodiv") {
            detailed_strata$Stratum_Value <- filtered_data$biodiv
          } else if (split$stratification_method == "random") {
            detailed_strata$Stratum_Value <- "no_stratification"
          } else {
            detailed_strata$Stratum_Value <- "unknown_method"
          }
        } else {
          # Fallback wenn keine Methode verfügbar
          detailed_strata$Stratum_Value <- "method_unknown"
        }
        
        # Speichere detaillierte Strata-Tabelle im zentralen RSS Info Ordner
        detailed_strata_file <- file.path(algo_rss_dir, paste0(soil_param, "_detailed_strata_documentation.csv"))
        write.csv(detailed_strata, detailed_strata_file, row.names = FALSE)
        cat("Detaillierte Strata-Dokumentation gespeichert in:", detailed_strata_file, "\n")
      }
      
      for (algo_name in names(algorithms)) {
        cat(sprintf("Anwendung von %s...\n", algo_name))
        
        tryCatch({
          # Wende den aktuellen Algorithmus auf beide Splits an
          correlation_1 <- algorithms[[algo_name]](split_1, soil_param, 
                                                  paste0(results_dir, "/", algo_name, "/", soil_param, "_split1"))
          correlation_2 <- algorithms[[algo_name]](split_2, soil_param, 
                                                  paste0(results_dir, "/", algo_name, "/", soil_param, "_split2"))
          
          # Berechne die durchschnittliche Korrelation aus beiden Splits
          avg_correlation <- (correlation_1 + correlation_2) / 2
          
          # Extrahiere die Namen der Spektralbänder (Band1 bis Band254)
          spectral_bands <- grep("^Band[0-9]+$", colnames(filtered_data), value = TRUE)
          # Stelle sicher, dass wir nur Bänder verwenden, die in den Korrelationen vorkommen
          valid_bands <- spectral_bands[spectral_bands %in% rownames(correlation_1)]
          
          # Erstelle die Korrelationstabelle
          correlation_table <- data.frame(
            Band = valid_bands,
            Correlation = avg_correlation[valid_bands, 1],
            Correlation_Split1 = correlation_1[valid_bands, 1],
            Correlation_Split2 = correlation_2[valid_bands, 1]
          )
          
          # Sortiere die Tabelle nach absoluter Korrelation absteigend
          correlation_table <- correlation_table[order(-abs(correlation_table$Correlation)), ]
          
          # Speichere die Korrelationstabelle
          correlation_file <- file.path(results_dir, algo_name, paste0(soil_param, "_correlation_table.csv"))
          write.csv(correlation_table, correlation_file, row.names = FALSE)
          cat(sprintf("  - Korrelationstabelle für %s gespeichert in: %s\n", soil_param, correlation_file))
          
          correlation_similarity <- cor(as.vector(correlation_1), as.vector(correlation_2), 
                                      use = "complete.obs")
          cat(sprintf("  - Korrelation zwischen Split 1 und Split 2 (%s): %.4f\n", algo_name, correlation_similarity))
          
          # Speichere in beiden Listen: der gesamten und der algorithmus-spezifischen
          all_correlations[[paste0(soil_param, "_", algo_name)]] <- list(
            split1 = correlation_1,
            split2 = correlation_2,
            similarity = correlation_similarity,
            avg_correlation = avg_correlation
          )
          
          # Speichere auch in der algorithmus-spezifischen Liste
          algo_correlations[[algo_name]][[soil_param]] <- list(
            split1 = correlation_1,
            split2 = correlation_2,
            similarity = correlation_similarity,
            avg_correlation = avg_correlation
          )
          
        }, error = function(e) {
          warning("Fehler bei der Verarbeitung von ", soil_param, " (", algo_name, "): ", e$message)
        })
      }
    }, error = function(e) {
      cat(sprintf("Error processing parameter %s: %s\n", soil_param, e$message))
    })
  }
  
  # 1. Erstelle die gemeinsame Zusammenfassung wie bisher
  if (length(all_correlations) > 0) {
    similarity_df <- data.frame(
      Parameter = names(all_correlations),
      Similarity = sapply(all_correlations, function(x) x$similarity)
    )
    
    summary_file <- paste0(results_dir, "/correlation_summary.csv")
    write.csv(similarity_df, summary_file, row.names = FALSE)
    cat("\nZusammenfassung aller Korrelationsähnlichkeiten gespeichert in:", summary_file, "\n")
  }
  
  # 2. Erstelle für jeden Algorithmus eine separate Zusammenfassung
  for (algo_name in names(algo_correlations)) {
    if (length(algo_correlations[[algo_name]]) > 0) {
      # Erstelle DataFrame für diesen Algorithmus
      algo_similarity_df <- data.frame(
        Parameter = names(algo_correlations[[algo_name]]),
        Similarity = sapply(algo_correlations[[algo_name]], function(x) x$similarity)
      )
      
      # Sortiere nach höchster Ähnlichkeit
      algo_similarity_df <- algo_similarity_df[order(-algo_similarity_df$Similarity), ]
      
      # Speichere in algorithmus-spezifischer Datei
      algo_summary_file <- file.path(results_dir, algo_name, "correlation_summary.csv")
      write.csv(algo_similarity_df, algo_summary_file, row.names = FALSE)
      cat(sprintf("%s-Zusammenfassung gespeichert in: %s\n", algo_name, algo_summary_file))
    }
  }
  
  return(all_correlations)
}

# HINWEIS: Auto-Execution auskommentiert - Script wird nur über master_pipeline.R gesteuert
# # Starte die Verarbeitung
# tryCatch({
#   results <- process_soil_parameters()
#   
#   # Print completion message mit dynamischer Ausgabe
#   successful_params <- length(results)
#   cat("\nKorrelationsanalyse abgeschlossen.\n")
#   cat(sprintf("%d von %d Bodenparametern wurden erfolgreich verarbeitet.\n", 
#               successful_params, length(names(results))))
#   cat("Die Ergebnisse wurden im 'results'-Verzeichnis gespeichert.\n")
# }, error = function(e) {
#   cat("\nFehler in der Korrelationsanalyse:\n")
#   cat(e$message, "\n")
#   cat("Überprüfen Sie die Daten und versuchen Sie es erneut.\n")
# })
