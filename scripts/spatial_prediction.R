library(raster)
library(terra)
library(dplyr)
library(readr)
library(caret)
library(pls)  # Für PLSR-Modelle
library(pracma)  # Für Savitzky-Golay Filter
library(parallel)
library(pbapply)
library(tmap)
library(RColorBrewer)
# Add explicit viridis dependency
if (!requireNamespace("viridis", quietly = TRUE)) {
  cat("Installing viridis package for turbo color palette...\n")
  install.packages("viridis", repos = "https://cloud.r-project.org")
}
library(viridis)  # Make sure viridis is explicitly loaded for turbo palette


# 🚀 PERFORMANCE OPTIMIZATIONS 🚀

# Global flag for quick testing (set to FALSE for production)
QUICK_TEST <- FALSE  # Keep at FALSE for full resolution processing
REDUCTION_FACTOR <- 1  # No reduction - keep original resolution
BATCH_SIZE <- 5000  # Process pixels in smaller batches

# Remove the quick_test_mode function - we want full resolution
# quick_test_mode function removed

# Enhanced batch processing with column alignment
process_in_batches <- function(pred_df, plsr_model, ncomp, batch_size = BATCH_SIZE) {
  n_rows <- nrow(pred_df)
  n_batches <- ceiling(n_rows / batch_size)
  
  cat("      📦 Verarbeite", n_rows, "Pixel in", n_batches, "Batches à", batch_size, "Pixel\n")
  
  all_predictions <- numeric(n_rows)
  
  # Progress bar for batches
  pb <- txtProgressBar(min = 0, max = n_batches, style = 3, char = "=")
  
  for (i in 1:n_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_rows)
    
    # Process batch - exclude coordinate columns
    batch_data <- pred_df[start_idx:end_idx, ]
    
    # Remove x,y coordinates if they exist
    coord_cols <- c("x", "y")
    spectral_cols <- setdiff(names(batch_data), coord_cols)
    batch_spectral <- batch_data[, spectral_cols, drop = FALSE]
    
    batch_pred <- predict(plsr_model, newdata = batch_spectral, ncomp = ncomp)
    all_predictions[start_idx:end_idx] <- as.numeric(batch_pred)
    
    # Update progress and cleanup memory every 10 batches
    if (i %% 10 == 0 || i == n_batches) {
      setTxtProgressBar(pb, i)
      gc(verbose = FALSE) # Silent memory cleanup
    }
  }
  
  close(pb)
  cat("\n")
  return(all_predictions)
}

# Test area cropping function
crop_to_test_area <- function(raster_obj, test_fraction = 0.25) {
  # Crop to center area for testing
  ext_orig <- ext(raster_obj)
  width <- ext_orig[2] - ext_orig[1]
  height <- ext_orig[4] - ext_orig[3]
  
  # Calculate center crop
  center_x <- (ext_orig[1] + ext_orig[2]) / 2
  center_y <- (ext_orig[3] + ext_orig[4]) / 2
  
  crop_width <- width * sqrt(test_fraction)
  crop_height <- height * sqrt(test_fraction)
  
  test_ext <- ext(center_x - crop_width/2, center_x + crop_width/2,
                  center_y - crop_height/2, center_y + crop_height/2)
  
  cropped <- crop(raster_obj, test_ext)
  
  cat("  ✂️ Test-Bereich:", ncell(cropped), "Pixel (", 
      round(test_fraction * 100, 1), "% der ursprünglichen Fläche)\n")
  
  return(cropped)
}

# Verbessertes Temporärdateien-Management
setup_temp_directory <- function() {
  # Erstelle ein temporäres Verzeichnis im Working Directory
  temp_dir <- file.path(getwd(), "temp_files")
  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
  }
  
  # Setze den temporären Verzeichnispfad für R
  old_temp_dir <- Sys.getenv("TMPDIR")
  Sys.setenv(TMPDIR = temp_dir)
  
  # Setze den temporären Verzeichnispfad für terra
  terraOptions(tempdir = temp_dir)
  
  cat("Temporäres Verzeichnis eingerichtet:", temp_dir, "\n")
  
  # Gib das alte temp dir zurück, um es später wiederherzustellen
  return(old_temp_dir)
}

# Funktion zum Aufräumen des temporären Verzeichnisses
cleanup_temp_files <- function(temp_dir, old_temp_dir) {
  # Lösche alle Dateien im temporären Verzeichnis
  temp_files <- list.files(temp_dir, full.names = TRUE)
  
  if (length(temp_files) > 0) {
    cat("Lösche", length(temp_files), "temporäre Dateien...\n")
    unlink(temp_files, recursive = TRUE, force = TRUE)
  }
  
  # Rufe den Garbage Collector auf, um Speicher freizugeben
  gc()
  
  # Setze den ursprünglichen temporären Verzeichnispfad wieder
  if (!is.null(old_temp_dir)) {
    Sys.setenv(TMPDIR = old_temp_dir)
  }
  
  cat("Temporäre Dateien aufgeräumt.\n")
}

# Funktion zur Überwachung des Speicherplatzes
check_disk_space <- function(min_space_gb = 5, path = "/") {
  if (.Platform$OS.type == "windows") {
    drive <- substr(path, 1, 2)
    df_info <- system(paste0("wmic logicaldisk where DeviceID='", drive, "' get FreeSpace"), intern = TRUE)
    free_space <- as.numeric(gsub("[^0-9]", "", df_info[2])) / (1024^3)  # Bytes zu GB
  } else {
    df_info <- system(paste("df -BG", path), intern = TRUE)
    free_space <- as.numeric(gsub("G", "", strsplit(df_info[2], "\\s+")[[1]][4]))
  }
  
  cat("Freier Speicherplatz:", round(free_space, 2), "GB\n")
  
  if (free_space < min_space_gb) {
    warning(paste0("Wenig Speicherplatz übrig (", round(free_space, 2), 
                   " GB). Mindestens ", min_space_gb, " GB empfohlen."))
    return(FALSE)
  }
  
  return(TRUE)
}

# Wrapper function to ensure devices get closed properly
safe_plot_wrapper <- function(plot_fn, ...) {
  # Make sure we start with a clean state
  while (dev.cur() > 1) { dev.off() }
  
  # Try the plot function
  result <- tryCatch({
    plot_fn(...)
    TRUE
  }, error = function(e) {
    # If an error occurs, make sure to close any open devices
    while (dev.cur() > 1) { dev.off() }
    cat("❌ Plot error:", e$message, "\n")
    FALSE
  })
  
  # Ensure devices are closed before returning
  while (dev.cur() > 1) { dev.off() }
  
  return(result)
}

# Function to write prediction statistics to a text file
write_stats_file <- function(param_dir, param, algorithm, plsr_type, stats_info) {
  stats_file <- file.path(param_dir, paste0(param, "_stats_", algorithm, "_", plsr_type, ".txt"))
  
  # Format the statistics into readable text
  stats_lines <- c(
    paste0("Parameter: ", param),
    paste0("Algorithm: ", algorithm),
    paste0("PLSR Type: ", plsr_type),
    "===========================",
    "Training Statistics:",
    paste0("R² (Training): ", if(!is.null(stats_info$R_Squared_Train) && !is.na(stats_info$R_Squared_Train)) round(stats_info$R_Squared_Train, 4) else "Not Available"),
    paste0("RMSE (Training): ", if(!is.null(stats_info$RMSE_Train) && !is.na(stats_info$RMSE_Train)) round(stats_info$RMSE_Train, 4) else "Not Available"),
    paste0("RPD (Training): ", if(!is.null(stats_info$RPD_Train) && !is.na(stats_info$RPD_Train)) round(stats_info$RPD_Train, 4) else "Not Available"),
    "===========================",
    "Validation Statistics:",
    paste0("R² (Validation): ", if(!is.null(stats_info$R_Squared_Valid) && !is.na(stats_info$R_Squared_Valid)) round(stats_info$R_Squared_Valid, 4) else "Not Available"),
    paste0("RMSE (Validation): ", if(!is.null(stats_info$RMSE_Valid) && !is.na(stats_info$RMSE_Valid)) round(stats_info$RMSE_Valid, 4) else "Not Available"),
    paste0("RPD (Validation): ", if(!is.null(stats_info$RPD_Valid) && !is.na(stats_info$RPD_Valid)) round(stats_info$RPD_Valid, 4) else "Not Available"),
    "===========================",
    paste0("Number of Components: ", if(!is.null(stats_info$N_Components) && !is.na(stats_info$N_Components)) stats_info$N_Components else "Not Available"),
    paste0("Generated on: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  )
  
  # Write statistics to file
  writeLines(stats_lines, stats_file)
  cat("  ✓ Statistics saved to", stats_file, "\n")
  
  return(invisible(NULL))
}

# Helper function to get consistent color schemes and ranges for parameters
get_param_visualization_settings <- function(param) {
  # Default settings
  settings <- list(
    scale_min = NULL,
    scale_max = NULL,
    color_palette = "turbo",    # Changed to turbo for all parameters
    reverse_colors = FALSE,     # No need to reverse with turbo
    n_breaks = 5,               # Number of breaks in legend
    legend_format = "general"   # Can be "general", "percent", "pH"
  )
  
  # Parameter-specific settings
  if (grepl("BS", param, ignore.case=TRUE)) {
    settings$scale_min <- 0
    settings$scale_max <- 100
    settings$color_palette <- "turbo"
    settings$reverse_colors <- FALSE
    settings$legend_format <- "percent"
    settings$n_colors <- 1000
  } 
  else if(grepl("pH", param, ignore.case=TRUE)) {
    settings$scale_min <- 3
    settings$scale_max <- 10
    settings$color_palette <- "turbo"
    settings$reverse_colors <- FALSE
    settings$legend_format <- "pH"
    settings$n_colors <- 1000
  } 
  else if(grepl("SWC|Corg_percentage|Norg_percentage", param, ignore.case=TRUE)) {
    # Nur diese Parameter als Prozent behandeln!
    settings$scale_min <- 0
    settings$scale_max <- 100
    settings$color_palette <- "turbo"
    settings$reverse_colors <- FALSE
    settings$legend_format <- "percent"
  }
  else if(grepl("X13C_perc", param, ignore.case=TRUE)) {
    # Optional: expliziter Bereich für X13C_perc, falls gewünscht
    settings$scale_min <- 1.07
    settings$scale_max <- 1.10
    settings$color_palette <- "turbo"
    settings$legend_format <- "general"
  }
  else if(grepl("X15N_perc", param, ignore.case=TRUE)) {
    # Optional: expliziter Bereich für X15N_perc, falls gewünscht
    settings$scale_min <- 0.95
    settings$scale_max <- 1.10
    settings$color_palette <- "turbo"
    settings$legend_format <- "general"
  }
  else if(grepl("CEC|_M3|_stock_|CN_ration", param, ignore.case=TRUE)) {
    settings$color_palette <- "turbo"
    settings$reverse_colors <- FALSE
    settings$legend_format <- "general"
    # For these parameters we'll use data-driven min/max
  }
  else {
    # For any other parameter use turbo
    settings$color_palette <- "turbo"
    settings$reverse_colors <- FALSE
  }
  
  return(settings)
}

# Funktion zum Erstellen einer PNG-Map mit Titel, Nordpfeil, Gitter, Legende und Skala
create_prediction_map <- function(raster_file, png_file, param, algorithm, plsr_type, scale_min = NULL, scale_max = NULL) {
  # Verify the output directory exists and is writable
  output_dir <- dirname(png_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("  Created output directory:", output_dir, "\n")
  }
  
  # Test file writing ability with a simple text file
  test_file <- file.path(output_dir, "write_test.txt")
  test_result <- tryCatch({
    writeLines("test", test_file)
    file.remove(test_file)
    TRUE
  }, error = function(e) {
    cat("  ⚠️ WARNING: Cannot write to output directory:", e$message, "\n")
    FALSE
  })
  
  if (!test_result) {
    cat("  ❌ ERROR: Directory", output_dir, "is not writable. PNG cannot be created.\n")
    return(FALSE)
  }
  
  cat("  Debug: Will create PNG at:", png_file, "\n")
  tryCatch({
    # Lade erforderliche Bibliotheken
    library(terra)
    library(tmap)
    library(RColorBrewer)
    suppressWarnings(library(cols4all))  # Für neue Farbpaletten - unterdrücke Warnungen falls nicht installiert
    library(viridis)  # Force load the viridis package
    
    # Prüfe tmap Version
    tmap_version <- packageVersion("tmap")
    is_tmap_v4 <- as.integer(strsplit(as.character(tmap_version), "\\.")[[1]][1]) >= 4
    
    cat("  Debug: Lade Raster", raster_file, "...\n")
    
    # Extrahiere Bildname für die Kartenunterschrift
    image_basename <- basename(dirname(raster_file))
    
    # Lade Raster und stelle sicher, dass es ein SpatRaster mit nur einem Layer ist
    r <- rast(raster_file)
    if (nlyr(r) > 1) {
      cat("  Debug: Raster hat mehr als einen Layer, reduziere auf ersten Layer\n")
      r <- r[[1]]
    }
    
    if (is.null(r) || !inherits(r, "SpatRaster")) {
      stop("Ungültiges Rasterobjekt")
    }
    
    # Debug-Info zu Raster-Werten
    r_summary <- summary(values(r))
    cat("  Debug: Raster-Werte - Min:", r_summary[1], "Max:", r_summary[6], "\n")
    
    # Get consistent visualization settings for this parameter
    viz_settings <- get_param_visualization_settings(param)
    
    # Override scale_min/max only if not provided explicitly
    if (is.null(scale_min)) scale_min <- viz_settings$scale_min
    if (is.null(scale_max)) scale_max <- viz_settings$scale_max
    
    # If still NULL, set based on data
    if (is.null(scale_min)) scale_min <- floor(min(values(r), na.rm=TRUE))
    if (is.null(scale_max)) scale_max <- ceiling(max(values(r), na.rm=TRUE))
    
    # Definiere einen einheitlichen Namen für die Legende
    map_title <- paste(param, "-", algorithm, plsr_type)
    if (grepl("BS", param, ignore.case=TRUE)) {
      map_title <- paste("Base Saturation -", algorithm, plsr_type)
    }
    
    # Berechne die Werte für die Legende basierend auf Min/Max und Parameter-Typ
    if (viz_settings$legend_format == "percent") {
      # Für Prozentparameter: Feste Beschriftung 0%, 25%, 50%, 75%, 100%
      legend_breaks <- c(0, 25, 50, 75, 100)
      legend_labels <- paste0(legend_breaks, "%")
      legend_title <- "Prozent (%)"
    } else if (viz_settings$legend_format == "pH") {
      # Für pH-Werte: Angepasste Beschriftung
      legend_breaks <- seq(scale_min, scale_max, length.out = 5)
      legend_labels <- as.character(round(legend_breaks, 1))
      legend_title <- "pH-Wert"
    } else {
      # Für allgemeine Parameter: Datenbasierte Beschriftung
      legend_breaks <- seq(scale_min, scale_max, length.out = 5)
      legend_labels <- as.character(round(legend_breaks, 1))
      legend_title <- param
    }
    
    cat("  Debug: Using turbo color palette from viridis package\n")
    
    # Generate and check turbo color palette
    turbo_pal <- viridis::turbo(1000)
    cat("  Debug: Turbo palette has", length(turbo_pal), "colors\n")
    
    # Erstelle Karte je nach tmap-Version
    if (is_tmap_v4) {
      cat("  Debug: Using tmap v4 for visualization\n")
      
      # tmap v4 Syntax with improved legend and north arrow
      tm <- tm_shape(r) +
            tm_raster(
              col.legend = tm_legend(
                title = legend_title,
                labels = legend_labels,
                width = 10,
                height = 18
              ),
              col.scale = tm_scale_continuous(  
                values = "turbo",
                n = 1000,
                breaks = legend_breaks,
                legend.is.portrait = TRUE,
                legend.reverse = FALSE
              ),
              style = NULL,
              interpolate = TRUE
            ) +
            tm_title(map_title, size = 1.2, fontface = "bold") +
            tm_layout(
              frame = TRUE,
              inner.margins = c(0.15, 0.05, 0.15, 0.15),
              outer.margins = c(0.02, 0.02, 0.02, 0.02),
              legend.outside = TRUE,
              legend.outside.position = "right",
              legend.frame = TRUE,
              legend.hist = FALSE,
              legend.width = 10,
              legend.text.size = 0.8,
              legend.title.size = 1.0
            ) +
            tm_scalebar(
              position = c("left", "bottom"),
              text.size = 0.8,
              breaks = c(0, 5, 10)
            ) +
            tm_compass(
              position = c("right", "top"), 
              type = "arrow",
              size = 3,
              text.size = 0.8
            ) +
            tm_graticules(
              lines = TRUE,
              labels.size = 0.6,
              alpha = 0.3
            )
      
      # Speichere die Karte
      cat("  Debug: Saving tmap v4 to:", png_file, "\n")
      tmap_save(tm, filename = png_file, width = 1400, height = 1000, dpi = 150)
      
      # Verify the file was created
      if (file.exists(png_file)) {
        cat("  ✓ PNG file successfully created at:", png_file, "\n")
        cat("  ✓ File size:", file.size(png_file) / 1024, "KB\n")
        return(TRUE)
      } else {
        cat("  ⚠️ tmap_save completed but file wasn't created!\n")
      }
    } else {
      cat("  Debug: Using tmap v3 for visualization\n")
      
      # tmap v3 Syntax with improved legend and north arrow
      tm <- tm_shape(r) +
            tm_raster(
              title = legend_title,
              style = "fixed",
              breaks = c(-Inf, legend_breaks, Inf),
              labels = c("", legend_labels, ""),
              palette = turbo_pal,
              n = 1000,
              legend.show = TRUE,
              legend.reverse = FALSE,
              legend.is.portrait = TRUE,
              legend.hist = FALSE,
              interpolate = TRUE,
              colorNA = "grey90",
              textNA = "Keine Daten"
            ) +
            tm_layout(
              title = map_title,
              title.size = 1.2,
              title.fontface = "bold",
              frame = TRUE,
              inner.margins = c(0.15, 0.05, 0.15, 0.15),
              outer.margins = c(0.02, 0.02, 0.02, 0.02),
              legend.outside = TRUE,
              legend.outside.position = "right",
              legend.frame = TRUE,
              legend.hist = FALSE,
              legend.width = 0.15,
              legend.text.size = 0.8,
              legend.title.size = 1.0
            ) +
            tm_scalebar(
              position = c("left", "bottom"),
              text.size = 0.8,
              breaks = c(0, 5, 10)
            ) +
            tm_compass(
              position = c("right", "top"), 
              type = "arrow",
              size = 3,
              text.size = 0.8
            ) +
            tm_grid(
              lines = TRUE,
              labels.size = 0.6,
              alpha = 0.3
            )
      
      # Speichere die Karte
      cat("  Debug: Saving tmap v3 to:", png_file, "\n")
      tmap_save(tm, filename = png_file, width = 1400, height = 1000, dpi = 150)
      
      # Verify the file was created
      if (file.exists(png_file)) {
        cat("  ✓ PNG file successfully created at:", png_file, "\n")
        cat("  ✓ File size:", file.size(png_file) / 1024, "KB\n")
        return(TRUE)
      } else {
        cat("  ⚠️ tmap_save completed but file wasn't created!\n")
      }
    }
    cat(" ✓\n")
    
  }, error = function(e) {
    cat(" ⚠️ FEHLER bei der tmap-Erstellung:", e$message, "\n")
    cat("  Versuche Fallback mit terra::plot...\n")
    
    # Fallback 1: Verbesserte terra::plot mit korrekter Beschriftung und Nordpfeil
    tryCatch({
      r <- rast(raster_file)
      if (nlyr(r) > 1) r <- r[[1]]
      
      # Get consistent visualization settings for this parameter
      viz_settings <- get_param_visualization_settings(param)
      
      # Override scale_min/max only if not provided explicitly
      if (is.null(scale_min)) scale_min <- viz_settings$scale_min
      if (is.null(scale_max)) scale_max <- viz_settings$scale_max
      
      # If still NULL, set based on data
      if (is.null(scale_min)) scale_min <- floor(min(values(r), na.rm=TRUE))
      if (is.null(scale_max)) scale_max <- ceiling(max(values(r), na.rm=TRUE))
      
      # Definiere Kartentitel und Legendenbeschriftung
      map_title <- paste(param, "-", algorithm, plsr_type)
      if (grepl("BS", param, ignore.case=TRUE)) {
        map_title <- paste("Base Saturation -", algorithm, plsr_type)
      }
      
      # Berechne korrekte Legendenwerte
      if (viz_settings$legend_format == "percent") {
        # Für Prozentparameter: Feste Beschriftung 0%, 25%, 50%, 75%, 100%
        legend_breaks <- c(0, 25, 50, 75, 100)
        legend_labels <- paste0(legend_breaks, "%")
        legend_title <- "Prozent (%)"
        plot_range <- c(0, 100)
      } else if (viz_settings$legend_format == "pH") {
        # Für pH-Werte
        legend_breaks <- seq(scale_min, scale_max, length.out = 5)
        legend_labels <- as.character(round(legend_breaks, 1))
        legend_title <- "pH-Wert"
        plot_range <- c(scale_min, scale_max)
      } else {
        # Für allgemeine Parameter
        legend_breaks <- seq(scale_min, scale_max, length.out = 5)
        legend_labels <- as.character(round(legend_breaks, 1))
        legend_title <- param
        plot_range <- c(scale_min, scale_max)
      }
      
      # Use turbo from viridis package
      cat("  Debug: Creating turbo color palette for terra::plot\n") 
      colors <- viridis::turbo(1000)
      cat("  Debug: Generated", length(colors), "colors for the plot\n")

      # Make sure output directory exists
      output_dir <- dirname(png_file)
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
      
      cat("  Debug: Creating PNG with terra::plot at:", png_file, "\n")
      
      # Verbesserte Kartenerstellung mit terra::plot
      png(png_file, width = 1600, height = 1200, res = 150)
      
      # Layout mit mehr Platz für Legende
      layout(matrix(c(1, 2), nrow = 1), widths = c(0.75, 0.25))
      
      # Hauptkarte
      par(mar = c(4, 4, 4, 1))
      plot(r, main = map_title, legend = FALSE, col = colors, 
           range = plot_range, axes = TRUE, 
           cex.main = 1.4, font.main = 2)
      
      # Koordinatengitter hinzufügen
      grid(col = "gray70", lty = 3, lwd = 0.5)
      
      # Maßstabsbalken unten links (verbessert)
      ext_r <- ext(r)
      x_range <- ext_r[2] - ext_r[1]
      y_range <- ext_r[4] - ext_r[3]
      
      # Berechne sinnvolle Maßstabslänge
      scale_length <- round(x_range * 0.15 / 1000) * 1000  # Gerundet auf km
      if (scale_length == 0) scale_length <- 1000  # Minimum 1 km
      
      x_sb <- ext_r[1] + x_range * 0.05
      y_sb <- ext_r[3] + y_range * 0.08
      
      # Maßstabsbalken zeichnen
      segments(x_sb, y_sb, x_sb + scale_length, y_sb, lwd = 6, col = "black")
      segments(x_sb, y_sb, x_sb + scale_length, y_sb, lwd = 4, col = "white")
      segments(x_sb, y_sb, x_sb + scale_length, y_sb, lwd = 2, col = "black")
      
      # Maßstabsbeschriftung
      text(x_sb + scale_length/2, y_sb - y_range * 0.03, 
           paste0(scale_length/1000, " km"), 
           cex = 1.2, font = 2, col = "black")
      
      # Nordpfeil oben rechts (verbessert)
      x_n <- ext_r[2] - x_range * 0.08
      y_n <- ext_r[4] - y_range * 0.12
      arrow_length <- y_range * 0.06
      
      # Nordpfeil-Schatten
      arrows(x_n + x_range * 0.002, y_n - y_range * 0.002, 
             x_n + x_range * 0.002, y_n + arrow_length - y_range * 0.002, 
             length = 0.15, lwd = 4, col = "gray50")
      
      # Nordpfeil
      arrows(x_n, y_n, x_n, y_n + arrow_length, 
             length = 0.15, lwd = 3, col = "red", angle = 20)
      
      # "N"-Beschriftung
      text(x_n, y_n + arrow_length + y_range * 0.02, "N", 
           cex = 1.8, col = "red", font = 2)
      
      # Separate Legende (rechts)
      par(mar = c(4, 1, 4, 4))
      
      # Erstelle Dummy-Plot für Legende
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", 
           xlim = c(0, 1), ylim = c(0, 1))
      
      # Farbbalken erstellen
      legend_y <- seq(0.15, 0.85, length.out = length(colors))
      legend_width <- 0.15
      
      for (i in 1:(length(colors)-1)) {
        rect(0.3, legend_y[i], 0.3 + legend_width, legend_y[i+1], 
             col = colors[i], border = NA)
      }
      
      # Legendenrahmen
      rect(0.3, 0.15, 0.3 + legend_width, 0.85, border = "black", lwd = 2)
      
      # Legendenbeschriftung
      legend_positions <- seq(0.15, 0.85, length.out = length(legend_breaks))
      
      for (i in 1:length(legend_breaks)) {
        # Markierungslinien
        segments(0.3 + legend_width, legend_positions[i], 
                 0.3 + legend_width + 0.05, legend_positions[i], 
                 col = "black", lwd = 1)
        
        # Beschriftung
        text(0.3 + legend_width + 0.08, legend_positions[i], 
             legend_labels[i], cex = 1.1, adj = 0)
      }
      
      # Legendentitel
      text(0.375, 0.95, legend_title, cex = 1.3, font = 2, adj = 0.5)
      
      dev.off()
      
      # Verify the file was created
      if (file.exists(png_file)) {
        cat("  ✓ PNG file successfully created with terra::plot at:", png_file, "\n")
        cat("  ✓ File size:", file.size(png_file) / 1024, "KB\n")
        return(TRUE)
      } else {
        cat("  ⚠️ terra::plot completed but file wasn't created!\n")
      }
    }, error = function(e2) {
      cat("  ⚠️ Auch terra::plot schlug fehl:", e2$message, "\n")
      cat("  Versuche letzten Fallback...\n")
      
      # Final fallback: Einfache base R Lösung mit korrekten Beschriftungen
      tryCatch({
        r <- rast(raster_file)
        if (nlyr(r) > 1) r <- r[[1]]
        
        # Get parameter settings
        viz_settings <- get_param_visualization_settings(param)
        
        # Create turbo color palette using viridis
        colors <- viridis::turbo(100)
        
        # Make sure output directory exists
        output_dir <- dirname(png_file)
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE)
        }
        
        cat("  Debug: Creating PNG with base R image() at:", png_file, "\n")
        
        # Einfache Karte mit base R
        png(png_file, width = 1200, height = 900, res = 100)
        par(mar = c(4, 4, 4, 6))
        
        # Hauptkarte
        plot(r, main = paste(param, "-", algorithm, plsr_type), 
             col = colors, axes = TRUE, cex.main = 1.2, font.main = 2)
        
        # Einfacher Nordpfeil
        ext_r <- ext(r)
        x_range <- ext_r[2] - ext_r[1]
        y_range <- ext_r[4] - ext_r[3]
        
        # Nordpfeil oben rechts
        x_n <- ext_r[2] - x_range * 0.1
        y_n <- ext_r[4] - y_range * 0.15
        arrows(x_n, y_n, x_n, y_n + y_range * 0.08, 
               length = 0.15, lwd = 3, col = "red")
        text(x_n, y_n + y_range * 0.1, "N", cex = 1.5, col = "red", font = 2)
        
        # Einfacher Maßstabsbalken
        x_sb <- ext_r[1] + x_range * 0.05
        y_sb <- ext_r[3] + y_range * 0.05
        scale_len <- x_range * 0.15
        segments(x_sb, y_sb, x_sb + scale_len, y_sb, lwd = 4, col = "black")
        text(x_sb + scale_len/2, y_sb - y_range * 0.03, 
             paste0(round(scale_len/1000, 1), " km"), cex = 1.0, font = 2)
        
        dev.off()
        
        # Verify the file was created
        if (file.exists(png_file)) {
          cat("  ✓ PNG file successfully created with base R image() at:", png_file, "\n")
          cat("  ✓ File size:", file.size(png_file) / 1024, "KB\n")
          return(TRUE)
        } else {
          cat("  ❌ All visualization methods failed for:", png_file, "\n")
          return(FALSE)
        }
      }, error = function(e3) {
        cat("  ❌ Keine Visualisierung möglich:", e3$message, "\n")
        return(FALSE)
      })
    })
  })
  
  # Final check if file exists
  if (file.exists(png_file)) {
    return(TRUE)
  } else {
    cat("  ❌ Failed to create PNG file after all attempts\n")
    return(FALSE)
  }
}

# Helper function to get parameter-specific value constraints
get_param_constraints <- function(param) {
  constraints <- list(
    min_value = NULL,
    max_value = NULL,
    apply_constraints = FALSE
  )
  
  # Parameter-specific constraints based on physical/chemical limits
  if (grepl("BS", param, ignore.case=TRUE)) {
    constraints$min_value <- 0
    constraints$max_value <- 100
    constraints$apply_constraints <- TRUE
  } 
  else if (grepl("pH", param, ignore.case=TRUE)) {
    constraints$min_value <- 3
    constraints$max_value <- 14
    constraints$apply_constraints <- TRUE
  } 
  else if (grepl("SWC", param, ignore.case=TRUE)) {
    # Soil Water Content should be between 0 and reasonable upper limit
    constraints$min_value <- 0
    constraints$max_value <- 50  # Reasonable upper limit for soil water content
    constraints$apply_constraints <- TRUE
  }
  else if (grepl("Corg_percentage|Norg_percentage", param, ignore.case=TRUE)) {
    constraints$min_value <- 0
    constraints$max_value <- 100
    constraints$apply_constraints <- TRUE
  }
  else if (grepl("CEC", param, ignore.case=TRUE)) {
    constraints$min_value <- 0
    constraints$max_value <- 200  # Reasonable upper limit for CEC
    constraints$apply_constraints <- TRUE
  }
  
  return(constraints)
}

# Function to write raster statistics to CSV
write_raster_stats_csv <- function(tif_file, param, algorithm, plsr_type, image_name) {
  # Create CSV filename in the same directory as the TIF
  param_dir <- dirname(tif_file)
  csv_file <- file.path(param_dir, paste0(param, "_raster_stats_", algorithm, "_", plsr_type, ".csv"))
  
  tryCatch({
    # Load raster and calculate statistics
    r <- rast(tif_file)
    raster_values <- values(r, na.rm = TRUE)
    
    if (length(raster_values) > 0) {
      # Calculate comprehensive statistics
      stats_data <- data.frame(
        Parameter = param,
        Algorithm = algorithm,
        PLSR_Type = plsr_type,
        Image_Name = image_name,
        TIF_File = basename(tif_file),
        Min_Value = min(raster_values, na.rm = TRUE),
        Max_Value = max(raster_values, na.rm = TRUE),
        Mean_Value = mean(raster_values, na.rm = TRUE),
        Median_Value = median(raster_values, na.rm = TRUE),
        Std_Dev = sd(raster_values, na.rm = TRUE),
        Q25 = quantile(raster_values, 0.25, na.rm = TRUE),
        Q75 = quantile(raster_values, 0.75, na.rm = TRUE),
        Valid_Pixels = length(raster_values),
        Total_Pixels = ncell(r),
        NA_Pixels = ncell(r) - length(raster_values),
        Coverage_Percent = round((length(raster_values) / ncell(r)) * 100, 2),
        Generated_On = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        stringsAsFactors = FALSE
      )
      
      # Write CSV file
      write.csv(stats_data, csv_file, row.names = FALSE)
      
      cat(sprintf("      📊 Raster-Statistik gespeichert: %s\n", basename(csv_file)))
      cat(sprintf("         Min: %.3f, Max: %.3f, Gültige Pixel: %d\n", 
                  stats_data$Min_Value, stats_data$Max_Value, stats_data$Valid_Pixels))
      
      return(TRUE)
    } else {
      cat(sprintf("      ⚠️ Keine gültigen Werte in Raster: %s\n", basename(tif_file)))
      return(FALSE)
    }
  }, error = function(e) {
    cat(sprintf("      ❌ Fehler beim Erstellen der Raster-Statistik: %s\n", e$message))
    return(FALSE)
  })
}

# Enhanced chunk processing for direct raster prediction with value constraints
process_raster_in_chunks <- function(smoothed_raster, plsr_model, ncomp, chunk_size = 50000, param_name = "unknown") {
  total_cells <- ncell(smoothed_raster)
  n_chunks <- ceiling(total_cells / chunk_size)
  
  cat("      📦 Verarbeite", total_cells, "Zellen in", n_chunks, "Chunks à max.", chunk_size, "Zellen\n")
  
  # Get parameter constraints
  constraints <- get_param_constraints(param_name)
  
  # Create result raster
  result <- rast(smoothed_raster[[1]])
  values(result) <- NA
  
  # Progress bar for chunks
  pb <- txtProgressBar(min = 0, max = n_chunks, style = 3, char = "=")
  
  # Track prediction statistics
  all_chunk_predictions <- c()
  
  for (i in 1:n_chunks) {
    start_cell <- (i - 1) * chunk_size + 1
    end_cell <- min(i * chunk_size, total_cells)
    
    # Extract values for this chunk
    chunk_cells <- start_cell:end_cell
    chunk_values <- values(smoothed_raster)[chunk_cells, , drop = FALSE]
    
    # Remove rows with too many NAs
    na_count <- rowSums(is.na(chunk_values))
    valid_rows <- na_count < (ncol(chunk_values) * 0.5)
    
    if (sum(valid_rows) > 0) {
      chunk_clean <- chunk_values[valid_rows, , drop = FALSE]
      valid_chunk_cells <- chunk_cells[valid_rows]
      
      # Interpolate remaining NAs
      for (j in 1:nrow(chunk_clean)) {
        row_data <- chunk_clean[j, ]
        na_indices <- which(is.na(row_data))
        
        if (length(na_indices) > 0) {
          for (na_idx in na_indices) {
            # Simple linear interpolation
            non_na_vals <- which(!is.na(row_data))
            if (length(non_na_vals) >= 2) {
              before_idx <- max(non_na_vals[non_na_vals < na_idx], default = 0)
              after_idx <- min(non_na_vals[non_na_vals > na_idx], default = 0)
              
              if (before_idx > 0 && after_idx > 0) {
                interpolated <- row_data[before_idx] + 
                               (row_data[after_idx] - row_data[before_idx]) * 
                               (na_idx - before_idx) / (after_idx - before_idx)
                chunk_clean[j, na_idx] <- interpolated
              } else if (before_idx > 0) {
                chunk_clean[j, na_idx] <- row_data[before_idx]
              } else if (after_idx > 0) {
                chunk_clean[j, na_idx] <- row_data[after_idx]
              } else {
                chunk_clean[j, na_idx] <- median(row_data, na.rm = TRUE)
              }
            }
          }
        }
      }
      
      # Create prediction dataframe with consistent column names
      chunk_df <- as.data.frame(chunk_clean)
      names(chunk_df) <- paste0("Band", 1:ncol(chunk_df))
      
      # Make predictions
      tryCatch({
        chunk_pred <- predict(plsr_model, newdata = chunk_df, ncomp = ncomp)
        chunk_pred <- as.numeric(chunk_pred)
        
        # Store raw predictions for statistics
        all_chunk_predictions <- c(all_chunk_predictions, chunk_pred[!is.na(chunk_pred)])
        
        # Apply parameter-specific constraints first
        if (constraints$apply_constraints) {
          # Constrain to parameter-specific limits
          if (!is.null(constraints$min_value)) {
            chunk_pred[chunk_pred < constraints$min_value] <- constraints$min_value
          }
          if (!is.null(constraints$max_value)) {
            chunk_pred[chunk_pred > constraints$max_value] <- constraints$max_value
          }
          
          cat(sprintf("      🔒 Parameter %s: Werte auf [%.1f, %.1f] begrenzt\n", 
                      param_name, 
                      ifelse(is.null(constraints$min_value), min(chunk_pred, na.rm=TRUE), constraints$min_value),
                      ifelse(is.null(constraints$max_value), max(chunk_pred, na.rm=TRUE), constraints$max_value)))
        } else {
          # For parameters without specific constraints, remove extreme outliers
          pred_median <- median(chunk_pred, na.rm = TRUE)
          pred_mad <- mad(chunk_pred, na.rm = TRUE)
          
          # Use robust outlier detection (median ± 3 * MAD)
          lower_bound <- pred_median - 3 * pred_mad
          upper_bound <- pred_median + 3 * pred_mad
          
          extremes <- chunk_pred < lower_bound | chunk_pred > upper_bound
          chunk_pred[extremes] <- NA
        }
        
        # Clean remaining invalid values
        chunk_pred[is.na(chunk_pred) | is.infinite(chunk_pred)] <- NA
        
        # Assign predictions to result raster
        result[valid_chunk_cells] <- chunk_pred
      }, error = function(e) {
        cat("      ⚠️ Chunk", i, "prediction failed:", e$message, "\n")
      })
    }
    
    # Update progress
    setTxtProgressBar(pb, i)
    
    # Memory cleanup every 10 chunks
    if (i %% 10 == 0) gc(verbose = FALSE)
  }
  
  close(pb)
  
  # Report prediction statistics
  if (length(all_chunk_predictions) > 0) {
    cat(sprintf("      📊 Rohe Vorhersage-Statistik für %s:\n", param_name))
    cat(sprintf("         Min: %.3f, Max: %.3f\n", 
                min(all_chunk_predictions, na.rm=TRUE), 
                max(all_chunk_predictions, na.rm=TRUE)))
    cat(sprintf("         Median: %.3f, MAD: %.3f\n", 
                median(all_chunk_predictions, na.rm=TRUE), 
                mad(all_chunk_predictions, na.rm=TRUE)))
  }
  
  cat("\n")
  return(result)
}

# Diese Funktion nutzt die CARS-PLSR Modelle für die Vorhersage
predict_from_rasters <- function(
  data_file = "data/final_before_correlation.csv",
  results_base_dir = "results",
  smoothed_dir = "data",
  smoothed_pattern = "_smoothed\\+adaptive_destriped_masked_snv_minmax\\.tif$",  # <-- updated pattern
  metadata_file = "data/metadata_mapping.csv",
  output_dir = NULL,
  algorithm = "Pearson",
  plsr_type = "filter+plsr",
  parameters_to_include = NULL,
  start_index = 1,
  skip_existing = TRUE,
  min_space_gb = 5,
  recursive_search = TRUE,
  create_detailed_viz = TRUE
) {
  # Performance mode info
  cat("\n⚡ VEREINFACHTER WORKFLOW ⚡\n")
  cat("- Direkte Raster-Vorhersage in Chunks\n")
  cat("- Anschließend PNG-Erstellung mit einheitlicher Skalierung\n\n")

  # Verbessertes Temporärdateien-Management
  old_temp_dir <- setup_temp_directory()
  
  # Prüfe verfügbaren Speicherplatz
  if (!check_disk_space(min_space_gb, getwd())) {
    cat("Warnung: Zu wenig Speicherplatz, Verarbeitung wird abgebrochen.\n")
    cleanup_temp_files(file.path(getwd(), "temp_files"), old_temp_dir)
    return(invisible(NULL))
  }
  
  # Wenn kein output_dir angegeben ist, erstelle einen standardmäßigen
  if (is.null(output_dir)) {
    output_dir <- file.path(results_base_dir, algorithm, plsr_type, "raster_predictions")
  }
  
  # Verzeichnis erstellen
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Prüfe ob bereits verarbeitete Bilder existieren und lade den letzten Index
  last_processed_file <- file.path(output_dir, "last_processed_image.txt")
  if (file.exists(last_processed_file)) {
    last_processed_index <- tryCatch({
      as.integer(readLines(last_processed_file, n = 1))
    }, error = function(e) {
      cat("⚠️ Konnte last_processed_image.txt nicht lesen, verwende start_index\n")
      start_index
    })
    
    # Verwende den höheren Wert zwischen start_index und last_processed + 1
    actual_start_index <- max(start_index, last_processed_index + 1)
    
    if (actual_start_index > start_index) {
      cat(sprintf("📄 Fortsetzen ab Bild %d (letzte verarbeitete: %d)\n", 
                  actual_start_index, last_processed_index))
    }
  } else {
    actual_start_index <- start_index
    cat(sprintf("📄 Beginne mit Bild %d (keine vorherige Verarbeitung gefunden)\n", 
                actual_start_index))
  }
  
  # Pfad zu PLSR-Modellen für den angegebenen Algorithmus und Typ
  plsr_dir <- file.path(results_base_dir, algorithm, plsr_type)
  
  if (!dir.exists(plsr_dir)) {
    cat("PLSR-Verzeichnis nicht gefunden für Algorithmus und Typ:", algorithm, plsr_type, "\n")
    cleanup_temp_files(file.path(getwd(), "temp_files"), old_temp_dir)
    return(invisible(NULL))
  }
  
  # Lade die PLSR-Zusammenfassung
  plsr_summary_file <- file.path(plsr_dir, "plsr_summary.csv")
  
  if (!file.exists(plsr_summary_file)) {
    cat("PLSR-Zusammenfassung nicht gefunden:", plsr_summary_file, "\n")
    cleanup_temp_files(file.path(getwd(), "temp_files"), old_temp_dir)
    return(invisible(NULL))  # Skip this algorithm
  }
  
  plsr_summary <- read_csv(plsr_summary_file, show_col_types = FALSE)
  
  # Remove R² threshold - use all available models
  cat(sprintf("\n=== Running spatial prediction without R² threshold for %s ===\n", plsr_type))
  
  # Use only models where R²_train > R²_valid > 0.4
  valid_models <- plsr_summary %>%
    dplyr::filter(R_Squared_Train > R_Squared_Valid & R_Squared_Valid > 0.4) %>%
    arrange(desc(R_Squared_Valid))

  if (nrow(valid_models) == 0) {
    cat("Keine PLSR-Modelle gefunden.\n")
    cleanup_temp_files(file.path(getwd(), "temp_files"), old_temp_dir)
    return(invisible(NULL))
  }
  
  # Nutze alle Parameter aus den valid_models (nicht mehr nur BS)
  parameters <- valid_models$Parameter

  if (length(parameters) == 0) {
    cat("Keine Parameter gefunden.\n")
    cleanup_temp_files(file.path(getwd(), "temp_files"), old_temp_dir)
    return(invisible(NULL))
  }
  
  cat(sprintf("\n✓ %d Parameter gefunden für die Vorhersage.\n", length(parameters)))
  
  # Lade die Metadaten
  metadata <- read_csv(metadata_file)

  # --- VERSCHOBEN: Lade die smoothed images mit rekursiver Suche, wenn aktiviert ---
  cat("\n=== Suche nach geglätteten Bildern ===\n")
  cat("Verzeichnis:     ", smoothed_dir, "\n")
  cat("Suchmuster:      ", smoothed_pattern, "\n")
  cat("Rekursive Suche: ", if(recursive_search) "Ja" else "Nein", "\n")
  
  smoothed_files <- list.files(
    smoothed_dir,
    pattern = smoothed_pattern,
    full.names = TRUE,
    recursive = recursive_search
  )
  
  # Fallback: If not found, try searching in 'data' directly (for legacy/flat structure)
  if (length(smoothed_files) == 0) {
    cat("Keine passenden Raster im angegebenen Verzeichnis gefunden. Suche im Hauptdatenverzeichnis...\n")
    smoothed_files <- list.files(
      "data",
      pattern = smoothed_pattern,
      full.names = TRUE,
      recursive = TRUE
    )
  }
  
  if (length(smoothed_files) == 0) {
    cat("\n⚠️ WARNUNG: Keine geglätteten Bilder gefunden!\n")
    cat("Bitte überprüfen Sie die Einstellungen:\n")
    cat("- Verzeichnis: ", smoothed_dir, "\n")
    cat("- Suchmuster:  ", smoothed_pattern, "\n")
    cat("- Rekursive Suche: ", if(recursive_search) "Ja" else "Nein", "\n")
    cat("Beispielhafte Dateinamen der originalen Bilder:\n")
    if(nrow(metadata) > 0) {
      cat(paste("  -", head(metadata$SpectralFile, 3)), sep="\n")
    }
    cleanup_temp_files(file.path(getwd(), "temp_files"), old_temp_dir)
    return(invisible(NULL))
  }
  
  cat(sprintf("\n✓ %d geglättete Bilder gefunden\n", length(smoothed_files)))
  cat("Beispiele: ", paste(basename(head(smoothed_files, 3)), collapse=", "), "\n")
  
  # Strip the correct suffix for matching
  smoothed_basenames <- gsub("_smoothed\\+adaptive_destriped_masked_snv_minmax\\.tif$", "", basename(smoothed_files))
  
  # Entferne ebenfalls das Suffix von den Einträgen in der Metadatendatei
  metadata$SpectralFile_base <- gsub("\\.(BSQ|TIF)$", "", basename(metadata$SpectralFile))
  
  # Verbesserte Matching-Logik: Prüfe auf teilweise Übereinstimmungen
  cat("\n=== Vergleiche Bildnamen mit Metadaten ===\n")
  cat("Prüfe Übereinstimmungen...\n")
  
  valid_indices <- c()
  pb <- txtProgressBar(min = 0, max = length(smoothed_basenames), style = 3)
  for (i in 1:length(smoothed_basenames)) {
    matches <- sapply(metadata$SpectralFile_base, function(x) {
      grepl(x, smoothed_basenames[i], fixed=TRUE) || 
      grepl(smoothed_basenames[i], x, fixed=TRUE)
    })
    if (any(matches)) {
      valid_indices <- c(valid_indices, i)
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  cat("\n")
  
  if (length(valid_indices) == 0) {
    cat("\n⚠️ WARNUNG: Keine Übereinstimmung zwischen geglätteten Bildern und Metadaten gefunden!\n")
    cat("Versuche alternative Matching-Methode...\n")
    clean_smoothed <- gsub("[0-9_]", "", smoothed_basenames)
    clean_metadata <- gsub("[0-9_]", "", metadata$SpectralFile_base)
    pb <- txtProgressBar(min = 0, max = length(clean_smoothed), style = 3)
    for (i in 1:length(clean_smoothed)) {
      for (j in 1:length(clean_metadata)) {
        if (grepl(clean_metadata[j], clean_smoothed[i], fixed=TRUE) || 
            grepl(clean_smoothed[i], clean_metadata[j], fixed=TRUE)) {
          valid_indices <- c(valid_indices, i)
          break
        }
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
    cat("\n")
    if (length(valid_indices) == 0) {
      cat("\n❌ FEHLER: Auch mit alternativer Matching-Methode keine Übereinstimmungen gefunden!\n")
      cleanup_temp_files(file.path(getwd(), "temp_files"), old_temp_dir)
      return(invisible(NULL))
    }
  }
  
  smoothed_files <- smoothed_files[valid_indices]
  smoothed_basenames <- smoothed_basenames[valid_indices]
  
  cat(sprintf("\n✓ %d geglättete Bilder mit Metadaten gefunden\n", length(smoothed_files)))
  # --- ENDE VERSCHOBEN ---

  # 3-PHASEN WORKFLOW: CSV → Min/Max → TIF/PNG
  cat(sprintf("\n=== 3-Phasen Workflow für alle Bilder ab Bild %d/%d ===\n", actual_start_index, length(smoothed_files)))
  
  # PHASE 1: Erstelle TIF-Vorhersagen für alle Bilder und Parameter
  cat(sprintf("\n=== PHASE 1: Erstelle Raster-Vorhersagen für alle Bilder ab Bild %d/%d ===\n", actual_start_index, length(smoothed_files)))
  
  tif_files_created <- 0
  
  # Prüfe ob alle Bilder bereits verarbeitet wurden
  if (actual_start_index > length(smoothed_files)) {
    cat("✅ Alle Bilder bereits verarbeitet. Springe zu Phase 2 (Min/Max Bestimmung).\n")
  } else {
    for (img_index in actual_start_index:length(smoothed_files)) {
      smoothed_file <- smoothed_files[img_index]
      image_name <- smoothed_basenames[img_index]
      
      cat(sprintf("  [%d/%d] Raster-Vorhersagen für Bild: %s\n", img_index, length(smoothed_files), image_name))
      
      # Prüfe ob die Raster-Datei existiert (immer mit dem gefundenen Pfad!)
      if (!file.exists(smoothed_file)) {
        cat(sprintf("    ❌ Raster-Datei nicht gefunden (Pfad aus smoothed_files): %s\n", smoothed_file))
        cat("    Debug: Beispiel gefundene Dateien:\n")
        print(head(smoothed_files, 3))
        next
      } else {
        cat(sprintf("    ✓ Raster-Datei gefunden: %s\n", smoothed_file))
      }
      
      # Erstelle Ausgabeordner für dieses Bild
      img_output_dir <- file.path(output_dir, image_name)
      if (!dir.exists(img_output_dir)) {
        dir.create(img_output_dir, recursive = TRUE)
        cat(sprintf("    Erstellt Bildordner: %s\n", img_output_dir))
      }
      
      # Prüfe ob bereits alle Parameter für dieses Bild verarbeitet wurden
      all_params_exist <- TRUE
      for (param in parameters) {
        param_dir <- file.path(img_output_dir, param)
        tif_file <- file.path(param_dir, paste0(param, "_predicted_", algorithm, "_", plsr_type, ".tif"))
        if (!file.exists(tif_file)) {
          all_params_exist <- FALSE
          break
        }
      }
      
      if (all_params_exist && skip_existing) {
        cat(sprintf("    ✅ Alle Parameter bereits für Bild %s verarbeitet, überspringe\n", image_name))
        # Aktualisiere trotzdem den Fortschritt
        writeLines(as.character(img_index), last_processed_file)
        next
      }
      
      # Lade Raster
      smoothed_raster <- tryCatch({
        start_time <- Sys.time()
        r <- rast(smoothed_file)
        load_time <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
        cat("    📊 Auflösung:", nrow(r), "x", ncol(r), "Pixel,", nlyr(r), "Bänder\n")
        cat("    ⏱️ Ladezeit:", load_time, "Sekunden\n")
        r
      }, error = function(e) {
        cat(sprintf("    ❌ Fehler beim Laden des Rasters: %s\n", e$message))
        return(NULL)
      })
      
      if (is.null(smoothed_raster)) next
      
      # Erstelle Vorhersage für jeden Parameter
      for (param_idx in 1:length(parameters)) {
        param <- parameters[param_idx]
        
        cat(sprintf("    🔄 Vorhersage für Parameter: %s\n", param))
        
        # Erstelle Parameter-Ordner
        param_dir <- file.path(img_output_dir, param)
        if (!dir.exists(param_dir)) {
          dir.create(param_dir, recursive = TRUE)
        }
        
        # TIF-Ausgabedatei
        tif_file <- file.path(param_dir, paste0(param, "_predicted_", algorithm, "_", plsr_type, ".tif"))
        
        # Prüfe ob TIF bereits existiert
        if (file.exists(tif_file) && skip_existing) {
          cat(sprintf("      TIF existiert bereits: %s\n", basename(tif_file)))
          tif_files_created <- tif_files_created + 1
          next
        }
        
        # Lade Modell
        model_file <- file.path(plsr_dir, param, paste0(param, "_plsr_model.rds"))
        if (!file.exists(model_file)) {
          cat(sprintf("      ❌ Modell nicht gefunden: %s\n", model_file))
          next
        }
        
        plsr_model <- tryCatch({
          readRDS(model_file)
        }, error = function(e) {
          cat(sprintf("      ❌ Fehler beim Laden des Modells: %s\n", e$message))
          return(NULL)
        })
        
        if (is.null(plsr_model)) next
        
        # Finde Anzahl der Komponenten
        ncomp_row <- valid_models[valid_models$Parameter == param, ]
        if (nrow(ncomp_row) == 0) {
          cat(sprintf("      ❌ Parameter nicht in valid_models: %s\n", param))
          next
        }
        ncomp <- ncomp_row$N_Components[1]
        
        cat(sprintf("      Verwende %d Komponenten\n", ncomp))
        
        # Erstelle Raster-Vorhersage in Chunks
        prediction_start <- Sys.time()
        result_raster <- process_raster_in_chunks(smoothed_raster, plsr_model, ncomp, chunk_size = 50000, param_name = param)
        prediction_time <- round(as.numeric(difftime(Sys.time(), prediction_start, units = "secs")), 2)
        
        cat(sprintf("      ⏱️ Vorhersagezeit: %s Sekunden\n", prediction_time))
        
        # Setze CRS und speichere TIF
        tryCatch({
          crs(result_raster) <- "EPSG:32719"
          writeRaster(result_raster, tif_file, overwrite = TRUE, filetype = "GTiff",
                     datatype = "FLT4S", gdal = c("COMPRESS=LZW", "TILED=YES"))
          tif_files_created <- tif_files_created + 1
          cat(sprintf("      ✓ TIF gespeichert: %s\n", basename(tif_file)))
          
          # Write raster statistics CSV
          write_raster_stats_csv(tif_file, param, algorithm, plsr_type, image_name)
          
        }, error = function(e) {
          cat(sprintf("      ❌ Fehler beim Speichern: %s\n", e$message))
        })
      }
      
      # Speichere Fortschritt nach jedem erfolgreich verarbeiteten Bild
      writeLines(as.character(img_index), last_processed_file)
      cat(sprintf("    📄 Fortschritt gespeichert: Bild %d/%d\n", img_index, length(smoothed_files)))
    }
  }
  
  cat(sprintf("✓ Phase 1 abgeschlossen. %d TIF-Dateien erstellt.\n", tif_files_created))

  # PHASE 2: Bestimme globale Min/Max für jeden Parameter
  cat("\n=== PHASE 2: Bestimme globale Min/Max für einheitliche Skalierung ===\n")
  
  global_ranges <- list()
  
  for (param in parameters) {
    cat(sprintf("  📊 Analysiere Parameter: %s", param))
    all_values <- c()
    tif_count <- 0
    
    # Sammle Werte aus allen TIF-Dateien für diesen Parameter (von Anfang an, nicht nur ab actual_start_index)
    for (img_index in 1:length(smoothed_files)) {
      image_name <- smoothed_basenames[img_index]
      param_dir <- file.path(output_dir, image_name, param)
      tif_file <- file.path(param_dir, paste0(param, "_predicted_", algorithm, "_", plsr_type, ".tif"))
      
      if (file.exists(tif_file)) {
        tif_count <- tif_count + 1
        tryCatch({
          # Warte kurz und versuche dann das Laden
          Sys.sleep(0.1)
          r <- rast(tif_file)
          
          # Für BS verwende feste Bereiche 0-100%
          if (grepl("BS", param, ignore.case=TRUE)) {
            # Für BS Parameter verwende immer 0-100%
            all_values <- c(all_values, 0, 100)
          } else {
            # Für andere Parameter sample aus den Daten
            sample_values <- sample(values(r), min(10000, ncell(r)), na.rm = TRUE)
            sample_values <- sample_values[!is.na(sample_values)]
            if (length(sample_values) > 0) {
              all_values <- c(all_values, sample_values)
            }
          }
        }, error = function(e) {
          cat(sprintf(" [Fehler bei %s: %s]", basename(tif_file), e$message))
        })
      }
    }
    
    if (length(all_values) > 0) {
      # Für BS verwende immer 0-100%
      if (grepl("BS", param, ignore.case=TRUE)) {
        global_ranges[[param]] <- list(
          min = 0,
          max = 100
        )
        cat(sprintf(" → %d TIFs, Range: 0%% - 100%% (fest für BS)\n", tif_count))
      } else {
        global_ranges[[param]] <- list(
          min = quantile(all_values, 0.02, na.rm = TRUE),  # 2% Quantil
          max = quantile(all_values, 0.98, na.rm = TRUE)   # 98% Quantil
        )
        cat(sprintf(" → %d TIFs, Range: %.3f - %.3f\n", 
                    tif_count, global_ranges[[param]]$min, global_ranges[[param]]$max))
      }
    } else {
      cat(sprintf(" → %d TIFs gefunden, aber keine Daten\n", tif_count))
    }
  }

  # PHASE 3: Erstelle PNG-Karten mit einheitlicher Skalierung
  cat("\n=== PHASE 3: Erstelle PNG-Karten mit einheitlicher Skalierung ===\n")
  
  png_files_created <- 0
  
  for (param in parameters) {
    if (!param %in% names(global_ranges)) {
      cat(sprintf("  ⚠️ Keine globalen Werte für Parameter %s - überspringe PNG-Erstellung\n", param))
      next
    }
    
    param_min <- global_ranges[[param]]$min
    param_max <- global_ranges[[param]]$max
    
    cat(sprintf("  🎨 Erstelle PNGs für Parameter %s (Skala: %.3f - %.3f)\n", param, param_min, param_max))
    
    # Erstelle PNGs für alle Bilder (nicht nur ab actual_start_index)
    for (img_index in 1:length(smoothed_files)) {
      image_name <- smoothed_basenames[img_index]
      param_dir <- file.path(output_dir, image_name, param)
      tif_file <- file.path(param_dir, paste0(param, "_predicted_", algorithm, "_", plsr_type, ".tif"))
      png_file <- file.path(param_dir, paste0(param, "_prediction_map_", algorithm, "_", plsr_type, ".png"))
      
      if (file.exists(tif_file) && (!file.exists(png_file) || !skip_existing)) {
        cat(sprintf("    [%d/%d] %s - %s\n", img_index, length(smoothed_files), image_name, param))
        
        png_success <- create_prediction_map(tif_file, png_file, param, algorithm, plsr_type, 
                                            scale_min = param_min, scale_max = param_max)
        if (png_success) {
          png_files_created <- png_files_created + 1
        }
      }
    }
  }
  
  cat(sprintf("✓ Phase 3 abgeschlossen. %d PNG-Dateien erstellt.\n", png_files_created))

  # Abschließendes Aufräumen
  cat("\n→ Abschließendes Aufräumen...\n")
  cleanup_temp_files(file.path(getwd(), "temp_files"), old_temp_dir)
  cat(sprintf("\n✓ Vereinfachter Workflow abgeschlossen: %d TIFs, %d PNGs erstellt.\n", 
              tif_files_created, png_files_created))
  
  return(invisible(NULL))
}

# Funktion zum Ausführen der Vorhersage für alle Algorithmen
run_all_predictions <- function(
  data_file = "data/final_before_correlation.csv",
  results_base_dir = "results",
  smoothed_dir = "data",
  metadata_file = "data/metadata_mapping.csv",
  normalization_method = "snv",  # "snv" oder "minmax"
  start_index = 1,
  skip_existing = TRUE,
  algorithms = c("Pearson", "Spearman"),
  plsr_types = c("filter+plsr"),
  recursive_search = TRUE,
  create_detailed_viz = TRUE
) {
  # Erstelle Suchmuster basierend auf Normalisierungsmethode
  if (normalization_method == "minmax") {
    smoothed_pattern <- "_smoothed\\+adaptive_destriped_masked_snv_minmax\\.tif$"
    method_name <- "MinMax"
  } else if (normalization_method == "snv") {
    smoothed_pattern <- "_smoothed\\+adaptive_destriped_masked_snv\\.tif$"
    method_name <- "SNV"
  } else {
    stop("❌ Invalid normalization_method: ", normalization_method, ". Must be 'snv' or 'minmax'.\n")
  }
  
  cat("Verwende", method_name, "normalisierte Bilder für räumliche Vorhersage.\n")
  # Nur für Pearson und Spearman, nur filter+plsr
  for (algorithm in algorithms) {
    predict_from_rasters(
      data_file = data_file,
      results_base_dir = results_base_dir,
      smoothed_dir = smoothed_dir,
      smoothed_pattern = smoothed_pattern,  # Uses the masked pattern
      metadata_file = metadata_file,
      algorithm = algorithm,
      plsr_type = "filter+plsr",
      start_index = start_index,
      skip_existing = skip_existing,
      recursive_search = recursive_search,
      create_detailed_viz = create_detailed_viz
    )
  }
  full_memory_cleanup()
}

# Aggressive memory cleanup function
full_memory_cleanup <- function() {
  # Store names of all large objects in environment
  large_objects <- ls(envir = .GlobalEnv, pattern = "enmap|raster|band|pred|result|smooth")
  
  # Remove these objects explicitly
  if (length(large_objects) > 0) {
    cat("Cleaning up", length(large_objects), "large objects from memory...\n")
    rm(list = large_objects, envir = .GlobalEnv)
  }
  
  # Force aggressive garbage collection multiple times
  for (i in 1:3) {
    gc(full = TRUE, verbose = TRUE)
  }
  
  # On Windows, try to force memory release back to the OS
  if (.Platform$OS.type == "windows") {
    invisible(tryCatch({
      memory.size(max = FALSE)
    }, error = function(e) NULL))
  }
  
  cat("Memory cleanup complete.\n")
}

# Hilfsfunktion: Liefert die Spaltennamen für die korrigierten Spektralbänder
get_corrected_band_names <- function(n_bands) {
  paste0("Band", 1:n_bands, "_corr")
}

# Beispielhafte Funktion, wie Prädiktorvariablen aus der Tabelle extrahiert werden sollten:
extract_predictors_from_table <- function(df, n_bands) {
  band_names <- get_corrected_band_names(n_bands)
  # Prüfe, ob alle Spalten existieren
  missing <- setdiff(band_names, names(df))
  if (length(missing) > 0) {
    stop(paste("Fehlende korrigierte Bands in Tabelle:", paste(missing, collapse = ", ")))
  }
  df[, band_names, drop = FALSE]
}

main <- function() {
  # Check if NORMALIZATION_METHOD is defined, otherwise default to "snv"
  if (!exists("NORMALIZATION_METHOD")) {
    NORMALIZATION_METHOD <- "snv"
    cat("⚠️  NORMALIZATION_METHOD not defined, defaulting to 'snv'\n")
  }
  
  cat("Starting spatial prediction processing...\n")
  cat("📊 Using normalization method:", toupper(NORMALIZATION_METHOD), "\n")
  
  run_all_predictions(
    data_file = "data/final_before_correlation.csv",
    results_base_dir = "results",
    smoothed_dir = "data",
    metadata_file = "data/metadata_mapping.csv",
    normalization_method = NORMALIZATION_METHOD,
    start_index = 1,
    skip_existing = TRUE,
    recursive_search = TRUE,
    create_detailed_viz = TRUE
  )
}

# Execute the main function if not sourced from master_pipeline
if (!exists("run_complete_pipeline")) {
  main()
}
