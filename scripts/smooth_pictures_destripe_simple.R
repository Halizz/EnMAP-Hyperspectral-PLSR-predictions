library(terra)
library(raster)
library(readr)
library(dplyr)    # Added for pipe operator and distinct function
library(pracma)  # Für Savitzky-Golay Filter
library(pbapply)
library(tools)   # Für Dateioperationen



# Check if needed functions exist in environment, define if not
if (!exists("configure_terra")) {
  configure_terra <- function(max_mem_percent = 50, cache_size = 100) {
    # Begrenzt Arbeitsspeicher auf einen Prozentsatz des verfügbaren RAMs
    total_mem <- NA
    if (.Platform$OS.type == "windows") {
      mem_lines <- system("wmic OS get TotalVisibleMemorySize /Value", intern = TRUE)
      mem_line <- mem_lines[grepl("^TotalVisibleMemorySize=", mem_lines)]
      if (length(mem_line) > 0) {
        total_mem <- as.numeric(sub("TotalVisibleMemorySize=", "", mem_line))
      }
    }
    # Fallback: 16GB, falls nicht ermittelbar
    if (is.na(total_mem) || total_mem <= 0) {
      warning("Konnte Gesamtspeicher nicht automatisch ermitteln, verwende 16GB als Standard.")
      total_mem <- 16 * 1024 * 1024 # in KB
    }
    max_mem <- floor(total_mem * max_mem_percent / 100) * 1024  # In Bytes

    # Konfiguriere terra Optionen
    terra::terraOptions(memfrac = max_mem_percent / 100,
                 memmax = max_mem,
                 progress = 10,
                 verbose = TRUE,
                 todisk = TRUE,
                 memmin = 0.1 * max_mem,
                 datatype = "FLT4S")

    cat("Terra konfiguriert: ", max_mem_percent, "% RAM (", round(max_mem/1e9, 1),
        "GB), GDAL Cache:", cache_size, "MB\n")
  }
}

# Konfiguriere terra für optimale Speichernutzung
configure_terra <- function(max_mem_percent = 50, cache_size = 100) {
  # Begrenzt Arbeitsspeicher auf einen Prozentsatz des verfügbaren RAMs
  total_mem <- NA
  if (.Platform$OS.type == "windows") {
    mem_lines <- system("wmic OS get TotalVisibleMemorySize /Value", intern = TRUE)
    mem_line <- mem_lines[grepl("^TotalVisibleMemorySize=", mem_lines)]
    if (length(mem_line) > 0) {
      total_mem <- as.numeric(sub("TotalVisibleMemorySize=", "", mem_line))
    }
  }
  # Fallback: 16GB, falls nicht ermittelbar
  if (is.na(total_mem) || total_mem <= 0) {
    warning("Konnte Gesamtspeicher nicht automatisch ermitteln, verwende 16GB als Standard.")
    total_mem <- 16 * 1024 * 1024 # in KB
  }
  max_mem <- floor(total_mem * max_mem_percent / 100) * 1024  # In Bytes

  # Konfiguriere terra Optionen
  terra::terraOptions(memfrac = max_mem_percent / 100,
               memmax = max_mem,
               progress = 10,
               verbose = TRUE,
               todisk = TRUE,
               memmin = 0.1 * max_mem,
               datatype = "FLT4S")

  # GDAL Cache-Größe (in MB)
  # terraOptions(gdal_cache_max = cache_size) # <-- Entfernt, da kein gültiger terra-Parameter

  cat("Terra konfiguriert: ", max_mem_percent, "% RAM (", round(max_mem/1e9, 1),
      "GB), GDAL Cache:", cache_size, "MB\n")
}

# Optimiere GDAL-Einstellungen für Speichereffizienz
optimize_gdal <- function() {
  # Besser für Speichernutzung bei großen Dateien
  Sys.setenv(GDAL_CACHEMAX = "512")
  Sys.setenv(GDAL_DISABLE_READDIR_ON_OPEN = "TRUE")
  Sys.setenv(CPL_VSIL_CURL_ALLOWED_EXTENSIONS = "TIF,tif,BSQ,bsq")
  Sys.setenv(VSI_CACHE = "TRUE")
  Sys.setenv(VSI_CACHE_SIZE = "50000000")  # 50MB Cache

  cat("GDAL für Speichereffizienz optimiert\n")
}



#' Savitzky-Golay Filter für Rasterdaten entlang der spektralen Dimension
#'
#' @param raster_obj Ein terra SpatRaster Objekt
#' @param window_size Fenstergröße für den Filter (ungerade Zahl)
#' @param poly_order Polynomordnung für den Filter
#' @param chunk_size Größe der Chunks für Blockverarbeitung (Zeilen)
#' @return Geglättetes SpatRaster Objekt
apply_savitzky_golay <- function(raster_obj, window_size = 11, poly_order = 3, chunk_size = 500) {
  cat("Wende Savitzky-Golay-Filterung auf spektrale Dimension an...\n")

  if (window_size %% 2 == 0) {
    window_size <- window_size + 1
    cat("Fenstergröße wurde auf", window_size, "angepasst (muss ungerade sein).\n")
  }

  n_bands <- terra::nlyr(raster_obj)
  if (n_bands <= window_size) {
    stop("Zu wenig spektrale Bänder (", n_bands, ") für Fenstergröße ", window_size)
  }
  cat("Verarbeite", n_bands, "Spektralbänder\n")

  # Erstelle Ausgaberaster mit gleicher Struktur
  smoothed_raster <- terra::rast(raster_obj)

  # Bestimme die Gesamtzahl der Zeilen und Chunks
  nrows <- terra::nrow(raster_obj)
  ncols <- terra::ncol(raster_obj) # Get number of columns
  n_chunks <- ceiling(nrows / chunk_size)
  cat("Räumliche Verarbeitung in", n_chunks, "Chunks (je", chunk_size, "Zeilen)\n")

  # Schreibe das geglättete Raster auf temporäre Datei
  temp_file <- tempfile(fileext = ".tif")
  terra::writeStart(smoothed_raster, temp_file, overwrite=TRUE,
             gdal=c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))

  pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)

  # Verarbeite jeden Chunk separat
  for (chunk in 1:n_chunks) {
    # Berechne die Zeilen für diesen Chunk
    start_row <- (chunk - 1) * chunk_size + 1
    end_row <- min(start_row + chunk_size - 1, nrows)
    rows_to_process <- end_row - start_row + 1

    cat("\nVerarbeite Chunk", chunk, "von", n_chunks, "(Zeilen", start_row, "bis", end_row, ")\n")

    chunk_data <- NULL # Initialize chunk_data to NULL
    start_cell <- (start_row - 1) * ncols + 1 # Calculate start_cell here for logging

    tryCatch({
      chunk_data <- tryCatch({
        # Try the current syntax first
        vals <- terra::values(raster_obj, row = start_row, nrows = rows_to_process)
        # Ensure correct shape: rows = rows_to_process * ncols, cols = bands
        if (!is.null(vals) && is.matrix(vals) && nrow(vals) != rows_to_process * ncols) {
          vals <- matrix(vals, nrow = rows_to_process * ncols, ncol = terra::nlyr(raster_obj))
        }
        vals
      }, error = function(e) {
        # Alternative approach if the above fails
        cat("Versuche alternative Methode für Datenzugriff...\n")
        tryCatch({
          cell_start <- (start_row - 1) * ncols + 1
          cell_end <- min(cell_start + rows_to_process * ncols - 1, terra::ncell(raster_obj))
          cells <- cell_start:cell_end
          vals <- terra::extract(raster_obj, cells)
          if (is.data.frame(vals) && "ID" %in% colnames(vals)) {
            vals <- vals[ , !(colnames(vals) %in% "ID"), drop=FALSE]
          }
          m <- as.matrix(vals)
          if (nrow(m) != length(cells) || ncol(m) != terra::nlyr(raster_obj)) {
            m <- matrix(m, nrow=length(cells), ncol=terra::nlyr(raster_obj))
          }
          m
        }, error = function(e2) {
          cat("FEHLER beim Datenzugriff:", e2$message, "\n")
          stop("Konnte keine Daten aus dem Raster extrahieren")
        })
      })

      # Diagnostische Information
      cat("Chunk Dimensionen:", nrow(chunk_data), "Zeilen x", ncol(chunk_data), "Bänder\n")
      cat("NA-Anteil:", round(100 * sum(is.na(chunk_data)) / length(chunk_data), 1), "%\n")

      # Robuste Prüfung, ob der Chunk leer oder nur NAs enthält
      if (is.null(chunk_data) || !is.matrix(chunk_data) || nrow(chunk_data) == 0) {
        cat("Warnung: Chunk", chunk, "ist NULL oder leer, wird übersprungen.\n")
        next
      }

      if (all(is.na(chunk_data))) {
        cat("Warnung: Chunk", chunk, "enthält nur NAs, wird übersprungen.\n")
        next
      }

      # Verarbeite jeden Pixel (Zeile der Matrix)
      for (pixel_idx in 1:nrow(chunk_data)) {
        pixel_spectrum <- chunk_data[pixel_idx, ]
        if (all(is.na(pixel_spectrum))) next

        # Interpoliere NA-Werte im Spektrum, wenn möglich
        if (any(is.na(pixel_spectrum))) {
          valid_bands <- which(!is.na(pixel_spectrum))
          if (length(valid_bands) > window_size) {
            for (b in 1:length(pixel_spectrum)) {
              if (is.na(pixel_spectrum[b])) {
                before <- max(which(!is.na(pixel_spectrum[1:b])), 0)
                after <- min(which(!is.na(pixel_spectrum[b:n_bands])) + b - 1, n_bands)
                if (before > 0 && after <= n_bands) {
                  pixel_spectrum[b] <- pixel_spectrum[before] +
                    (pixel_spectrum[after] - pixel_spectrum[before]) *
                    ((b - before) / (after - before))
                }
              }
            }
          }
        }

        # Wende Savitzky-Golay-Filter an, wenn genügend gültige Werte vorhanden
        valid_values <- which(!is.na(pixel_spectrum))
        if (length(valid_values) > window_size) {
          smoothed_spectrum <- try(
            savgol(pixel_spectrum[valid_values], fl = (window_size-1)/2, forder = poly_order),
            silent = TRUE
          )
          if (!inherits(smoothed_spectrum, "try-error")) {
            pixel_spectrum[valid_values] <- smoothed_spectrum
            chunk_data[pixel_idx, ] <- pixel_spectrum
          }
        }
      }

      # Write values - Fix dimensional mismatch issue
      if (!is.null(chunk_data) && is.matrix(chunk_data)) {
        # cell_start <- (start_row - 1) * ncols + 1 # This was redundant here and not for writeValues' start

        # Debug information
        cat("Matrix dimensions:", nrow(chunk_data), "x", ncol(chunk_data), "\n")
        cat("Expected dimensions:", rows_to_process * ncols, "x", terra::nlyr(raster_obj), "\n")

        # Check if matrix needs reshaping
        if (nrow(chunk_data) != rows_to_process * ncols) {
          cat("Reshaping matrix to match expected dimensions...\n")
          # Only reshape if dimensions are compatible
          if (nrow(chunk_data) * ncol(chunk_data) == rows_to_process * ncols * terra::nlyr(raster_obj)) {
            chunk_data <- matrix(as.vector(chunk_data),
                               nrow = rows_to_process * ncols,
                               ncol = terra::nlyr(raster_obj))
          } else {
            cat("WARNUNG: Nicht kompatible Dimensionen für Reshaping!\n")
            # If reshaping fails or is not possible, this might lead to an error in writeValues
            # or incorrect data being written. Consider stopping or handling this case.
          }
        }

        # Final check before writing
        if (nrow(chunk_data) == rows_to_process * ncols) {
          cat("Schreibe Werte: start_row =", start_row, ", nrows =", rows_to_process, "...\n")
          # Pass start_row as start, and rows_to_process as nrows
          terra::writeValues(smoothed_raster, chunk_data, start=start_row, nrows=rows_to_process)
          cat("writeValues erfolgreich für Chunk", chunk, "\n")
        } else {
          cat("WARNUNG: Unerwartete Chunk-Dimensionen, überspringe Schreiben.\n")
          cat("Erwartete Anzahl Elemente:", rows_to_process * ncols,
              "Tatsächliche Anzahl:", nrow(chunk_data), "\n")
        }
      } else {
         if (is.null(chunk_data)) {
            cat("WARNUNG: chunk_data ist NULL. Überspringe Schreiben.\n")
         } else if (!is.matrix(chunk_data)) {
            cat("WARNUNG: chunk_data ist keine Matrix. Überspringe Schreiben.\n")
            cat("Typ von chunk_data:", class(chunk_data), "\n")
         }
      }

    }, error = function(e) {
      # Enhanced error reporting
      cat("-------------------- ERROR --------------------\n")
      cat("FEHLER bei Chunk", chunk, ":", conditionMessage(e), "\n")
      # Attempt to capture the specific call if possible
      call_info <- try(deparse(conditionCall(e)), silent = TRUE)
      if (!inherits(call_info, "try-error")) {
        cat("Function Call:", call_info, "\n")
      }
      cat("Parameters passed to writeValues: start_row =", start_row, ", nrows =", rows_to_process, "\n")
      cat("Calculated cell_start for this chunk (for reference/extract):", (start_row - 1) * ncols + 1, "\n")
      cat("Raster Info: nrow =", terra::nrow(smoothed_raster), ", ncol =", terra::ncol(smoothed_raster), ", ncell =", terra::ncell(smoothed_raster), "\n")
      cat("Chunk Info: start_row =", start_row, ", rows_to_process =", rows_to_process, "\n")
      cat("Chunk Data Dims:", if(!is.null(chunk_data) && is.matrix(chunk_data)) paste(dim(chunk_data), collapse="x") else "NULL or not matrix", "\n")
      cat("-------------------------------------------------\n")
      # Stop execution to prevent cascading errors
      stop("Chunk processing failed due to writeValues error. See details above.")
    })

    # Speicher freigeben
    gc(verbose = FALSE)
    setTxtProgressBar(pb, chunk)
    if (chunk %% 5 == 0) {
      memory_cleanup()
    }
  }

  # Schließe das Schreiben ab
  terra::writeStop(smoothed_raster)
  close(pb)

  # Create a fresh connection to the temporary file to break any file locks
  # This is the key fix - we reload the raster from the file instead of using the original object
  result_raster <- try(terra::rast(temp_file), silent = TRUE)

  if(!inherits(result_raster, "try-error")) {
    smoothed_raster <- result_raster
    cat("\nTemporär gespeichertes Raster erfolgreich neu geladen.\n")
  } else {
    cat("\nWarnung: Konnte temporäres Raster nicht neu laden. Versuche mit originalem Objekt.\n")
    # Keep using original smoothed_raster
  }

  cat("\nSpektrale Glättung abgeschlossen.\n")
  # Return both raster and temp_file for fallback
  return(list(raster = smoothed_raster, temp_file = temp_file))
}

#' Verbesserte Speicherbereinigung
#'
memory_cleanup <- function() {
  # Setze nicht genutzte Objekte auf NULL
  rm(list = ls(envir = .GlobalEnv, all.names = TRUE, pattern = "^temp_|^_tmp_"),
     envir = .GlobalEnv)

  # Force aggressive garbage collection multiple times
  for (i in 1:3) {
    gc(full = TRUE, verbose = FALSE)
  }

  # Auf Windows auch versuchen, Speicher ans OS zurückzugeben
  if (.Platform$OS.type == "windows") {
    invisible(tryCatch({
      memory.size(max = FALSE)
    }, error = function(e) NULL))
  }

  # Cache leeren in terra
  terra::tmpFiles(remove = TRUE)
}

#' Verbessertes Einrichten temporäres Verzeichnis
#'
setup_temp_directory <- function() {
  # Erstelle ein temporäres Verzeichnis auf einer schnellen Festplatte wenn möglich
  # SSD bevorzugen, falls verfügbar
  system_drives <- c("C:", "D:", "E:", "F:")
  available_drives <- system_drives[sapply(paste0(system_drives, "/"), dir.exists)]

  # Prüfe Speicherplatz auf verfügbaren Laufwerken
  drive_space <- numeric(length(available_drives))
  for (i in seq_along(available_drives)) {
    drive_info <- system(paste0("wmic logicaldisk where DeviceID='",
                               available_drives[i], "' get FreeSpace"), intern = TRUE)
    space_line <- drive_info[grep("^[0-9]", drive_info)]
    if (length(space_line) > 0) {
      drive_space[i] <- as.numeric(space_line)
    }
  }

  # Wähle Laufwerk mit meistem freien Speicher
  best_drive <- available_drives[which.max(drive_space)]
  temp_dir <- file.path(best_drive, "R_temp", "enmap_processing")

  if (!dir.exists(temp_dir)) {
    dir.create(temp_dir, recursive = TRUE)
  }

  # Setze den temporären Verzeichnispfad für R
  old_temp_dir <- Sys.getenv("TMPDIR")
  Sys.setenv(TMPDIR = temp_dir)

  # Setze den temporären Verzeichnispfad für terra
  terra::terraOptions(tempdir = temp_dir)

  cat("Temporäres Verzeichnis eingerichtet:", temp_dir, "\n")
  cat("Freier Speicherplatz:", round(max(drive_space)/1e9, 1), "GB\n")

  # Gib das alte temp dir zurück,
  return(old_temp_dir)
}

# Downsampling-Funktion für Entwicklungszwecke
#' @param raster_obj Ein terra SpatRaster Objekt
#' @param factor Faktor für das Downsampling (z.B. 2 = halbe Auflösung)
downsample_raster <- function(raster_obj, factor = 2) {
  if (factor <= 1) return(raster_obj)

  # Berechne neue Dimensionen
  new_nrow <- ceiling(nrow(raster_obj) / factor)
  new_ncol <- ceiling(ncol(raster_obj) / factor)

  # Erstelle ein neues Raster mit reduzierter Auflösung
  aggregated <- aggregate(raster_obj, fact=factor, fun="mean", na.rm=TRUE)

  return(aggregated)
}

# Robuste Funktion zum Schreiben eines Rasters mit Fallbacks
write_raster_with_fallbacks <- function(raster_obj, input_file, output_suffix, temp_file = NULL) {
  # Schreibe das geglättete Raster in das gleiche Verzeichnis wie das Eingabebild
  input_basename <- tools::file_path_sans_ext(basename(input_file))
  output_dir <- dirname(input_file)
  out_file <- file.path(output_dir, paste0(input_basename, output_suffix, ".tif"))
  ok <- FALSE
  tryCatch({
    terra::writeRaster(raster_obj, filename=out_file, overwrite=TRUE,
                       gdal=c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
    ok <- file.exists(out_file) && file.size(out_file) > 0
  }, error=function(e) {
    cat("FEHLER: Konnte Bild nicht speichern:", out_file, "\n")
    ok <- FALSE
    # Fallback: try copying temp_file if available
    if (!is.null(temp_file) && file.exists(temp_file)) {
      cat("Versuche temporäre Datei als Fallback zu kopieren...\n")
      # Ensure all connections to temp_file are closed
      gc()
      Sys.sleep(1) # Give OS time to release file lock
      # Try to create output directory if it does not exist
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      }
      # Retry copy up to 3 times
      for (attempt in 1:3) {
        success <- tryCatch({
          file.copy(temp_file, out_file, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
        }, error = function(e) FALSE)
        if (success && file.exists(out_file) && file.size(out_file) > 0) {
          ok <- TRUE
          cat("Fallback erfolgreich: Temporäre Datei kopiert.\n")
          break
        } else {
          cat(sprintf("Fallback Versuch %d fehlgeschlagen.\n", attempt))
          Sys.sleep(1)
        }
      }
      if (!ok) {
        # Try to copy to a different location for debugging
        alt_out_file <- file.path(getwd(), basename(out_file))
        alt_success <- tryCatch({
          file.copy(temp_file, alt_out_file, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
        }, error = function(e) FALSE)
        if (alt_success && file.exists(alt_out_file) && file.size(alt_out_file) > 0) {
          cat("Konnte temporäre Datei an alternativen Ort kopieren:", alt_out_file, "\n")
        }
        cat("Fallback fehlgeschlagen: Konnte temporäre Datei nicht kopieren.\n")
        cat("Temp-Dateipfad:", temp_file, "\n")
        cat("Zielpfad:", out_file, "\n")
        if (file.exists(temp_file)) {
          cat("Temp-Dateigröße:", file.size(temp_file), "Bytes\n")
        } else {
          cat("Temp-Datei existiert nicht mehr.\n")
        }
      }
    }
  })
  return(ok)
}

# # Funktion zur Korrektur der Spektralbänder nach dem Sonnenhöhenwinkel
# correct_sun_elevation <- function(raster_obj, sun_elevation_angle) {
#   cat("Korrigiere Spektralbänder nach Sonnenhöhenwinkel:", sun_elevation_angle, "Grad\n")
  
#   # Konvertiere Winkel von Grad zu Radian
#   sun_elevation_rad <- sun_elevation_angle * pi / 180
  
#   # Berechne den Korrekturfaktor (1 / sin(Winkel))
#   correction_factor <- 1 / sin(sun_elevation_rad)
  
#   cat("Korrekturfaktor:", round(correction_factor, 4), "\n")
  
#   # Wende Korrektur auf alle Bänder an
#   corrected_raster <- raster_obj * correction_factor
  
#   return(corrected_raster)
# }

# Calculate striping ratio (robust: column vs row mean stddev)
calculate_striping_ratio <- function(band_matrix) {
  # Skip if matrix contains all NAs
  if(all(is.na(band_matrix))) {
    return(NA)
  }
  # Remove rows/cols with all NA
  band_matrix <- band_matrix[rowSums(is.na(band_matrix)) < ncol(band_matrix), , drop=FALSE]
  band_matrix <- band_matrix[, colSums(is.na(band_matrix)) < nrow(band_matrix), drop=FALSE]
  if (nrow(band_matrix) < 2 || ncol(band_matrix) < 2) return(NA)
  # Calculate stddev of column means (vertical stripes) and row means (horizontal stripes)
  col_means <- apply(band_matrix, 2, function(x) mean(x, na.rm=TRUE))
  row_means <- apply(band_matrix, 1, function(x) mean(x, na.rm=TRUE))
  col_sd <- sd(col_means, na.rm=TRUE)
  row_sd <- sd(row_means, na.rm=TRUE)
  if (row_sd == 0) return(NA)
  ratio <- col_sd / row_sd
  return(ratio)
}

# Determine optimal strength based on column/row mean stddev ratio
determine_optimal_strength <- function(band_matrix) {
  # Remove rows/cols with all NA
  band_matrix <- band_matrix[rowSums(is.na(band_matrix)) < ncol(band_matrix), , drop=FALSE]
  band_matrix <- band_matrix[, colSums(is.na(band_matrix)) < nrow(band_matrix), drop=FALSE]
  if (nrow(band_matrix) < 2 || ncol(band_matrix) < 2) return(0.8)
  col_means <- apply(band_matrix, 2, function(x) mean(x, na.rm=TRUE))
  row_means <- apply(band_matrix, 1, function(x) mean(x, na.rm=TRUE))
  col_sd <- sd(col_means, na.rm=TRUE)
  row_sd <- sd(row_means, na.rm=TRUE)
  if (is.na(col_sd) || is.na(row_sd) || row_sd == 0) return(0.8)
  ratio <- col_sd / row_sd
  # Map ratio to strength: ratio 1 → 0.8, ratio 10 → 0.99, clamp
  strength <- min(0.99, max(0.8, 0.8 + 0.19 * (min(ratio, 10) - 1) / 9))
  return(strength)
}

# ================================================================================
# HORN & WOODHAM DESTRIPING IMPLEMENTATION
# ================================================================================
# 
# Reference Paper: "DESTRIPING SATELLITE IMAGES" 
# Authors: B. K. P. Horn and R. J. Woodham
# Institution: Massachusetts Institute of Technology - Artificial Intelligence Laboratory
#
# OVERVIEW:
# This implementation addresses systematic striping artifacts in satellite imagery
# caused by sensor calibration differences. The method assumes linear, time-invariant
# sensors with individual gain and offset characteristics.
#
# MATHEMATICAL FOUNDATION:
# Each sensor has a linear transfer function: x' = a + b*x
# Where:
#   - x' = observed sensor output (what we measure)
#   - x = true scene radiance (what we want to recover)
#   - a = sensor offset (bias)
#   - b = sensor gain
#
# CORRECTION FORMULA:
# To recover scene radiance: x = (x' - a)/b
#
# STATISTICAL ESTIMATION (when ground truth is unavailable):
# - Gain estimation: b = σ'/σ (ratio of standard deviations)
# - Offset estimation: a = (μ'σ - μσ')/σ
# Where μ', σ' are sensor statistics and μ, σ are reference statistics
#
# ADVANTAGES OVER SIMPLE MEAN SUBTRACTION:
# 1. Preserves contrast through gain correction
# 2. Physically motivated (based on actual sensor characteristics)
# 3. Handles varying noise levels across sensors
# 4. Mathematically rigorous approach
#
# IMPLEMENTATION DETAILS:
# - Each column represents one sensor's output
# - Whole-image statistics serve as reference
# - Robust estimation prevents extreme corrections
# - Quality metrics quantify improvement
# ================================================================================

# Horn & Woodham Destriping Implementation
# Reference: "DESTRIPING SATELLITE IMAGES" by B. K. P. Horn and R. J. Woodham
# Massachusetts Institute of Technology - Artificial Intelligence Laboratory
#
# This implementation follows the linear sensor model approach where each sensor
# has individual gain and offset characteristics that can be estimated statistically.
destripe_raster_horn_woodham <- function(raster_obj, preserve_statistics = TRUE, 
                                        robust_estimation = TRUE, min_valid_pixels = 100,
                                        verbose = TRUE, quality_assessment = TRUE,
                                        use_enmap_specs = TRUE, adaptive_destriping = TRUE,
                                        ratio_threshold_high = 1.15,     # Oberer Schwellenwert
                                        ratio_threshold_low = 0.85) {    # Unterer Schwellenwert
  if(verbose) cat("  Applying Horn & Woodham destriping method...\n")
  if(verbose) cat("    Reference: Horn & Woodham, MIT AI Lab\n")
  
  # Load official EnMAP sensor configuration
  sensor_config <- NULL
  if (use_enmap_specs) {
    sensor_config <- get_enmap_official_sensor_config()
    if(verbose) cat("    Using official EnMAP HSI sensor specifications\n")
  }
  
  if(adaptive_destriping) {
    if(verbose) cat("    Adaptive destriping enabled - only skip bands with ratio between", ratio_threshold_low, "and", ratio_threshold_high, "\n")
  }
  
  n_bands <- nlyr(raster_obj)
  
  # Initialize quality metrics for assessment
  quality_metrics <- data.frame(
    band = 1:n_bands,
    striping_ratio_before = numeric(n_bands),
    striping_ratio_after = numeric(n_bands),
    mean_gain_correction = numeric(n_bands),
    mean_offset_correction = numeric(n_bands),
    valid_columns = numeric(n_bands),
    destriping_applied = logical(n_bands),  # Angewendet oder nicht
    skip_reason = character(n_bands)        # Grund für Überspringen
  )
  
  if(verbose) cat("  Processing", n_bands, "spectral bands using Horn & Woodham method...\n")
  
  # Process bands in chunks for memory efficiency
  if(n_bands > 100) {
    if(verbose) cat("  Large number of bands detected. Processing in chunks...\n")
    chunk_size <- 50
    for(chunk_start in seq(1, n_bands, by=chunk_size)) {
      chunk_end <- min(chunk_start + chunk_size - 1, n_bands)
      if(verbose) cat("  Processing bands", chunk_start, "to", chunk_end, "of", n_bands, "\n")
      
      for(i in chunk_start:chunk_end) {
        # Calculate striping ratio before correction
        band_mat <- as.matrix(raster_obj[[i]], wide=TRUE)
        ratio_before <- calculate_striping_ratio(band_mat)
        quality_metrics$striping_ratio_before[i] <- ratio_before
        
        # KORRIGIERTE Entscheidungslogik: Nur "gute" Bänder überspringen
        should_destripe <- TRUE
        skip_reason <- ""
        
        if(adaptive_destriping && !is.na(ratio_before)) {
          # Nur Bänder im "guten" Bereich überspringen (zwischen den Schwellenwerten)
          if(ratio_before >= ratio_threshold_low && ratio_before <= ratio_threshold_high) {
            # Ratio ist im guten Bereich - NICHT bearbeiten
            should_destripe <- FALSE
            skip_reason <- "good_quality"
          } else if(ratio_before < ratio_threshold_low) {
            # Ratio ist zu niedrig - BEARBEITEN (mögliche horizontale Streifen)
            should_destripe <- TRUE
            skip_reason <- "horizontal_stripes_suspected"
          } else {
            # Ratio > ratio_threshold_high - BEARBEITEN (vertikale Streifen)
            should_destripe <- TRUE
            skip_reason <- "vertical_stripes_detected"
          }
        } else {
          # Fallback: immer bearbeiten wenn adaptive destriping aus ist
          skip_reason <- "processed_default"
        }
        
        quality_metrics$destriping_applied[i] <- should_destripe
        quality_metrics$skip_reason[i] <- skip_reason
        
        if(should_destripe) {
          # Apply Horn & Woodham destriping with EnMAP specifications
          result <- destripe_single_band_horn_woodham(
            raster_obj[[i]], 
            preserve_statistics = preserve_statistics,
            robust_estimation = robust_estimation,
            min_valid_pixels = min_valid_pixels,
            verbose = FALSE,
            band_number = i,
            sensor_config = sensor_config
          )
          
          # Calculate striping ratio after correction
          band_mat_after <- as.matrix(result$raster, wide=TRUE)
          ratio_after <- calculate_striping_ratio(band_mat_after)
          
          # QUALITÄTSKONTROLLE: Prüfe ob Korrektur tatsächlich verbessert hat
          if(!is.na(ratio_after) && !is.na(ratio_before)) {
            improvement_factor <- ratio_after / ratio_before
            
            # Für horizontale Streifen: Akzeptiere auch wenn Ratio steigt (aber näher zu 1.0 kommt)
            is_horizontal_case <- ratio_before < ratio_threshold_low
            gets_closer_to_one <- abs(ratio_after - 1.0) < abs(ratio_before - 1.0)
            
            # Verwerfe Korrektur nur wenn sie deutlich verschlechtert UND nicht näher zu 1.0 kommt
            if(improvement_factor > 1.5 && !gets_closer_to_one) {
              if(verbose) cat(sprintf("    WARNUNG Band %d: Korrektur verschlechtert (%.3f → %.3f), wird verworfen\n", 
                                      i, ratio_before, ratio_after))
              # Behalte das ursprüngliche Band
              quality_metrics$striping_ratio_after[i] <- ratio_before
              quality_metrics$mean_gain_correction[i] <- 1.0
              quality_metrics$mean_offset_correction[i] <- 0.0
              quality_metrics$valid_columns[i] <- 0
              quality_metrics$skip_reason[i] <- "correction_failed"
              quality_metrics$destriping_applied[i] <- FALSE
            } else {
              # Korrektur wird akzeptiert
              raster_obj[[i]] <- result$raster
              quality_metrics$mean_gain_correction[i] <- result$mean_gain
              quality_metrics$mean_offset_correction[i] <- result$mean_offset
              quality_metrics$valid_columns[i] <- result$valid_columns
              quality_metrics$striping_ratio_after[i] <- ratio_after
              
              # Spezielle Markierung für horizontale Streifen
              if(is_horizontal_case && gets_closer_to_one) {
                quality_metrics$skip_reason[i] <- "horizontal_stripes_improved"
              }
            }
          } else {
            # Falls Ratio-Berechnung fehlschlägt, behalte Original
            quality_metrics$striping_ratio_after[i] <- ratio_before
            quality_metrics$mean_gain_correction[i] <- 1.0
            quality_metrics$mean_offset_correction[i] <- 0.0
            quality_metrics$valid_columns[i] <- 0
            quality_metrics$skip_reason[i] <- "ratio_calculation_failed"
            quality_metrics$destriping_applied[i] <- FALSE
          }
        } else {
          # Keine Korrektur angewendet - Werte bleiben unverändert
          quality_metrics$striping_ratio_after[i] <- ratio_before
          quality_metrics$mean_gain_correction[i] <- 1.0
          quality_metrics$mean_offset_correction[i] <- 0.0
          quality_metrics$valid_columns[i] <- 0
        }
        
        if(verbose && i %% 25 == 0) {
          status_map <- c(
            "good_quality" = "SKIPPED (good)",
            "vertical_stripes_detected" = "PROCESSED (vertical)",
            "horizontal_stripes_suspected" = "PROCESSED (horizontal)",
            "horizontal_stripes_improved" = "PROCESSED (horizontal improved)",
            "correction_failed" = "FAILED (discarded)",
            "processed_default" = "PROCESSED"
          )
          
          status <- status_map[quality_metrics$skip_reason[i]]
          if(is.na(status)) status <- "UNKNOWN"
          
          cat(sprintf("    Band %d: ratio %.3f → %.3f (gain: %.3f, offset: %.3f) [%s]\n", 
                      i, quality_metrics$striping_ratio_before[i], 
                      quality_metrics$striping_ratio_after[i], 
                      quality_metrics$mean_gain_correction[i],
                      quality_metrics$mean_offset_correction[i], status))
        }
      }
      gc() # Force garbage collection after each chunk
    }
  } else {
    # Process all bands at once for smaller rasters - GLEICHE LOGIK
    pb <- txtProgressBar(min = 1, max = n_bands, style = 3)
    for(i in 1:n_bands) {
      # Calculate striping ratio before correction
      band_mat <- as.matrix(raster_obj[[i]], wide=TRUE)
      ratio_before <- calculate_striping_ratio(band_mat)
      quality_metrics$striping_ratio_before[i] <- ratio_before
      
      # KORRIGIERTE Entscheidungslogik
      should_destripe <- TRUE
      skip_reason <- ""
      
      if(adaptive_destriping && !is.na(ratio_before)) {
        if(ratio_before >= ratio_threshold_low && ratio_before <= ratio_threshold_high) {
          should_destripe <- FALSE
          skip_reason <- "good_quality"
        } else if(ratio_before < ratio_threshold_low) {
          should_destripe <- TRUE
          skip_reason <- "horizontal_stripes_suspected"
        } else {
          should_destripe <- TRUE
          skip_reason <- "vertical_stripes_detected"
        }
      } else {
        skip_reason <- "processed_default"
      }
      
      quality_metrics$destriping_applied[i] <- should_destripe
      quality_metrics$skip_reason[i] <- skip_reason
      
      if(should_destripe) {
        # Apply Horn & Woodham destriping
        result <- destripe_single_band_horn_woodham(
          raster_obj[[i]], 
          preserve_statistics = preserve_statistics,
          robust_estimation = robust_estimation,
          min_valid_pixels = min_valid_pixels,
          verbose = FALSE,
          band_number = i,
          sensor_config = sensor_config
        )
        
        # Calculate striping ratio after correction
        band_mat_after <- as.matrix(result$raster, wide=TRUE)
        ratio_after <- calculate_striping_ratio(band_mat_after)
        
        # Quality control (same as above)
        if(!is.na(ratio_after) && !is.na(ratio_before)) {
          improvement_factor <- ratio_after / ratio_before
          is_horizontal_case <- ratio_before < ratio_threshold_low
          gets_closer_to_one <- abs(ratio_after - 1.0) < abs(ratio_before - 1.0)
          
          if(improvement_factor > 1.5 && !gets_closer_to_one) {
            # Discard correction
            quality_metrics$striping_ratio_after[i] <- ratio_before
            quality_metrics$mean_gain_correction[i] <- 1.0
            quality_metrics$mean_offset_correction[i] <- 0.0
            quality_metrics$valid_columns[i] <- 0
            quality_metrics$skip_reason[i] <- "correction_failed"
            quality_metrics$destriping_applied[i] <- FALSE
          } else {
            # Accept correction
            raster_obj[[i]] <- result$raster
            quality_metrics$mean_gain_correction[i] <- result$mean_gain
            quality_metrics$mean_offset_correction[i] <- result$mean_offset
            quality_metrics$valid_columns[i] <- result$valid_columns
            quality_metrics$striping_ratio_after[i] <- ratio_after
            
            if(is_horizontal_case && gets_closer_to_one) {
              quality_metrics$skip_reason[i] <- "horizontal_stripes_improved"
            }
          }
        } else {
          # Fallback
          quality_metrics$striping_ratio_after[i] <- ratio_before
          quality_metrics$mean_gain_correction[i] <- 1.0
          quality_metrics$mean_offset_correction[i] <- 0.0
          quality_metrics$valid_columns[i] <- 0
          quality_metrics$skip_reason[i] <- "ratio_calculation_failed"
          quality_metrics$destriping_applied[i] <- FALSE
        }
      } else {
        # No correction applied
        quality_metrics$striping_ratio_after[i] <- ratio_before
        quality_metrics$mean_gain_correction[i] <- 1.0
        quality_metrics$mean_offset_correction[i] <- 0.0
        quality_metrics$valid_columns[i] <- 0
      }
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  
  # Report improvement metrics
  if(verbose && quality_assessment) {
    # Separate metrics for different categories
    good_quality_bands <- quality_metrics$skip_reason == "good_quality"
    vertical_stripes_bands <- quality_metrics$skip_reason %in% c("vertical_stripes_detected")
    horizontal_stripes_bands <- quality_metrics$skip_reason %in% c("horizontal_stripes_suspected", "horizontal_stripes_improved")
    failed_bands <- quality_metrics$skip_reason == "correction_failed"
    
    avg_before <- mean(quality_metrics$striping_ratio_before, na.rm=TRUE)
    avg_after <- mean(quality_metrics$striping_ratio_after, na.rm=TRUE)
    improvement <- (1 - avg_after/avg_before) * 100
    
    cat("  ✓ Horn & Woodham destriping complete\n")
    cat("    Bands skipped (good quality 0.85-1.15):", sum(good_quality_bands, na.rm=TRUE), "\n")
    cat("    Bands processed (vertical stripes >1.15):", sum(vertical_stripes_bands, na.rm=TRUE), "\n")
    cat("    Bands processed (horizontal stripes <0.85):", sum(horizontal_stripes_bands, na.rm=TRUE), "\n")
    cat("    Bands failed (correction discarded):", sum(failed_bands, na.rm=TRUE), "\n")
    cat("    Average striping ratio before:", round(avg_before, 4), "\n")
    cat("    Average striping ratio after:", round(avg_after, 4), "\n")
    cat("    Overall improvement:", round(improvement, 1), "%\n")
  }
  
  return(list(raster = raster_obj, metrics = quality_metrics))
}

# Horn & Woodham single band destriping implementation
# 
# Mathematical foundation from Horn & Woodham (1979):
# - Sensor model: observed = gain × true + offset
# - Correction: true = (observed - offset) / gain
# - Statistical estimation:
#   * Gain: b = s'/s (sensor_std / reference_std)
#   * Offset: a = (x̄'·s - x̄·s')/s
#
# This method assumes that each column should have the same statistical distribution
# as the entire image, and corrects for systematic gain and offset errors.
destripe_single_band_horn_woodham <- function(band_rast, preserve_statistics = FALSE, 
                                             robust_estimation = TRUE, min_valid_pixels = 50,
                                             verbose = FALSE, band_number = NULL, sensor_config = NULL) {
  # Convert to matrix for processing
  mat <- as.matrix(band_rast, wide=TRUE)
  
  # Skip if too many NAs
  if(sum(is.na(mat)) / length(mat) > 0.9) {
    if(verbose) cat("    Skipping band with >90% NAs\n")
    return(list(raster = band_rast, mean_gain = 1.0, mean_offset = 0.0, valid_columns = 0))
  }
  
  # Get EnMAP-specific parameters for this band
  correction_params <- NULL
  if (!is.null(sensor_config) && !is.null(band_number)) {
    correction_params <- get_enmap_band_correction_params(band_number, sensor_config)
    if(verbose) cat("    Using", correction_params$sensor_type, "parameters for band", band_number, "\n")
    
    # Override min_valid_pixels with EnMAP-specific value
    min_valid_pixels <- correction_params$min_valid_pixels
  }
  
  # Create NA mask
  na_mask <- is.na(mat)
  
  # Calculate reference (target) statistics
  # Horn & Woodham: Each column should have the same statistical distribution
  # We use the median column statistics as the target (more robust than mean)
  col_means <- apply(mat, 2, mean, na.rm=TRUE)
  col_sds <- apply(mat, 2, sd, na.rm=TRUE)
  
  # Use median statistics as reference (target that all columns should match)
  reference_mean <- median(col_means, na.rm=TRUE)
  reference_sd <- median(col_sds, na.rm=TRUE)
  
  # Alternative: Use overall image statistics (original approach)
  # reference_mean <- mean(mat, na.rm=TRUE)
  # reference_sd <- sd(mat, na.rm=TRUE)
  
  if(is.na(reference_mean) || is.na(reference_sd) || reference_sd == 0) {
    if(verbose) cat("    Invalid reference statistics, skipping\n")
    return(list(raster = band_rast, mean_gain = 1.0, mean_offset = 0.0, valid_columns = 0))
  }
  
  # Initialize correction tracking
  gains <- numeric(ncol(mat))
  offsets <- numeric(ncol(mat))
  valid_columns <- 0
  
  # Apply Horn & Woodham correction column by column
  # Each column represents one detector/sensor element
  for(col in 1:ncol(mat)) {
    col_data <- mat[, col]
    
    # Skip columns with insufficient valid data
    valid_pixels <- sum(!is.na(col_data))
    if(valid_pixels < min_valid_pixels) {
      gains[col] <- 1.0
      offsets[col] <- 0.0
      next
    }
    
    valid_columns <- valid_columns + 1
    
    # Calculate column statistics
    col_mean <- mean(col_data, na.rm=TRUE)
    col_sd <- sd(col_data, na.rm=TRUE)
    
    # Handle edge cases
    if(is.na(col_mean) || is.na(col_sd) || col_sd == 0) {
      gains[col] <- 1.0
      offsets[col] <- 0.0
      next
    }
    
    # Horn & Woodham statistical estimation
    # Gain: b = s'/s (sensor_std / reference_std)
    gain <- col_sd / reference_sd
    
    # Offset: a = (x̄'·s - x̄·s')/s
    # where x̄' = col_mean, s' = col_sd, x̄ = reference_mean, s = reference_sd
    offset <- (col_mean * reference_sd - reference_mean * col_sd) / reference_sd
    
    # Apply EXTREMELY CONSERVATIVE robust estimation
    if(robust_estimation) {
      # Use EnMAP-specific gain and offset ranges if available
      if (!is.null(correction_params)) {
        # EnMAP-specific clamping ranges - EXTREMELY CONSERVATIVE
        gain_min <- correction_params$gain_range[1]
        gain_max <- correction_params$gain_range[2]
        offset_range <- correction_params$offset_range
        
        # Apply extremely tight gain clamping
        gain <- max(gain_min, min(gain_max, gain))
        
        # Apply extremely tight offset clamping
        max_offset <- max(abs(offset_range)) * 0.05  # Extremely conservative factor
        offset <- max(-max_offset, min(max_offset, offset))
        
        # Apply very weak correction strength
        correction_factor <- correction_params$correction_strength
        gain <- correction_factor * gain + (1 - correction_factor) * 1.0
        offset <- correction_factor * offset
        
        # Special handling for atmospheric absorption bands
        if (correction_params$is_atmospheric) {
          # Almost no correction for atmospheric features
          atmospheric_factor <- correction_params$atmospheric_strength * 0.02  # Extremely weak
          gain <- atmospheric_factor * gain + (1 - atmospheric_factor) * 1.0
          offset <- atmospheric_factor * offset
        }
        
        # Special handling for overlap region (900-1000 nm)
        if (correction_params$is_overlap) {
          # Virtually no correction in detector overlap region
          gain <- 0.02 * gain + 0.98 * 1.0  # Almost no correction
          offset <- 0.02 * offset           # Almost no offset correction
        }
      } else {
        # Default robust estimation - EXTREMELY CONSERVATIVE
        gain <- max(0.99, min(1.01, gain))  # Extremely narrow range
        data_range <- max(mat, na.rm=TRUE) - min(mat, na.rm=TRUE)
        max_offset <- 0.005 * data_range  # Extremely small offset range
        offset <- max(-max_offset, min(max_offset, offset))
        
        # Apply extremely weak overall correction strength
        gain <- 0.1 * gain + 0.9 * 1.0  # Extremely conservative
        offset <- 0.1 * offset          # Extremely conservative
      }
    }
    
    # Store correction values
    gains[col] <- gain
    offsets[col] <- offset
    
    # Apply Horn & Woodham correction: true = (observed - offset) / gain
    # Only correct non-NA values
    valid_mask <- !is.na(col_data)
    if(sum(valid_mask) > 0) {
      corrected_values <- (col_data[valid_mask] - offset) / gain
      mat[valid_mask, col] <- corrected_values
    }
  }
  
  # Restore NAs
  mat[na_mask] <- NA
  
  # Create output raster with corrected values
  corrected_rast <- rast(mat, extent=ext(band_rast), crs=crs(band_rast))
  names(corrected_rast) <- names(band_rast)
  
  # Calculate average correction statistics
  mean_gain <- mean(gains[gains != 1.0], na.rm=TRUE)
  mean_offset <- mean(offsets[offsets != 0.0], na.rm=TRUE)
  
  if(is.na(mean_gain)) mean_gain <- 1.0
  if(is.na(mean_offset)) mean_offset <- 0.0
  
  return(list(
    raster = corrected_rast,
    mean_gain = mean_gain,
    mean_offset = mean_offset,
    valid_columns = valid_columns,
    gains = gains,
    offsets = offsets
  ))
}


# ================================================================================
# EnMAP HSI OFFICIAL SENSOR SPECIFICATIONS
# ================================================================================
# 
# Based on official EnMAP HSI Instrument Specification:
# - Spectral range: 420 nm - 2450 nm
# - VNIR: 420-1000 nm, 6.5 nm sampling, SNR > 500 @ 495 nm
# - SWIR: 900-2450 nm, 10 nm sampling, SNR > 150 @ 2200 nm
# - Radiometric resolution: ≥ 14 bits
# - Radiometric accuracy: 5% / stability: 2.5%
# ================================================================================

#' EnMAP HSI Sensor Configuration based on Official Specifications
get_enmap_official_sensor_config <- function() {
  
  vnir_bands <- 1:88    # Based on typical EnMAP band numbering
  swir_bands <- 89:224  # Total 224 bands typical for EnMAP
  
  list(
    # VNIR Detector Array (Silicon-based)
    vnir = list(
      bands = vnir_bands,
      wavelength_range = c(420, 1000),          # nm (official spec)
      spectral_sampling = 6.5,                  # nm (official spec)
      snr_reference = 500,                      # @ 495 nm (official spec)
      radiometric_accuracy = 5.0,               # % (official spec)
      radiometric_stability = 2.5,              # % (official spec)
      
      # Destriping parameters - EVEN MORE CONSERVATIVE
      typical_gain_range = c(0.98, 1.02),       # Even narrower range
      typical_offset_range = c(-10, 10),        # Even smaller offset range
      correction_strength = 0.15,               # Much weaker correction
      min_valid_pixels = 200,                   # Higher threshold
      robust_clamp_factor = 0.6                 # Even tighter clamping
    ),
    
    # SWIR Detector Array (InGaAs-based)
    swir = list(
      bands = swir_bands,
      wavelength_range = c(900, 2450),          # nm (official spec)
      spectral_sampling = 10.0,                 # nm (official spec)
      snr_reference = 150,                      # @ 2200 nm (official spec)
      radiometric_accuracy = 5.0,               # % (official spec)
      radiometric_stability = 2.5,              # % (official spec)
      
      # Destriping parameters - EVEN MORE CONSERVATIVE
      typical_gain_range = c(0.985, 1.015),     # Even narrower range
      typical_offset_range = c(-8, 8),          # Even smaller offset range
      correction_strength = 0.12,               # Even weaker correction
      min_valid_pixels = 250,                   # Higher threshold
      robust_clamp_factor = 0.5                 # Even tighter clamping
    ),
    
    # Overlap region (900-1000 nm) - both detectors active
    overlap_region = list(
      bands = 85:92,                            # Approximate overlap
      correction_strength = 0.08,               # Minimal correction
      robust_clamp_factor = 0.3                 # Very tight clamping
    ),
    
    # Known atmospheric absorption features
    atmospheric_absorption = list(
      water_940nm = list(bands = 78:82, strength = 0.05),    # Very weak
      water_1130nm = list(bands = 108:112, strength = 0.08), # Very weak
      water_1380nm = list(bands = 137:142, strength = 0.03), # Very weak
      water_1900nm = list(bands = 187:192, strength = 0.05), # Very weak
      co2_2000nm = list(bands = 200:205, strength = 0.1)     # Very weak
    )
  )
}

#' Get band-specific correction parameters based on official EnMAP specifications
get_enmap_band_correction_params <- function(band_number, sensor_config) {
  
  # Determine base sensor parameters
  if (band_number %in% sensor_config$vnir$bands) {
    base_params <- sensor_config$vnir
    sensor_type <- "VNIR"
  } else if (band_number %in% sensor_config$swir$bands) {
    base_params <- sensor_config$swir
    sensor_type <- "SWIR"
  } else {
    # Fallback
    base_params <- sensor_config$swir
    sensor_type <- "UNKNOWN"
  }
  
  # Check for special regions
  is_overlap <- band_number %in% sensor_config$overlap_region$bands
  
  # Check for atmospheric absorption features
  is_atmospheric <- FALSE
  atmospheric_strength <- 1.0
  
  for (feature in sensor_config$atmospheric_absorption) {
    if (band_number %in% feature$bands) {
      is_atmospheric <- TRUE
      atmospheric_strength <- feature$strength
      break
    }
  }
  
  # Calculate correction strength based on specifications
  correction_strength <- base_params$correction_strength
  
  # Modify for special cases
  if (is_overlap) {
    correction_strength <- sensor_config$overlap_region$correction_strength
  }
  
  if (is_atmospheric) {
    correction_strength <- correction_strength * atmospheric_strength
  }
  
  # Calculate expected gain/offset ranges based on radiometric specifications
  radiometric_factor <- base_params$radiometric_accuracy / 100
  gain_spread <- radiometric_factor * 2  # Convert to gain range
  
  # Adjust gain range based on SNR
  snr_factor <- base_params$snr_reference / 500.0  # Normalize to VNIR reference
  gain_range <- c(
    1.0 - gain_spread / snr_factor,
    1.0 + gain_spread / snr_factor
  )
  
  # Ensure gain range stays within reasonable bounds
  gain_range[1] <- max(0.3, gain_range[1])
  gain_range[2] <- min(3.0, gain_range[2])
  
  # Calculate offset range based on radiometric stability
  stability_factor <- base_params$radiometric_stability / 100
  offset_range <- c(-200 * stability_factor, 200 * stability_factor)
  
  return(list(
    sensor_type = sensor_type,
    band_number = band_number,
    gain_range = gain_range,
    offset_range = offset_range,
    correction_strength = correction_strength,
    min_valid_pixels = base_params$min_valid_pixels,
    robust_clamp_factor = base_params$robust_clamp_factor,
    snr_reference = base_params$snr_reference,
    is_overlap = is_overlap,
    is_atmospheric = is_atmospheric,
    atmospheric_strength = atmospheric_strength
  ))
}

# Hauptfunktion zum Glätten der Bilder
smooth_spectral_images <- function(
  spectral_dir = "data/extracted_files",
  output_suffix = "_smoothed+destripedx",
  window_size = 11,
  poly_order = 3,
  chunk_size = 500,
  dev_mode = FALSE,
  dev_factor = 4,
  overwrite = FALSE,
  start_index = 1,
  metadata_file = "data/metadata_mapping.csv",
  apply_destriping = TRUE,
  destriping_preserve_statistics = FALSE,
  destriping_robust_estimation = TRUE,
  destriping_min_valid_pixels = 50,
  adaptive_destriping = TRUE,        # Neuer Parameter
  ratio_threshold_high = 1.15,       # Oberer Schwellenwert
  ratio_threshold_low = 0.85,        # Unterer Schwellenwert
  shapefile_path = NULL,             # Shapefile mit Koordinaten zum Filtern
  apply_normalization = TRUE,
  normalization_method = "rayleigh"
) {
  # Prüfe, ob die Daten existieren
  if (!file.exists(metadata_file)) {
    stop("Metadaten-Datei nicht gefunden: ", metadata_file)
  }
  if (!dir.exists(spectral_dir)) {
    stop("Spektralbildverzeichnis nicht gefunden: ", spectral_dir)
  }

  # Konfiguriere terra und GDAL für optimale Speichernutzung
  configure_terra(max_mem_percent = 60, cache_size = 512)
  optimize_gdal()
  old_temp_dir <- setup_temp_directory()
  
  # Metadaten mit Sonnenhöhenwinkeln laden
  cat("Lade Metadaten mit Sonnenhöhenwinkeln aus", metadata_file, "...\n")
  metadata <- tryCatch({
    readr::read_csv(metadata_file, show_col_types = FALSE)
  }, error = function(e) {
    stop("Fehler beim Laden der Metadatendatei: ", e$message)
  })
  
  # Überprüfe, ob die erforderlichen Spalten in der Metadatendatei vorhanden sind
  if (!all(c("SpectralFile", "sun_elevation") %in% colnames(metadata))) {
    stop("Metadatendatei muss die Spalten 'SpectralFile' und 'sun_elevation' enthalten.")
  }
  
  # Bereinige Duplikate in den Metadaten und erstelle ein Lookup-Dictionary
  metadata <- metadata %>% distinct(SpectralFile, .keep_all = TRUE)
  
  # Load shapefile to filter raster files based on coordinates
  raster_files_to_process <- NULL
  if (!is.null(shapefile_path) && file.exists(shapefile_path)) {
    cat("Lade Shapefile aus", shapefile_path, "zum Filtern der Raster-Dateien...\n")
    
    tryCatch({
      # Lade Shapefile mit sf
      shp <- sf::st_read(shapefile_path, quiet = TRUE)
      cat("Shapefile geladen mit", nrow(shp), "Punkten\n")
      
      # Extrahiere die Koordinaten
      coords <- sf::st_coordinates(shp)
      
      if (nrow(coords) > 0) {
        cat("Extrahierte", nrow(coords), "Koordinaten aus Shapefile\n")
        
        # Finde alle verfügbaren Raster-Dateien (TIF oder BSQ)
        all_spectral_files <- list.files(
          spectral_dir,
          pattern = "SPECTRAL_IMAGE\\.(TIF|BSQ)$",
          full.names = FALSE,
          recursive = TRUE,
          ignore.case = TRUE
        )
        
        cat("Gefundene Spektralbilder:", length(all_spectral_files), "\n")
        
        # Bestimme, welche Raster die Koordinaten enthalten
        raster_files_to_process <- c()
        
        cat("Prüfe, welche Raster die Shapefile-Koordinaten enthalten...\n")
        for (raster_file in all_spectral_files) {
          raster_path <- file.path(spectral_dir, raster_file)
          
          # Versuche, das Raster zu laden und die Ausdehnung zu prüfen
          tryCatch({
            r <- terra::rast(raster_path)
            ext <- terra::ext(r)
            
            # Prüfe, ob mindestens eine Koordinate innerhalb der Raster-Ausdehnung liegt
            coords_in_extent <- coords[, 1] >= ext[1] & coords[, 1] <= ext[2] &
                                coords[, 2] >= ext[3] & coords[, 2] <= ext[4]
            
            if (any(coords_in_extent)) {
              raster_files_to_process <- c(raster_files_to_process, basename(raster_file))
              cat("  ✓", basename(raster_file), "enthält", sum(coords_in_extent), "Koordinaten\n")
            }
          }, error = function(e) {
            cat("  ⚠ Warnung: Konnte", basename(raster_file), "nicht prüfen:", e$message, "\n")
          })
        }
        
        if (length(raster_files_to_process) > 0) {
          cat("Filter aktiv: Nur", length(raster_files_to_process), "Raster mit Shapefile-Koordinaten werden verarbeitet\n")
        } else {
          cat("⚠ Warnung: Keine Raster enthalten die Shapefile-Koordinaten. Alle Raster werden verarbeitet.\n")
          raster_files_to_process <- NULL
        }
      }
    }, error = function(e) {
      cat("⚠ Warnung: Konnte Shapefile nicht verarbeiten:", e$message, "\n")
      cat("Alle Raster werden verarbeitet.\n")
    })
  } else if (!is.null(shapefile_path)) {
    cat("⚠ Warnung: Shapefile nicht gefunden:", shapefile_path, "\n")
    cat("Alle Raster werden verarbeitet.\n")
  }

  # Extrahiere alle Spektralbilddateien
  spectral_files <- list.files(
    spectral_dir,
    pattern = "SPECTRAL_IMAGE\\.(BSQ|TIF)$",
    full.names = TRUE,
    recursive = TRUE
  )
  if (length(spectral_files) == 0) {
    stop("Keine Spektralbilder gefunden in: ", spectral_dir)
  }
  cat("Gefundene Spektralbilder insgesamt:", length(spectral_files), "\n")

  # Improved matching logic
  if (!is.null(raster_files_to_process)) {
    # Use the filter list - match spectral files to the ones in filter CSV
    cat("Filtere Spektralbilder basierend auf Filter-CSV...\n")
    
    # Extract basenames from filter list for matching
    # Remove the -SPECTRAL_IMAGE suffix from filter basenames if present
    filter_basenames <- gsub("-SPECTRAL_IMAGE\\.(BSQ|TIF)$", "", basename(raster_files_to_process))
    filter_basenames <- gsub("-SPECTRAL_IMAGE$", "", filter_basenames)  # Remove -SPECTRAL_IMAGE suffix
    filter_basenames <- gsub("\\.(BSQ|TIF)$", "", filter_basenames)     # Also remove .BSQ/.TIF if present
    
    # Extract basenames from found spectral files (these already have -SPECTRAL_IMAGE removed)
    spectral_basenames <- gsub("-SPECTRAL_IMAGE\\.(BSQ|TIF)$", "", basename(spectral_files))
    
    # Find matches using the corrected basenames
    matching_indices <- which(spectral_basenames %in% filter_basenames)
    
    if (length(matching_indices) == 0) {
      cat("DEBUG: Keine Übereinstimmungen gefunden nach erster Bereinigung.\n")
      cat("Filter-Basenames (erste 5):", head(filter_basenames, 5), "\n")
      cat("Spektral-Basenames (erste 5):", head(spectral_basenames, 5), "\n")
      
      # Try alternative matching approach - more aggressive cleaning
      cat("Versuche erweiterte Matching-Strategien...\n")
      
      # Strategy 1: Try even more aggressive basename extraction
      # Remove common prefixes and suffixes to get core identifiers
      filter_core <- gsub("^ENMAP01-____L2A-", "", filter_basenames)
      spectral_core <- gsub("^ENMAP01-____L2A-", "", spectral_basenames)
      
      matching_indices <- which(spectral_core %in% filter_core)
      
      if (length(matching_indices) == 0) {
        # Strategy 2: Try partial matching on the DT identifier
        cat("Versuche Matching über DT-Identifikatoren...\n")
        
        # Extract DT identifiers (e.g., DT0000028792_20230720T151638Z_003)
        filter_dt_pattern <- regmatches(filter_basenames, regexpr("DT[0-9]+_[0-9T]+Z_[0-9]+", filter_basenames))
        spectral_dt_pattern <- regmatches(spectral_basenames, regexpr("DT[0-9]+_[0-9T]+Z_[0-9]+", spectral_basenames))
        
        # Find matches based on DT patterns
        for (i in seq_along(filter_dt_pattern)) {
          if (length(filter_dt_pattern[i]) > 0) {
            matches <- which(spectral_dt_pattern == filter_dt_pattern[i])
            if (length(matches) > 0) {
              matching_indices <- c(matching_indices, matches)
            }
          }
        }
        matching_indices <- unique(matching_indices)
        
        if (length(matching_indices) == 0) {
          # Strategy 3: Most aggressive - match any overlapping substring
          cat("Versuche substring Matching...\n")
          for (i in seq_along(filter_basenames)) {
            pattern <- filter_basenames[i]
            # Try to find spectral files that contain parts of this pattern
            for (j in seq_along(spectral_basenames)) {
              # Check if they share a significant common substring
              if (nchar(pattern) > 10 && nchar(spectral_basenames[j]) > 10) {
                # Extract date-time parts for comparison
                filter_datetime <- regmatches(pattern, regexpr("20[0-9]{6}T[0-9]{6}Z", pattern))
                spectral_datetime <- regmatches(spectral_basenames[j], regexpr("20[0-9]{6}T[0-9]{6}Z", spectral_basenames[j]))
                
                if (length(filter_datetime) > 0 && length(spectral_datetime) > 0) {
                  # Check if datetimes are within a few seconds of each other
                  if (abs(as.numeric(difftime(
                    as.POSIXct(gsub("T|Z", "", filter_datetime), format="%Y%m%d%H%M%S"),
                    as.POSIXct(gsub("T|Z", "", spectral_datetime), format="%Y%m%d%H%M%S"),
                    units="secs"))) < 30) {  # Within 30 seconds
                    matching_indices <- c(matching_indices, j)
                  }
                }
              }
            }
          }
          matching_indices <- unique(matching_indices)
        }
      }
      
      if (length(matching_indices) == 0) {
        # Final debug output
        cat("FINAL DEBUG OUTPUT:\n")
        cat("Beispiel Filter-Filename:", raster_files_to_process[1], "\n")
        cat("Beispiel Spektral-Filename:", basename(spectral_files[1]), "\n")
        cat("Filter-Core (erste 3):", head(filter_core, 3), "\n")
        cat("Spektral-Core (erste 3):", head(spectral_core, 3), "\n")
        
        stop("Keine Spektralbilder gefunden, die mit den Dateien im Filter-CSV übereinstimmen.\n",
             "Die Dateinamen-Formate scheinen sich zu unterscheiden.\n",
             "Filter-Format: ", basename(raster_files_to_process[1]), "\n",
             "Spektral-Format: ", basename(spectral_files[1]), "\n",
             "Prüfen Sie die Dateinamen in: ", filter_csv_file)
      } else {
        cat("Erfolgreiches Matching mit erweiterten Strategien gefunden!\n")
      }
    }
    
    spectral_files <- spectral_files[matching_indices]
    cat("Nach Filterung verbleibende Bilder:", length(spectral_files), "\n")
    
    # Debug: Show which files were matched
    if (length(spectral_files) > 0) {
      cat("Erste 3 gematchte Spektral-Dateien:\n")
      for (i in 1:min(3, length(spectral_files))) {
        cat("  ", basename(spectral_files[i]), "\n")
      }
    }
    
  } else {
    # Kein Shapefile-Filter: Verwende alle verfügbaren Raster
    # Optional: Filtere basierend auf Metadaten (wenn RasterFile-Spalte vorhanden)
    if ("RasterFile" %in% colnames(metadata)) {
      cat("Filtere Spektralbilder basierend auf Metadaten (RasterFile-Spalte)...\n")
      raster_files_in_metadata <- unique(metadata$RasterFile)
      raster_files_in_metadata <- raster_files_in_metadata[!is.na(raster_files_in_metadata)]
      
      if (length(raster_files_in_metadata) > 0) {
        # Normalisiere Pfade für Vergleich (verwende basename und entferne Backslashes)
        metadata_basenames <- basename(raster_files_in_metadata)
        spectral_basenames <- basename(spectral_files)
        
        # Finde Übereinstimmungen basierend auf Dateinamen
        matching_indices <- which(spectral_basenames %in% metadata_basenames)
        
        if (length(matching_indices) > 0) {
          spectral_files <- spectral_files[matching_indices]
          cat("Nach Metadaten-Filterung:", length(spectral_files), "Bilder\n")
        } else {
          cat("⚠ Warnung: Keine Übereinstimmung mit Metadaten. Verarbeite alle gefundenen Raster.\n")
          cat("   Metadata enthält:", length(metadata_basenames), "Einträge\n")
          cat("   Gefunden wurden:", length(spectral_basenames), "Spektralbilder\n")
          if (length(metadata_basenames) > 0) cat("   Beispiel Metadata:", metadata_basenames[1], "\n")
          if (length(spectral_basenames) > 0) cat("   Beispiel Spektral:", spectral_basenames[1], "\n")
        }
      }
    } else {
      cat("Keine RasterFile-Spalte in Metadaten. Verarbeite alle gefundenen Raster.\n")
    }
  }

  cat("Zu verarbeitende Bilder:", length(spectral_files), "\n")

  # Prüfe und passe den Startindex an
  if (start_index < 1) start_index <- 1
  if (start_index > length(spectral_files)) {
    cat("Startindex überschreitet Anzahl der Bilder. Beende.\n")
    return(invisible(NULL))
  }

  cat("Beginne Verarbeitung bei Bild", start_index, "von", length(spectral_files), "\n")
  cat("Neue Ausgabedateien werden mit Suffix '", output_suffix, "' erstellt\n")
  
  for (i in start_index:length(spectral_files)) {
    file_path <- spectral_files[i]
    
    # Bestimme die zugehörige Zeile in den Metadaten
    # Entferne den Suffix "-SPECTRAL_IMAGE.BSQ" oder "-SPECTRAL_IMAGE.TIF" vom Dateinamen
    file_base <- basename(file_path)
    file_base <- gsub("-SPECTRAL_IMAGE\\.(BSQ|TIF)$", "", file_base, ignore.case = TRUE)
    
    meta_row <- metadata[metadata$SpectralFile == file_base, ]
    
    if (nrow(meta_row) == 0) {
      cat("WARNUNG: Keine Metadaten gefunden für Datei:", file_base, "\n")
      cat("  Gesuchter Dateiname:", file_base, "\n")
      cat("  Verfügbare SpectralFiles:", paste(head(metadata$SpectralFile, 3), collapse = ", "), "\n")
      next
    }
    
    # Extrahiere Sonnenhöhenwinkel
    sun_elevation_angle <- meta_row$sun_elevation
    
    # Korrigiere Spektralbänder nach Sonnenhöhenwinkel (optional)
    # raster_obj <- correct_sun_elevation(raster_obj, sun_elevation_angle)
    
    # Lade das Rasterbild
    cat("Lade Rasterbild:", basename(file_path), "...\n")
    raster_obj <- terra::rast(file_path)
    
    # Downsampling für Entwicklungszwecke (optional)
    if (dev_mode) {
      cat("Entwicklungsmodus aktiv: Downsampling des Rasters um Faktor", dev_factor, "...\n")
      raster_obj <- downsample_raster(raster_obj, factor = dev_factor)
    }
    
    # Speichere temporär das heruntergestufte Bild (falls Downsampling durchgeführt)
    temp_file <- NULL
    if (dev_mode) {
      temp_file <- tempfile(fileext = ".tif")
      cat("Speichere temporäre Datei für Entwicklungszwecke:", temp_file, "\n")
      terra::writeRaster(raster_obj, filename=temp_file, overwrite=TRUE,
                         gdal=c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
    }
    
    # Wende Savitzky-Golay-Filter an
    cat("Wende Savitzky-Golay-Filter an...\n")
    savitzky_result <- apply_savitzky_golay(raster_obj, window_size = window_size, 
                                            poly_order = poly_order, chunk_size = chunk_size)
    smoothed_raster <- savitzky_result$raster
    
    # Speicherort für das geglättete Bild
    output_file <- file.path(dirname(file_path), 
                             paste0(tools::file_path_sans_ext(basename(file_path)), output_suffix, ".tif"))
    
    # Schreibe das geglättete Bild
    cat("Speichere geglättetes Bild:", output_file, "\n")
    success <- write_raster_with_fallbacks(smoothed_raster, file_path, output_suffix, temp_file = temp_file)
    
    if (!success) {
      cat("FEHLER: Konnte das geglättete Bild nicht speichern. Verarbeitung abgebrochen.\n")
      return(invisible(NULL))
    }
    
    # Speicher freigeben
    gc(verbose = FALSE)
    
    # Destriping anwenden (falls gewünscht)
    if (apply_destriping) {
      cat("Wende adaptives Horn & Woodham Destriping an...\n")
      
      # Get sensor configuration
      sensor_config <- get_enmap_official_sensor_config()
      
      cat("  Verwende EnMAP-Sensorkonfiguration für alle Bänder\n")
      
      destripe_result <- destripe_raster_horn_woodham(
        smoothed_raster, 
        preserve_statistics = destriping_preserve_statistics,
        robust_estimation = destriping_robust_estimation,
        min_valid_pixels = destriping_min_valid_pixels,
        verbose = TRUE,
        quality_assessment = TRUE,
        adaptive_destriping = adaptive_destriping,     
        ratio_threshold_high = ratio_threshold_high,   # Oberer Schwellenwert
        ratio_threshold_low = ratio_threshold_low      # Unterer Schwellenwert
      )
      
      destriped_raster <- destripe_result$raster
      
      # Speicherort für das destripte Bild
      destriped_output_file <- file.path(dirname(file_path), 
                                         paste0(tools::file_path_sans_ext(basename(file_path)), 
                                                output_suffix, "_destriped.tif"))
      
      # Schreibe das destripte Bild
      cat("Speichere destriptes Bild:", destriped_output_file, "\n")
      success <- write_raster_with_fallbacks(destriped_raster, file_path, 
                                             paste0(output_suffix, "_destriped"), temp_file = temp_file)
      
      if (!success) {
        cat("FEHLER: Konnte das destripte Bild nicht speichern. Verarbeitung abgebrochen.\n")
        return(invisible(NULL))
      }
    }
    
    # Normalisierung anwenden (falls gewünscht)
    if (apply_normalization) {
      cat("Wende Normalisierung an...\n")
      
      # Hier könnte die Normalisierungslogik implementiert werden
      # Als Platzhalter: Schreibe einfach das geglättete Bild als normalisiertes Bild
      normalized_raster <- smoothed_raster  # Placeholder for actual normalization logic
      
      # Speicherort für das normalisierte Bild
      normalized_output_file <- file.path(dirname(file_path), 
                                          paste0(tools::file_path_sans_ext(basename(file_path)), 
                                                 output_suffix, "_normalized.tif"))
      
      # Schreibe das normalisierte Bild
      cat("Speichere normalisiertes Bild:", normalized_output_file, "\n")
      success <- write_raster_with_fallbacks(normalized_raster, file_path, 
                                              paste0(output_suffix, "_normalized"), temp_file = temp_file)
      
      if (!success) {
        cat("FEHLER: Konnte das normalisierte Bild nicht speichern. Verarbeitung abgebrochen.\n")
        return(invisible(NULL))
      }
    }
  }
  
  cat("Verarbeitung abgeschlossen.\n")
  return(invisible(NULL))
}

# ================================================================================
# SCRIPT EXECUTION - MAIN PROCESSING CALL
# ================================================================================

# HINWEIS: Diese Sektion ist auskommentiert, damit das Script nur über die
# master_pipeline.R gesteuert wird und nicht automatisch beim Laden ausgeführt wird.
# Falls Sie dieses Script einzeln ausführen möchten, entfernen Sie die Kommentare.

# # Starte die Verarbeitung, falls Skript direkt ausgeführt wird
# last_index <- 1
# 
# # Use progress file for Horn & Woodham + EnMAP specs version
# if (file.exists("last_enmap_horn_woodham_destriped_image_index.txt")) {
#   last_index <- as.integer(readLines("last_enmap_horn_woodham_destriped_image_index.txt")) + 1
#   cat("Fortsetzen der Verarbeitung ab Bild", last_index, "\n")
# }
# 
# # Hauptverarbeitung mit adaptivem Horn & Woodham Destriping
# # Nutzt Shapefile anstatt CSV zum Filtern der zu verarbeitenden Raster
# cat("Starte Verarbeitung mit adaptivem Horn & Woodham Destriping und EnMAP-Spezifikationen...\n")
# 
# smooth_spectral_images(
#   spectral_dir = "data/extracted_files",
#   output_suffix = "_smoothed+adaptive_destriped",
#   window_size = 11,
#   poly_order = 3,
#   chunk_size = 500,
#   dev_mode = FALSE,
#   dev_factor = 4,
#   overwrite = TRUE,
#   start_index = last_index,
#   metadata_file = "data/metadata_mapping.csv",
#   apply_destriping = TRUE,
#   destriping_preserve_statistics = FALSE,
#   destriping_robust_estimation = TRUE,
#   destriping_min_valid_pixels = 200,
#   adaptive_destriping = TRUE,           
#   ratio_threshold_high = 1.15,         # Nur Bänder mit Ratio > 1.15 werden bearbeitet
#   ratio_threshold_low = 0.85,          # Bänder mit Ratio < 0.85 werden auch übersprungen
#   shapefile_path = "data/Peru_Bolivien_reprojiziert.shp",  # Shapefile statt CSV
#   apply_normalization = FALSE,
#   normalization_method = "rayleigh"
# )
# 
# cat("Verarbeitung abgeschlossen!\n")
# cat("Ergebnisse mit Suffix '_smoothed+adaptive_destriped' wurden erstellt.\n")