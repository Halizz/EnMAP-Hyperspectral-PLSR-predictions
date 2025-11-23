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

  # Gib das alte temp dir zurück, um es später wiederherzustellen
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

#' Scientific Destriping based on Adaptive Anisotropic Total Variation and Truncated Nuclear Norm
#' 
#' Implementation of the method described in:
#' Hu, T., Liu, N., Li, W., Tao, R., Zhang, F., & Scheunders, P.
#' "Destriping hyperspectral imagery by adaptive anisotropic total variation and truncated nuclear norm"
#' 
#' The model: Y = I + S, where Y is striped image, I is clean image, S is stripes
#' Objective: min_{I,S} {1/2||Y-I-S||_F^2 + ||I||_AATV + sum_b ||S_b||_r}

# AATT Destriping Implementation following Hu et al. (2021)
# Adaptive Anisotropic Total Variation and Truncated Nuclear Norm

# Soft thresholding function for ADMM algorithm
soft_threshold <- function(theta, phi) {
  sign(theta) * pmax(abs(theta) - phi, 0)
}

# Horizontal gradient operator (convolution with [1, -1])
horizontal_gradient <- function(mat) {
  grad <- array(0, dim = dim(mat))
  if (ncol(mat) > 1) {
    grad[, 1:(ncol(mat)-1)] <- mat[, 2:ncol(mat)] - mat[, 1:(ncol(mat)-1)]
  }
  return(grad)
}

# Divergence operator (adjoint of gradient)
divergence_horizontal <- function(grad) {
  div <- array(0, dim = dim(grad))
  if (ncol(grad) > 1) {
    div[, 1] <- -grad[, 1]
    div[, 2:(ncol(grad)-1)] <- grad[, 1:(ncol(grad)-2)] - grad[, 2:(ncol(grad)-1)]
    div[, ncol(grad)] <- grad[, ncol(grad)-1]
  }
  return(div)
}

# Truncated Nuclear Norm (TNM) for stripe modeling
truncated_nuclear_norm <- function(mat, r = 2) {
  if (any(dim(mat) == 0) || all(is.na(mat))) {
    return(mat)
  }
  
  # Handle NAs with interpolation
  na_mask <- is.na(mat)
  if (any(na_mask)) {
    # Row-wise interpolation
    for (i in 1:nrow(mat)) {
      row_na <- is.na(mat[i, ])
      if (any(row_na) && !all(row_na)) {
        mat[i, row_na] <- mean(mat[i, !row_na], na.rm = TRUE)
      }
    }
    # Column-wise interpolation
    for (j in 1:ncol(mat)) {
      col_na <- is.na(mat[, j])
      if (any(col_na) && !all(col_na)) {
        mat[col_na, j] <- mean(mat[!col_na, j], na.rm = TRUE)
      }
    }
    # Final fallback
    if (any(is.na(mat))) {
      overall_mean <- mean(mat, na.rm = TRUE)
      mat[is.na(mat)] <- if (!is.na(overall_mean)) overall_mean else 0
    }
  }
  
  # SVD decomposition for TNM
  svd_result <- tryCatch({
    svd(mat)
  }, error = function(e) {
    # Fallback for singular matrices
    return(mat)
  })
  
  # Apply truncated nuclear norm: keep only r largest singular values
  if (length(svd_result$d) > r) {
    svd_result$d[(r+1):length(svd_result$d)] <- 0
  }
  
  # Reconstruct matrix
  result <- svd_result$u %*% diag(svd_result$d, nrow = length(svd_result$d)) %*% t(svd_result$v)
  
  # Restore NAs
  result[na_mask] <- NA
  
  return(result)
}

# AATT Destriping for single band following Hu et al. (2021)
# Model: Y = I + S (striped image = clean image + stripe component)
# Objective: min_{I,S} {1/2||Y-I-S||_F^2 + ||I||_AATV + sum_b ||S_b||_r}
aatt_destripe_single_band <- function(Y_band, max_iter = 30, alpha = 0.1, beta_x = 0.1, beta_z = 0.1, 
                                     tol = 1e-5, verbose = FALSE, band_id = NULL) {
  # Ensure Y_band is a matrix
  if (is.vector(Y_band)) {
    Y_band <- matrix(Y_band, ncol = 1)
  }
  
  M <- nrow(Y_band)
  N <- ncol(Y_band)
  
  # Handle NAs
  na_mask <- is.na(Y_band)
  Y_band_filled <- Y_band
  if (any(na_mask)) {
    # Simple interpolation for NAs
    for (i in 1:M) {
      row_na <- is.na(Y_band_filled[i, ])
      if (any(row_na) && !all(row_na)) {
        Y_band_filled[i, row_na] <- mean(Y_band_filled[i, !row_na], na.rm = TRUE)
      }
    }
    for (j in 1:N) {
      col_na <- is.na(Y_band_filled[, j])
      if (any(col_na) && !all(col_na)) {
        Y_band_filled[col_na, j] <- mean(Y_band_filled[!col_na, j], na.rm = TRUE)
      }
    }
    if (any(is.na(Y_band_filled))) {
      overall_mean <- mean(Y_band_filled, na.rm = TRUE)
      Y_band_filled[is.na(Y_band_filled)] <- if (!is.na(overall_mean)) overall_mean else 0
    }
  }
  
  # ADMM Algorithm following Hu et al. (2021)
  # Initialization
  I <- Y_band_filled  # Clean image estimate
  S <- matrix(0, M, N)  # Stripe component
  
  # Auxiliary variables for ADMM
  W1 <- matrix(0, M, N)  # Auxiliary variable for horizontal gradient
  W2 <- matrix(0, M, N)  # Auxiliary variable for spectral gradient (simplified for 2D)
  
  # Lagrangian multipliers
  A1 <- matrix(0, M, N)  # Multiplier for horizontal gradient constraint
  A2 <- matrix(0, M, N)  # Multiplier for spectral gradient constraint
  
  prev_energy <- Inf
  
  for (k in 1:max_iter) {
    # Step 1: Update I (clean image) using closed-form solution
    # This involves solving: I = argmin ||Y - I - S||_F^2 + penalty terms
    
    # Data fidelity term
    data_term <- Y_band_filled - S
    
    # Regularization terms from ADMM constraints
    h_grad_I <- horizontal_gradient(I)
    reg_term1 <- divergence_horizontal(W1 - h_grad_I + A1/alpha)
    reg_term2 <- (W2 + A2/alpha)  # Simplified for 2D case
    
    # Update I with adaptive step size
    step_size <- 0.5 / (1 + alpha)  # Adaptive step size for stability
    I_new <- data_term + alpha * (reg_term1 + reg_term2)
    I <- (1 - step_size) * I + step_size * I_new
    
    # Step 2: Update S (stripe component) using Truncated Nuclear Norm
    residual <- Y_band_filled - I
    S <- truncated_nuclear_norm(residual, r = 2)
    
    # Step 3: Update auxiliary variables W1, W2 using soft thresholding
    h_grad_I <- horizontal_gradient(I)
    
    # Adaptive regularization parameters (key feature of AATT)
    # Lambda1 and Lambda2 adapt based on local gradient magnitudes
    Lambda1 <- pmin(abs(h_grad_I), beta_x)  # Horizontal direction
    Lambda2 <- pmin(abs(I), beta_z)         # Spectral direction (simplified)
    
    # Soft thresholding updates
    W1 <- soft_threshold(h_grad_I - A1/alpha, Lambda1/alpha)
    W2 <- soft_threshold(-A2/alpha, Lambda2/alpha)
    
    # Step 4: Update Lagrangian multipliers
    A1 <- A1 + alpha * (W1 - h_grad_I)
    A2 <- A2 + alpha * W2
    
    # Convergence check
    energy <- 0.5 * sum((Y_band_filled - I - S)^2, na.rm = TRUE)
    
    if (verbose && k %% 10 == 0) {
      band_info <- if(!is.null(band_id)) paste0("Band ", band_id, " ") else ""
      cat(sprintf("    AATT %sIteration %d: Energy = %.6f\n", band_info, k, energy))
    }
    
    # Check for convergence
    if (k > 1 && abs(prev_energy - energy) < tol) {
      if (verbose) cat(sprintf("    AATT converged at iteration %d\n", k))
      break
    }
    
    # Stability check - prevent divergence
    if (energy > 10 * sum(Y_band_filled^2, na.rm = TRUE)) {
      if (verbose) cat("    AATT: Energy diverging, using fallback method\n")
      # Fallback to simple column-wise correction
      col_means <- apply(Y_band_filled, 2, mean, na.rm = TRUE)
      overall_mean <- mean(col_means, na.rm = TRUE)
      col_deviations <- col_means - overall_mean
      for (j in 1:N) {
        if (abs(col_deviations[j]) > 2 * sd(col_deviations, na.rm = TRUE)) {
          I[, j] <- Y_band_filled[, j] - col_deviations[j] * 0.5
        }
      }
      break
    }
    
    prev_energy <- energy
  }
  
  # Restore NAs in the final result
  I[na_mask] <- NA
  
  return(I)
}

# AATT Destriping for hyperspectral raster following Hu et al. (2021)
# Uses standard parameters: alpha=0.1, beta_x=0.1, beta_z=0.1, max_iter=30
destripe_raster_aatt <- function(raster_obj, max_iter = 30, alpha = 0.1, beta_x = 0.1, beta_z = 0.1,
                                verbose = TRUE, preserve_stats = TRUE) {
  if(verbose) cat("  Applying AATT destriping (Hu et al. 2021 method)...\n")
  if(verbose) cat("  Parameters: ADMM penalty α=", alpha, ", gradient thresholds βx=", beta_x, ", βz=", beta_z, "\n")
  
  n_bands <- nlyr(raster_obj)
  
  # Store quality metrics
  quality_metrics <- data.frame(
    band = 1:n_bands,
    energy_before = numeric(n_bands),
    energy_after = numeric(n_bands),
    stripe_variance_reduction = numeric(n_bands),
    iterations_used = numeric(n_bands)
  )
  
  # Process bands in moderate chunks to balance memory and stability
  chunk_size <- min(15, n_bands)
  n_chunks <- ceiling(n_bands / chunk_size)
  
  if(verbose) {
    cat("  Processing", n_bands, "bands in", n_chunks, "chunks (chunk size:", chunk_size, ")...\n")
  }
  
  for(chunk in 1:n_chunks) {
    chunk_start <- (chunk - 1) * chunk_size + 1
    chunk_end <- min(chunk_start + chunk_size - 1, n_bands)
    
    if(verbose) {
      cat("  Processing chunk", chunk, "of", n_chunks, "(bands", chunk_start, "to", chunk_end, ")\n")
    }
    
    # Process each band in the current chunk
    for(i in chunk_start:chunk_end) {
      band_mat <- as.matrix(raster_obj[[i]], wide=TRUE)
      
      # Skip if band is empty or all NAs
      if (all(is.na(band_mat)) || nrow(band_mat) == 0 || ncol(band_mat) == 0) {
        if (verbose) cat("    Skipping empty band", i, "\n")
        next
      }
      
      # Calculate energy before destriping
      quality_metrics$energy_before[i] <- sum(band_mat^2, na.rm = TRUE)
      
      # Store original statistics for preservation
      original_mean <- mean(band_mat, na.rm = TRUE)
      original_sd <- sd(as.vector(band_mat), na.rm = TRUE)
      
      # Apply AATT destriping with standard parameters
      destriped_band <- aatt_destripe_single_band(
        band_mat, 
        max_iter = max_iter,    # Standard: 30 iterations
        alpha = alpha,          # Standard: 0.1 ADMM penalty
        beta_x = beta_x,        # Standard: 0.1 horizontal gradient threshold
        beta_z = beta_z,        # Standard: 0.1 spectral gradient threshold
        verbose = verbose && (i == chunk_start || i %% 10 == 0),
        band_id = i             # Pass band ID for clearer output
      )
      
      # Preserve original statistics if requested
      if (preserve_stats && !is.na(original_mean) && !is.na(original_sd) && original_sd > 0) {
        current_mean <- mean(destriped_band, na.rm = TRUE)
        current_sd <- sd(as.vector(destriped_band), na.rm = TRUE)
        
        if (!is.na(current_mean) && !is.na(current_sd) && current_sd > 0) {
          # Rescale to match original statistics
          destriped_band <- (destriped_band - current_mean) * (original_sd / current_sd) + original_mean
        }
      }
      
      # Calculate energy after destriping
      quality_metrics$energy_after[i] <- sum(destriped_band^2, na.rm = TRUE)
      
      # Estimate stripe variance reduction
      if (ncol(band_mat) > 1) {
        column_means_before <- apply(band_mat, 2, mean, na.rm = TRUE)
        column_means_after <- apply(destriped_band, 2, mean, na.rm = TRUE)
        stripe_var_before <- var(column_means_before, na.rm = TRUE)
        stripe_var_after <- var(column_means_after, na.rm = TRUE)
        
        if (!is.na(stripe_var_before) && stripe_var_before > 0) {
          quality_metrics$stripe_variance_reduction[i] <- 
            (stripe_var_before - stripe_var_after) / stripe_var_before
        }
      }
      
      # Update the raster band
      raster_obj[[i]] <- rast(destriped_band, extent=ext(raster_obj[[i]]), crs=crs(raster_obj[[i]]))
      names(raster_obj[[i]]) <- names(raster_obj)[i]
    }
    
    # Force garbage collection after each chunk
    gc()
  }
  
  # Report overall performance
  if(verbose) {
    avg_reduction <- mean(quality_metrics$stripe_variance_reduction, na.rm=TRUE) * 100
    energy_change <- mean((quality_metrics$energy_after - quality_metrics$energy_before) / 
                         quality_metrics$energy_before, na.rm=TRUE) * 100
    
    cat("  ✓ AATT destriping complete (Hu et al. 2021 method)\n")
    cat("    Average stripe variance reduction:", round(avg_reduction, 1), "%\n")
    cat("    Average energy change:", round(energy_change, 1), "%\n")
    cat("    Method: Adaptive Anisotropic Total Variation + Truncated Nuclear Norm\n")
    cat("    Applied smoothing only in striped regions, preserving details elsewhere\n")
  }
  
  return(list(raster = raster_obj, metrics = quality_metrics))
}

# Hauptfunktion zum Glätten der Bilder
smooth_spectral_images <- function(
  data_file = "data/final_before_correlation.csv",
  spectral_dir = "data/extracted_files/ich_everzwlfe",
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
  destriping_adaptive = TRUE,
  destriping_edge_tapering = TRUE,
  destriping_visualize = FALSE,
  filter_csv_file = "data_Sac/final_before_correlation.csv",  # New parameter for filtering
  apply_normalization = TRUE,
  normalization_method = "rayleigh"
) {
  # Prüfe, ob die Daten existieren
  if (!file.exists(data_file)) {
    stop("Datendatei nicht gefunden: ", data_file)
  }
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

  # Datendatei laden
  cat("Lade Daten aus", data_file, "...\n")
  data <- tryCatch({
    readr::read_csv(data_file, show_col_types = FALSE)
  }, error = function(e) {
    stop("Fehler beim Laden der Datendatei: ", e$message)
  })
  
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
  
  # Load filter CSV if provided to limit processing to specific raster files
  raster_files_to_process <- NULL
  if (!is.null(filter_csv_file) && file.exists(filter_csv_file)) {
    cat("Lade Filter-CSV aus", filter_csv_file, "...\n")
    filter_data <- tryCatch({
      readr::read_csv(filter_csv_file, show_col_types = FALSE)
    }, error = function(e) {
      cat("Warnung: Konnte Filter-CSV nicht laden:", e$message, "\n")
      NULL
    })
    
    if (!is.null(filter_data)) {
      # Find the raster file column
      filter_raster_col <- NULL
      for (col in c("RasterFile", "SpectralFile", "ImageFile")) {
        if (col %in% colnames(filter_data)) {
          filter_raster_col <- col
          break
        }
      }
      
      if (!is.null(filter_raster_col)) {
        raster_files_to_process <- unique(filter_data[[filter_raster_col]])
        cat("Filter aktiv: Nur", length(raster_files_to_process), "spezifische Raster-Dateien werden verarbeitet\n")
      }
    }
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
    # Original logic: match with data file
    spectral_basenames <- gsub("-SPECTRAL_IMAGE\\.(BSQ|TIF)$", "", basename(spectral_files))
    possible_cols <- c("RasterFile", "SpectralFile", "ImageFile")
    raster_col <- NULL
    for (col in possible_cols) {
      if (col %in% colnames(data)) {
        raster_col <- col
        break
      }
    }
    if (is.null(raster_col)) {
      stop("Keine Spalte mit Bildpfaden in den Daten gefunden.")
    }
    raster_files_in_data <- unique(data[[raster_col]])
    raster_basenames_in_data <- gsub("-SPECTRAL_IMAGE\\.(BSQ|TIF)$", "", raster_files_in_data)
    matching_indices <- which(spectral_basenames %in% raster_basenames_in_data)
    if (length(matching_indices) == 0) {
      stop("Keine Bilder aus den Daten im Spektralverzeichnis gefunden.")
    }
    spectral_files <- spectral_files[matching_indices]
  }

  cat("Zu verarbeitende Bilder:", length(spectral_files), "\n")

  # Prüfe und passe den Startindex an
  if (start_index < 1) start_index <- 1
  if (start_index > length(spectral_files)) {
    cat("Startindex überschreitet Anzahl der Bilder. Beende.\n")
    return(invisible(NULL))
  }

  cat("Beginne Verarbeitung bei Bild", start_index, "von", length(spectral_files), "\n")
  cat("Neue Ausgabedateien werden mit Suffix '", output_suffix, "' erstellt (AATT-Methode)\n")
  
  for (i in start_index:length(spectral_files)) {
    input_file <- spectral_files[i]
    input_basename <- tools::file_path_sans_ext(basename(input_file))
    output_dir <- dirname(input_file)
    output_file <- file.path(output_dir, paste0(input_basename, output_suffix, ".tif"))

    cat("\n====== Verarbeite Bild", i, "von", length(spectral_files), "======\n")
    cat("Eingabe:", input_file, "\n")
    cat("Ausgabe:", output_file, "\n")

    # Prüfe, ob Ausgabedatei bereits existiert
    if (file.exists(output_file)) {
      if (overwrite) {
        cat("WARNUNG: Ausgabedatei existiert bereits und wird gelöscht:", output_file, "\n")
        file.remove(output_file)
      } else {
        # NEW: Always overwrite aatt_destriped files since we want to replace them with scientific method
        if (grepl("aatt_destriped", output_suffix)) {
          cat("INFO: AATT-destriped Datei existiert bereits und wird mit wissenschaftlicher Methode ersetzt:", basename(output_file), "\n")
          file.remove(output_file)
        } else {
          cat("Ausgabedatei existiert bereits. Wird übersprungen:", basename(output_file), "\n")
          cat("(Verwende overwrite=TRUE zum Überschreiben)\n")
          next
        }
      }
    }

    # Lade das Bild
    cat("Lade Spektralbild...\n")
    img <- tryCatch({
      terra::rast(input_file)
    }, error = function(e) {
      cat("FEHLER beim Laden des Bildes:", e$message, "\n")
      return(NULL)
    })
    if (is.null(img)) {
      cat("Überspringe fehlerhafte Datei.\n")
      next
    }

    # Im Entwicklungsmodus Auflösung reduzieren
    if (dev_mode && dev_factor > 1) {
      cat("ENTWICKLUNGSMODUS: Reduziere Auflösung um Faktor", dev_factor, "\n")
      img <- downsample_raster(img, dev_factor)
    }
    
    # Extrahiere den Dateinamen für die Metadaten-Suche
    file_basename <- basename(input_file)
    
    # Suche den Sonnenhöhenwinkel in den Metadaten
    sun_elevation_angle <- NA
    metadata_match <- metadata[metadata$SpectralFile == file_basename, ]
    
    if (nrow(metadata_match) > 0) {
      sun_elevation_angle <- metadata_match$sun_elevation[1]
      cat("Gefundener Sonnenhöhenwinkel:", sun_elevation_angle, "Grad\n")
    } else {
      cat("WARNUNG: Kein Sonnenhöhenwinkel für", file_basename, "gefunden. Überspringe Korrektur.\n")
    }
    
    # # Wende Sonnenhöhenwinkel-Korrektur an, wenn verfügbar
    # if (!is.na(sun_elevation_angle)) {
    #   cat("Wende Sonnenhöhenwinkel-Korrektur an...\n")
    #   img <- correct_sun_elevation(img, sun_elevation_angle)
    # }
    
    # AATT-Destriping nach Hu et al. (2021) mit bewährten Standardparametern
    if (apply_destriping) {
      cat("Wende AATT-Destriping an (Hu et al. 2021 - Adaptive Anisotropic Total Variation)...\n")
      destriping_result <- destripe_raster_aatt(
        img,
        max_iter = 30,      # Standard: 30 Iterationen für ausreichende Konvergenz
        alpha = 0.1,        # Standard: ADMM-Strafparameter für stabile Konvergenz
        beta_x = 0.1,       # Standard: Gradientenschwellwert horizontal
        beta_z = 0.1,       # Standard: Gradientenschwellwert spektral
        verbose = TRUE,
        preserve_stats = TRUE
      )
      img <- destriping_result$raster
      
      # Speichere AATT-Destriping-Metriken
      if (!is.null(destriping_result$metrics)) {
        metrics_file <- file.path(output_dir, paste0(input_basename, "_aatt_destriping_metrics.csv"))
        if (!dir.exists(output_dir)) {
          dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        }
        write.csv(destriping_result$metrics, metrics_file, row.names = FALSE)
        cat("AATT-Destriping-Metriken gespeichert in:", metrics_file, "\n")
        cat("Methode: Adaptive Anisotropic Total Variation + Truncated Nuclear Norm (Hu et al. 2021)\n")
      }
    }

    # Wende Savitzky-Golay-Filter mit Chunking an
    smoothed_img_result <- tryCatch({
      apply_savitzky_golay(img, window_size, poly_order, chunk_size)
    }, error = function(e) {
      cat("FEHLER während apply_savitzky_golay:", conditionMessage(e), "\n")
      return(NULL)
    })

    # Check if smoothing was successful before writing
    if (is.null(smoothed_img_result) || !is.list(smoothed_img_result) ||
        is.null(smoothed_img_result$raster) || !inherits(smoothed_img_result$raster, "SpatRaster")) {
      cat("FEHLER: Savitzky-Golay-Filterung fehlgeschlagen oder abgebrochen für Bild", i, ". Überspringe Speichern.\n")
      next
    }
    smoothed_img <- smoothed_img_result$raster
    temp_file <- smoothed_img_result$temp_file

    # Speichere das geglättete Bild
    cat("Speichere verbessertes geglättetes und destriped Bild...\n")
    write_success <- write_raster_with_fallbacks(smoothed_img, input_file, output_suffix, temp_file = temp_file)
    if (write_success) {
    cat("✓ AATT-destriped und SG-geglättetes Bild gespeichert:", basename(output_file), "\n")
    } else {
      cat("FEHLER: Konnte Bild nicht speichern. Bitte prüfen Sie Pfad und Rechte.\n")
    }

    # Speicherbereinigung
    img <- smoothed_img <- smoothed_img_result <- NULL
    memory_cleanup()

    # Fortschritt speichern
    writeLines(as.character(i), "last_aatt_destriped_image_index.txt")
    cat("Fortschritt:", i, "von", length(spectral_files), "Bildern verarbeitet\n")
  }
  cat("\nAATT-Destriping-Verarbeitung abgeschlossen (Hu et al. 2021). Alle", length(spectral_files) - start_index + 1, "Bilder wurden verarbeitet.\n")
  cat("Neue Dateien mit AATT-Destriping-Methode wurden mit Suffix '", output_suffix, "' erstellt.\n")
  cat("Methode: Adaptive Anisotropic Total Variation + Truncated Nuclear Norm\n")
  cat("Parameter: ADMM α=0.1, Gradientenschwellwerte βx=βz=0.1, max. 30 Iterationen\n")
  # Temporäre Dateien erst jetzt löschen!
  terra::tmpFiles(remove=TRUE)
  # Optional: altes tempdir wiederherstellen
  # Sys.setenv(TMPDIR = old_temp_dir)
  # terra::terraOptions(tempdir = old_temp_dir)
}

# Starte die Verarbeitung, falls Skript direkt ausgeführt wird
last_index <- 1
# Progress file for AATT destriping (Hu et al. 2021)
if (file.exists("last_aatt_destriped_image_index.txt")) {
  last_index <- as.integer(readLines("last_aatt_destriped_image_index.txt")) + 1
  cat("Fortsetzen der AATT-Destriping-Verarbeitung (Hu et al. 2021) ab Bild", last_index, "\n")
}

# Schritt 1: CSV einlesen und RasterFile-Liste extrahieren
library(data.table)
csv_path <- "f:/Uni/Bachelorarbeit/R working directory/enmap-soil-correlation/data_Sac/final_before_correlation.csv"
df <- fread(csv_path)
raster_files_to_process <- unique(df$RasterFile)

# Schritt 2: Nur Raster aus dieser Liste bearbeiten
# Beispiel: Wenn du eine Schleife über Raster hast, dann so filtern:
# for (raster_file in all_raster_files) {
#   if (!(raster_file %in% raster_files_to_process)) next
#   ... smoothing code ...
# }

smooth_spectral_images(
  data_file = "data/final_before_correlation.csv",
  spectral_dir = "data/extracted_files/ich_everzwlfe",
  output_suffix = "_smoothed+aat_destriped",
  window_size = 11,
  poly_order = 3,
  chunk_size = 500,
  dev_mode = FALSE,
  dev_factor = 4,
  overwrite = TRUE,  # Set to TRUE to ensure replacement of existing files
  start_index = last_index,
  metadata_file = "data/metadata_mapping.csv",
  apply_destriping = TRUE,
  destriping_adaptive = TRUE,
  destriping_edge_tapering = TRUE,
  destriping_visualize = FALSE,
  filter_csv_file = "data_Sac/final_before_correlation.csv"  # Use the filter CSV
)
