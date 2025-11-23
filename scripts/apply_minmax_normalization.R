# Min/Max-Normalisierung für alle ENMAP Raster im data/-Ordner

library(terra)

# --- Progress-Tracking-System ---
progress_file <- "minmax_normalization_progress.txt"

get_last_processed_index <- function() {
  if (file.exists(progress_file)) {
    last_index <- as.integer(readLines(progress_file, warn = FALSE))
    return(if (is.na(last_index)) 0 else last_index)
  }
  return(0)
}

save_progress <- function(current_index) {
  writeLines(as.character(current_index), progress_file)
}

# --- Hilfsfunktion für Min/Max-Normalisierung eines Vektors ---
minmax_norm <- function(x) {
  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)
  rng <- max_x - min_x
  if (is.na(rng) || rng == 0) {
    out <- rep(0, length(x))
    out[is.na(x)] <- NA
    return(out)
  }
  out <- (x - min_x) / rng
  out[is.na(x)] <- NA
  return(out)
}

# --- Suche nach geglätteten, entstreiften und maskierten ENMAP-Bildern ---
input_dir <- "data/extracted_files/ich_everzwlfe"  # Korrektes Verzeichnis
output_suffix <- "_minmax.tif"

# Pattern für die neuen maskierten Bilder
raster_files <- list.files(
  input_dir,
  pattern = "_smoothed\\+adaptive_destriped_masked_snv\\.tif$",
  full.names = TRUE,
  recursive = TRUE
)

cat("Gefundene smoothed+adaptive_destriped_masked_snv Raster:", length(raster_files), "\n")

if (length(raster_files) == 0) {
  cat("Keine maskierten Raster gefunden. Suche in alternativen Verzeichnissen...\n")
  
  # Suche auch im Hauptdatenverzeichnis
  raster_files_alt <- list.files(
    "data",
    pattern = "_smoothed\\+adaptive_destriped_masked_snv\\.tif$",
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (length(raster_files_alt) > 0) {
    raster_files <- raster_files_alt
    cat("Gefunden in alternativen Verzeichnissen:", length(raster_files), "\n")
  } else {
    stop("Keine maskierten Raster-Dateien gefunden! Prüfen Sie, ob die Maskierung abgeschlossen ist.")
  }
}

# Progress-Tracking: Hole letzten verarbeiteten Index
last_processed <- get_last_processed_index()
start_index <- last_processed + 1

if (start_index > length(raster_files)) {
  cat("Alle Dateien bereits verarbeitet!\n")
  cat("Zum Neustart löschen Sie die Datei:", progress_file, "\n")
  quit()
}

if (last_processed > 0) {
  cat("Fortsetzen der Verarbeitung ab Datei", start_index, "von", length(raster_files), "\n")
  cat("Bereits verarbeitet:", last_processed, "Dateien\n")
} else {
  cat("Beginne neue Verarbeitung mit", length(raster_files), "Dateien\n")
}

# Verarbeite jedes gefundene Raster ab dem letzten Index
for (i in start_index:length(raster_files)) {
  input_raster <- raster_files[i]
  cat("\nBearbeite [", i, "/", length(raster_files), "]:", basename(input_raster), "\n")
  
  # Prüfe, ob Ausgabedatei bereits existiert
  # GEÄNDERT: Verwende das korrekte Suffix für spatial_prediction.R
  output_raster <- gsub("_smoothed\\+adaptive_destriped_masked_snv\\.tif$", 
                        "_smoothed+adaptive_destriped_masked_snv_minmax.tif", 
                        input_raster)
  
  if (file.exists(output_raster) && file.size(output_raster) > 0) {
    cat("  -> Ausgabedatei existiert bereits, überspringe\n")
    save_progress(i)
    next
  }
  
  # Lade das Raster
  r <- tryCatch({
    rast(input_raster)
  }, error = function(e) {
    cat("  FEHLER beim Laden:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(r)) {
    cat("  Überspringe Datei aufgrund von Ladefehlern.\n")
    save_progress(i)
    next
  }
  
  n_bands <- nlyr(r)
  cat("  Verarbeite", n_bands, "Bänder...\n")
  
  # Erstelle Kopie für normalisierte Werte
  r_norm <- r
  
  # Verarbeite jedes Band einzeln
  valid_bands <- 0
  for (band_i in 1:n_bands) {
    tryCatch({
      vals <- values(r[[band_i]], mat = FALSE)
      
      # Prüfe ob genügend gültige Werte vorhanden sind
      valid_count <- sum(!is.na(vals))
      if (valid_count < 10) {
        cat("    Band", band_i, ": Zu wenige gültige Werte (", valid_count, "), überspringe\n")
        next
      }
      
      # Wende Min-Max-Normalisierung an
      vals_norm <- minmax_norm(vals)
      values(r_norm[[band_i]]) <- vals_norm
      valid_bands <- valid_bands + 1
      
      if (band_i %% 50 == 0) {
        cat("    Verarbeitet:", band_i, "von", n_bands, "Bändern\n")
      }
      
    }, error = function(e) {
      cat("    FEHLER bei Band", band_i, ":", e$message, "\n")
    })
  }
  
  cat("  Erfolgreich normalisierte Bänder:", valid_bands, "von", n_bands, "\n")
  
  if (valid_bands == 0) {
    cat("  Keine Bänder erfolgreich normalisiert, überspringe Speichern.\n")
    save_progress(i)
    next
  }
  
  # Erstelle Ausgabedatei-Namen
  # ...existing code...
  
  # Speichere das normalisierte Raster
  tryCatch({
    writeRaster(r_norm, output_raster, overwrite = TRUE, 
               gdal = c("COMPRESS=LZW", "TILED=YES"))
    cat("  -> Erfolgreich gespeichert als:", basename(output_raster), "\n")
    
    # Validierung: Prüfe ob Datei korrekt geschrieben wurde
    if (file.exists(output_raster) && file.size(output_raster) > 0) {
      cat("  -> Dateivalidierung: OK (", round(file.size(output_raster)/1e6, 1), "MB)\n")
    } else {
      cat("  -> WARNUNG: Ausgabedatei scheint beschädigt zu sein!\n")
    }
    
  }, error = function(e) {
    cat("  FEHLER beim Speichern:", e$message, "\n")
  })
  
  # Speicher freigeben
  rm(r, r_norm)
  gc()
  
  # Speichere Progress nach jeder erfolgreich verarbeiteten Datei
  save_progress(i)
}

cat("\nMin-Max-Normalisierung abgeschlossen.\n")
cat("Alle normalisierten Dateien haben das Suffix '_smoothed+adaptive_destriped_masked_snv_minmax.tif'\n")
cat("Progress-Datei:", progress_file, "kann gelöscht werden.\n")
