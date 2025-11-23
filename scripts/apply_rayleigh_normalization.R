# Script zur Anwendung der Rayleigh-Normalisierung auf maskierte Hyperspektralbilder
# Anforderung: spectral_normalization.R und terra Paket müssen installiert sein

library(terra)
library(tools)
library(pbapply)

# Optimiere Speichernutzung für terra
optimize_terra_memory <- function() {
  # Ermittelt verfügbaren RAM und konfiguriert terra entsprechend
  max_mem_percent <- 70  # Prozent des verfügbaren RAMs für terra
  
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
  
  # Berechne max. Speicher in Bytes (70% des RAM)
  max_mem <- floor(total_mem * max_mem_percent / 100) * 1024
  
  # Konfiguriere terra
  terra::terraOptions(memfrac = max_mem_percent / 100,
                     memmax = max_mem,
                     progress = 10,
                     verbose = TRUE,
                     todisk = TRUE,
                     memmin = 0.1 * max_mem,
                     datatype = "FLT4S")
  
  # GDAL-Cache optimieren
  Sys.setenv(GDAL_CACHEMAX = "512")
  Sys.setenv(GDAL_DISABLE_READDIR_ON_OPEN = "TRUE")
  Sys.setenv(VSI_CACHE = "TRUE")
  Sys.setenv(VSI_CACHE_SIZE = "50000000")  # 50MB Cache
  
  cat("Terra für Speichereffizienz optimiert:", round(max_mem/1e9, 1), "GB verfügbar\n")
}

# Anwenden der Rayleigh-Normalisierung auf maskierte Bilder
apply_rayleigh_normalization <- function(
  masked_dir = "data/extracted_files/ich_everzwlfe",  # Verzeichnis mit maskierten Bildern
  input_suffix = "_smoothed\\+destripedx_masked",     # Suffix der zu normalisierenden Dateien
  output_suffix = "_smoothed+destripedx_masked+rayleigh",  # Suffix der Ausgabedateien
  wavelength_file = "data/wavelength_reference.csv",  # Pfad zur Wellenlängenreferenzdatei
  overwrite = FALSE,                                  # Überschreibe bestehende Dateien?
  max_batch_size = 10000,                             # Max. Anzahl an Pixeln pro Batch für Normalisierung
  start_index = 1                                     # Startindex für die Verarbeitung
) {
  cat("\n=== ANWENDUNG DER RAYLEIGH-NORMALISIERUNG AUF MASKIERTE BILDER ===\n\n")
  
  # Lade die spectral_normalization Funktionen
  if (!file.exists("scripts/spectral_normalization.R")) {
    stop("spectral_normalization.R nicht gefunden. Bitte erst dieses Script erstellen.")
  }
  source("scripts/spectral_normalization.R")
  
  # Lade Wellenlängenreferenzdatei
  if (!file.exists(wavelength_file)) {
    stop("Wellenlängenreferenzdatei nicht gefunden: ", wavelength_file)
  }
  
  wavelength_ref <- read.csv(wavelength_file)
  if (!"wavelength_nm" %in% colnames(wavelength_ref)) {
    stop("Wellenlängenreferenzdatei muss eine Spalte 'wavelength_nm' enthalten")
  }
  
  wavelengths <- wavelength_ref$wavelength_nm
  cat("✓ Wellenlängenreferenz geladen:", length(wavelengths), "Bänder\n")
  
  # Finde alle maskierten Bilder
  masked_files <- list.files(
    masked_dir,
    pattern = paste0(input_suffix, "\\.tif$"),
    full.names = TRUE,
    recursive = TRUE
  )
  
  if (length(masked_files) == 0) {
    cat("❌ Keine maskierten Bilder gefunden in:", masked_dir, "\n")
    cat("   Gesuchtes Pattern:", paste0(input_suffix, "\\.tif$"), "\n")
    return(invisible(NULL))
  }
  
  cat("✓ Gefundene maskierte Bilder:", length(masked_files), "\n")
  
  # Prüfe, ob der Startindex gültig ist
  if (start_index < 1 || start_index > length(masked_files)) {
    cat("⚠️ Ungültiger Startindex:", start_index, "\n")
    cat("   Der Index wurde auf 1 zurückgesetzt.\n")
    start_index <- 1
  } else if (start_index > 1) {
    cat("ℹ️ Starte Verarbeitung ab Bild", start_index, "von", length(masked_files), "\n")
  }
  
  # Erstelle einen Fortschrittbalken
  pb <- txtProgressBar(min = 0, max = length(masked_files) - start_index + 1, style = 3)
  progress_counter <- 0
  
  # Verarbeite jedes maskierte Bild
  for (i in start_index:length(masked_files)) {
    img_file <- masked_files[i]
    
    # Generiere Ausgabedateiname
    output_file <- gsub(input_suffix, output_suffix, img_file, fixed = FALSE)
    
    # Aktualisiere Fortschrittsbalken und zeige Details
    progress_counter <- progress_counter + 1
    setTxtProgressBar(pb, progress_counter)
    
    cat("\n\n[", i, "/", length(masked_files), "] Verarbeite:", basename(img_file), "\n")
    cat("   Ausgabe:", basename(output_file), "\n")
    
    # Prüfe, ob Ausgabedatei bereits existiert
    if (file.exists(output_file) && !overwrite) {
      cat("   ⏩ Ausgabedatei existiert bereits. Überspringe. (Verwende overwrite=TRUE zum Überschreiben)\n")
      next
    }
    
    # Lade das Bild
    img <- NULL
    tryCatch({
      img <- rast(img_file)
      cat("   ✓ Bild geladen:", nlyr(img), "Bänder,", nrow(img), "×", ncol(img), "Pixel\n")
    }, error = function(e) {
      cat("   ❌ Fehler beim Laden:", e$message, "\n")
    })
    
    if (is.null(img)) {
      cat("   ⚠️ Überspringe Bild wegen Ladefehler\n")
      next
    }
    
    # Prüfe Anzahl der Bänder gegen Wellenlängen
    if (nlyr(img) != length(wavelengths)) {
      cat("   ⚠️ Anzahl der Bänder (", nlyr(img), ") stimmt nicht mit Wellenlängenreferenz (", 
          length(wavelengths), ") überein\n")
      
      # Versuche anzupassen
      if (nlyr(img) < length(wavelengths)) {
        cat("   ⚙️ Kürze Wellenlängenliste auf", nlyr(img), "Einträge\n")
        wavelengths <- wavelengths[1:nlyr(img)]
      } else {
        cat("   ❌ Mehr Bänder als Wellenlängen, Bild wird übersprungen!\n")
        next
      }
    }
    
    # Starte Timer für die Verarbeitungszeit
    start_time <- Sys.time()
    
    # Verarbeite das Bild in Batches
    cat("   ⚙️ Wende Rayleigh-Normalisierung an...\n")
    
    # Hole Pixelmatrix und bestimme Anzahl gültiger Pixel
    pix_matrix <- values(img)
    valid_pixels <- sum(!rowSums(is.na(pix_matrix)) == ncol(pix_matrix))
    cat("   ℹ️ Gültige Pixel:", valid_pixels, "von", nrow(pix_matrix), 
        "(", round(valid_pixels/nrow(pix_matrix)*100, 1), "%)\n")
    
    if (valid_pixels == 0) {
      cat("   ⚠️ Keine gültigen Pixel für Normalisierung. Überspringe.\n")
      next
    }
    
    # Initialisiere normalisierte Matrix
    norm_matrix <- matrix(NA, nrow = nrow(pix_matrix), ncol = ncol(pix_matrix))
    
    # Verarbeite Daten in Batches
    batch_size <- min(max_batch_size, 10000)  # Anpassbare Batch-Größe
    num_batches <- ceiling(nrow(pix_matrix) / batch_size)
    
    cat("   ⚙️ Verarbeite in", num_batches, "Batches mit je bis zu", batch_size, "Pixeln\n")
    
    # Inner progress bar for batches
    batch_pb <- txtProgressBar(min = 0, max = num_batches, style = 3)
    
    for (b in 1:num_batches) {
      # Bestimme Batch-Grenzen
      start_idx <- (b - 1) * batch_size + 1
      end_idx <- min(b * batch_size, nrow(pix_matrix))
      
      # Extrahiere aktuellen Batch
      batch_data <- pix_matrix[start_idx:end_idx, , drop = FALSE]
      
      # Finde gültige Zeilen (nicht komplett NA)
      valid_indices <- which(!rowSums(is.na(batch_data)) == ncol(batch_data))
      
      if (length(valid_indices) > 0) {
        # Extrahiere nur die gültigen Pixel für die Verarbeitung
        valid_batch <- batch_data[valid_indices, , drop = FALSE]
        
        # Anwendung der Rayleigh-Normalisierung
        # Das zweite Argument ist die Wellenlänge für den Rayleigh-Faktor
        normalized_batch <- tryCatch({
          rayleigh_normalize(valid_batch, wavelengths, 
                            fit_range = 1:ceiling(length(wavelengths) * 0.1),
                            subtraction_mode = "divide")
        }, error = function(e) {
          cat("\n   ❌ Fehler bei Rayleigh-Normalisierung (Batch", b, "):", e$message, "\n")
          return(NULL)
        })
        
        # Patch: Wenn Input gültig, aber Output NA, dann setze Output = Input
        if (!is.null(normalized_batch)) {
          for (j in seq_len(nrow(valid_batch))) {
            input_row <- valid_batch[j, ]
            output_row <- normalized_batch[j, ]
            if (all(!is.na(input_row)) && any(is.na(output_row))) {
              normalized_batch[j, is.na(output_row)] <- input_row[is.na(output_row)]
            }
          }
        }
        
        if (!is.null(normalized_batch)) {
          # Setze normalisierte Werte zurück in die Gesamtmatrix
          norm_matrix[start_idx + valid_indices - 1, ] <- normalized_batch
        }
      }
      
      # Update batch progress bar
      setTxtProgressBar(batch_pb, b)
    }
    
    close(batch_pb)
    
    # Erstelle neues Raster mit normalisierten Werten
    normalized_raster <- rast(img)  # Behalte Geo-Informationen
    values(normalized_raster) <- norm_matrix
    
    # Speichere normalisiertes Raster
    cat("   ⚙️ Speichere normalisiertes Raster...\n")
    
    tryCatch({
      writeRaster(normalized_raster, output_file, overwrite = TRUE, 
                 gdal = c("COMPRESS=LZW", "TILED=YES"))
      
      end_time <- Sys.time()
      processing_time <- difftime(end_time, start_time, units = "mins")
      
      cat("   ✅ Rayleigh-normalisiertes Bild gespeichert:", basename(output_file), "\n")
      cat("   ⏱️ Verarbeitungszeit:", round(processing_time, 2), "Minuten\n")
      
      # Speichere Fortschritt
      writeLines(as.character(i), "last_rayleigh_normalized_image_index.txt")
      
    }, error = function(e) {
      cat("   ❌ Fehler beim Speichern:", e$message, "\n")
    })
    
    # Speicherbereinigung
    rm(img, pix_matrix, norm_matrix, normalized_raster)
    invisible(gc(full = TRUE))
  }
  
  close(pb)
  cat("\n\n✅ Rayleigh-Normalisierung auf", progress_counter, "Bilder angewendet\n")
  cat("   Ausgabedateien wurden mit Suffix '", output_suffix, "' gespeichert\n")
}

# Funktion zum Erstellen eines Wellenlängenreferenzfiles, falls noch nicht vorhanden
create_wavelength_reference <- function(
  output_file = "data/wavelength_reference.csv",
  num_bands = 224  # EnMAP Standard
) {
  if (file.exists(output_file)) {
    cat("✓ Wellenlängenreferenzdatei existiert bereits:", output_file, "\n")
    return(invisible(NULL))
  }
  
  cat("⚙️ Erstelle Wellenlängenreferenzdatei...\n")
  
  # EnMAP Wellenlängen (approximiert)
  # VNIR: ~420-1000 nm (Band 1-88)
  # SWIR: ~900-2450 nm (Band 89-224)
  
  vnir_start <- 420
  vnir_end <- 1000
  swir_start <- 900
  swir_end <- 2450
  
  vnir_bands <- 88
  swir_bands <- 136
  
  # Erzeuge gleichmäßig verteilte Wellenlängen
  vnir_wavelengths <- seq(vnir_start, vnir_end, length.out = vnir_bands)
  swir_wavelengths <- seq(swir_start, swir_end, length.out = swir_bands)
  
  # Kombiniere zu einer Referenz
  wavelength_ref <- data.frame(
    band_number = 1:num_bands,
    band_name = paste0("Band", 1:num_bands),
    wavelength_nm = c(vnir_wavelengths, swir_wavelengths)[1:num_bands]
  )
  
  # Speichere als CSV
  write.csv(wavelength_ref, output_file, row.names = FALSE)
  
  cat("✅ Wellenlängenreferenz erstellt und gespeichert:", output_file, "\n")
  cat("   Enthält", nrow(wavelength_ref), "Bänder\n")
  
  return(invisible(wavelength_ref))
}

# Hauptprogramm
main <- function() {
  # Optimiere terra für Speichereffizienz
  optimize_terra_memory()
  
  # Erstelle/prüfe Wellenlängenreferenz
  wavelength_file <- "data/wavelength_reference.csv"
  create_wavelength_reference(wavelength_file)
  
  # Lese letzten Verarbeitungsindex, falls vorhanden
  start_index <- 1
  if (file.exists("last_rayleigh_normalized_image_index.txt")) {
    last_index <- as.integer(readLines("last_rayleigh_normalized_image_index.txt"))
    if (!is.na(last_index) && last_index > 0) {
      start_index <- last_index + 1
      cat("ℹ️ Fortsetzen ab Bild", start_index, "(letzter verarbeiteter Index:", last_index, ")\n")
    }
  }
  
  # Anwenden der Rayleigh-Normalisierung auf maskierte Bilder
  apply_rayleigh_normalization(
    masked_dir = "data/extracted_files/ich_everzwlfe",
    input_suffix = "_smoothed\\+destripedx_masked",
    output_suffix = "_smoothed+destripedx_masked+rayleigh",
    wavelength_file = wavelength_file,
    overwrite = TRUE,  # Überschreibe bestehende normalisierte Dateien
    max_batch_size = 10000,  # Anpassen je nach verfügbarem Arbeitsspeicher
    start_index = start_index
  )
}

# HINWEIS: Auto-Execution auskommentiert - Script wird nur über master_pipeline.R gesteuert
# # Starte das Hauptprogramm
# main()
