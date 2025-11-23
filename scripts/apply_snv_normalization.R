# Script zur Anwendung von SNV und Kontinuum-Removal Normalisierung auf maskierte Hyperspektralbilder
# Alternative zur Rayleigh-Normalisierung

# Überprüfen und installieren benötigter Pakete
required_packages <- c("terra", "tools", "pbapply", "matrixStats")
for(pkg in required_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    cat(paste0("Installiere benötigtes Paket: ", pkg, "\n"))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Pakete laden mit expliziter Fehlermeldung
load_package <- function(package_name) {
  tryCatch({
    library(package_name, character.only = TRUE)
    cat(paste0("✓ Paket geladen: ", package_name, "\n"))
  }, error = function(e) {
    stop(paste0("Fehler beim Laden des Pakets ", package_name, ": ", e$message))
  })
}

# Lade alle benötigten Pakete
for(pkg in required_packages) {
  load_package(pkg)
}

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

# Standard Normal Variate (SNV) Normalisierung
# Diese Methode zentriert und skaliert jedes Spektrum individuell
snv_normalize <- function(spectra_matrix) {
  # Berechne Zeilenmittelwerte und Standardabweichungen
  row_means <- rowMeans(spectra_matrix, na.rm = TRUE)
  row_sds <- matrixStats::rowSds(spectra_matrix, na.rm = TRUE)
  
  # Ersetze Nullwerte in den Standardabweichungen, um Division durch Null zu vermeiden
  row_sds[row_sds == 0 | is.na(row_sds)] <- 1
  
  # Normalisiere jede Zeile: (x - mean) / sd
  normalized <- matrix(NA, nrow = nrow(spectra_matrix), ncol = ncol(spectra_matrix))
  for (i in seq_len(nrow(spectra_matrix))) {
    normalized[i,] <- (spectra_matrix[i,] - row_means[i]) / row_sds[i]
  }
  
  normalized
}

# Kontinuumsremoval-Normalisierung
# Normalisiert gegen eine Hüllkurve (Kontinuum) über die Spektralpeaks
continuum_removal <- function(spectra_matrix, wavelengths) {
  normalized <- matrix(NA, nrow = nrow(spectra_matrix), ncol = ncol(spectra_matrix))
  
  for (i in seq_len(nrow(spectra_matrix))) {
    spectrum <- spectra_matrix[i,]
    valid_idx <- !is.na(spectrum)
    
    if (sum(valid_idx) > 2) {
      # Erstelle die konvexe Hülle über dem Spektrum
      valid_spectrum <- spectrum[valid_idx]
      valid_wavelengths <- wavelengths[valid_idx]
      
      # Finde lokale Maxima (vereinfachte Version)
      hull_points <- c(1, which(diff(sign(diff(c(-Inf, valid_spectrum, -Inf)))) == -2))
      
      if (length(hull_points) > 1) {
        # Interpoliere die Hüllkurve
        continuum_values <- approx(valid_wavelengths[hull_points], 
                                  valid_spectrum[hull_points], 
                                  xout = valid_wavelengths, 
                                  method = "linear")$y
        
        # Dividiere das Spektrum durch die Hüllkurve
        normalized_spectrum <- valid_spectrum / continuum_values
        
        # Setze die normalisierten Werte zurück
        normalized[i, valid_idx] <- normalized_spectrum
      } else {
        # Fallback bei Problem mit der Hülle
        normalized[i, valid_idx] <- valid_spectrum
      }
    } else {
      # Nicht genügend gültige Datenpunkte
      normalized[i, valid_idx] <- spectrum[valid_idx]
    }
  }
  
  normalized
}

# Anwenden der SNV oder Kontinuum-Removal Normalisierung auf maskierte Bilder
apply_spectral_normalization <- function(
  masked_dir = "data/extracted_files/ich_everzwlfe",  # Verzeichnis mit maskierten Bildern
  input_suffix = "_smoothed\\+destripedx_masked",     # Suffix der zu normalisierenden Dateien
  output_suffix = "_normalized",                      # Suffix der Ausgabedateien
  method = "snv",                                    # Normalisierungsmethode: "snv" oder "continuum"
  wavelength_file = "data/wavelength_reference.csv",  # Pfad zur Wellenlängenreferenzdatei (für continuum)
  overwrite = TRUE,                                 # Überschreibe bestehende Dateien?
  max_batch_size = 10000,                            # Max. Anzahl an Pixeln pro Batch für Normalisierung
  start_index = 1                                    # Startindex für die Verarbeitung
) {
  # Bestimme das korrekte Ausgabesuffix basierend auf der Methode
  method_suffix <- switch(method,
                         "snv" = "_snv",
                         "continuum" = "_cont",
                         "_norm") # Fallback
  
  output_suffix <- paste0(output_suffix, method_suffix)
  
  cat(sprintf("\n=== ANWENDUNG DER %s-NORMALISIERUNG AUF MASKIERTE BILDER ===\n\n", toupper(method)))
  
  # Lade Wellenlängenreferenzdatei (vor allem für Kontinuum-Methode wichtig)
  if (method == "continuum" && !file.exists(wavelength_file)) {
    stop("Wellenlängenreferenzdatei nicht gefunden: ", wavelength_file)
  }
  
  wavelength_ref <- NULL
  wavelengths <- NULL
  
  if (file.exists(wavelength_file)) {
    wavelength_ref <- read.csv(wavelength_file)
    if (!"wavelength_nm" %in% colnames(wavelength_ref)) {
      warning("Wellenlängenreferenzdatei enthält keine Spalte 'wavelength_nm', generische Wellenlängen werden verwendet")
    } else {
      wavelengths <- wavelength_ref$wavelength_nm
      cat("✓ Wellenlängenreferenz geladen:", length(wavelengths), "Bänder\n")
    }
  }
  
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
    cat("   Methode:", method, "\n")
    
    # Prüfe, ob Ausgabedatei bereits existiert
    if (file.exists(output_file) && !overwrite) {
      cat("   ⏩ Ausgabedatei existiert bereits. Überspringe. (Verwende overwrite=TRUE zum Überschreiben)\n")
      next
    }
    
    # Lade das Bild
    img <- NULL
    tryCatch({
      # Verwende expliziten Namensraum für rast-Funktion
      img <- terra::rast(img_file)
      cat("   ✓ Bild geladen:", terra::nlyr(img), "Bänder,", terra::nrow(img), "×", terra::ncol(img), "Pixel\n")
    }, error = function(e) {
      cat("   ❌ Fehler beim Laden:", e$message, "\n")
      # Versuche zusätzliche Diagnoseinformationen
      cat("   🔍 Terra-Paket geladen:", isNamespaceLoaded("terra"), "\n")
      if(isNamespaceLoaded("terra")) {
        cat("   🔍 Terra-Version:", packageVersion("terra"), "\n")
      }
    })
    
    if (is.null(img)) {
      cat("   ⚠️ Überspringe Bild wegen Ladefehler\n")
      next
    }
    
    # Prüfe Anzahl der Bänder gegen Wellenlängen (nur für Kontinuum-Methode wichtig)
    if (!is.null(wavelengths) && terra::nlyr(img) != length(wavelengths)) {
      cat("   ⚠️ Anzahl der Bänder (", terra::nlyr(img), ") stimmt nicht mit Wellenlängenreferenz (", 
          length(wavelengths), ") überein\n")
      
      # Versuche anzupassen
      if (terra::nlyr(img) < length(wavelengths)) {
        cat("   ⚙️ Kürze Wellenlängenliste auf", terra::nlyr(img), "Einträge\n")
        wavelengths <- wavelengths[seq_len(terra::nlyr(img))]
      } else if (terra::nlyr(img) > length(wavelengths) && method == "continuum") {
        cat("   ⚙️ Erstelle synthetische Wellenlängenliste\n")
        # Erstelle synthetische Wellenlängenliste bei fehlenden Daten
        wavelengths <- seq(400, 2500, length.out = terra::nlyr(img))
      }
    } else if (is.null(wavelengths) && method == "continuum") {
      # Erstelle synthetische Wellenlängenliste für Kontinuum-Methode
      cat("   ⚙️ Erstelle synthetische Wellenlängenliste\n")
      wavelengths <- seq(400, 2500, length.out = terra::nlyr(img))
    }
    
    # Starte Timer für die Verarbeitungszeit
    start_time <- Sys.time()
    
    # Verarbeite das Bild in Batches
    cat(sprintf("   ⚙️ Wende %s-Normalisierung an...\n", toupper(method)))
    
    # Hole Pixelmatrix und bestimme Anzahl gültiger Pixel
    pix_matrix <- terra::values(img)
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
        
        # Anwendung der gewählten Normalisierungsmethode
        normalized_batch <- tryCatch({
          switch(method,
                "snv" = snv_normalize(valid_batch),
                "continuum" = continuum_removal(valid_batch, wavelengths),
                # Fallback auf SNV, falls unbekannte Methode
                snv_normalize(valid_batch))
        }, error = function(e) {
          cat("\n   ❌ Fehler bei", toupper(method), "-Normalisierung (Batch", b, "):", e$message, "\n")
          return(NULL)
        })
        
        if (!is.null(normalized_batch)) {
          # Setze normalisierte Werte zurück in die Gesamtmatrix
          norm_matrix[start_idx + valid_indices - 1, ] <- normalized_batch
        }
      }
      
      # Update batch progress bar
      setTxtProgressBar(batch_pb, b)
    }
    
    close(batch_pb)
    
    # --- NA-Handling nach der SNV-Normalisierung ---
    # Ziel: Jeder NA-Input bleibt NA im Output, jeder nicht-NA-Input darf kein NA im Output sein
    # 1. Erzeuge eine Maske der Input-NAs
    input_na_mask <- is.na(pix_matrix)
    # 2. Setze Output an denselben Stellen auf NA
    norm_matrix[input_na_mask] <- NA
    # 3. Alle verbleibenden NAs im Output (wo Input nicht NA war) auf 0 setzen
    non_na_input_mask <- !input_na_mask
    na_in_output <- is.na(norm_matrix) & non_na_input_mask
    if (any(na_in_output)) {
      norm_matrix[na_in_output] <- 0
      cat("   ⚠️ NA-Werte im Output (bei nicht-NA Input) wurden auf 0 gesetzt.\n")
    }
    
    # Erstelle neues Raster mit normalisierten Werten
    normalized_raster <- terra::rast(img)  # Behalte Geo-Informationen
    terra::values(normalized_raster) <- norm_matrix
    
    # Speichere normalisiertes Raster
    cat("   ⚙️ Speichere normalisiertes Raster...\n")
    
    tryCatch({
      terra::writeRaster(normalized_raster, output_file, overwrite = TRUE, 
                      gdal = c("COMPRESS=LZW", "TILED=YES"))
      
      end_time <- Sys.time()
      processing_time <- difftime(end_time, start_time, units = "mins")
      
      cat(sprintf("   ✅ %s-normalisiertes Bild gespeichert:", toupper(method)), basename(output_file), "\n")
      cat("   ⏱️ Verarbeitungszeit:", round(processing_time, 2), "Minuten\n")
      
      # Speichere Fortschritt
      progress_file <- paste0("last_", method, "_normalized_image_index.txt")
      writeLines(as.character(i), progress_file)
      
    }, error = function(e) {
      cat("   ❌ Fehler beim Speichern:", e$message, "\n")
    })
    
    # Speicherbereinigung
    rm(img, pix_matrix, norm_matrix, normalized_raster)
    invisible(gc(full = TRUE))
  }
  
  close(pb)
  cat(sprintf("\n\n✅ %s-Normalisierung auf %d Bilder angewendet\n", toupper(method), progress_counter))
  cat("   Ausgabedateien wurden mit Suffix '", output_suffix, "' gespeichert\n")
}

# Hauptprogramm
main <- function() {
  # Überprüfe ob terra korrekt geladen ist
  if(!isNamespaceLoaded("terra")) {
    cat("⚠️ Terra Paket nicht geladen. Versuche erneut zu laden...\n")
    load_package("terra")
  }
  
  # Optimiere terra für Speichereffizienz
  optimize_terra_memory()
  
  # Erstelle/prüfe Wellenlängenreferenz für Kontinuum-Methode
  wavelength_file <- "data/wavelength_reference.csv"
  
  # Prüfe ob Wellenlängenreferenz existiert
  if (!file.exists(wavelength_file)) {
    cat("ℹ️ Wellenlängenreferenz nicht gefunden.\n")
    cat("  Falls Kontinuum-Normalisierung verwendet wird, bitte Wellenlängenreferenz manuell erstellen.\n")
  }
  
  # Bestimme die zu verwendende Normalisierungsmethode
  # Standardeinstellung: "snv", Alternative: "continuum"
  normalization_method <- "snv"
  
  # Konfigurieren Sie hier die Methode oder fügen Sie einen Parameterdialog hinzu
  cat("ℹ️ Verwendete Normalisierungsmethode:", toupper(normalization_method), "\n")
  
  # Lese letzten Verarbeitungsindex, falls vorhanden
  start_index <- 1
  progress_file <- paste0("last_", normalization_method, "_normalized_image_index.txt")
  if (file.exists(progress_file)) {
    last_index <- as.integer(readLines(progress_file))
    if (!is.na(last_index) && last_index > 0) {
      start_index <- last_index + 1
      cat("ℹ️ Fortsetzen ab Bild", start_index, "(letzter verarbeiteter Index:", last_index, ")\n")
    }
  }
  
  # Anwenden der gewählten Normalisierung auf maskierte Bilder
  apply_spectral_normalization(
    masked_dir = "data/extracted_files/ich_everzwlfe",
    input_suffix = "_smoothed\\+destripedx_masked",
    output_suffix = "_smoothed+destripedx_masked",
    method = normalization_method,
    wavelength_file = wavelength_file,
    overwrite = TRUE,  # Überschreibe bestehende normalisierte Dateien
    max_batch_size = 10000,  # Anpassen je nach verfügbarem Arbeitsspeicher
    start_index = start_index
  )
}

# HINWEIS: Auto-Execution auskommentiert - Script wird nur über master_pipeline.R gesteuert
# # Starte das Hauptprogramm
# main()
