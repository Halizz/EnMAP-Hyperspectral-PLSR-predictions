# Script: extract_metadata.R
# Zweck: Extrahiert Metadaten aus EnMAP XML-Dateien und erstellt metadata_mapping.csv
# PLUS: Extrahiert Wellenlängen-Referenz für Spektralbänder
# Wird benötigt VOR dem Smoothing-Schritt, um Sonnenhöhenwinkel und andere Parameter verfügbar zu machen

library(xml2)
library(sf)
library(terra)
library(data.table)
library(readr)

#' Extrahiert Wellenlängen-Referenz aus EnMAP XML-Datei
#'
#' @param metadata_file Pfad zur ersten gefundenen Metadata-XML-Datei
#' @param output_file Ausgabedatei für wavelength_reference.csv
#' @return Data frame mit Wellenlängen-Referenz
extract_wavelength_reference <- function(
  metadata_file,
  output_file = "data/wavelength_reference.csv"
) {
  cat("\n--- Wellenlängen-Extraktion ---\n")
  
  # Prüfe, ob Datei bereits existiert
  if (file.exists(output_file)) {
    cat("✓ Wellenlängen-Referenz existiert bereits:", output_file, "\n")
    return(invisible(read_csv(output_file, show_col_types = FALSE)))
  }
  
  cat("Extrahiere Wellenlängen aus:", basename(metadata_file), "\n")
  
  # Lade XML
  metadata <- tryCatch({
    read_xml(metadata_file)
  }, error = function(e) {
    stop("❌ Fehler beim Lesen der Metadaten-Datei: ", e$message)
  })
  
  # Extrahiere zentrale Wellenlängen
  wavelengths <- tryCatch({
    xml_find_all(metadata, ".//bandID/wavelengthCenterOfBand") %>%
      xml_double()
  }, error = function(e) {
    stop("❌ Fehler beim Extrahieren der Wellenlängen: ", e$message)
  })
  
  if (length(wavelengths) == 0) {
    stop("❌ Keine Wellenlängen in XML-Datei gefunden")
  }
  
  # Erstelle Data Frame mit Band-Namen und Wellenlängen
  band_names <- paste0("Band", seq_along(wavelengths), "_corr")
  wavelength_df <- data.frame(
    band_name = band_names,
    band_number = seq_along(wavelengths),
    wavelength_nm = wavelengths
  )
  
  # Erstelle Ausgabeverzeichnis falls nötig
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Speichere Wellenlängen-Referenz
  write_csv(wavelength_df, output_file)
  
  cat("✓ Wellenlängen-Referenz gespeichert:", output_file, "\n")
  cat("  → Anzahl Bänder:", length(wavelengths), "\n")
  cat("  → Wellenlängen-Bereich:", round(min(wavelengths)), "-", round(max(wavelengths)), "nm\n\n")
  
  return(invisible(wavelength_df))
}

#' Extrahiert Metadaten aus allen EnMAP XML-Dateien
#'
#' @param extract_dir Verzeichnis mit extrahierten EnMAP-Daten
#' @param shapefile_path Pfad zum Shapefile mit Koordinaten
#' @param output_file Ausgabedatei für metadata_mapping.csv
#' @param wavelength_output Ausgabedatei für wavelength_reference.csv
#' @return Data frame mit Metadaten
extract_metadata_mapping <- function(
  extract_dir = "data/extracted_files",
  shapefile_path = NULL,
  output_file = "data/metadata_mapping.csv",
  wavelength_output = "data/wavelength_reference.csv"
) {
  cat("\n════════════════════════════════════════════════════════════════════════════════\n")
  cat("                    METADATEN-EXTRAKTION & KOORDINATENABGLEICH\n")
  cat("════════════════════════════════════════════════════════════════════════════════\n\n")
  
  # Prüfe, ob Output bereits existiert
  if (file.exists(output_file)) {
    cat("⚠️  Metadaten-Datei existiert bereits:", output_file, "\n")
    cat("Soll sie neu erstellt werden? (j/n): ")
    response <- tolower(trimws(readline()))
    if (response != "j" && response != "y") {
      cat("✓ Verwende vorhandene Metadaten-Datei\n")
      return(invisible(fread(output_file)))
    }
    cat("→ Erstelle Metadaten-Datei neu...\n\n")
  }
  
  cat("Lade Metadaten aus Verzeichnis:", extract_dir, "\n")
  
  # Finde alle Metadata XML-Dateien
  metadata_files <- list.files(
    extract_dir,
    pattern = "-METADATA\\.XML$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  if (length(metadata_files) == 0) {
    stop("❌ Keine Metadaten-XML-Dateien gefunden in: ", extract_dir)
  }
  
  cat("✓ Gefundene Metadaten-Dateien:", length(metadata_files), "\n\n")
  
  # Extrahiere Wellenlängen-Referenz aus erster Metadaten-Datei
  extract_wavelength_reference(
    metadata_file = metadata_files[1],
    output_file = wavelength_output
  )
  
  # Optional: Lade Shapefile um relevante Raster zu filtern
  relevant_rasters <- NULL
  if (!is.null(shapefile_path) && file.exists(shapefile_path)) {
    cat("Lade Shapefile für Koordinatenabgleich:", shapefile_path, "\n")
    tryCatch({
      shp <- st_read(shapefile_path, quiet = TRUE)
      coords <- st_coordinates(shp)
      cat("✓ Shapefile geladen mit", nrow(coords), "Koordinaten\n")
      
      # Finde relevante Spektralbilder
      cat("\nPrüfe, welche Raster die Shapefile-Koordinaten enthalten...\n")
      spectral_files <- list.files(
        extract_dir,
        pattern = "SPECTRAL_IMAGE\\.TIF$",
        full.names = TRUE,
        recursive = TRUE,
        ignore.case = TRUE
      )
      
      relevant_rasters <- c()
      for (raster_file in spectral_files) {
        tryCatch({
          r <- rast(raster_file)
          ext <- ext(r)
          coords_in_extent <- coords[, 1] >= ext[1] & coords[, 1] <= ext[2] &
                              coords[, 2] >= ext[3] & coords[, 2] <= ext[4]
          
          if (any(coords_in_extent)) {
            # Extrahiere Basisnamen
            base_name <- gsub("-SPECTRAL_IMAGE\\.TIF$", "", basename(raster_file), ignore.case = TRUE)
            relevant_rasters <- c(relevant_rasters, base_name)
            cat("  ✓", base_name, "-", sum(coords_in_extent), "Koordinaten\n")
          }
        }, error = function(e) {
          # Fehler ignorieren
        })
      }
      
      if (length(relevant_rasters) > 0) {
        cat("\n→ Nur", length(relevant_rasters), "relevante Raster werden in Metadaten aufgenommen\n\n")
      } else {
        cat("\n⚠️  Keine Raster enthalten Shapefile-Koordinaten. Alle Metadaten werden extrahiert.\n\n")
        relevant_rasters <- NULL
      }
    }, error = function(e) {
      cat("⚠️  Warnung: Shapefile-Abgleich fehlgeschlagen:", e$message, "\n")
      cat("→ Alle Metadaten werden extrahiert.\n\n")
    })
  }
  
  # Initialisiere Liste für Metadaten
  metadata_list <- list()
  processed_count <- 0
  
  # Verarbeite jede Metadaten-Datei
  for (i in seq_along(metadata_files)) {
    metadata_file <- metadata_files[i]
    
    # Extrahiere Basisnamen für Matching mit Spektraldateien
    base_name <- gsub("-METADATA\\.XML$", "", basename(metadata_file), ignore.case = TRUE)
    
    # Filtere, falls nur relevante Raster gewünscht
    if (!is.null(relevant_rasters) && !(base_name %in% relevant_rasters)) {
      next
    }
    
    processed_count <- processed_count + 1
    cat(sprintf("[%d/%d] Verarbeite: %s\n", processed_count, length(metadata_files), base_name))
    
    # Lade XML
    metadata <- tryCatch({
      read_xml(metadata_file)
    }, error = function(e) {
      cat("  ❌ Fehler beim Lesen:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(metadata)) {
      next
    }
    
    # Finde zugehöriges Spektralbild im gleichen Verzeichnis wie die Metadaten-Datei
    metadata_dir <- dirname(metadata_file)
    spectral_image_pattern <- paste0(base_name, "-SPECTRAL_IMAGE\\.(TIF|BSQ)$")
    
    spectral_image_path <- list.files(
      metadata_dir,
      pattern = spectral_image_pattern,
      full.names = TRUE,
      recursive = FALSE,
      ignore.case = TRUE
    )
    
    if (length(spectral_image_path) == 0) {
      cat("  ⚠️  Warnung: Kein Spektralbild gefunden für", base_name, "\n")
      cat("     Gesucht in:", metadata_dir, "\n")
      cat("     Pattern:", spectral_image_pattern, "\n")
      spectral_image_path <- NA_character_
    } else {
      spectral_image_path <- spectral_image_path[1]
      cat("  ✓ Spektralbild gefunden:", basename(spectral_image_path), "\n")
    }
    
    # Extrahiere Metadaten-Parameter
    tryCatch({
      metadata_params <- list(
        SpectralFile = base_name,
        RasterFile = spectral_image_path,
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
      
      metadata_list[[length(metadata_list) + 1]] <- metadata_params
      cat("  ✓ Sonnenhöhenwinkel:", round(metadata_params$sun_elevation, 2), "°\n")
      
    }, error = function(e) {
      cat("  ❌ Fehler beim Extrahieren:", e$message, "\n")
    })
  }
  
  # Kombiniere alle Metadaten in einen Data Frame
  if (length(metadata_list) == 0) {
    stop("❌ Keine Metadaten konnten erfolgreich extrahiert werden")
  }
  
  metadata_mapping <- do.call(rbind, lapply(metadata_list, as.data.frame))
  
  # Erstelle Ausgabeverzeichnis falls nötig
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Speichere Metadaten-Mapping
  write.csv(metadata_mapping, output_file, row.names = FALSE)
  
  cat("\n════════════════════════════════════════════════════════════════════════════════\n")
  cat("✅ METADATEN-EXTRAKTION ABGESCHLOSSEN\n")
  cat("════════════════════════════════════════════════════════════════════════════════\n")
  cat("→ Gespeichert:", output_file, "\n")
  cat("→ Anzahl Datensätze:", nrow(metadata_mapping), "\n")
  cat("→ Spalten:", paste(colnames(metadata_mapping), collapse = ", "), "\n\n")
  
  return(invisible(metadata_mapping))
}

# Hauptfunktion zum direkten Ausführen
main <- function() {
  extract_metadata_mapping(
    extract_dir = "data/extracted_files",
    shapefile_path = "data/Peru_Bolivien_reprojiziert.shp",
    output_file = "data/metadata_mapping.csv",
    wavelength_output = "data/wavelength_reference.csv"
  )
}

# HINWEIS: Auto-Execution auskommentiert - Script wird nur über master_pipeline.R gesteuert
# # Ausführen wenn direkt aufgerufen
# if (!interactive()) {
#   main()
# }
