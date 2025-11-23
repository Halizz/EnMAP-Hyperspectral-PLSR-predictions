library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# ============================================================================
# Skript: plot_spectra_from_shapefile.R
# -----------------------------------------------------------------------------
# Beschreibung:
# Dieses Skript extrahiert Spektren aus einem EnMAP-Raster an den Zentroiden
# von Polygonen eines Shapefiles und erstellt einen Plot der Spektralkurven.
#
# Voraussetzungen:
# - R-Pakete: terra, sf, ggplot2, dplyr, readr, tidyr
# - Korrekte Pfade zu den Eingabedateien (Shapefile, Raster, CSV mit Wellenlängen)
#
# Nutzung:
# 1. Funktion plot_spectra_from_shapefile() mit gewünschten Parametern aufrufen
# 2. Beispiel: plot_spectra_from_shapefile(wanted_ids = 5:8)
#
# Autor: [Dein Name]
# Datum: [Aktuelles Datum]
# ============================================================================

# -----------------------------
# Funktion: Spektren aus Shapefile extrahieren und plotten
# -----------------------------
plot_spectra_from_shapefile <- function(
  shapefile_path = "data/Spektren_fluss2.shp", 
  wanted_ids = 5:8,
  raster_dir = "F:/Uni/Bachelorarbeit/R working directory/enmap-soil-correlation/data/extracted_files/Ich_everzwlfe",
  raster_suffix = "_smoothed\\+adaptive_destriped_masked_snv\\.tif$",
  wavelength_csv = "data/wavelength_reference.csv",
  output_dir = "F:/Uni/Bachelorarbeit/R working directory/enmap-soil-correlation/Spektren",
  output_filename = NULL,
  predict_bs = TRUE  # Neuer Parameter zum Aktivieren/Deaktivieren der BS-Vorhersagen
) {
  
  # Ausgabeverzeichnis erstellen, falls nicht vorhanden
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Standardwert für output_filename, wenn nicht angegeben
  if (is.null(output_filename)) {
    id_range <- paste0("id", min(wanted_ids), "-", max(wanted_ids))
    output_filename <- paste0("spektren_plot_", id_range, ".png")
  }
  
  # Lade Wellenlängen aus CSV
  wavelengths_df <- read_csv(wavelength_csv, show_col_types = FALSE)
  wavelength_vector <- wavelengths_df$wavelength_nm
  
  # Hauptlogik in tryCatch für robuste Fehlerbehandlung
  tryCatch({
    # 1. Shapefile laden und prüfen
    if (!file.exists(shapefile_path)) {
      stop(paste("Shapefile nicht gefunden:", shapefile_path))
    }
    shp <- st_read(shapefile_path, quiet = TRUE)
    print("Shapefile erfolgreich geladen mit Anzahl Objekte:")
    print(nrow(shp))
    
    # Prüfe die verfügbaren IDs im Shapefile
    print("Verfügbare IDs im Shapefile:")
    print(shp)
    
    # Behandle ID-Spalte (Fall-insensitiv)
    if ("id" %in% names(shp) && !"ID" %in% names(shp)) {
      print("Kleingeschriebene 'id'-Spalte gefunden, wird in 'ID' umbenannt")
      shp <- shp %>% rename(ID = id)
    } else if (!"ID" %in% names(shp) && !"id" %in% names(shp)) {
      print("Keine ID-Spalte gefunden, erzeuge IDs 1 bis n:")
      shp$ID <- 1:nrow(shp)
    }
    
    # 2. Rasterbilder mit dem gewünschten Suffix suchen
    raster_files <- list.files(raster_dir, pattern = raster_suffix, full.names = TRUE, recursive = TRUE)
    print("Gefundene Rasterbilder mit passendem Suffix:")
    print(raster_files)
    
    # Prüfe alle Rasterbilder und wähle das mit den meisten gefundenen IDs
    best_raster_idx <- NULL
    best_raster_valid_coords <- NULL
    best_raster_ids <- NULL
    best_raster_img <- NULL
    max_found_ids <- 0
    
    coords <- st_coordinates(st_centroid(shp))
    coords_ids <- coords[shp$ID %in% wanted_ids, , drop = FALSE]
    ids_valid <- shp$ID[shp$ID %in% wanted_ids]
    
    if (nrow(coords_ids) == 0) stop("Keine passenden IDs im Shapefile gefunden.")
    
    print(paste("Suche nach Rastern für", length(wanted_ids), "gewünschte IDs:", paste(wanted_ids, collapse=", ")))
    
    for (raster_idx in 1:length(raster_files)) {
      raster_file <- raster_files[raster_idx]
      print(paste("Prüfe Raster:", basename(raster_file)))
      
      img <- rast(raster_file)
      ext_img <- ext(img)
      
      # Prüfe, ob Koordinaten im Rasterbereich liegen
      valid_coords <- coords_ids[,1] >= ext_img[1] & coords_ids[,1] <= ext_img[2] &
                      coords_ids[,2] >= ext_img[3] & coords_ids[,2] <= ext_img[4]
      
      coords_valid <- coords_ids[valid_coords, , drop = FALSE]
      ids_valid_in_raster <- ids_valid[valid_coords]
      
      if (nrow(coords_valid) > 0) {
        print(paste("  Gefunden:", length(ids_valid_in_raster), "IDs:", paste(ids_valid_in_raster, collapse=", ")))
        
        # Wenn dieses Raster mehr IDs enthält als das bisher beste, wähle dieses
        if (length(ids_valid_in_raster) > max_found_ids) {
          max_found_ids <- length(ids_valid_in_raster)
          best_raster_idx <- raster_idx
          best_raster_valid_coords <- coords_valid
          best_raster_ids <- ids_valid_in_raster
          best_raster_img <- img
          
          # Wenn wir alle gewünschten IDs gefunden haben, können wir abbrechen
          if (max_found_ids == length(wanted_ids)) {
            print("  Alle gewünschten IDs gefunden!")
            break
          }
        }
      } else {
        print("  Keine gültigen Koordinaten im Rasterbereich für dieses Raster")
      }
    }
    
    if (is.null(best_raster_idx)) {
      stop("Keine gültigen Koordinaten in irgendeinem Rasterbereich gefunden.")
    }
    
    # Verwende das beste gefundene Raster
    raster_file <- raster_files[best_raster_idx]
    img <- best_raster_img
    coords_valid <- best_raster_valid_coords
    ids_valid_in_raster <- best_raster_ids
    
    print(paste("Verwende Raster mit den meisten gefundenen IDs:", basename(raster_file)))
    print(paste("Gefundene IDs:", paste(ids_valid_in_raster, collapse=", ")))
    print(paste("Fehlende IDs:", paste(setdiff(wanted_ids, ids_valid_in_raster), collapse=", ")))
    
    # 3. Extrahiere Spektren
    spectra <- terra::extract(img, coords_valid)
    spectra$ID <- ids_valid_in_raster
    band_names <- names(img)
    
    # Extrahiere Basennamen des Rasterordners für die Vorhersagepfade
    raster_basename <- basename(dirname(raster_file))
    
    # Extrahiere BS-Vorhersagewerte, wenn gewünscht
    bs_predictions <- NULL
    if (predict_bs) {
      # Definiere Pfade zu den BS-Vorhersageraster
      pearson_bs_path <- file.path("results", "Pearson", "filter+plsr", "raster_predictions", 
                                  raster_basename, "BS", "BS_predicted_Pearson_filter+plsr.tif")
      spearman_bs_path <- file.path("results", "Spearman", "filter+plsr", "raster_predictions", 
                                   raster_basename, "BS", "BS_predicted_Spearman_filter+plsr.tif")
      
      # Überprüfe die Existenz der Dateien
      pearson_exists <- file.exists(pearson_bs_path)
      spearman_exists <- file.exists(spearman_bs_path)
      
      if (pearson_exists || spearman_exists) {
        print("Extrahiere BS-Vorhersagewerte...")
        
        # Initialisiere Dataframe für BS-Vorhersagen
        bs_predictions <- data.frame(
          ID = ids_valid_in_raster,
          Pearson_BS = NA,
          Spearman_BS = NA
        )
        
        # Extrahiere Pearson BS-Werte, wenn verfügbar
        if (pearson_exists) {
          print(paste("  Pearson BS-Vorhersageraster gefunden:", pearson_bs_path))
          pearson_bs_raster <- rast(pearson_bs_path)
          pearson_bs_values <- terra::extract(pearson_bs_raster, coords_valid)
          if (ncol(pearson_bs_values) > 0) {
            bs_predictions$Pearson_BS <- pearson_bs_values[,1]
          } else {
            print("  Keine Werte aus Pearson BS-Raster extrahiert.")
          }
        } else {
          print("  Pearson BS-Vorhersageraster nicht gefunden.")
        }
        
        # Extrahiere Spearman BS-Werte, wenn verfügbar
        if (spearman_exists) {
          print(paste("  Spearman BS-Vorhersageraster gefunden:", spearman_bs_path))
          spearman_bs_raster <- rast(spearman_bs_path)
          spearman_bs_values <- terra::extract(spearman_bs_raster, coords_valid)
          if (ncol(spearman_bs_values) > 0) {
            bs_predictions$Spearman_BS <- spearman_bs_values[,1]
          } else {
            print("  Keine Werte aus Spearman BS-Raster extrahiert.")
          }
        } else {
          print("  Spearman BS-Vorhersageraster nicht gefunden.")
        }
        
        # Ausgabe der extrahierten BS-Werte
        print("Extrahierte BS-Vorhersagewerte:")
        print(bs_predictions)
      } else {
        print("Keine BS-Vorhersageraster gefunden für diesen Rasterdatensatz.")
      }
    }
    
    # 4. Spektraldaten für alle Punkte aufbereiten
    all_spectra <- data.frame()
    for (i in 1:nrow(spectra)) {
      point_id <- spectra$ID[i]
      if (is.null(point_id) || is.na(point_id)) next  # Überspringe Punkte ohne ID
      
      spectral_values <- as.numeric(unlist(spectra[i, band_names]))
      n_bands <- length(spectral_values)
      
      # Verwende entweder die Wellenlängen aus der CSV oder eine Sequenz
      if (length(wavelength_vector) >= n_bands) {
        wavelengths <- wavelength_vector[1:n_bands]
      } else {
        wavelengths <- seq_len(n_bands)
      }
      
      # Erstelle Datenrahmen für diesen Punkt
      point_df <- data.frame(
        ID = rep(point_id, n_bands),
        band_number = seq_len(n_bands),
        reflectance = spectral_values,
        wavelength_nm = wavelengths
      )
      all_spectra <- rbind(all_spectra, point_df)
    }
    
    # Prüfe, ob Spektraldaten gefunden wurden
    if (nrow(all_spectra) == 0) stop("Keine Spektraldaten für die gewünschten IDs gefunden.")
    
    # Entferne NA-Werte, die das Plotten stören würden
    all_spectra <- all_spectra[!is.na(all_spectra$reflectance), ]
    
    print("Struktur der gefilterten Spektraldaten:")
    print(str(all_spectra))
    
    # 5. Plotten und Speichern als PNG
    output_file <- file.path(output_dir, output_filename)
    
    # Erhalte Liste der tatsächlich gefundenen IDs für den Titel
    found_ids <- sort(unique(all_spectra$ID))
    found_ids_str <- paste(found_ids, collapse=", ")
    
    plot <- ggplot(all_spectra, aes(x = wavelength_nm, y = reflectance, color = as.factor(ID), group = ID)) +
      
      geom_point(size = 0.5, alpha = 0.5) +  # Kleine Punkte zur Markierung der Datenpunkte
      labs(
        x = "Wellenlänge [nm]",
        y = "Reflektanz",
        color = "Polygon ID",
        title = paste("Spektren der IDs", found_ids_str, "aus", basename(raster_file))
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 10)
      )
    
    # Speichere Plot als PNG
    ggsave(
      filename = output_file,
      plot = plot,
      width = 10,
      height = 6,
      dpi = 300
    )
    
    print(plot)
    
    # Erweiterte Rückgabe mit BS-Vorhersagen
    return_list <- list(
      plot = plot,
      data = all_spectra,
      raster_file = raster_file,
      output_file = output_file,
      found_ids = found_ids
    )
    
    if (!is.null(bs_predictions)) {
      return_list$bs_predictions <- bs_predictions
    }
    
    return(return_list)
    
  }, error = function(e) {
    stop(paste("Fehler beim Verarbeiten:", e$message))
  })
}

# -----------------------------
# Ausführung mit Standardparametern
# -----------------------------
# HINWEIS: Auto-Execution auskommentiert - Script wird nur über master_pipeline.R gesteuert
# # Funktion aufrufen und Ergebnis speichern
# result <- plot_spectra_from_shapefile()
# 
# # BS-Vorhersagen anzeigen, falls verfügbar
# if (!is.null(result$bs_predictions)) {
#   print("Zusammenfassung der BS-Vorhersagen:")
#   print(result$bs_predictions)
# }
