#
# Wissenschaftliche RGB-Visualisierung von geglätteten Raster-Daten
# Dieses Script erstellt hochwertige RGB-Karten mit wissenschaftlichem Anspruch

# Verbesserte Package-Verwaltung
load_required_packages <- function() {
  required_packages <- c("terra", "raster", "RColorBrewer", "dplyr", "readr", 
                         "pbapply", "tmap", "mapview", "viridis", "rgdal",
                         "sf", "ggplot2", "gridExtra")
  
  # Prüfe, ob Pakete installiert sind
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    cat("WARNUNG: Folgende Pakete fehlen und werden installiert:", paste(missing_packages, collapse=", "), "\n")
    install.packages(missing_packages)
  }
  
  # Lade Pakete
  invisible(sapply(required_packages, function(pkg) {
    success <- suppressWarnings(suppressPackageStartupMessages(
      require(pkg, character.only = TRUE)
    ))
    
    if (!success) {
      cat("FEHLER: Konnte Paket", pkg, "nicht laden!\n")
    }
  }))
  
  # Prüfe, ob terra verfügbar ist
  terra_available <<- suppressWarnings(requireNamespace("terra", quietly = TRUE))
  if (!terra_available) {
    cat("WARNUNG: terra Paket konnte nicht geladen werden. Fallback zu raster wird verwendet.\n")
    requireNamespace("raster", quietly = TRUE)
  }
  
  # Prüfe, ob tmap verfügbar ist
  tmap_available <<- suppressWarnings(requireNamespace("tmap", quietly = TRUE))
  if (!tmap_available) {
    cat("WARNUNG: tmap Paket konnte nicht geladen werden. Fallback zu base plot wird verwendet.\n")
  } else {
    # Aktiviere tmap interaktiven Modus für Vorschau
    tmap::tmap_mode("plot")
  }
}

# Lade die benötigten Pakete
load_required_packages()

# Verbesserte Wrapper-Funktion für Raster-Operationen
read_raster <- function(file_path, max_attempts = 3, verbose = TRUE) {
  if(verbose) cat(sprintf("Lade Raster: %s\n", basename(file_path)))
  
  for(attempt in 1:max_attempts) {
    if(attempt > 1 && verbose) cat(sprintf("  Versuch %d von %d...\n", attempt, max_attempts))
    
    # Terra hat Priorität, wenn verfügbar
    if(exists("terra_available") && terra_available) {
      tryCatch({
        if(verbose) cat("  Verwende terra::rast()\n")
        r <- terra::rast(file_path)
        if(verbose) cat(sprintf("  ✓ Geladen: %d Zeilen, %d Spalten, %d Bänder\n", 
                               nrow(r), ncol(r), terra::nlyr(r)))
        return(r)
      }, error = function(e) {
        if(verbose) cat("  ⚠️ terra::rast() fehlgeschlagen:", e$message, "\n")
      })
    }
    
    # Fallback zu raster, wenn terra nicht funktioniert hat
    tryCatch({
      if(verbose) cat("  Verwende raster::brick()\n")
      r <- raster::brick(file_path)
      if(verbose) cat(sprintf("  ✓ Geladen: %d Zeilen, %d Spalten, %d Bänder\n", 
                             nrow(r), ncol(r), raster::nlayers(r)))
      return(r)
    }, error = function(e) {
      if(verbose) cat("  ⚠️ raster::brick() fehlgeschlagen:", e$message, "\n")
    })
    
    # Weitere Fallbacks probieren
    tryCatch({
      if(verbose) cat("  Verwende raster::stack()\n")
      r <- raster::stack(file_path)
      if(verbose) cat(sprintf("  ✓ Geladen: %d Zeilen, %d Spalten, %d Bänder\n", 
                             nrow(r), ncol(r), raster::nlayers(r)))
      return(r)
    }, error = function(e) {
      if(verbose) cat("  ⚠️ raster::stack() fehlgeschlagen:", e$message, "\n")
    })
    
    # Pause zwischen Versuchen
    if(attempt < max_attempts) {
      Sys.sleep(1)
      gc()  # Speicher bereinigen
    }
  }
  
  stop(sprintf("Konnte Rasterdatei nicht laden nach %d Versuchen: %s", max_attempts, file_path))
}

# ======================= WISSENSCHAFTLICHE KARTENFUNKTIONEN =======================

# Professionelle Nordpfeil-Funktion für wissenschaftliche Karten
create_scientific_north_arrow <- function(width = 1.5, height = 2, style = "standard") {
  # Erstelle eine leere Zeichenfläche
  north_arrow <- function() {
    par(mai = c(0, 0, 0, 0))
    plot(0, 0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1), 
         axes = FALSE, xlab = "", ylab = "", asp = 1)
    
    # Verschiedene Nordpfeil-Stile
    if(style == "compass") {
      # Kompass-Stil mit Kreisscheibe
      # Zeichne Hintergrundkreis
      symbols(0, 0, circles = 0.8, inches = FALSE, add = TRUE, 
              fg = "black", bg = "white", lwd = 2)
      
      # Zeichne die Hauptrichtungen
      directions <- c("N", "NO", "O", "SO", "S", "SW", "W", "NW")
      angles <- seq(90, -270, by = -45) * pi / 180
      
      for(i in 1:length(directions)) {
        # Pfeil
        x_end <- 0.7 * cos(angles[i])
        y_end <- 0.7 * sin(angles[i])
        
        # N und S hervorheben
        if(directions[i] %in% c("N", "S")) {
          arrows(0, 0, x_end, y_end, length = 0.1, angle = 25, 
                 code = 2, col = "red", lwd = 2.5)
        } else {
          arrows(0, 0, x_end, y_end, length = 0.1, angle = 20, 
                 code = 2, col = "black", lwd = 1.2)
        }
        
        # Beschriftung
        x_text <- 0.9 * cos(angles[i])
        y_text <- 0.9 * sin(angles[i])
        text(x_text, y_text, directions[i], cex = 1.0, font = ifelse(directions[i] == "N", 2, 1),
             col = ifelse(directions[i] == "N", "red", "black"))
      }
      
    } else {
      # Standard-Nordpfeil (eleganter und platzsparend)
      # Pfeil nach Norden
      polygon(c(-0.15, 0, 0.15), c(0, 0.8, 0), col = "black")
      polygon(c(-0.15, 0, 0.15), c(0, -0.8, 0), col = "white", border = "black")
      
      # N-Beschriftung
      text(0, 0.9, "N", cex = 1.5, font = 2)
    }
  }
  
  # Erstelle den Nordpfeil als Grafikobjekt
  tf <- tempfile(fileext = ".png")
  png(tf, width = width*100, height = height*100, res = 100, bg = "transparent")
  north_arrow()
  dev.off()
  
  # Gib den Pfad zum Nordpfeil zurück
  return(tf)
}

# Professionelle Maßstabsleiste für wissenschaftliche Karten
create_scientific_scale_bar <- function(extent, n_divs = 3, 
                                        unit = "m", style = "alternating") {
  # Berechne eine sinnvolle Maßstabslänge (ca. 1/5 der Kartenbreite)
  width_map <- extent[2] - extent[1]
  
  # Finde "schöne" Maßstabswerte (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, etc.)
  log10_width <- floor(log10(width_map/5))
  base_value <- 10^log10_width
  
  candidates <- c(base_value, 2*base_value, 5*base_value)
  scale_length <- candidates[which.min(abs(candidates - width_map/5))]
  
  # Einheiten konvertieren wenn nötig
  if(unit == "km" && scale_length >= 1000) {
    scale_length_display <- scale_length / 1000
    unit_display <- "km"
  } else {
    scale_length_display <- scale_length
    unit_display <- "m"
  }
  
  segment_length <- scale_length / n_divs
  
  # Funktion zum Zeichnen der Skala
  draw_scale_bar <- function() {
    # Einstellungen
    bar_height <- 0.2
    
    # Leeren Plot erstellen
    par(mai = c(0.3, 0.1, 0.1, 0.1))
    plot(0, 0, type = "n", xlim = c(0, scale_length), ylim = c(0, 1), 
         axes = FALSE, xlab = "", ylab = "")
    
    # Segmente zeichnen
    for(i in 1:n_divs) {
      x_start <- (i-1) * segment_length
      x_end <- i * segment_length
      
      # Alternierende Farben
      if(style == "alternating") {
        col <- ifelse(i %% 2 == 1, "black", "white")
      } else {
        col <- "white"
      }
      
      # Segment zeichnen
      rect(x_start, 0.3, x_end, 0.3 + bar_height, 
           col = col, border = "black", lwd = 1.2)
      
      # Beschriftung für jeden Abschnitt
      if(i == 1 || i == n_divs) {
        text(x_start, 0.15, paste0(round((i-1)*segment_length_display, 2)), 
             adj = c(0.5, 0.5), cex = 0.9)
      }
    }
    
    # Beschriftung rechts
    text(scale_length, 0.15, 
         paste0(round(scale_length_display, 1), " ", unit_display), 
         adj = c(0.5, 0.5), cex = 0.9)
    
    # Optionale Umrissbox
    rect(0, 0, scale_length, 0.7, border = "black", lwd = 1, lty = 3)
  }
  
  # Erstelle die Skala als Grafikobjekt
  tf <- tempfile(fileext = ".png")
  png(tf, width = 300, height = 80, res = 100, bg = "transparent")
  draw_scale_bar()
  dev.off()
  
  return(list(
    image = tf,
    length = scale_length,
    unit = unit_display,
    display_value = scale_length_display
  ))
}

# ======================= VERBESSERTE TRUE-COLOR-VISUALISIERUNG =======================

# Vereinfachte und robuste Funktion zur True-Color-Visualisierung
create_true_color_visualization <- function(raster_data, wavelength_file = "data/wavelength_reference.csv") {
  # Status-Update
  cat("Erstelle True-Color-Darstellung...\n")
  
  # Prüfe, ob die Wellenlängen-Datei existiert
  if (!file.exists(wavelength_file)) {
    warning("Wellenlängen-Datei nicht gefunden: ", wavelength_file)
    
    # Fallback zu Standardbändern für EnMAP
    cat("Fallback zu Standardbändern: R=45 (650nm), G=28 (550nm), B=7 (450nm)\n")
    
    # Verwende direkt die bekannten Bandnummern für EnMAP
    red_band <- 45    # ca. 650 nm
    green_band <- 28  # ca. 550 nm
    blue_band <- 7    # ca. 450 nm
    
    # Extrahiere die Bänder
    if (class(raster_data)[1] == "SpatRaster") {
      if (terra::nlyr(raster_data) >= max(c(red_band, green_band, blue_band))) {
        rgb_stack <- raster_data[[c(red_band, green_band, blue_band)]]
      } else {
        # Wenn nicht genug Bänder vorhanden sind, verwende die ersten drei
        warning("Nicht genügend Bänder für True-Color-Darstellung. Verwende die ersten drei Bänder.")
        rgb_stack <- raster_data[[1:3]]
      }
    } else {
      if (raster::nlayers(raster_data) >= max(c(red_band, green_band, blue_band))) {
        rgb_stack <- raster::subset(raster_data, c(red_band, green_band, blue_band))
      } else {
        warning("Nicht genügend Bänder für True-Color-Darstellung. Verwende die ersten drei Bänder.")
        rgb_stack <- raster::subset(raster_data, 1:3)
      }
    }
    
    # Setze Bandnamen explizit für weitere Verarbeitung
    if (class(rgb_stack)[1] == "SpatRaster") {
      names(rgb_stack) <- c("Red", "Green", "Blue")
    } else {
      names(rgb_stack) <- c("Red", "Green", "Blue")
    }
    
    return(rgb_stack)
  }
  
  # Lade Wellenlängen-Referenzdaten
  cat("Lade Wellenlängen-Referenz aus:", wavelength_file, "\n")
  tryCatch({
    wavelengths <- read.csv(wavelength_file)
    
    # Validiere Format
    if (!all(c("wavelength_nm", "band_number") %in% colnames(wavelengths))) {
      stop("Ungültiges Format der Wellenlängen-Datei. Benötigte Spalten: wavelength_nm, band_number")
    }
    
    # Finde Bänder mit den besten Wellenlängen für True Color
    target_red <- 650    # nm
    target_green <- 550  # nm
    target_blue <- 450   # nm
    
    red_idx <- which.min(abs(wavelengths$wavelength_nm - target_red))
    green_idx <- which.min(abs(wavelengths$wavelength_nm - target_green))
    blue_idx <- which.min(abs(wavelengths$wavelength_nm - target_blue))
    
    red_band <- wavelengths$band_number[red_idx]
    green_band <- wavelengths$band_number[green_idx]
    blue_band <- wavelengths$band_number[blue_idx]
    
    cat("Gefundene True-Color-Bänder:\n")
    cat(sprintf("  Rot:   Band %3d (%.1f nm)\n", red_band, wavelengths$wavelength_nm[red_idx]))
    cat(sprintf("  Grün:  Band %3d (%.1f nm)\n", green_band, wavelengths$wavelength_nm[green_idx]))
    cat(sprintf("  Blau:  Band %3d (%.1f nm)\n", blue_band, wavelengths$wavelength_nm[blue_idx]))
    
    # Erstelle RGB-Stack basierend auf der Bibliothek
    if (class(raster_data)[1] == "SpatRaster") {
      if (terra::nlyr(raster_data) >= max(c(red_band, green_band, blue_band))) {
        rgb_stack <- raster_data[[c(red_band, green_band, blue_band)]]
        names(rgb_stack) <- c("Red", "Green", "Blue")
      } else {
        stop(sprintf("Raster hat zu wenige Bänder (%d). Benötigt: Band %d", 
                     terra::nlyr(raster_data), max(c(red_band, green_band, blue_band))))
      }
    } else {
      if (raster::nlayers(raster_data) >= max(c(red_band, green_band, blue_band))) {
        rgb_stack <- raster::subset(raster_data, c(red_band, green_band, blue_band))
        names(rgb_stack) <- c("Red", "Green", "Blue")
      } else {
        stop(sprintf("Raster hat zu wenige Bänder (%d). Benötigt: Band %d", 
                     raster::nlayers(raster_data), max(c(red_band, green_band, blue_band))))
      }
    }
    
    cat("True-Color-Stack erfolgreich erstellt.\n")
    return(rgb_stack)
    
  }, error = function(e) {
    warning(sprintf("Fehler bei True-Color-Erstellung: %s", e$message))
    
    # Standardbänder verwenden als Fallback
    cat("Fallback zu Standardbändern: R=45 (650nm), G=28 (550nm), B=7 (450nm)\n")
    
    red_band <- 45    # ca. 650 nm
    green_band <- 28  # ca. 550 nm
    blue_band <- 7    # ca. 450 nm
    
    # Prüfe, ob genügend Bänder vorhanden sind
    if (class(raster_data)[1] == "SpatRaster") {
      n_bands <- terra::nlyr(raster_data)
    } else {
      n_bands <- raster::nlayers(raster_data)
    }
    
    if (n_bands < max(c(red_band, green_band, blue_band))) {
      # Verwende die ersten drei Bänder als Notlösung
      warning("Nicht genügend Bänder für True-Color-Darstellung. Verwende die ersten drei Bänder.")
      
      if (class(raster_data)[1] == "SpatRaster") {
        rgb_stack <- raster_data[[1:min(3, n_bands)]]
      } else {
        rgb_stack <- raster::subset(raster_data, 1:min(3, n_bands))
      }
    } else {
      # Verwende die Standardbänder
      if (class(raster_data)[1] == "SpatRaster") {
        rgb_stack <- raster_data[[c(red_band, green_band, blue_band)]]
      } else {
        rgb_stack <- raster::subset(raster_data, c(red_band, green_band, blue_band))
      }
    }
    
    # Setze Namen
    if (class(rgb_stack)[1] == "SpatRaster") {
      names(rgb_stack) <- c("Red", "Green", "Blue")
    } else {
      names(rgb_stack) <- c("Red", "Green", "Blue")
    }
    
    return(rgb_stack)
  })
}

# Erstellt eine wissenschaftliche RGB-Karte mit vollständigen kartografischen Elementen
create_scientific_rgb_map <- function(
  raster_file,                 # Pfad zum Raster
  output_file,                 # Ausgabe-PNG-Datei
  rgb_bands = NULL,            # Bänder für RGB, wenn NULL -> automatische True-Color
  stretch = "hist",            # Stretch-Methode: "hist" oder "lin"
  width = 2400,                # Bildbreite in Pixel (erhöht für bessere Qualität)
  height = 1800,               # Bildhöhe in Pixel
  dpi = 300,                   # DPI für wissenschaftliche Druckqualität
  title = NULL,                # Haupttitel
  subtitle = NULL,             # Untertitel
  add_grid = TRUE,             # Koordinatengitter hinzufügen?
  add_scale = TRUE,            # Maßstabsleiste hinzufügen?
  add_north = TRUE,            # Nordpfeil hinzufügen?
  add_points = FALSE,          # Probenahmepunkte hinzufügen?
  points_file = "data/sampling_points.csv"  # Datei mit Punktkoordinaten
) {
  # Start-Nachricht
  cat(sprintf("\n==== Erstelle wissenschaftliche RGB-Karte ====\n"))
  cat(sprintf("Eingabe: %s\n", basename(raster_file)))
  cat(sprintf("Ausgabe: %s\n", basename(output_file)))
  
  # Starte Timer
  start_time <- Sys.time()
  
  # Lade das Raster
  raster_data <- read_raster(raster_file)
  
  # Prüfe Raster-Dimensionen
  if (class(raster_data)[1] == "SpatRaster") {
    n_bands <- terra::nlyr(raster_data)
    n_rows <- nrow(raster_data)
    n_cols <- ncol(raster_data)
  } else {
    n_bands <- raster::nlayers(raster_data)
    n_rows <- nrow(raster_data)
    n_cols <- ncol(raster_data)
  }
  
  cat(sprintf("Raster: %d×%d Pixel, %d Bänder\n", n_rows, n_cols, n_bands))
  
  # Erstelle RGB-Kombination
  if (is.null(rgb_bands)) {
    # True-Color-Darstellung
    cat("Erstelle automatische True-Color-Darstellung...\n")
    rgb_subset <- create_true_color_visualization(raster_data)
  } else {
    # Benutzerdefinierte RGB-Kombination
    cat(sprintf("Erstelle benutzerdefinierte RGB-Darstellung mit Bändern %s\n", 
               paste(rgb_bands, collapse=", ")))
    
    if (max(rgb_bands) > n_bands) {
      stop(sprintf("Band %d für RGB angefordert, aber Raster hat nur %d Bänder", 
                  max(rgb_bands), n_bands))
    }
    
    # Extrahiere die gewünschten Bänder
    if (class(raster_data)[1] == "SpatRaster") {
      rgb_subset <- raster_data[[rgb_bands]]
      names(rgb_subset) <- c("Red", "Green", "Blue")
    } else {
      rgb_subset <- raster::subset(raster_data, rgb_bands)
      names(rgb_subset) <- c("Red", "Green", "Blue")
    }
  }
  
  # Prüfe, ob alle RGB-Bänder gültige Werte haben
  cat("Validiere RGB-Bänder...\n")
  band_stats <- list()
  
  for (i in 1:3) {
    # Extrahiere Werte
    if (class(rgb_subset)[1] == "SpatRaster") {
      band_values <- terra::values(rgb_subset[[i]])
    } else {
      band_values <- raster::values(rgb_subset[[i]])
    }
    
    # Berechne Statistiken
    valid_values <- band_values[!is.na(band_values)]
    
    band_stats[[i]] <- list(
      min = if(length(valid_values) > 0) min(valid_values, na.rm = TRUE) else NA,
      max = if(length(valid_values) > 0) max(valid_values, na.rm = TRUE) else NA,
      na_percent = sum(is.na(band_values)) / length(band_values) * 100
    )
    
    cat(sprintf("  Band %d: Min=%.4f, Max=%.4f, NA=%.1f%%\n", 
               i, band_stats[[i]]$min, band_stats[[i]]$max, band_stats[[i]]$na_percent))
    
    # Warne bei leeren Bändern
    if (is.na(band_stats[[i]]$min) || length(valid_values) == 0) {
      warning(sprintf("Band %d enthält keine gültigen Daten!", i))
    }
  }
  
  # Setze das Koordinatensystem (UTM Zone 19S für die Chile-Daten)
  if (class(rgb_subset)[1] == "SpatRaster") {
    terra::crs(rgb_subset) <- "EPSG:32719"
  } else {
    raster::crs(rgb_subset) <- sp::CRS("+init=epsg:32719")
  }
  cat("Koordinatensystem auf EPSG:32719 (UTM Zone 19S) gesetzt\n")
  
  # Extrahiere Bildname für den Titel
  image_name <- gsub("_smoothed\\+.*$|_smoothed$", "", basename(raster_file))
  
  # Setze Standardtitel, wenn keiner angegeben
  if (is.null(title)) {
    if (is.null(rgb_bands)) {
      title <- "True-Color-Darstellung"
    } else {
      title <- sprintf("RGB-Darstellung (R=%d, G=%d, B=%d)", 
                      rgb_bands[1], rgb_bands[2], rgb_bands[3])
    }
  }
  
  # Setze Standarduntertitel, wenn keiner angegeben
  if (is.null(subtitle)) {
    subtitle <- sprintf("Bild: %s", image_name)
  }
  
  # Lade Probenahmepunkte, wenn gewünscht
  sample_points <- NULL
  if (add_points && file.exists(points_file)) {
    cat("Lade Probenahmepunkte aus:", points_file, "\n")
    tryCatch({
      sample_points <- read.csv(points_file)
      
      # Prüfe, ob erforderliche Spalten vorhanden sind
      if (!all(c("x", "y") %in% colnames(sample_points))) {
        warning("Punktdaten benötigen 'x' und 'y' Spalten")
        sample_points <- NULL
      } else {
        cat(sprintf("  %d Punkte geladen\n", nrow(sample_points)))
      }
    }, error = function(e) {
      cat("Fehler beim Laden der Punkte:", e$message, "\n")
    })
  }
  
  # Erstelle Verzeichnis für die Ausgabe
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("Verzeichnis erstellt: %s\n", output_dir))
  }
  
  # ----- Kartenerstellung mit tmap (wenn verfügbar) -----
  if (exists("tmap_available") && tmap_available) {
    cat("Erstelle Karte mit tmap...\n")
    
    # Benutze tmap für hochwertige Karten
    tryCatch({
      # Vorbereitung für Koordinatengitter
      if (class(rgb_subset)[1] == "SpatRaster") {
        ext <- terra::ext(rgb_subset)
      } else {
        ext <- raster::extent(rgb_subset)
      }
      
      # Berechne sinnvolle Grid-Intervalle
      x_range <- ext[2] - ext[1]
      y_range <- ext[4] - ext[3]
      
      x_magnitude <- 10^floor(log10(x_range/5))
      y_magnitude <- 10^floor(log10(y_range/5))
      
      x_candidates <- c(1, 2, 5) * x_magnitude
      y_candidates <- c(1, 2, 5) * y_magnitude
      
      grid_x <- x_candidates[which.min(abs(x_candidates - x_range/5))]
      grid_y <- y_candidates[which.min(abs(y_candidates - y_range/5))]
      
      # Angepasstes tmap-Layout
      tm <- tmap::tm_shape(rgb_subset) +
        tmap::tm_rgb(stretch.palette = stretch == "hist") +
        tmap::tm_layout(
          title = title,
          title.position = c("center", "top"),
          title.size = 1.5,
          title.fontface = "bold",
          frame = TRUE,
          frame.double.line = TRUE,
          outer.margins = 0.02,
          inner.margins = c(0.04, 0.04, 0.08, 0.04),
          asp = 0  # Automatisches Seitenverhältnis beibehalten
        ) +
        tmap::tm_credits(subtitle, position = c("center", "bottom"), size = 1.0)
      
      # Füge Gitter hinzu
      if (add_grid) {
        tm <- tm + tmap::tm_grid(
          x = seq(ext[1], ext[2], by = grid_x),
          y = seq(ext[3], ext[4], by = grid_y),
          projection = "EPSG:32719",
          labels.size = 0.7,
          labels.format = list(fun = function(x) format(x, scientific = FALSE, big.mark = ",")),
          labels.inside.frame = FALSE,
          lines = TRUE
        )
      }
      
      # Füge Maßstab hinzu
      if (add_scale) {
        tm <- tm + tmap::tm_scale_bar(
          position = c("left", "bottom"),
          width = 0.25,
          text.size = 0.8,
          breaks = c(0, 0.5, 1) * 5 * grid_x
        )
      }
      
      # Füge Nordpfeil hinzu
      if (add_north) {
        north_arrow_file <- create_scientific_north_arrow(style = "compass")
        
        tm <- tm + tmap::tm_compass(
          type = "arrow",
          position = c("right", "top"),
          size = 2.0
        )
      }
      
      # Füge Probenahmepunkte hinzu
      if (!is.null(sample_points) && nrow(sample_points) > 0) {
        # Konvertiere zu sf-Punkten
        if (requireNamespace("sf", quietly = TRUE)) {
          sf_points <- sf::st_as_sf(sample_points, coords = c("x", "y"), crs = 32719)
          
          tm <- tm + tmap::tm_shape(sf_points) +
            tmap::tm_dots(size = 0.1, col = "red", border.col = "black", 
                         border.lwd = 0.8, alpha = 0.7)
        }
      }
      
      # Speichere die Karte
      cat(sprintf("Speichere tmap-Karte: %s\n", basename(output_file)))
      tmap::tmap_save(tm, output_file, width = width, height = height, dpi = dpi)
      
      # Erfolgsmeldung
      if (file.exists(output_file)) {
        end_time <- Sys.time()
        duration <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
        cat(sprintf("✓ tmap-Karte erstellt: %s (%.2f MB, %.2f Sek.)\n", 
                   basename(output_file), file.size(output_file)/1024/1024, duration))
        return(TRUE)
      } else {
        warning("tmap_save scheint erfolgreich, aber die Datei wurde nicht erstellt")
      }
    }, error = function(e) {
      cat("⚠️ tmap fehlgeschlagen:", e$message, "\n")
      cat("Fallback zu terra/raster plot...\n")
      # Weiter mit Fallback-Methode
    })
  }
  
  # ----- Fallback: Kartenerstellung mit base R und terra/raster -----
  cat("Erstelle Karte mit terra/raster...\n")
  
  tryCatch({
    # Öffne Ausgabedatei
    png(output_file, width = width, height = height, res = dpi, bg = "white")
    
    # Layout mit Platz für Titel und Legende
    layout(matrix(c(1, 2), nrow = 1), widths = c(0.8, 0.2))
    
    # Hauptkartenfenster mit Rand
    par(mar = c(4, 4, 6, 1))
    
    # Plotte RGB-Bild
    if (class(rgb_subset)[1] == "SpatRaster") {
      terra::plotRGB(rgb_subset, r = 1, g = 2, b = 3, 
                    stretch = ifelse(stretch == "hist", "hist", "lin"), 
                    axes = TRUE, main = "", las = 1)
    } else {
      raster::plotRGB(rgb_subset, r = 1, g = 2, b = 3, 
                     stretch = ifelse(stretch == "hist", "hist", "lin"), 
                     axes = TRUE, main = "", las = 1)
    }
    
    # Titel mit mehreren Zeilen
    mtext(title, side = 3, line = 4, cex = 1.5, font = 2)
    mtext(subtitle, side = 3, line = 2, cex = 1.0)
    
    # Koordinatengitter
    if (add_grid) {
      if (class(rgb_subset)[1] == "SpatRaster") {
        ext <- terra::ext(rgb_subset)
      } else {
        ext <- raster::extent(rgb_subset)
      }
      
      # Berechne sinnvolle Grid-Intervalle
      x_range <- ext[2] - ext[1]
      y_range <- ext[4] - ext[3]
      
      x_magnitude <- 10^floor(log10(x_range/5))
      y_magnitude <- 10^floor(log10(y_range/5))
      
      x_candidates <- c(1, 2, 5) * x_magnitude
      y_candidates <- c(1, 2, 5) * y_magnitude
      
      grid_x <- x_candidates[which.min(abs(x_candidates - x_range/5))]
      grid_y <- y_candidates[which.min(abs(y_candidates - y_range/5))]
      
      # Erzeuge Gitterlinien
      x_seq <- seq(ext[1], ext[2], by = grid_x)
      y_seq <- seq(ext[3], ext[4], by = grid_y)
      
      # Zeichne Gitterlinien
      abline(v = x_seq, col = "gray70", lty = 3, lwd = 0.8)
      abline(h = y_seq, col = "gray70", lty = 3, lwd = 0.8)
      
      # Koordinatenbeschriftungen verbessern
      axis(1, at = x_seq, labels = format(x_seq, scientific = FALSE, big.mark = ","), 
           las = 1, cex.axis = 0.8, tck = -0.01, col.tick = "gray50")
      axis(2, at = y_seq, labels = format(y_seq, scientific = FALSE, big.mark = ","), 
           las = 1, cex.axis = 0.8, tck = -0.01, col.tick = "gray50")
      
      # Koordinatensystem-Info
      mtext("UTM Zone 19S (EPSG:32719)", side = 1, line = 2.5, cex = 0.7, adj = 1)
    }
    
    # Probenahmepunkte
    if (!is.null(sample_points) && nrow(sample_points) > 0) {
      points(sample_points$x, sample_points$y, pch = 21, bg = "red", 
             col = "black", cex = 0.8, lwd = 1.2)
    }
    
    # Maßstabsleiste
    if (add_scale) {
      if (class(rgb_subset)[1] == "SpatRaster") {
        ext <- terra::ext(rgb_subset)
      } else {
        ext <- raster::extent(rgb_subset)
      }
      
      # Bessere Maßstabsberechnung
      width_meters <- ext[2] - ext[1]
      scale_length <- round(width_meters/5, -floor(log10(width_meters/5)))
      
      # Position links unten mit Abstand zum Rand
      x_start <- ext[1] + (ext[2] - ext[1]) * 0.05
      y_pos <- ext[3] + (ext[4] - ext[3]) * 0.05
      
      # Zeichne Maßstabsleiste mit Hintergrund
      rect(x_start - width_meters*0.005, y_pos - (ext[4] - ext[3])*0.008, 
           x_start + scale_length + width_meters*0.005, y_pos + (ext[4] - ext[3])*0.015, 
           col = "white", border = NA)
      
      # Maßstabsleiste
      segments(x_start, y_pos, x_start + scale_length, y_pos, lwd = 3, col = "black")
      segments(x_start, y_pos - width_meters*0.003, x_start, y_pos + width_meters*0.003, lwd = 2)
      segments(x_start + scale_length, y_pos - width_meters*0.003, 
               x_start + scale_length, y_pos + width_meters*0.003, lwd = 2)
      
      # Beschriftung
      if (scale_length >= 1000) {
        scale_text <- paste(scale_length/1000, "km")
      } else {
        scale_text <- paste(scale_length, "m")
      }
      
      text(x_start + scale_length/2, y_pos + width_meters*0.01, scale_text, cex = 0.9)
    }
    
    # Nordpfeil
    if (add_north) {
      if (class(rgb_subset)[1] == "SpatRaster") {
        ext <- terra::ext(rgb_subset)
      } else {
        ext <- raster::extent(rgb_subset)
      }
      
      # Position rechts oben mit Abstand zum Rand
      x_arrow <- ext[2] - (ext[2] - ext[1]) * 0.08
      y_arrow <- ext[4] - (ext[4] - ext[3]) * 0.12
      arrow_length <- (ext[4] - ext[3]) * 0.06
      
      # Hintergrund für bessere Sichtbarkeit
      symbols(x_arrow, y_arrow, circles = arrow_length * 1.2, 
              inches = FALSE, add = TRUE, bg = "#FFFFFF99", fg = NA)
      
      # Nordpfeil mit Schatten
      arrows(x_arrow, y_arrow, x_arrow, y_arrow + arrow_length, 
             length = 0.12, angle = 20, code = 2, lwd = 3, col = "red")
      
      # N-Beschriftung
      text(x_arrow, y_arrow + arrow_length * 1.3, "N", cex = 1.3, font = 2, col = "red")
    }
    
    # Legendenfenster
    par(mar = c(3, 0, 6, 3))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    
    # Legendenüberschrift
    if (is.null(rgb_bands)) {
      legend_title <- "True Color"
      red_text <- "Rot (~650nm)"
      green_text <- "Grün (~550nm)"
      blue_text <- "Blau (~450nm)"
    } else {
      legend_title <- "RGB-Bänder"
      red_text <- sprintf("Rot: Band %d", rgb_bands[1])
      green_text <- sprintf("Grün: Band %d", rgb_bands[2])
      blue_text <- sprintf("Blau: Band %d", rgb_bands[3])
    }
    
    text(0.5, 0.95, legend_title, font = 2, cex = 1.2)
    
    # Farbboxen
    rect(0.15, 0.80, 0.35, 0.88, col = "red", border = "black", lwd = 1.5)
    rect(0.15, 0.65, 0.35, 0.73, col = "green", border = "black", lwd = 1.5)
    rect(0.15, 0.50, 0.35, 0.58, col = "blue", border = "black", lwd = 1.5)
    
    # Beschriftungen
    text(0.4, 0.84, red_text, adj = 0, cex = 1.0)
    text(0.4, 0.69, green_text, adj = 0, cex = 1.0)
    text(0.4, 0.54, blue_text, adj = 0, cex = 1.0)
    
    # Bildinfo
    text(0.5, 0.35, "Bildinfo:", font = 3, cex = 1.0)
    
    img_info <- sprintf(
      "Raster: %d×%d Pixel\nBänder: %d\nAuflösung: %.1fm",
      n_rows, n_cols, n_bands,
      if(class(rgb_subset)[1] == "SpatRaster") terra::res(rgb_subset)[1] else raster::res(rgb_subset)[1]
    )
    
    text(0.5, 0.25, img_info, cex = 0.8, adj = 0.5)
    
    # Metadaten am unteren Rand
    text(0.5, 0.06, paste("Erstellt am:", format(Sys.time(), "%d.%m.%Y")), cex = 0.7)
    
    # Schließe die Ausgabedatei
    dev.off()
    
    # Erfolgsmeldung
    if (file.exists(output_file)) {
      end_time <- Sys.time()
      duration <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
      cat(sprintf("✓ terra/raster-Karte erstellt: %s (%.2f MB, %.2f Sek.)\n", 
                 basename(output_file), file.size(output_file)/1024/1024, duration))
      return(TRUE)
    } else {
      warning("Kartenerstellung scheint erfolgreich, aber die Datei wurde nicht erstellt")
      return(FALSE)
    }
  }, error = function(e) {
    cat("❌ Fehler bei Kartenerstellung:", e$message, "\n")
    return(FALSE)
  })
}

# Erstellt alle RGB-Visualisierungen für ein Raster mit wissenschaftlichem Anspruch
create_scientific_rgb_visualizations <- function(
  raster_file,                 # Pfad zum geglätteten Raster
  output_dir = NULL,           # Ausgabeverzeichnis (wenn NULL, wird automatisch erstellt)
  rgb_combinations = list(     # Verschiedene RGB-Kombinationen zur Erstellung
    "TRUE-COLOR" = NULL,       # True-Color (automatisch)
    "NIR-G-B" = c(40, 20, 10),
    "SWIR-NIR-R" = c(100, 40, 30)
  ),
  width = 2400,                # Bildbreite für wissenschaftliche Qualität
  height = 1800,               # Bildhöhe
  dpi = 300,                   # DPI für Druckqualität
  add_points = FALSE,          # Probenahmepunkte hinzufügen?
  points_file = "data/sampling_points.csv"
) {
  # Start-Nachricht
  cat(sprintf("\n===== Wissenschaftliche RGB-Visualisierungen für %s =====\n", 
             basename(raster_file)))
  
  # Starte Timer
  start_time <- Sys.time()
  
  # Erstelle Ausgabeverzeichnis
  image_name <- gsub("_smoothed\\+.*$|_smoothed$", "", basename(raster_file))
  
  if (is.null(output_dir)) {
    output_dir <- file.path("results/rgb_visualizations", image_name)
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("Ausgabeverzeichnis erstellt: %s\n", output_dir))
  }
  
  # Erstelle alle RGB-Kombinationen
  results <- list()
  
  for (combo_name in names(rgb_combinations)) {
    rgb_bands <- rgb_combinations[[combo_name]]
    
    # Ausgabedateiname
    output_file <- file.path(output_dir, paste0(image_name, "_", combo_name, ".png"))
    
    # Setze passenden Titel
    if (is.null(rgb_bands)) {
      # True-Color
      title <- "True-Color-Darstellung"
    } else {
      # RGB-Kombination
      title <- sprintf("RGB-Kombination: %s", combo_name)
    }
    
    # Erstelle die Karte
    cat(sprintf("\nErstelle %s Visualisierung...\n", combo_name))
    
    result <- create_scientific_rgb_map(
      raster_file = raster_file,
      output_file = output_file,
      rgb_bands = rgb_bands,
      title = title,
      subtitle = paste("Bild:", image_name),
      width = width,
      height = height,
      dpi = dpi,
      add_points = add_points,
      points_file = points_file
    )
    
    results[[combo_name]] <- list(
      output_file = output_file,
      success = result
    )
  }
  
  # Erfolgsmeldung
  end_time <- Sys.time()
  duration <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
  
  successful <- sum(sapply(results, function(r) r$success))
  cat(sprintf("\n✓ %d von %d RGB-Visualisierungen erstellt (%.2f Sek.)\n", 
             successful, length(results), duration))
  
  return(invisible(results))
}

# Batch-Verarbeitung aller Raster-Dateien
process_all_rasters <- function(
  input_dir = "data",
  file_pattern = "_smoothed\\+adaptive_destriped_masked_snv_minmax\\.tif$",
  output_dir = "results/rgb_visualizations",
  rgb_combinations = list(
    "TRUE-COLOR" = NULL,       # True-Color (automatisch)
    "NIR-G-B" = c(40, 20, 10)  # NIR-Grün-Blau
  ),
  recursive_search = TRUE
) {
  # Start-Nachricht
  cat("\n===== BATCH-VERARBEITUNG ALLER RASTER-DATEIEN =====\n")
  
  # Starte Timer
  start_time <- Sys.time()
  
  # Suche nach Rasterdateien
  cat(sprintf("Suche nach Rasterdateien mit Muster '%s'...\n", file_pattern))
  
  raster_files <- list.files(
    input_dir,
    pattern = file_pattern,
    full.names = TRUE,
    recursive = recursive_search
  )
  
  if (length(raster_files) == 0) {
    cat("⚠️ Keine Rasterdateien gefunden! Alternative Muster probieren...\n")
    
    # Alternative Muster probieren
    alt_pattern <- "_smoothed\\.tif$"
    cat(sprintf("Versuche alternatives Muster '%s'...\n", alt_pattern))
    
    raster_files <- list.files(
      input_dir,
      pattern = alt_pattern,
      full.names = TRUE,
      recursive = recursive_search
    )
    
    if (length(raster_files) == 0) {
      cat("❌ Keine Rasterdaten gefunden. Bitte Verzeichnis und Dateinamen überprüfen.\n")
      return(invisible(NULL))
    }
  }
  
  cat(sprintf("✓ %d Rasterdateien gefunden\n", length(raster_files)))
  for (i in 1:min(5, length(raster_files))) {
    cat(sprintf("  %d. %s\n", i, basename(raster_files[i])))
  }
  if (length(raster_files) > 5) {
    cat(sprintf("  ... und %d weitere\n", length(raster_files) - 5))
  }
  
  # Frage nach Bestätigung, wenn viele Dateien
  if (length(raster_files) > 10) {
    cat(sprintf("\n⚠️ Es wurden %d Dateien gefunden. Verarbeitung könnte lange dauern.\n", 
               length(raster_files)))
    cat("Fortfahren? (y/n): ")
    answer <- readLines(n=1)
    
    if (tolower(substr(answer, 1, 1)) != "y") {
      cat("Abgebrochen.\n")
      return(invisible(NULL))
    }
  }
  
  # Erstelle Ausgabeverzeichnis
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("Ausgabeverzeichnis erstellt: %s\n", output_dir))
  }
  
  # Verarbeite alle Rasterdateien
  results <- list()
  
  for (i in 1:length(raster_files)) {
    raster_file <- raster_files[i]
    image_name <- gsub("_smoothed\\+.*$|_smoothed$", "", basename(raster_file))
    
    cat(sprintf("\n[%d/%d] Verarbeite: %s\n", i, length(raster_files), image_name))
    
    # Bildspezifisches Ausgabeverzeichnis
    img_output_dir <- file.path(output_dir, image_name)
    
    # Erstelle Visualisierungen
    result <- create_scientific_rgb_visualizations(
      raster_file = raster_file,
      output_dir = img_output_dir,
      rgb_combinations = rgb_combinations,
      add_points = FALSE  # In Batch-Modus keine Punkte hinzufügen
    )
    
    results[[image_name]] <- result
    
    # Kurze Pause
    Sys.sleep(0.5)
    gc()  # Speicher freigeben
  }
  
  # Erfolgsmeldung
  end_time <- Sys.time()
  total_duration <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)
  
  cat(sprintf("\n✓ Batch-Verarbeitung abgeschlossen!\n"))
  cat(sprintf("  %d Bilder verarbeitet in %.2f Minuten\n", 
             length(raster_files), total_duration))
  cat(sprintf("  Ergebnisse in: %s\n\n", output_dir))
  
  return(invisible(results))
}

# ========================= KOPIEREN DER BILDER IN ERGEBNISVERZEICHNISSE =========================

# Funktion zum Kopieren der RGB-Bilder in die raster_prediction Verzeichnisse
copy_rgb_to_prediction_dirs <- function(
  source_dir = "results/rgb_visualizations",
  results_dir = "results",
  algorithms = c("Pearson", "Spearman"),
  plsr_types = c("filter+plsr")
) {
  cat("\n===== KOPIERE RGB-BILDER IN ERGEBNISVERZEICHNISSE =====\n")
  
  # Finde alle RGB-Bilder
  rgb_files <- list.files(source_dir, pattern = "\\.png$", 
                         full.names = TRUE, recursive = TRUE)
  
  if (length(rgb_files) == 0) {
    cat("Keine RGB-Bilder gefunden in:", source_dir, "\n")
    return(invisible(NULL))
  }
  
  cat(sprintf("Gefunden: %d RGB-Bilder\n", length(rgb_files)))
  
  # Finde alle raster_prediction Verzeichnisse
  target_dirs <- c()
  for (algo in algorithms) {
    for (plsr in plsr_types) {
      pred_dir <- file.path(results_dir, algo, plsr, "raster_predictions")
      if (dir.exists(pred_dir)) {
        # Finde alle Bildverzeichnisse darin
        image_dirs <- list.dirs(pred_dir, recursive = FALSE)
        target_dirs <- c(target_dirs, image_dirs)
      }
    }
  }
  
  if (length(target_dirs) == 0) {
    cat("Keine Zielverzeichnisse für raster_predictions gefunden.\n")
    return(invisible(NULL))
  }
  
  cat(sprintf("Gefunden: %d Zielverzeichnisse\n", length(target_dirs)))
  
  # Zuordnungstabelle für Bilder und Zielverzeichnisse
  copy_count <- 0
  
  # Für jedes RGB-Bild
  for (rgb_file in rgb_files) {
    # Extrahiere Bildname aus Pfad
    # Annahme: Format ist results/rgb_visualizations/BILDNAME/BILDNAME_RGB-TYP.png
    parts <- strsplit(rgb_file, "/")[[1]]
    if (length(parts) >= 3) {
      image_name <- parts[length(parts) - 1]  # Vorletzte Komponente ist Bildname
      
      # Finde passende Zielverzeichnisse für dieses Bild
      matching_dirs <- target_dirs[grep(image_name, target_dirs, fixed = TRUE)]
      
      if (length(matching_dirs) > 0) {
        cat(sprintf("\nKopiere %s nach %d Verzeichnisse(n)\n", 
                   basename(rgb_file), length(matching_dirs)))
        
        for (target_dir in matching_dirs) {
          # Ziel-RGB-Datei
          target_file <- file.path(target_dir, paste0("RGB_", basename(rgb_file)))
          
          # Kopiere die Datei
          tryCatch({
            file.copy(rgb_file, target_file, overwrite = TRUE)
            cat(sprintf("  ✓ Kopiert nach: %s\n", gsub(results_dir, "results/...", target_file)))
            copy_count <- copy_count + 1
            
            # Kopiere auch in Parameter-Unterverzeichnisse
            param_dirs <- list.dirs(target_dir, recursive = FALSE)
            if (length(param_dirs) > 0) {
              for (param_dir in param_dirs) {
                param_target <- file.path(param_dir, paste0("RGB_", basename(rgb_file)))
                file.copy(rgb_file, param_target, overwrite = TRUE)
                copy_count <- copy_count + 1
              }
            }
          }, error = function(e) {
            cat(sprintf("  ❌ Fehler beim Kopieren: %s\n", e$message))
          })
        }
      }
    }
  }
  
  cat(sprintf("\n✓ %d Dateien insgesamt kopiert\n", copy_count))
  return(invisible(copy_count))
}

# ========================= HAUPTFUNKTION =========================

# Hauptfunktion
main <- function() {
  cat("\n===== WISSENSCHAFTLICHE RGB-VISUALISIERUNG =====\n")
  
  # Finde und verwende die TRUE Color Bänder für EnMAP
  tryCatch({
    # Berechne die True-Color-Bänder automatisch aus der Wellenlängentabelle
    wavelength_file <- "data/wavelength_reference.csv"
    
    # Finde die passenden Bänder direkt im Verarbeitungscode
    cat("Automatische RGB-Band-Auswahl aktiviert\n")
    
    # Erstelle RGB-Visualisierungen mit verbesserten wissenschaftlichen Standards
    process_all_rasters(
      input_dir = "data",
      file_pattern = "_smoothed\\+adaptive_destriped_masked_snv_minmax\\.tif$",
      output_dir = "results/rgb_visualizations",
      rgb_combinations = list(
        "TRUE-COLOR" = NULL,          # Automatische True-Color-Erkennung
        "NIR-G-B" = c(40, 20, 10)     # NIR-Grün-Blau für Vegetation
      ),
      recursive_search = TRUE
    )
    
    # Kopiere die Ergebnisse in die entsprechenden Vorhersageverzeichnisse
    copy_rgb_to_prediction_dirs(
      source_dir = "results/rgb_visualizations",
      results_dir = "results",
      algorithms = c("Pearson", "Spearman"),
      plsr_types = c("filter+plsr")
    )
    
  }, error = function(e) {
    cat("⚠️ Fehler in der Hauptfunktion:", e$message, "\n")
  })
  
  cat("\nWissenschaftliche RGB-Visualisierung abgeschlossen.\n")
}

# HINWEIS: Auto-Execution auskommentiert - Script wird nur über master_pipeline.R gesteuert
# # Führe das Hauptprogramm aus
# main()
