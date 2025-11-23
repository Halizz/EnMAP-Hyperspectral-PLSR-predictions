################################################################################
#                                                                              #
#                    ENMAP HYPERSPEKTRAL-PIPELINE                              #
#                    Vollständige Datenverarbeitung                            #
#                                                                              #
################################################################################
#
# Beschreibung:
# Dieses Master-Script führt die komplette Pipeline der EnMAP-Hyperspektral-
# Datenverarbeitung durch - von der Extraktion der Satellitenbilder bis zur
# Erstellung von Vorhersagekarten für Bodenparameter.
#
# Pipeline-Schritte:
#  1. Extraktion der EnMAP .tar.gz Archive
#  2. Metadaten-Extraktion & Koordinatenabgleich
#  3. Savitzky-Golay Glättung (SG-Filter) & Destriping (Horn & Woodham)
#  4. Maskierung (Wolken, Schatten, Schnee, etc.)
#  5. Min-Max-Normalisierung oder SNV-Normalisierung
#  6. Koordinaten angleichen & Spektraldaten extrahieren
#  7. Korrelationsanalyse (Pearson & Spearman mit RSS 50:50)
#  8. Signifikanztest (Multiple Test Correction - Benjamini-Hochberg)
#  9. Visualisierung (optional)
# 10. Filter & PLSR-Training (RSS 70:30, 5-fold Kreuzvalidierung)
# 11. Räumliche Vorhersagekarten erstellen
#
# Autor: [Automatisch generiert]
# Datum: 2025-11-22
################################################################################

# ══════════════════════════════════════════════════════════════════════════════
# 📋 KONFIGURATION - HIER PFADE ANPASSEN
# ══════════════════════════════════════════════════════════════════════════════

# Hauptverzeichnis mit allen Daten
WORKING_DIR <- "C:/Users/leoha/Desktop/scripts - Kopie"

# Input-Pfade
ENMAP_TAR_DIR <- file.path(WORKING_DIR, "data/ENMAP")              # Ordner mit .tar.gz Dateien
SHAPEFILE_PATH <- file.path(WORKING_DIR, "data/Peru_Bolivien_reprojiziert.shp")   # Shapefile mit Bodenprobenpunkten
SOIL_DATA_DIR <- file.path(WORKING_DIR, "data/soil_measurements")  # Bodenmessungsdaten

# Output-Pfade
EXTRACT_DIR <- file.path(WORKING_DIR, "data/extracted_files")      # Extrahierte Bilder
RESULTS_DIR <- file.path(WORKING_DIR, "results")                   # Endergebnisse
SCRIPTS_DIR <- file.path(WORKING_DIR, "scripts")                   # R-Scripts

# Normalisierungsmethode: "minmax" oder "snv"
NORMALIZATION_METHOD <- "snv"

# Seeds für Reproduzierbarkeit (gemäß Methodik)
SEED_CORRELATION_SPLIT <- 358  # RSS 50:50 für Korrelationsanalyse
SEED_PLSR_SPLIT <- 696         # RSS 70:30 für PLSR Training
SEED_PLSR_CV <- 263            # 5-fold Kreuzvalidierung für PLSR

# Savitzky-Golay Parameter
SG_WINDOW_SIZE <- 11
SG_POLY_ORDER <- 3

# Pipeline-Optionen
SKIP_EXTRACTION <- FALSE       # TRUE falls Bilder bereits extrahiert sind
SKIP_VISUALIZATION <- FALSE    # TRUE um Visualisierungen zu überspringen
CREATE_DETAILED_VIZ <- TRUE    # TRUE um detaillierte Visualisierungen erstellen

# ══════════════════════════════════════════════════════════════════════════════
# 📦 PAKETE LADEN
# ══════════════════════════════════════════════════════════════════════════════
# Setze Working Directory
setwd(WORKING_DIR)
cat("📂 Working Directory:", getwd(), "\n\n")

cat("\n")
cat("════════════════════════════════════════════════════════════════════════════════\n")
cat("                        ENMAP HYPERSPEKTRAL-PIPELINE                            \n")
cat("════════════════════════════════════════════════════════════════════════════════\n")
cat("\n")

cat("🔧 Lade erforderliche R-Pakete...\n")

required_packages <- c(
  "terra", "raster", "sf", "dplyr", "readr", "readxl", "tidyr", "xml2",
  "data.table", "tools", "pracma", "pbapply", "caret", "pls", "glmnet",
  "ggplot2", "scales", "reshape2", "tmap", "RColorBrewer", "viridis",
  "vegan", "matrixStats", "parallel"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("   📥 Installiere Paket:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

cat("✅ Alle Pakete erfolgreich geladen\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# 🛠️ HILFSFUNKTIONEN
# ══════════════════════════════════════════════════════════════════════════════

#' Berechnet die geschätzte Verarbeitungszeit basierend auf Datengröße
#'
#' @param data_size_gb Größe der zu verarbeitenden Daten in GB
#' @param processing_type Art der Verarbeitung ("extract", "smooth", "normalize", etc.)
#' @return Geschätzte Zeit als lesbarer String
estimate_processing_time <- function(data_size_gb, processing_type) {
  # Geschätzte Verarbeitungszeiten pro GB (in Minuten)
  time_per_gb <- list(
    extract = 2,
    smooth = 8,
    normalize = 3,
    mask = 2,
    data_prep = 5,
    correlation = 1,
    significance = 0.5,
    visualization = 2,
    plsr = 3,
    prediction = 10
  )
  
  base_time <- time_per_gb[[processing_type]] * data_size_gb
  
  if (base_time < 1) {
    return(sprintf("< 1 Minute"))
  } else if (base_time < 60) {
    return(sprintf("ca. %.0f Minuten", base_time))
  } else {
    hours <- floor(base_time / 60)
    minutes <- round(base_time %% 60)
    return(sprintf("ca. %d Std. %d Min.", hours, minutes))
  }
}

#' Berechnet die Größe eines Verzeichnisses in GB
#'
#' @param dir_path Pfad zum Verzeichnis
#' @return Größe in GB
get_directory_size_gb <- function(dir_path) {
  if (!dir.exists(dir_path)) return(0)
  
  files <- list.files(dir_path, recursive = TRUE, full.names = TRUE)
  total_bytes <- sum(file.info(files)$size, na.rm = TRUE)
  return(total_bytes / (1024^3))
}

#' Prüft ob Ausgabedateien bereits existieren
#'
#' @param check_paths Vektor mit Pfaden oder Patterns die geprüft werden sollen
#' @param pattern_search Falls TRUE, werden Patterns als Dateisuche verwendet
#' @return TRUE falls Dateien gefunden wurden, FALSE sonst
check_output_exists <- function(check_paths, pattern_search = FALSE) {
  if (pattern_search) {
    for (path_pattern in check_paths) {
      # Extrahiere Verzeichnis und Pattern
      # Bei Patterns wie "data/extracted_files/*_smoothed.tif"
      # oder "data/extracted_files/subdir/*_smoothed.tif"
      
      # Finde das letzte Verzeichnis ohne Wildcards
      path_parts <- strsplit(path_pattern, "/")[[1]]
      base_dir_parts <- character(0)
      pattern_part <- NULL
      
      for (i in seq_along(path_parts)) {
        if (grepl("\\*|\\?", path_parts[i])) {
          # Erstes Teil mit Wildcard gefunden
          pattern_part <- paste(path_parts[i:length(path_parts)], collapse = "/")
          break
        } else {
          base_dir_parts <- c(base_dir_parts, path_parts[i])
        }
      }
      
      if (is.null(pattern_part)) {
        # Kein Wildcard gefunden, nutze als normalen Pfad
        pattern_part <- basename(path_pattern)
        base_dir <- dirname(path_pattern)
      } else {
        base_dir <- paste(base_dir_parts, collapse = "/")
        if (base_dir == "") base_dir <- "."
      }
      
      # Prüfe ob Verzeichnis existiert
      if (!dir.exists(base_dir)) next
      
      # Konvertiere Pattern für regulären Ausdruck
      regex_pattern <- gsub("\\.", "\\\\.", pattern_part)
      regex_pattern <- gsub("\\*", ".*", regex_pattern)
      regex_pattern <- gsub("\\?", ".", regex_pattern)
      
      # Suche Dateien
      files <- list.files(base_dir, 
                         pattern = regex_pattern, 
                         full.names = TRUE,
                         recursive = TRUE,
                         ignore.case = TRUE)
      
      if (length(files) > 0) return(TRUE)
    }
    return(FALSE)
  } else {
    for (path in check_paths) {
      if (file.exists(path) || dir.exists(path)) return(TRUE)
    }
    return(FALSE)
  }
}

#' Fordert Benutzerbestätigung an bevor ein Verarbeitungsschritt ausgeführt wird
#'
#' @param step_number Schrittnummer
#' @param step_name Name des Verarbeitungsschritts
#' @param description Detaillierte Beschreibung
#' @param estimated_time Geschätzte Verarbeitungszeit
#' @param check_existing_output Pfade oder Patterns die auf existierende Ausgaben prüfen
#' @return Liste mit skip (TRUE/FALSE) und overwrite (TRUE/FALSE)
ask_user_confirmation <- function(step_number, step_name, description, estimated_time, 
                                 check_existing_output = NULL) {
  cat("\n")
  cat("────────────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("📍 SCHRITT %d: %s\n", step_number, toupper(step_name)))
  cat("────────────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("   %s\n", description))
  cat(sprintf("   ⏱️  Geschätzte Dauer: %s\n", estimated_time))
  
  # Prüfe ob Ausgaben bereits existieren
  outputs_exist <- FALSE
  if (!is.null(check_existing_output)) {
    pattern_search <- attr(check_existing_output, "pattern_search") %||% FALSE
    outputs_exist <- check_output_exists(check_existing_output, pattern_search)
  }
  
  if (outputs_exist) {
    cat("   ⚠️  Ergebnisse bereits vorhanden!\n")
    cat("────────────────────────────────────────────────────────────────────────────────\n")
    response <- readline(prompt = "   🔄 Schritt wiederholen? (j=Ja/überschreiben, n=Nein/vorhandene nutzen): ")
    
    if (tolower(response) == "n") {
      cat("   ⏭️  Überspringe Schritt, nutze vorhandene Ergebnisse.\n\n")
      return(list(skip = TRUE, overwrite = FALSE))
    } else if (tolower(response) == "j") {
      cat("   ♻️  Überschreibe vorhandene Ergebnisse...\n\n")
      return(list(skip = FALSE, overwrite = TRUE))
    } else {
      cat("   ⏸️  Ungültige Eingabe. Pipeline angehalten.\n\n")
      return(list(skip = TRUE, overwrite = FALSE, abort = TRUE))
    }
  } else {
    cat("────────────────────────────────────────────────────────────────────────────────\n")
    response <- readline(prompt = "   ▶️  Fortfahren? (j/n): ")
    
    if (tolower(response) != "j") {
      cat("   ⏸️  Pipeline angehalten. Starte erneut um fortzufahren.\n\n")
      return(list(skip = TRUE, overwrite = FALSE, abort = TRUE))
    }
    
    cat("\n")
    return(list(skip = FALSE, overwrite = FALSE))
  }
}

#' Logging-Funktion mit Zeitstempel
#'
#' @param message Nachricht die geloggt werden soll
#' @param level Log-Level ("INFO", "WARNING", "ERROR", "SUCCESS")
log_message <- function(message, level = "INFO") {
  timestamp <- format(Sys.time(), "%H:%M:%S")
  
  prefix <- switch(level,
                   "INFO" = "ℹ️ ",
                   "WARNING" = "⚠️ ",
                   "ERROR" = "❌",
                   "SUCCESS" = "✅",
                   "  ")
  
  cat(sprintf("[%s] %s%s\n", timestamp, prefix, message))
}

# ══════════════════════════════════════════════════════════════════════════════
# 📥 SCHRITT 1: EXTRAKTION DER ENMAP-BILDER
# ══════════════════════════════════════════════════════════════════════════════

step1_extract_images <- function() {
  if (SKIP_EXTRACTION) {
    log_message("Überspringe Extraktion (SKIP_EXTRACTION = TRUE)", "INFO")
    return(TRUE)
  }
  
  # Berechne Datengröße
  tar_size_gb <- get_directory_size_gb(ENMAP_TAR_DIR)
  estimated_time <- estimate_processing_time(tar_size_gb, "extract")
  
  # Prüfe auf existierende Ausgaben
  check_outputs <- structure(
    c(file.path(EXTRACT_DIR, "*SPECTRAL_IMAGE*")),
    pattern_search = TRUE
  )
  
  # Benutzerbestätigung
  confirmation <- ask_user_confirmation(
    step_number = 1,
    step_name = "Extraktion der EnMAP .tar.gz Archive",
    description = "Extrahiert alle .tar.gz Dateien in den Zielordner.",
    estimated_time = estimated_time,
    check_existing_output = check_outputs
  )
  
  if (confirmation$skip) {
    if (!is.null(confirmation$abort) && confirmation$abort) return(FALSE)
    return(TRUE)
  }
  
  # Quelle das extract_data.r Script
  log_message("Starte Extraktion der EnMAP-Archive...", "INFO")
  
  source(file.path(SCRIPTS_DIR, "extract_data.r"), local = TRUE)
  
  # Führe Extraktion durch
  tryCatch({
    run_extraction(
      tar_dir = ENMAP_TAR_DIR,
      extract_dir = EXTRACT_DIR,
      pause_duration = 5
    )
    log_message("Extraktion erfolgreich abgeschlossen", "SUCCESS")
    return(TRUE)
  }, error = function(e) {
    log_message(sprintf("Fehler bei Extraktion: %s", e$message), "ERROR")
    return(FALSE)
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# 📊 SCHRITT 2: METADATEN-EXTRAKTION & KOORDINATENABGLEICH
# ══════════════════════════════════════════════════════════════════════════════

step2_extract_metadata <- function() {
  # Berechne Datengröße
  data_size_gb <- get_directory_size_gb(EXTRACT_DIR)
  estimated_time <- estimate_processing_time(data_size_gb, "metadata")
  
  # Benutzerbestätigung
  if (!ask_user_confirmation(
    step_number = 2,
    step_name = "Metadaten-Extraktion & Koordinatenabgleich",
    description = "Extrahiert Metadaten aus XML-Dateien (Sonnenhöhenwinkel, etc.) und gleicht mit Shapefile-Koordinaten ab.",
    estimated_time = estimated_time
  )) return(FALSE)
  
  log_message("Starte Metadaten-Extraktion...", "INFO")
  
  # Lade das extract_metadata Script
  source(file.path(SCRIPTS_DIR, "extract_metadata.R"), local = TRUE)
  
  tryCatch({
    extract_metadata_mapping(
      extract_dir = EXTRACT_DIR,
      shapefile_path = SHAPEFILE_PATH,
      output_file = "data/metadata_mapping.csv"
    )
    log_message("Metadaten-Extraktion erfolgreich abgeschlossen", "SUCCESS")
    return(TRUE)
  }, error = function(e) {
    log_message(sprintf("Fehler bei Metadaten-Extraktion: %s", e$message), "ERROR")
    return(FALSE)
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# 📊 SCHRITT 2: METADATEN-EXTRAKTION & KOORDINATENABGLEICH
# ══════════════════════════════════════════════════════════════════════════════

step2_extract_metadata <- function() {
  # Berechne Datengröße
  data_size_gb <- get_directory_size_gb(EXTRACT_DIR)
  estimated_time <- "< 1 Minute"
  
  # Prüfe auf existierende Ausgaben
  check_outputs <- c(
    "data/metadata_mapping.csv",
    "data/wavelength_reference.csv"
  )
  
  # Benutzerbestätigung
  confirmation <- ask_user_confirmation(
    step_number = 2,
    step_name = "Metadaten-Extraktion & Koordinatenabgleich",
    description = "Extrahiert Metadaten aus XML-Dateien (Sonnenhöhenwinkel, etc.) und gleicht mit Shapefile-Koordinaten ab.",
    estimated_time = estimated_time,
    check_existing_output = check_outputs
  )
  
  if (confirmation$skip) {
    if (!is.null(confirmation$abort) && confirmation$abort) return(FALSE)
    return(TRUE)
  }
  
  log_message("Starte Metadaten-Extraktion...", "INFO")
  
  # Lade das extract_metadata Script
  source(file.path(SCRIPTS_DIR, "extract_metadata.R"), local = TRUE)
  
  tryCatch({
    extract_metadata_mapping(
      extract_dir = EXTRACT_DIR,
      shapefile_path = SHAPEFILE_PATH,
      output_file = "data/metadata_mapping.csv"
    )
    log_message("Metadaten-Extraktion erfolgreich abgeschlossen", "SUCCESS")
    return(TRUE)
  }, error = function(e) {
    log_message(sprintf("Fehler bei Metadaten-Extraktion: %s", e$message), "ERROR")
    return(FALSE)
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# 🎛️ SCHRITT 3: SAVITZKY-GOLAY GLÄTTUNG & DESTRIPING
# ══════════════════════════════════════════════════════════════════════════════

step3_smooth_and_destripe <- function() {
  # Berechne Datengröße
  data_size_gb <- get_directory_size_gb(EXTRACT_DIR)
  estimated_time <- estimate_processing_time(data_size_gb, "smooth")
  
  # Prüfe auf existierende Ausgaben - suche nach geglätteten UND destriped Dateien
  # Pattern muss eindeutig smoothed+adaptive_destriped im Dateinamen enthalten
  smoothed_files <- list.files(EXTRACT_DIR, 
                               pattern = ".*_smoothed\\+adaptive_destriped\\.tif$",
                               full.names = TRUE, 
                               recursive = TRUE,
                               ignore.case = TRUE)
  
  outputs_exist <- length(smoothed_files) > 0
  
  # Benutzerbestätigung
  cat("\n")
  cat("────────────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("📍 SCHRITT 3: %s\n", toupper("Savitzky-Golay Glättung & Destriping")))
  cat("────────────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("   %s\n", sprintf(
    "Wendet SG-Filter (Fenster=%d, Ordnung=%d) an und führt Destriping durch (Horn & Woodham Methode).",
    SG_WINDOW_SIZE, SG_POLY_ORDER
  )))
  cat(sprintf("   ⏱️  Geschätzte Dauer: %s\n", estimated_time))
  
  if (outputs_exist) {
    cat(sprintf("   ⚠️  %d geglättete Datei(en) bereits vorhanden!\n", length(smoothed_files)))
    cat("────────────────────────────────────────────────────────────────────────────────\n")
    response <- readline(prompt = "   🔄 Schritt wiederholen? (j=Ja/überschreiben, n=Nein/vorhandene nutzen): ")
    
    if (tolower(response) == "n") {
      cat("   ⏭️  Überspringe Schritt, nutze vorhandene Ergebnisse.\n\n")
      return(TRUE)
    } else if (tolower(response) != "j") {
      cat("   ⏸️  Ungültige Eingabe. Pipeline angehalten.\n\n")
      return(FALSE)
    }
    cat("   ♻️  Überschreibe vorhandene Ergebnisse...\n\n")
  } else {
    cat("────────────────────────────────────────────────────────────────────────────────\n")
    response <- readline(prompt = "   ▶️  Fortfahren? (j/n): ")
    
    if (tolower(response) != "j") {
      cat("   ⏸️  Pipeline angehalten. Starte erneut um fortzufahren.\n\n")
      return(FALSE)
    }
    cat("\n")
  }
  
  log_message("Starte SG-Glättung und Destriping...", "INFO")
  
  # Lade das smooth_pictures Script
  source(file.path(SCRIPTS_DIR, "smooth_pictures_destripe_simple.R"), local = TRUE)
  
  tryCatch({
    smooth_spectral_images(
      spectral_dir = EXTRACT_DIR,
      output_suffix = "_smoothed+destriped",
      window_size = SG_WINDOW_SIZE,
      poly_order = SG_POLY_ORDER,
      chunk_size = 500,
      dev_mode = FALSE,
      overwrite = TRUE,
      apply_destriping = TRUE,
      destriping_preserve_statistics = FALSE,
      adaptive_destriping = TRUE,
      shapefile_path = SHAPEFILE_PATH  # Shapefile zum Filtern der Raster
    )
    log_message("Glättung und Destriping erfolgreich abgeschlossen", "SUCCESS")
    return(TRUE)
  }, error = function(e) {
    log_message(sprintf("Fehler bei Glättung/Destriping: %s", e$message), "ERROR")
    return(FALSE)
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# 🎭 SCHRITT 4: MASKIERUNG
# ══════════════════════════════════════════════════════════════════════════════

step4_mask_images <- function() {
  # Berechne Datengröße
  data_size_gb <- get_directory_size_gb(EXTRACT_DIR)
  estimated_time <- estimate_processing_time(data_size_gb, "mask")
  
  # Prüfe auf existierende Ausgaben
  check_outputs <- structure(
    c(file.path(EXTRACT_DIR, "*_masked*.tif")),
    pattern_search = TRUE
  )
  
  # Benutzerbestätigung
  confirmation <- ask_user_confirmation(
    step_number = 4,
    step_name = "Maskierung",
    description = "Maskiert Wolken, Wolkenschatten, Dunst, Schnee und Sensorfehler.",
    estimated_time = estimated_time,
    check_existing_output = check_outputs
  )
  
  if (confirmation$skip) {
    if (!is.null(confirmation$abort) && confirmation$abort) return(FALSE)
    return(TRUE)
  }
  
  log_message("Starte Maskierung...", "INFO")
  
  # Lade Maskierungs-Script
  source(file.path(SCRIPTS_DIR, "masking_pictures.R"), local = TRUE)
  
  tryCatch({
    mask_smoothed_images(
      smoothed_dir = EXTRACT_DIR,
      mask_dir = EXTRACT_DIR,
      mask_types = c("QL_QUALITY_CLOUD", "QL_QUALITY_CLOUDSHADOW", 
                    "QL_QUALITY_HAZE", "QL_QUALITY_CIRRUS", 
                    "QL_QUALITY_SNOW", "QL_PIXELMASK")
    )
    log_message("Maskierung erfolgreich abgeschlossen", "SUCCESS")
    return(TRUE)
  }, error = function(e) {
    log_message(sprintf("Fehler bei Maskierung: %s", e$message), "ERROR")
    return(FALSE)
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# 🔢 SCHRITT 5: NORMALISIERUNG
# ══════════════════════════════════════════════════════════════════════════════

step5_normalize <- function() {
  # Berechne Datengröße
  data_size_gb <- get_directory_size_gb(EXTRACT_DIR)
  estimated_time <- estimate_processing_time(data_size_gb, "normalize")
  
  method_name <- ifelse(NORMALIZATION_METHOD == "minmax", 
                       "Min-Max-Normalisierung",
                       "SNV-Normalisierung")
  
  # Prüfe auf existierende Ausgaben (je nach Methode)
  if (NORMALIZATION_METHOD == "minmax") {
    check_outputs <- structure(
      c(file.path(EXTRACT_DIR, "*_minmax.tif")),
      pattern_search = TRUE
    )
  } else {
    check_outputs <- structure(
      c(file.path(EXTRACT_DIR, "*_masked_snv.tif")),
      pattern_search = TRUE
    )
  }
  
  # Benutzerbestätigung
  confirmation <- ask_user_confirmation(
    step_number = 5,
    step_name = method_name,
    description = sprintf(
      "Normalisiert alle Spektralbänder mit %s Methode auf einheitlichen Wertebereich.",
      method_name
    ),
    estimated_time = estimated_time,
    check_existing_output = check_outputs
  )
  
  if (confirmation$skip) {
    if (!is.null(confirmation$abort) && confirmation$abort) return(FALSE)
    return(TRUE)
  }
  
  log_message(sprintf("Starte %s...", method_name), "INFO")
  
  tryCatch({
    if (NORMALIZATION_METHOD == "minmax") {
      # Lade Min-Max Normalisierungs-Script
      source(file.path(SCRIPTS_DIR, "apply_minmax_normalization.R"), local = TRUE)
      
      # Führe Min-Max Normalisierung durch
      # (Annahme: Funktion heißt apply_minmax_normalization)
      if (exists("apply_minmax_normalization")) {
        apply_minmax_normalization(
          smoothed_dir = EXTRACT_DIR,
          input_suffix = "_smoothed\\+destriped",
          output_suffix = "_normalized",
          overwrite = TRUE
        )
      } else {
        stop("Funktion apply_minmax_normalization nicht gefunden")
      }
    } else {
      # Lade SNV Normalisierungs-Script
      source(file.path(SCRIPTS_DIR, "apply_snv_normalization.R"), local = TRUE)
      
      apply_spectral_normalization(
        masked_dir = EXTRACT_DIR,
        input_suffix = "_smoothed\\+adaptive_destriped_masked",
        output_suffix = "_smoothed+adaptive_destriped_masked",
        method = "snv",
        overwrite = TRUE
      )
    }
    
    log_message("Normalisierung erfolgreich abgeschlossen", "SUCCESS")
    return(TRUE)
  }, error = function(e) {
    log_message(sprintf("Fehler bei Normalisierung: %s", e$message), "ERROR")
    return(FALSE)
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# 🗺️ SCHRITT 6: KOORDINATEN ANGLEICHEN & DATENEXTRAKTION
# ══════════════════════════════════════════════════════════════════════════════

step6_prepare_data <- function() {
  # Berechne Datengröße
  data_size_gb <- get_directory_size_gb(EXTRACT_DIR)
  estimated_time <- estimate_processing_time(data_size_gb, "data_prep")
  
  # Prüfe auf existierende Ausgaben (basierend auf gewählter Normalisierungsmethode)
  check_outputs <- c(
    paste0("data/combined_data_", NORMALIZATION_METHOD, ".csv"),
    "data/final_before_correlation.csv"
  )
  
  # Benutzerbestätigung
  confirmation <- ask_user_confirmation(
    step_number = 6,
    step_name = "Koordinaten angleichen & Spektraldaten extrahieren",
    description = sprintf("Projiziert Koordinaten auf EPSG:32719 und extrahiert Spektraldaten für alle Bodenprobenpunkte (Normalisierung: %s).", toupper(NORMALIZATION_METHOD)),
    estimated_time = estimated_time,
    check_existing_output = check_outputs
  )
  
  if (confirmation$skip) {
    if (!is.null(confirmation$abort) && confirmation$abort) return(FALSE)
    return(TRUE)
  }
  
  log_message("Starte Datenaufbereitung...", "INFO")
  
  # Lade Datenaufbereitungs-Script
  source(file.path(SCRIPTS_DIR, "data_preperation.R"), local = TRUE)
  
  tryCatch({
    data_preparation(
      tar_dir = ENMAP_TAR_DIR,
      extract_dir = EXTRACT_DIR,
      shapefile_path = SHAPEFILE_PATH,
      soil_data_dir = SOIL_DATA_DIR,
      normalization_method = NORMALIZATION_METHOD
    )
    output_file <- paste0("data/combined_data_", NORMALIZATION_METHOD, ".csv")
    log_message("Datenaufbereitung erfolgreich abgeschlossen", "SUCCESS")
    log_message(sprintf("Ausgabe: %s", output_file), "INFO")
    return(TRUE)
  }, error = function(e) {
    log_message(sprintf("Fehler bei Datenaufbereitung: %s", e$message), "ERROR")
    return(TRUE) # Rückgabe TRUE, um Pipeline fortzusetzen
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# 📊 SCHRITT 7: KORRELATIONSANALYSE
# ══════════════════════════════════════════════════════════════════════════════

step7_correlation_analysis <- function() {
  # Geschätzte Zeit basierend auf Anzahl Datenpunkte
  estimated_time <- estimate_processing_time(0.1, "correlation")
  
  # Prüfe auf existierende Ausgaben
  check_outputs <- c(
    file.path(RESULTS_DIR, "Pearson"),
    file.path(RESULTS_DIR, "Spearman")
  )
  
  # Benutzerbestätigung
  confirmation <- ask_user_confirmation(
    step_number = 7,
    step_name = "Korrelationsanalyse",
    description = sprintf(
      "Berechnet Pearson- und Spearman-Korrelationen mit RSS 50:50 Split (Seed=%d).",
      SEED_CORRELATION_SPLIT
    ),
    estimated_time = estimated_time,
    check_existing_output = check_outputs
  )
  
  if (confirmation$skip) {
    if (!is.null(confirmation$abort) && confirmation$abort) return(FALSE)
    return(TRUE)
  }
  
  log_message("Starte Korrelationsanalyse...", "INFO")
  
  # Lade Korrelations-Script
  source(file.path(SCRIPTS_DIR, "correlation_algorithms.R"), local = TRUE)
  
  # Verwende die richtige Datendatei basierend auf Normalisierungsmethode
  data_file <- paste0("data/combined_data_", NORMALIZATION_METHOD, ".csv")
  
  tryCatch({
    process_soil_parameters(
      data_file = data_file
    )
    log_message("Korrelationsanalyse erfolgreich abgeschlossen", "SUCCESS")
    log_message("Ausgabe: results/[Algorithm]/", "INFO")
    return(TRUE)
  }, error = function(e) {
    log_message(sprintf("Fehler bei Korrelationsanalyse: %s", e$message), "ERROR")
    return(FALSE)
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# 🎯 SCHRITT 8: SIGNIFIKANZTEST (MULTIPLE TEST CORRECTION)
# ══════════════════════════════════════════════════════════════════════════════

step8_significance_test <- function() {
  # Geschätzte Zeit
  estimated_time <- estimate_processing_time(0.05, "significance")
  
  # Prüfe auf existierende Ausgaben - suche direkt nach korrigierten Dateien
  corrected_files <- c()
  for (algo in c("Pearson", "Spearman", "BarlowTwin")) {
    algo_dir <- file.path(RESULTS_DIR, algo)
    if (dir.exists(algo_dir)) {
      algo_files <- list.files(algo_dir, 
                               pattern = "_correlation_table_corrected\\.csv$",
                               full.names = TRUE)
      corrected_files <- c(corrected_files, algo_files)
    }
  }
  
  outputs_exist <- length(corrected_files) > 0
  
  # Manuelle Benutzerbestätigung
  cat("\n")
  cat("────────────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("📍 SCHRITT 8: %s\n", toupper("Signifikanztest (Multiple Test Correction)")))
  cat("────────────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("   %s\n", "Wendet Benjamini-Hochberg FDR-Korrektur auf alle Korrelationen an (α = 0.05)."))
  cat(sprintf("   ⏱️  Geschätzte Dauer: %s\n", estimated_time))
  
  if (outputs_exist) {
    cat(sprintf("   ⚠️  %d korrigierte Datei(en) bereits vorhanden!\n", length(corrected_files)))
    cat("────────────────────────────────────────────────────────────────────────────────\n")
    response <- readline(prompt = "   🔄 Schritt wiederholen? (j=Ja/überschreiben, n=Nein/vorhandene nutzen): ")
    
    if (tolower(response) == "n") {
      cat("   ⏭️  Überspringe Schritt, nutze vorhandene Ergebnisse.\n\n")
      return(TRUE)
    } else if (tolower(response) != "j") {
      cat("   ⏸️  Ungültige Eingabe. Pipeline angehalten.\n\n")
      return(FALSE)
    }
    cat("   ♻️  Überschreibe vorhandene Ergebnisse...\n\n")
  } else {
    cat("────────────────────────────────────────────────────────────────────────────────\n")
    response <- readline(prompt = "   ▶️  Fortfahren? (j/n): ")
    
    if (tolower(response) != "j") {
      cat("   ⏸️  Pipeline angehalten. Starte erneut um fortzufahren.\n\n")
      return(FALSE)
    }
    cat("\n")
  }
  
  # Wenn wir hier ankommen, wurde bestätigt
  confirmation <- list(skip = FALSE, overwrite = outputs_exist)
  
  if (confirmation$skip) {
    if (!is.null(confirmation$abort) && confirmation$abort) return(FALSE)
    return(TRUE)
  }
  
  log_message("Starte Signifikanztest...", "INFO")
  
  # Lade Multiple Test Correction Script
  source(file.path(SCRIPTS_DIR, "multiple_test_correction.R"), local = TRUE)
  
  tryCatch({
    # Das Script führt sich automatisch aus
    log_message("Signifikanztest erfolgreich abgeschlossen", "SUCCESS")
    log_message("Ausgabe: results/[Algorithm]/*_corrected.csv", "INFO")
    return(TRUE)
  }, error = function(e) {
    log_message(sprintf("Fehler bei Signifikanztest: %s", e$message), "ERROR")
    return(FALSE)
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# 📊 SCHRITT 9: VISUALISIERUNG (OPTIONAL)
# ══════════════════════════════════════════════════════════════════════════════

step9_visualize <- function() {
  if (SKIP_VISUALIZATION) {
    log_message("Überspringe Visualisierung (SKIP_VISUALIZATION = TRUE)", "INFO")
    return(TRUE)
  }
  
  # Geschätzte Zeit
  estimated_time <- estimate_processing_time(0.2, "visualization")
  
  # Prüfe auf existierende Ausgaben - direkte Suche nach Boxplot-Dateien
  boxplot_files <- c()
  for (algo in c("Pearson", "Spearman", "BarlowTwin")) {
    algo_dir <- file.path(RESULTS_DIR, algo, "boxplots")
    if (dir.exists(algo_dir)) {
      algo_plots <- list.files(algo_dir, 
                               pattern = "\\.png$",
                               full.names = TRUE)
      boxplot_files <- c(boxplot_files, algo_plots)
    }
  }
  
  outputs_exist <- length(boxplot_files) > 0
  
  # Manuelle Benutzerbestätigung
  cat("\n")
  cat("────────────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("📍 SCHRITT 9: %s\n", toupper("Visualisierung (optional)")))
  cat("────────────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("   %s\n", "Erstellt Boxplots, Balkendiagramme und Liniendiagramme der Korrelationen."))
  cat(sprintf("   ⏱️  Geschätzte Dauer: %s\n", estimated_time))
  
  if (outputs_exist) {
    cat(sprintf("   ⚠️  %d Visualisierung(en) bereits vorhanden!\n", length(boxplot_files)))
    cat("────────────────────────────────────────────────────────────────────────────────\n")
    response <- readline(prompt = "   🔄 Schritt wiederholen? (j=Ja/überschreiben, n=Nein/vorhandene nutzen): ")
    
    if (tolower(response) == "n") {
      cat("   ⏭️  Überspringe Schritt, nutze vorhandene Ergebnisse.\n\n")
      return(TRUE)
    } else if (tolower(response) != "j") {
      cat("   ⏸️  Ungültige Eingabe. Pipeline angehalten.\n\n")
      return(FALSE)
    }
    cat("   ♻️  Überschreibe vorhandene Ergebnisse...\n\n")
  } else {
    cat("────────────────────────────────────────────────────────────────────────────────\n")
    response <- readline(prompt = "   ▶️  Fortfahren? (j/n): ")
    
    if (tolower(response) != "j") {
      cat("   ⏸️  Pipeline angehalten. Starte erneut um fortzufahren.\n\n")
      return(FALSE)
    }
    cat("\n")
  }
  
  # Wenn wir hier ankommen, wurde bestätigt
  confirmation <- list(skip = FALSE, overwrite = outputs_exist)
  
  if (confirmation$skip) {
    if (!is.null(confirmation$abort) && confirmation$abort) return(FALSE)
    return(TRUE)
  }
  
  log_message("Starte Visualisierung...", "INFO")
  
  # Lade Visualisierungs-Script
  source(file.path(SCRIPTS_DIR, "boxplot_visualization.R"), local = TRUE)
  
  tryCatch({
    visualize_correlations(
      results_base_dir = RESULTS_DIR,
      significance_threshold = 0.3
    )
    log_message("Visualisierung erfolgreich abgeschlossen", "SUCCESS")
    log_message("Ausgabe: results/[Algorithm]/boxplots/", "INFO")
    return(TRUE)
  }, error = function(e) {
    log_message(sprintf("Fehler bei Visualisierung: %s", e$message), "ERROR")
    return(FALSE)
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# 🧮 SCHRITT 10: FILTER & PLSR-TRAINING
# ══════════════════════════════════════════════════════════════════════════════

step10_filter_and_plsr <- function() {
  # Geschätzte Zeit
  estimated_time <- estimate_processing_time(0.15, "plsr")
  
  # Prüfe auf existierende Ausgaben - suche direkt nach PLSR Summary-Dateien
  plsr_summary_files <- c()
  for (algo in c("Pearson", "Spearman", "BarlowTwin")) {
    algo_dir <- file.path(RESULTS_DIR, algo, "filter+plsr")
    if (dir.exists(algo_dir)) {
      summary_file <- file.path(algo_dir, "plsr_summary.csv")
      if (file.exists(summary_file)) {
        plsr_summary_files <- c(plsr_summary_files, summary_file)
      }
    }
  }
  
  outputs_exist <- length(plsr_summary_files) > 0
  
  # Manuelle Benutzerbestätigung
  cat("\n")
  cat("────────────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("📍 SCHRITT 10: %s\n", toupper("Filter & PLSR-Training")))
  cat("────────────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("   %s\n", sprintf(
    "Filtert NA-Bänder und trainiert PLSR-Modelle mit RSS 70:30 Split (Seed=%d) und 5-fold CV (Seed=%d).",
    SEED_PLSR_SPLIT, SEED_PLSR_CV
  )))
  cat(sprintf("   ⏱️  Geschätzte Dauer: %s\n", estimated_time))
  
  if (outputs_exist) {
    cat(sprintf("   ⚠️  %d PLSR-Zusammenfassung(en) bereits vorhanden!\n", length(plsr_summary_files)))
    cat("────────────────────────────────────────────────────────────────────────────────\n")
    response <- readline(prompt = "   🔄 Schritt wiederholen? (j=Ja/überschreiben, n=Nein/vorhandene nutzen): ")
    
    if (tolower(response) == "n") {
      cat("   ⏭️  Überspringe Schritt, nutze vorhandene Ergebnisse.\n\n")
      return(TRUE)
    } else if (tolower(response) != "j") {
      cat("   ⏸️  Ungültige Eingabe. Pipeline angehalten.\n\n")
      return(FALSE)
    }
    cat("   ♻️  Überschreibe vorhandene Ergebnisse...\n\n")
  } else {
    cat("────────────────────────────────────────────────────────────────────────────────\n")
    response <- readline(prompt = "   ▶️  Fortfahren? (j/n): ")
    
    if (tolower(response) != "j") {
      cat("   ⏸️  Pipeline angehalten. Starte erneut um fortzufahren.\n\n")
      return(FALSE)
    }
    cat("\n")
  }
  
  # Wenn wir hier ankommen, wurde bestätigt
  log_message("Starte Filter und PLSR-Training...", "INFO")
  
  # Lade Filter+PLSR Script
  source(file.path(SCRIPTS_DIR, "filter+plsr.R"), local = TRUE)
  
  tryCatch({
    perform_filter_plsr_by_algorithm(
      data_file = "data/final_before_correlation.csv",
      results_base_dir = RESULTS_DIR,
      min_bands = 4,
      min_correlation = 0.3,
      use_validation_split = TRUE,
      validation_ratio = 0.3,
      max_components = 5
    )
    log_message("Filter und PLSR-Training erfolgreich abgeschlossen", "SUCCESS")
    log_message("Ausgabe: results/[Algorithm]/PLSR_models/", "INFO")
    return(TRUE)
  }, error = function(e) {
    log_message(sprintf("Fehler bei PLSR-Training: %s", e$message), "ERROR")
    return(FALSE)
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# 🗺️ SCHRITT 11: RÄUMLICHE VORHERSAGEKARTEN
# ══════════════════════════════════════════════════════════════════════════════

step11_spatial_prediction <- function() {
  # Berechne Datengröße
  data_size_gb <- get_directory_size_gb(EXTRACT_DIR)
  estimated_time <- estimate_processing_time(data_size_gb, "prediction")
  
  # Prüfe auf existierende Ausgaben - suche direkt nach Vorhersage-TIF-Dateien
  prediction_tif_files <- c()
  for (algo in c("Pearson", "Spearman")) {
    algo_dir <- file.path(RESULTS_DIR, algo)
    if (dir.exists(algo_dir)) {
      # Suche nach allen *_predicted_*.tif Dateien
      tif_files <- list.files(algo_dir, 
                               pattern = "_predicted_.*\\.tif$",
                               full.names = TRUE,
                               recursive = TRUE)
      prediction_tif_files <- c(prediction_tif_files, tif_files)
    }
  }
  
  outputs_exist <- length(prediction_tif_files) > 0
  
  # Manuelle Benutzerbestätigung
  cat("\n")
  cat("────────────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("📍 SCHRITT 11: %s\n", toupper("Räumliche Vorhersagekarten erstellen")))
  cat("────────────────────────────────────────────────────────────────────────────────\n")
  cat(sprintf("   %s\n", "Wendet trainierte PLSR-Modelle auf alle Pixel an und erstellt Vorhersagekarten."))
  cat(sprintf("   ⏱️  Geschätzte Dauer: %s\n", estimated_time))
  
  if (outputs_exist) {
    cat(sprintf("   ⚠️  %d Vorhersagekarte(n) bereits vorhanden!\n", length(prediction_tif_files)))
    cat("────────────────────────────────────────────────────────────────────────────────\n")
    response <- readline(prompt = "   🔄 Schritt wiederholen? (j=Ja/überschreiben, n=Nein/vorhandene nutzen): ")
    
    if (tolower(response) == "n") {
      cat("   ⏭️  Überspringe Schritt, nutze vorhandene Ergebnisse.\n\n")
      return(TRUE)
    } else if (tolower(response) != "j") {
      cat("   ⏸️  Ungültige Eingabe. Pipeline angehalten.\n\n")
      return(FALSE)
    }
    cat("   ♻️  Überschreibe vorhandene Ergebnisse...\n\n")
  } else {
    cat("────────────────────────────────────────────────────────────────────────────────\n")
    response <- readline(prompt = "   ▶️  Fortfahren? (j/n): ")
    
    if (tolower(response) != "j") {
      cat("   ⏸️  Pipeline angehalten. Starte erneut um fortzufahren.\n\n")
      return(FALSE)
    }
    cat("\n")
  }
  
  # Wenn wir hier ankommen, wurde bestätigt
  log_message("Starte räumliche Vorhersage...", "INFO")
  
  # Lade Spatial Prediction Script
  source(file.path(SCRIPTS_DIR, "spatial_prediction.R"), local = TRUE)
  
  tryCatch({
    run_all_predictions(
      data_file = "data/final_before_correlation.csv",
      results_base_dir = RESULTS_DIR,
      smoothed_dir = EXTRACT_DIR,
      metadata_file = "data/metadata_mapping.csv",
      normalization_method = NORMALIZATION_METHOD,
      start_index = 1,
      skip_existing = TRUE,
      algorithms = c("Pearson", "Spearman"),
      plsr_types = c("filter+plsr"),
      recursive_search = TRUE,
      create_detailed_viz = CREATE_DETAILED_VIZ
    )
    log_message("Räumliche Vorhersage erfolgreich abgeschlossen", "SUCCESS")
    log_message("Ausgabe: results/prediction_maps/", "INFO")
    return(TRUE)
  }, error = function(e) {
    log_message(sprintf("Fehler bei räumlicher Vorhersage: %s", e$message), "ERROR")
    return(FALSE)
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# 🚀 HAUPTPROGRAMM - PIPELINE AUSFÜHREN
# ══════════════════════════════════════════════════════════════════════════════

run_complete_pipeline <- function() {
  start_time <- Sys.time()
  
  cat("\n")
  cat("════════════════════════════════════════════════════════════════════════════════\n")
  cat("                            PIPELINE GESTARTET                                  \n")
  cat("════════════════════════════════════════════════════════════════════════════════\n")
  cat(sprintf("   Startzeit: %s\n", format(start_time, "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("   Working Directory: %s\n", WORKING_DIR))
  cat(sprintf("   Normalisierungsmethode: %s\n", toupper(NORMALIZATION_METHOD)))
  cat("════════════════════════════════════════════════════════════════════════════════\n")
  cat("\n")
  
  # Pipeline-Schritte ausführen
  steps <- list(
    step1_extract_images,
    step2_extract_metadata,
    step3_smooth_and_destripe,
    step4_mask_images,
    step5_normalize,
    step6_prepare_data,
    step7_correlation_analysis,
    step8_significance_test,
    step9_visualize,
    step10_filter_and_plsr,
    step11_spatial_prediction
  )
  
  # Führe jeden Schritt aus
  for (i in seq_along(steps)) {
    result <- steps[[i]]()
    
    if (!result) {
      log_message("Pipeline abgebrochen durch Benutzer oder Fehler", "WARNING")
      return(FALSE)
    }
  }
  
  # Abschluss
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  
  cat("\n")
  cat("════════════════════════════════════════════════════════════════════════════════\n")
  cat("                       PIPELINE ERFOLGREICH ABGESCHLOSSEN                       \n")
  cat("════════════════════════════════════════════════════════════════════════════════\n")
  cat(sprintf("   Endzeit: %s\n", format(end_time, "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("   Gesamtdauer: %.1f Minuten\n", as.numeric(duration)))
  cat("════════════════════════════════════════════════════════════════════════════════\n")
  cat("\n")
  cat("📁 Ergebnisse finden Sie in:\n")
  cat(sprintf("   - Verarbeitete Bilder: %s\n", EXTRACT_DIR))
  cat(sprintf("   - Analyseergebnisse: %s\n", RESULTS_DIR))
  cat(sprintf("   - Vorhersagekarten: %s\n", file.path(RESULTS_DIR, "prediction_maps")))
  cat("\n")
  
  return(TRUE)
}

# ══════════════════════════════════════════════════════════════════════════════
# ▶️ PIPELINE STARTEN
# ══════════════════════════════════════════════════════════════════════════════

# Prüfe ob alle erforderlichen Pfade existieren
if (!dir.exists(ENMAP_TAR_DIR) && !SKIP_EXTRACTION) {
  stop(sprintf("❌ EnMAP Verzeichnis nicht gefunden: %s", ENMAP_TAR_DIR))
}

if (!file.exists(SHAPEFILE_PATH) && !SKIP_EXTRACTION) {
  warning(sprintf("⚠️  Shapefile nicht gefunden: %s", SHAPEFILE_PATH))
}

# Erstelle Output-Verzeichnisse
dir.create(EXTRACT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Führe Pipeline aus
tryCatch({
  run_complete_pipeline()
}, error = function(e) {
  cat("\n")
  log_message(sprintf("KRITISCHER FEHLER: %s", e$message), "ERROR")
  cat("\n")
  stop("Pipeline abgebrochen")
})




































