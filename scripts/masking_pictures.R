library(terra)
library(tools)

# Funktion: Maskiere alle smoothed+adaptive_destriped-Bilder mit den relevanten Masken
mask_smoothed_images <- function(
  smoothed_dir = "data/extracted_files",
  mask_dir = "data/extracted_files",
  output_suffix = "_smoothed+adaptive_destriped_masked",
  mask_types = c("QL_QUALITY_CLOUD", "QL_QUALITY_CLOUDSHADOW", "QL_QUALITY_HAZE", "QL_QUALITY_CIRRUS", "QL_QUALITY_SNOW", "QL_PIXELMASK")
) {
  # Finde alle smoothed+adaptive_destriped-Bilder
  smoothed_files <- list.files(
    smoothed_dir,
    pattern = "_smoothed\\+adaptive_destriped\\.tif$",
    full.names = TRUE,
    recursive = TRUE
  )
  cat("Gefundene smoothed+adaptive_destriped-Bilder:", length(smoothed_files), "\n")

  for (img_file in smoothed_files) {
    cat("\nBearbeite:", basename(img_file), "\n")
    base_name <- file_path_sans_ext(basename(img_file))
    
    # Verbesserte Extraktion des ursprünglichen Namens
    # Entferne "_smoothed+adaptive_destriped" vom Ende
    original_base <- sub("_smoothed\\+adaptive_destriped$", "", base_name)
    cat("  Original basename:", original_base, "\n")
    
    # Verschiedene mögliche Masken-Basename Varianten testen
    possible_mask_bases <- c(
      original_base,                             # Direkt ohne Suffix
      gsub("_", "-", original_base),            # Unterstriche durch Bindestriche ersetzen
      sub("-SPECTRAL_IMAGE$", "", original_base) # Entferne -SPECTRAL_IMAGE am Ende falls vorhanden
    )
    
    cat("  Teste Masken-Basenamen:\n")
    for (i in seq_along(possible_mask_bases)) {
      cat("   ", i, ":", possible_mask_bases[i], "\n")
    }

    # Zusätzliche Debug-Info: Liste alle verfügbaren Mask-Dateien im Verzeichnis auf
    cat("  Verfügbare Mask-Dateien im Verzeichnis:\n")
    all_mask_files <- list.files(mask_dir, pattern = "\\.(TIF|tif)$", recursive = TRUE, ignore.case = TRUE)
    
    # Suche nach Dateien mit ähnlichem Namen (ohne _smoothed+adaptive_destriped)
    base_pattern <- sub("-SPECTRAL_IMAGE$", "", original_base)
    matching_files <- grep(base_pattern, all_mask_files, value = TRUE, ignore.case = TRUE)
    
    if (length(matching_files) > 0) {
      cat("  Passende Dateien gefunden:\n")
      for (mf in matching_files[1:min(20, length(matching_files))]) {
        cat("   ", mf, "\n")
      }
      if (length(matching_files) > 20) cat("   ... und", length(matching_files) - 20, "weitere\n")
    } else {
      cat("   Keine passenden Dateien gefunden für Pattern:", base_pattern, "\n")
      
      # Zeige alle verfügbaren TIF-Dateien zur Analyse
      cat("  Alle verfügbaren TIF-Dateien (erste 20):\n")
      for (af in all_mask_files[1:min(20, length(all_mask_files))]) {
        cat("   ", af, "\n")
      }
      if (length(all_mask_files) > 20) cat("   ... und", length(all_mask_files) - 20, "weitere\n")
    }

    # Masken suchen mit flexiblerem Pattern
    mask_files <- sapply(mask_types, function(masktype) {
      # Teste verschiedene Patterns
      patterns_to_test <- c()
      
      for (mask_base in possible_mask_bases) {
        # Direkte Pattern mit dem exakten Masktype
        patterns_to_test <- c(patterns_to_test,
                             paste0(mask_base, "-", masktype, "\\.TIF$"),
                             paste0(mask_base, "_", masktype, "\\.TIF$"))
      }
      
      # Teste auch ohne SPECTRAL_IMAGE
      base_without_spectral <- sub("-SPECTRAL_IMAGE$", "", original_base)
      patterns_to_test <- c(patterns_to_test,
                           paste0(base_without_spectral, "-", masktype, "\\.TIF$"),
                           paste0(base_without_spectral, "_", masktype, "\\.TIF$"))
      
      for (pattern in patterns_to_test) {
        f <- list.files(mask_dir, pattern = pattern, full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
        if (length(f) > 0) {
          cat("  Gefunden mit Pattern '", pattern, "': ", basename(f[1]), "\n")
          return(f[1])
        }
      }
      return(NA)
    }, USE.NAMES = TRUE)

    # Debug: Liste alle gefundenen Masken auf
    cat("  Gefundene Masken:\n")
    for (masktype in names(mask_files)) {
      cat("   ", masktype, ": ", ifelse(is.na(mask_files[masktype]), "NICHT gefunden", basename(mask_files[masktype])), "\n")
    }
    
    # Prüfe, ob mindestens eine Maske gefunden wurde
    if (all(is.na(mask_files))) {
      cat("  Keine Masken gefunden, überspringe.\n")
      next
    }

    # Lade das Bild
    img <- tryCatch(rast(img_file), error = function(e) {cat("  Fehler beim Laden:", e$message, "\n"); return(NULL)})
    if (is.null(img)) next

    # Lade und kombiniere Masken
    mask_stack <- list()
    for (masktype in mask_types) {
      mask_path <- mask_files[masktype]
      if (!is.na(mask_path) && file.exists(mask_path)) {
        m <- tryCatch(rast(mask_path), error = function(e) NULL)
        if (!is.null(m)) mask_stack[[masktype]] <- m
      }
    }
    if (length(mask_stack) == 0) {
      cat("  Keine gültigen Masken geladen, überspringe.\n")
      next
    }

    # Resample Masken auf Bildgröße (falls nötig)
    for (i in seq_along(mask_stack)) {
      if (!compareGeom(img, mask_stack[[i]], stopOnError = FALSE)) {
        mask_stack[[i]] <- resample(mask_stack[[i]], img, method = "near")
      }
    }

    # Kombiniere alle Masken zu einer Gesamtmaske (TRUE = maskiert)
    mask_total <- mask_stack[[1]] > 0
    if (length(mask_stack) > 1) {
      for (i in 2:length(mask_stack)) {
        mask_total <- mask_total | (mask_stack[[i]] > 0)
      }
    }

    # Setze maskierte Pixel auf NA
    img_masked <- mask(img, mask_total, maskvalue = TRUE, updatevalue = NA)

    # Speichere das maskierte Bild
    # Behalte den vollständigen Namen (inkl. SPECTRAL_IMAGE und smoothed+adaptive_destriped)
    output_file <- file.path(dirname(img_file), paste0(original_base, output_suffix, ".tif"))
    writeRaster(img_masked, output_file, overwrite = TRUE, gdal = c("COMPRESS=LZW", "TILED=YES"))
    cat("  Maskiertes Bild gespeichert:", basename(output_file), "\n")
  }
  cat("\nAlle Bilder wurden maskiert.\n")
}

# HINWEIS: Auto-Execution auskommentiert - Script wird nur über master_pipeline.R gesteuert
# # Aufruf der Funktion
# mask_smoothed_images(
#   smoothed_dir = "data/extracted_files/ich_everzwlfe",
#   mask_dir = "data/extracted_files",
#   output_suffix = "_smoothed+adaptive_destriped_masked_snv"
# )
