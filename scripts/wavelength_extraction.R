library(xml2)
library(dplyr)
library(readr)

# Funktion zum Extrahieren und Speichern der Wellenlängen-Referenz
extract_and_save_wavelengths <- function(metadata_file, output_file = "data/wavelength_reference.csv") {
  # Lade die XML-Datei
  if (!file.exists(metadata_file)) {
    stop("Die Metadatendatei wurde nicht gefunden: ", metadata_file)
  }

  metadata <- read_xml(metadata_file)

  # Extrahiere die zentralen Wellenlängen
  wavelengths <- xml_find_all(metadata, ".//bandID/wavelengthCenterOfBand") %>%
    xml_double()

  # Erstelle einen Dataframe mit Band-IDs und Wellenlängen
  band_names <- paste0("Band", seq_along(wavelengths), "_corr")
  wavelength_df <- data.frame(
    band_name = band_names,
    band_number = seq_along(wavelengths),
    wavelength_nm = wavelengths
  )
  
  # Speichere die Referenzdaten
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  write_csv(wavelength_df, output_file)
  
  cat("Wellenlängen-Referenz gespeichert in:", output_file, "\n")
  cat("Anzahl der gespeicherten Bänder:", length(wavelengths), "\n")
  
  return(wavelength_df)
}

# Funktion zum Laden der Wellenlängen-Referenz
load_wavelength_reference <- function(reference_file = "data/wavelength_reference.csv") {
  if (!file.exists(reference_file)) {
    stop("Die Wellenlängen-Referenzdatei wurde nicht gefunden: ", reference_file)
  }
  
  wavelength_df <- read_csv(reference_file, show_col_types = FALSE)
  return(wavelength_df)
}

# HINWEIS: Auto-Execution auskommentiert - Script wird nur über master_pipeline.R gesteuert
# # Extrahiere und speichere die Wellenlängen-Referenz (nur einmal ausführen)
# if (!file.exists("data/wavelength_reference.csv")) {
#   metadata_file <- "data/extracted_files/Ich_everzwlfe/ENMAP01-____L2A-DT0000005022_20221104T152642Z_001_V010502_20250327T170118Z/ENMAP01-____L2A-DT0000005022_20221104T152642Z_001_V010502_20250327T170118Z-METADATA.XML"
#   wavelength_reference <- extract_and_save_wavelengths(metadata_file)
# } else {
#   cat("Wellenlängen-Referenz bereits vorhanden.\n")
# }