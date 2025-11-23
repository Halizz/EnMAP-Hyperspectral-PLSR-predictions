# Script: csv_to_shapefile.R
# Zweck: Wandelt die Tabelle final_before_correlation.csv in ein ESRI Shapefile um
# Voraussetzung: Die CSV muss Spalten für Koordinaten (z.B. lon, lat oder x, y) enthalten

library(sf)
library(readr)
library(dplyr)

# Pfad zur CSV-Datei und zum Ausgabeverzeichnis
csv_path <- "data/final_before_correlation.csv"
shapefile_path <- "data/final_before_correlation.shp"

# Prüfe, ob Datei existiert und zeige die ersten Zeilen zur Kontrolle
if (!file.exists(csv_path)) {
  cat("Arbeitsverzeichnis: ", getwd(), "\n")
  stop("Datei existiert nicht: ", csv_path)
}
cat("Erste Zeilen der CSV:\n")
print(readLines(csv_path, n = 5))

# CSV einlesen mit explizitem Encoding
df <- read_csv(csv_path, col_types = cols(), locale = locale(encoding = "UTF-8"))

# Prüfe, ob Spaltennamen korrekt erkannt wurden
cat("Spaltennamen in der CSV:\n", paste(names(df), collapse=", "), "\n")

# Spaltennamen bereinigen und ausgeben
coord_names <- names(df)
coord_names_trim <- trimws(coord_names)
cat("Spaltennamen in der CSV:\n"); print(coord_names_trim)

# Suche nach X und Y (exakt oder getrimmt, unabhängig von Groß-/Kleinschreibung)
x_idx <- which(tolower(coord_names_trim) == "x")
y_idx <- which(tolower(coord_names_trim) == "y")
if (length(x_idx) == 1 && length(y_idx) == 1) {
  lon_col <- coord_names[x_idx]
  lat_col <- coord_names[y_idx]
} else {
  # Suche nach Koordinatenspalten (lon/lat/x/y)
  lon_candidates <- grep("lon|x|longitude", coord_names_trim, ignore.case = TRUE, value = TRUE)
  lat_candidates <- grep("lat|y|latitude", coord_names_trim, ignore.case = TRUE, value = TRUE)
  if (length(lon_candidates) == 0 || length(lat_candidates) == 0) {
    stop("Konnte keine Koordinatenspalten (lon/lat/x/y oder X/Y) erkennen!\nSpaltennamen: ", paste(coord_names_trim, collapse=", "))
  }
  lon_col <- lon_candidates[1]
  lat_col <- lat_candidates[1]
  # Zusätzliche Prüfung auf closures (Namenskonflikt mit base::X)
  if (is.function(df[[lon_col]]) || is.function(df[[lat_col]])) {
    stop("Fehler: Koordinatenspalten wurden als Funktion erkannt. Prüfe die Spaltennamen!")
  }
}

# DataFrame zu sf-Objekt konvertieren
# Beispiel: UTM Zone 18S = EPSG:32718, ggf. anpassen!
sf_obj <- st_as_sf(df, coords = c(lon_col, lat_col), crs = 32718, remove = FALSE)

# Als Shapefile speichern
st_write(sf_obj, shapefile_path, delete_layer = TRUE)

cat("Shapefile erfolgreich erstellt:", shapefile_path, "\n")
cat("Shapefile erfolgreich erstellt:", shapefile_path, "\n")
