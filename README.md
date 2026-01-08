# EnMAP Hyperspektral-PLSR Pipeline

## Übersicht

Bei diesem Projekt handelt es sich um eine Verallgemeinerung des Codes, welcher für meine Bachelorarbeit verwendet wurde.

Es enthält eine vollautomatisierte Pipeline zur Verarbeitung von EnMAP-Hyperspektraldaten und zur Vorhersage von Bodenparametern mittels Partial Least Squares Regression (PLSR). Die Pipeline ermöglicht es, aus Feldmessungen und EnMAP-Satellitenbildern räumliche .TIF Vorhersagekarten für jeden Pixel zu erstellen.

## Anwendungsbereiche

### Primäre Anwendungen
- **Bodenparameter-Mapping**: Vorhersage von organischem Kohlenstoff, Stickstoff, Basensättigung und anderen Bodeneigenschaften !!!Achtung: Logarithmische Werte, wie der pH-Wert sollten reskaliert werden, bevor sie ins script gefüttert werden!!!
- **Präzisionslandwirtschaft**: Erstellung von Bewirtschaftungskarten basierend auf Bodeneigenschaften
- **Umweltmonitoring**: Überwachung von Bodenqualität und -degradation über große Flächen
- **Forschung**: Korrelationsanalysen zwischen Spektralsignaturen und Bodenparametern

### Geeignete Szenarien
- Sie haben Feldmessungen von Bodenparametern an spezifischen Koordinaten
- Für Ihr Untersuchungsgebiet liegen EnMAP-Satellitenbilder vor
- Sie möchten räumliche Vorhersagen für das gesamte Untersuchungsgebiet erstellen
- Sie benötigen statistisch validierte Korrelationen zwischen Spektraldaten und Bodenparametern

## Methodische Grundlagen

### Pipeline-Architektur
Die Pipeline implementiert einen wissenschaftlich fundierten Workflow mit folgenden Kernelementen:

1. **Datenvorverarbeitung**:
   - Savitzky-Golay-Filter zur spektralen Glättung
   - Horn & Woodham Destriping zur Streifenreduzierung
   - Maskierung von Wolken, Schatten, Schnee und Sensorfehlerern
   - Min-Max oder SNV-Normalisierung

2. **Statistische Analyse**:
   - Korrelationsanalyse (Pearson & Spearman) mit Random Sample Split (RSS 50:50)
   - Multiple Test Correction nach Benjamini-Hochberg (FDR)
   - Identifikation signifikanter Wellenlängenbereiche

3. **Modellierung**:
   - PLSR-Training mit RSS 70:30 Split
   - 5-fache Kreuzvalidierung
   - Optimierung der Komponentenanzahl

4. **Räumliche Vorhersage**:
   - Pixelweise Anwendung der trainierten Modelle
   - Erstellung der Vorhersagekarten

### Reproduzierbarkeit
Alle Zufallsprozesse verwenden festgelegte Seeds:
- `SEED_CORRELATION_SPLIT = 358` (RSS 50:50 für Korrelation)
- `SEED_PLSR_SPLIT = 696` (RSS 70:30 für PLSR)
- `SEED_PLSR_CV = 263` (5-fold Kreuzvalidierung)

## Projektstruktur

```
Enmap_PLSR_prediction/
├── master_pipeline.R          # Haupt-Pipeline-Script
├── data/
│   ├── ENMAP/                 # EnMAP .tar.gz Archive (Input)
│   ├── extracted_files/       # Extrahierte und verarbeitete Bilder
│   ├── soil_measurements/     # Bodenmessungsdaten (CSV/Excel)
│   └── *.shp                  # Shapefile mit Probenpunkt-Koordinaten
├── scripts/
│   ├── extract_data.r         # Extraktion der .tar.gz Archive
│   ├── extract_metadata.R     # Metadaten-Extraktion
│   ├── smooth_pictures_destripe_simple.R  # SG-Filter & Destriping
│   ├── masking_pictures.R     # Maskierung
│   ├── apply_minmax_normalization.R       # Min-Max-Normalisierung
│   ├── apply_snv_normalization.R          # SNV-Normalisierung
│   ├── data_preperation.R     # Datenaufbereitung & Extraktion
│   ├── correlation_algorithms.R           # Korrelationsanalyse
│   ├── multiple_test_correction.R         # Signifikanztests
│   ├── boxplot_visualization.R            # Visualisierungen
│   ├── filter+plsr.R          # PLSR-Modellierung
│   └── spatial_prediction.R   # Räumliche Vorhersagekarten
└── results/                   # Ausgabeverzeichnis
    ├── Pearson/               # Pearson-Korrelationsergebnisse
    ├── Spearman/              # Spearman-Korrelationsergebnisse
    └── prediction_maps/       # Vorhersagekarten (.tif)
```

## Installation & Systemanforderungen

### Systemanforderungen
- **R Version**: ≥ 4.0.0
- **RAM**: Mindestens 16 GB (empfohlen: 32 GB für große Datensätze)
- **Festplattenspeicher**: Ca. 50 GB frei (abhängig von Anzahl der EnMAP-Bilder)
- **Prozessor**: Multi-Core empfohlen für parallele Verarbeitung

### Erforderliche R-Pakete
Die Pipeline installiert automatisch alle benötigten Pakete:

```r
# Geodatenverarbeitung
terra, raster, sf

# Datenmanipulation
dplyr, readr, readxl, tidyr, data.table, tools

# Numerische Verarbeitung
pracma, matrixStats, parallel

# Statistische Modellierung
caret, pls, glmnet, vegan

# Visualisierung
ggplot2, scales, reshape2, tmap, RColorBrewer, viridis

# XML-Verarbeitung
xml2

# Fortschrittsanzeige
pbapply
```

### Manuelle Installation (falls erforderlich)
```r
install.packages(c("terra", "raster", "sf", "dplyr", "readr", "readxl", 
                   "tidyr", "xml2", "data.table", "tools", "pracma", 
                   "pbapply", "caret", "pls", "glmnet", "ggplot2", 
                   "scales", "reshape2", "tmap", "RColorBrewer", 
                   "viridis", "vegan", "matrixStats", "parallel"))
```

## Gebrauchsanweisung

### Schritt 1: Daten vorbereiten

#### 1.1 EnMAP-Satellitenbilder
- Laden Sie EnMAP-L2A-Produkte (atmosphärisch korrigiert) herunter
- Platzieren Sie die `.tar.gz`-Archive im Ordner `data/ENMAP/`
- Keine manuelle Extraktion erforderlich (wird automatisch durchgeführt)

#### 1.2 Bodenmessungsdaten
- Erstellen Sie eine Excel- oder CSV-Datei mit Ihren Messungen
- Erforderliche Spalten:
  - `ID`: Eindeutige Proben-ID
  - Bodenparameter (z.B. `pH`, `OC`, `N`, `P`, `K`, etc.)
  - Weitere Metadaten nach Bedarf
- Speichern Sie die Datei in `data/soil_measurements/`

#### 1.3 Koordinaten-Shapefile
- Erstellen Sie ein Shapefile mit den Koordinaten Ihrer Probenpunkte
- Erforderliche Attribute:
  - `ID`: Muss mit IDs in Bodenmessungsdaten übereinstimmen
  - Geometrie: Punktgeometrie (Point)
- Koordinatensystem: Achten Sie darauf, dass das Koordinatensystem Ihres Shapefiles mit dem der Satellitenbilder übereinstimmt
- Speichern Sie das Shapefile im `data/`-Ordner

**Beispiel-Datenstruktur:**

Bodenmessungsdatei (`soil_data.xlsx`):
| ID | pH | OC | N | Latitude | Longitude |
|----|----|----|---|----------|-----------|
| P01 | 6.5 | 2.3 | 0.18 | -12.0456 | -76.2341 |
| P02 | 7.1 | 1.8 | 0.15 | -12.0523 | -76.2398 |

### Schritt 2: Konfiguration anpassen

Öffnen Sie `master_pipeline.R` und passen Sie die Konfiguration an:

```r
# ══════════════════════════════════════════════════════════════════════════════
# KONFIGURATION - HIER PFADE ANPASSEN
# ══════════════════════════════════════════════════════════════════════════════

# Hauptverzeichnis
WORKING_DIR <- "C:/Users/IhrName/Desktop/Enmap_PLSR_prediction"

# Input-Pfade
ENMAP_TAR_DIR <- file.path(WORKING_DIR, "data/ENMAP")
SHAPEFILE_PATH <- file.path(WORKING_DIR, "data/IhrShapefile.shp")
SOIL_DATA_DIR <- file.path(WORKING_DIR, "data/soil_measurements")

# Normalisierungsmethode wählen
NORMALIZATION_METHOD <- "snv"  # "minmax" oder "snv"

# Seeds für Reproduzierbarkeit
SEED_CORRELATION_SPLIT <- 358
SEED_PLSR_SPLIT <- 696
SEED_PLSR_CV <- 263

# Savitzky-Golay Parameter
SG_WINDOW_SIZE <- 11  # Fenstergröße (ungerade Zahl)
SG_POLY_ORDER <- 3    # Polynomordnung

# Pipeline-Optionen
SKIP_EXTRACTION <- FALSE       # TRUE falls bereits extrahiert
SKIP_VISUALIZATION <- FALSE    # TRUE um Visualisierungen zu überspringen
```

### Schritt 3: Pipeline ausführen

#### Option A: Vollständige automatische Ausführung
```r
# In R oder RStudio
source("master_pipeline.R")
```

Die Pipeline führt Sie interaktiv durch alle Schritte:
- Vor jedem Schritt werden Sie um Bestätigung gebeten
- Sie sehen eine Zeitschätzung für jeden Schritt
- Bei bereits vorhandenen Ergebnissen werden Sie gefragt, ob diese überschrieben werden sollen

#### Option B: Einzelne Schritte ausführen
```r
# Master-Pipeline laden
source("master_pipeline.R")

# Einzelne Schritte manuell ausführen
step1_extract_images()           # Extraktion
step2_extract_metadata()         # Metadaten
step3_smooth_and_destripe()      # Glättung & Destriping
step4_mask_images()              # Maskierung
step5_normalize()                # Normalisierung
step6_prepare_data()             # Datenextraktion
step7_correlation_analysis()     # Korrelationsanalyse
step8_significance_test()        # Signifikanztest
step9_visualize()                # Visualisierung
step10_filter_and_plsr()         # PLSR-Training
step11_spatial_prediction()      # Vorhersagekarten
```

### Schritt 4: Ergebnisse interpretieren

Nach erfolgreicher Ausführung finden Sie folgende Ergebnisse:

#### 4.1 Korrelationsergebnisse (`results/Pearson/` & `results/Spearman/`)
- `*_correlation_results.csv`: R- und p-Werte für jede Wellenlänge
- `*_correlation_results_BH_corrected.csv`: FDR-korrigierte p-Werte
- `*_significant_correlations.csv`: Nur signifikante Korrelationen (α = 0.05)

**Interpretation:**
- **R > 0.7**: Starke Korrelation (gut für Vorhersagen)
- **0.5 < R < 0.7**: Moderate Korrelation
- **R < 0.5**: Schwache Korrelation
- **p-value < 0.05**: Statistisch signifikant nach FDR-Korrektur

#### 4.2 PLSR-Modelle (`results/Pearson/PLSR/` & `results/Spearman/PLSR/`)
- `*_PLSR_summary.txt`: Modellperformanz (R², RMSE, Bias)
- `*_PLSR_model.rds`: Gespeichertes PLSR-Modell (für spätere Verwendung)
- `*_predictions_vs_observed.csv`: Vorhersagen vs. Beobachtungen

**Bewertung der Modellgüte:**
- **R² > 0.8**: Exzellente Vorhersage
- **0.6 < R² < 0.8**: Gute Vorhersage
- **0.4 < R² < 0.6**: Moderate Vorhersage
- **R² < 0.4**: Schwache Vorhersage (nicht empfohlen)

#### 4.3 Vorhersagekarten (`results/prediction_maps/`)
- `*_prediction.tif`: Räumliche Vorhersagekarte (GeoTIFF)
- Georeferenziert (EPSG:32719)
- Visualisierbar in QGIS, ArcGIS, R (terra/raster)

**Verwendung der Karten:**
```r
# Karte laden und visualisieren
library(terra)
prediction <- rast("results/prediction_maps/pH_Pearson_prediction.tif")
plot(prediction, main = "pH-Vorhersagekarte")
```

#### 4.4 Visualisierungen (`results/Pearson/Plots/` & `results/Spearman/Plots/`)
- Boxplots der Korrelationen
- Balkendiagramme signifikanter Wellenlängen
- Liniendiagramme des Korrelationsverlaufs

## Erweiterte Konfiguration

### Normalisierungsmethoden

#### Min-Max-Normalisierung
```r
NORMALIZATION_METHOD <- "minmax"
```
- Skaliert Werte auf Bereich [0, 1]
- Empfohlen für Vergleiche zwischen verschiedenen Aufnahmezeitpunkten
- Formel: `(x - min) / (max - min)`

#### SNV-Normalisierung (Standard Normal Variate)
```r
NORMALIZATION_METHOD <- "snv"
```
- Entfernt Streulichteffekte
- Empfohlen für Bodenspektroskopie
- Formel: `(x - mean) / sd`

### Savitzky-Golay-Parameter anpassen

```r
# Stärkere Glättung (für verrauschte Daten)
SG_WINDOW_SIZE <- 15
SG_POLY_ORDER <- 3

# Schwächere Glättung (für saubere Daten)
SG_WINDOW_SIZE <- 7
SG_POLY_ORDER <- 2
```

**Richtlinien:**
- `WINDOW_SIZE` muss ungerade sein
- Größeres Fenster = stärkere Glättung, aber Verlust von Details
- `POLY_ORDER` typischerweise 2-4

### Seeds ändern (für neue Experimente)

```r
# Neue Seeds für unterschiedliche Trainings-/Test-Splits
SEED_CORRELATION_SPLIT <- sample(1:1000, 1)
SEED_PLSR_SPLIT <- sample(1:1000, 1)
SEED_PLSR_CV <- sample(1:1000, 1)
```

## Fehlerbehebung

### Problem: "Keine EnMAP-Archive gefunden"
**Lösung:**
- Prüfen Sie, ob `.tar.gz`-Dateien in `data/ENMAP/` liegen
- Stellen Sie sicher, dass der Pfad `ENMAP_TAR_DIR` korrekt ist

### Problem: "Shapefile-Koordinaten stimmen nicht überein"
**Lösung:**
- Stellen Sie sicher, dass IDs im Shapefile und in Bodenmessungsdaten identisch sind
- Prüfen Sie, ob Koordinaten im gültigen Bereich der EnMAP-Bilder liegen
- -Prüfen Sie ob die Koordinatenreferenzsysteme Ihres Shapefiles mit denen der EnMAP Bilder übereinstimmt.

### Problem: "Zu wenig Speicher (Memory)"
**Lösung:**
```r
# Speicherlimit erhöhen (Windows)
memory.limit(size = 56000)  # 56 GB

# Parallele Verarbeitung reduzieren
options(mc.cores = 2)  # Weniger CPU-Kerne verwenden
```

### Problem: "PLSR-Modell konvergiert nicht"
**Lösung:**
- Prüfen Sie auf Ausreißer in Bodenmessungsdaten
- Reduzieren Sie Anzahl der Komponenten
- Erhöhen Sie Stichprobengröße (mehr Feldproben)

### Problem: "Keine signifikanten Korrelationen gefunden"
**Mögliche Ursachen:**
- Zu wenige Probenpunkte (< 30)
- Große räumliche Variabilität
- Ungeeigneter Bodenparameter für spektrale Vorhersage
- Fehler in Datenaufbereitung

## Datenqualitätskriterien

### Minimale Anforderungen
- **Probenpunkte**: ≥ 30 (besser: ≥ 50)
- **Räumliche Verteilung**: Gleichmäßig über Untersuchungsgebiet
- **Variabilität**: Ausreichende Variation im Zielparameter
- **EnMAP-Bildqualität**: Wolkenbedeckung < 20%

### Empfohlene Praktiken
- **Stratifizierte Probenahme**: Verschiedene Bodentypen/Landnutzungen
- **Zeitliche Übereinstimmung**: Feldproben ±30 Tage zur EnMAP-Aufnahme
- **Mehrfachmessungen**: 3 Wiederholungen pro Probenpunkt
- **Qualitätskontrolle**: Laboranalysen nach ISO-Standards

## Typische Verarbeitungszeiten

| Schritt | Kleine Datensätze* | Große Datensätze** |
|---------|-------------------|-------------------|
| Extraktion | 5 Min. | 30 Min. |
| Metadaten | < 1 Min. | 2 Min. |
| Glättung & Destriping | 20 Min. | 3 Std. |
| Maskierung | 5 Min. | 30 Min. |
| Normalisierung | 10 Min. | 1 Std. |
| Datenextraktion | 5 Min. | 20 Min. |
| Korrelationsanalyse | 2 Min. | 5 Min. |
| Signifikanztest | < 1 Min. | 1 Min. |
| Visualisierung | 3 Min. | 10 Min. |
| PLSR-Training | 5 Min. | 15 Min. |
| Räumliche Vorhersage | 30 Min. | 4 Std. |
| **GESAMT** | **~1.5 Std.** | **~10 Std.** |

*Kleine Datensätze: 1-2 EnMAP-Szenen, 30-50 Probenpunkte  
**Große Datensätze: 5-10 EnMAP-Szenen, 100+ Probenpunkte

## Wissenschaftliche Grundlagen

### Zitierbare Methoden
Dieses Script implementiert etablierte Methoden aus der Fernerkundung und Chemometrie:

1. **Savitzky-Golay-Filter**: Savitzky, A., & Golay, M. J. (1964). Smoothing and differentiation of data by simplified least squares procedures. Analytical Chemistry, 36(8), 1627-1639.

2. **SNV-Normalisierung**: Barnes, R. J., Dhanoa, M. S., & Lister, S. J. (1989). Standard normal variate transformation and de-trending of near-infrared diffuse reflectance spectra. Applied Spectroscopy, 43(5), 772-777.

3. **PLSR**: Wold, S., Sjöström, M., & Eriksson, L. (2001). PLS-regression: a basic tool of chemometrics. Chemometrics and Intelligent Laboratory Systems, 58(2), 109-130.

4. **Multiple Test Correction**: Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society: Series B, 57(1), 289-300.

### Validierung
Die Pipeline verwendet wissenschaftlich etablierte Validierungsstrategien:
- **Hold-out-Validierung**: 70/30 Train-Test-Split
- **Kreuzvalidierung**: 5-fach für Parameteroptimierung
- **Random Sample Split**: Vermeidung von Bias durch randomisierte Aufteilung

## Support & Beitrag

### Bei Problemen
1. Prüfen Sie die Fehlerbehebung in dieser Anleitung
2. Erstellen Sie ein Issue auf GitHub mit:
   - Fehlermeldung (vollständig)
   - R-Version und Betriebssystem
   - Schritte zur Reproduktion

### Verbesserungsvorschläge
- Fork des Repositories
- Implementierung der Verbesserung
- Pull Request mit Beschreibung

## Lizenz

Dieses Projekt steht unter einer Open-Source-Lizenz. Bei Verwendung in wissenschaftlichen Publikationen bitten wir um entsprechende Zitation.

## Danksagung

Entwickelt wurde dieses Projekt mit der Beihilfe von Prof. Dr. Sabine Chabrillat und Dr. Jens Boy, beide aus dem Institut für Erdsystemwissenschaften (IESW) an der Leibniz Universität Hannover.

Außerdem wurde für die Erstellung vieler Scriptpassagen der GitHub Copilot verwendet.

---

**Version**: 1.0  
**Letztes Update**: November 2025  
**Autor**: Leo Hampe
