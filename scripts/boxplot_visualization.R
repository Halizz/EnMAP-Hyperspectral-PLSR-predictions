###############################################################################
# Skript: boxplot_visualization.R
# -----------------------------------------------------------------------------
# Beschreibung:
# Dieses Skript verarbeitet die Ergebnisse von Korrelationsanalysen zwischen
# Bodenparametern und Spektralbändern (verschiedene Algorithmen) und erstellt
# für jeden Algorithmus und Parameter verschiedene Visualisierungen:
# - Boxplots der Korrelationskoeffizienten
# - Balkendiagramme und Liniendiagramme der Korrelationen über das Spektrum
# Die Visualisierungen werden als PNG-Dateien gespeichert.
#
# Voraussetzungen:
# - R-Pakete: dplyr, readr, ggplot2, tidyr, scales, reshape2
# - Korrekte Pfade zu den Ergebnis- und Referenzdateien
#
# Nutzung:
# 1. Pfade zu den Eingabedateien ggf. anpassen.
# 2. Skript ausführen. Die Visualisierungen werden in den jeweiligen
#    Ergebnisordnern gespeichert.
#
# Autor: [Leo Hampe]
# Datum: [22.11.2025]
###############################################################################

# Bibliotheken laden
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(scales)
library(reshape2)

# Deklariere globale Variablen für NSE (Non-Standard Evaluation)
utils::globalVariables(c(
  "Correlation_corrected",
  "Band",
  "Spectral_Band",
  "Parameter",
  "wavelength_nm",
  "Correlation_Value",
  "Correlation_Strength"
))

# Lade die Wellenlängen-Referenz (Zuordnung Bandname <-> Wellenlänge)
wavelength_reference_file <- "data/wavelength_reference.csv"
if (!file.exists(wavelength_reference_file)) {
  stop("Wellenlängen-Referenzdatei nicht gefunden:", wavelength_reference_file)
}
wavelength_reference <- read_csv(wavelength_reference_file, show_col_types = FALSE)

# Funktion zur Erstellung von Visualisierungen für jeden Algorithmus
# Laden der verifizierten Korrelationen mit Debugging-Informationen
visualize_correlations <- function(
  results_base_dir = "results",
  significance_threshold = 0.5
) {
  # Liste aller Algorithmus-Verzeichnisse
  algorithms <- c("BarlowTwin", "Pearson", "Spearman")

  # Lade die gemeinsame Korrelations-Zusammenfassung
  summary_file <- file.path(results_base_dir, "correlation_summary.csv")
  if (!file.exists(summary_file)) {
    stop("Keine globale Korrelationszusammenfassung gefunden:", summary_file)
  }
  
  correlation_summary <- read_csv(summary_file, show_col_types = FALSE)
  cat("Globale Korrelationszusammenfassung geladen mit", nrow(correlation_summary), "Einträgen\n")

  for (algo_name in algorithms) {
    algo_dir <- file.path(results_base_dir, algo_name)

    if (!dir.exists(algo_dir)) {
      warning("Algorithmus-Verzeichnis nicht gefunden: ", algo_dir)
      next
    }

    cat("\n===== Visualisiere Korrelationen für", algo_name, "=====\n")

    # Erstelle boxplot-Verzeichnis für diesen Algorithmus
    boxplot_dir <- file.path(algo_dir, "boxplots")
    if (!dir.exists(boxplot_dir)) {
      dir.create(boxplot_dir, recursive = TRUE)
      cat("Boxplot-Verzeichnis erstellt:", boxplot_dir, "\n")
    }

    # Lade die Zusammenfassung der Korrelationen (falls vorhanden)
    verified_file <- file.path(algo_dir, "verified_significant_correlations.csv")

    if (!file.exists(verified_file)) {
      # Wenn keine verifizierten Korrelationen existieren, erstelle sie
      cat("Erstelle verifizierte Korrelationen für", algo_name, "...\n")
      verified_correlations <- process_verified_correlations(
        algo_dir, algo_name, correlation_summary, significance_threshold
      )
    } else {
      # Sonst lade bestehende verifizierte Korrelationen
      verified_correlations <- read_csv(verified_file, show_col_types = FALSE)
      cat("Verifizierte Korrelationen geladen:", nrow(verified_correlations), "Einträge\n")
      
      # Debug: Überprüfe die Struktur der Spearman-Korrelationen
      if (algo_name == "Spearman") {
        cat("Debugging Spearman-Korrelationen:\n")
        # Zeige die ersten paar Zeilen
        if (nrow(verified_correlations) > 0) {
          cat("Erste Zeilen der verifizierten Korrelationen:\n")
          print(head(verified_correlations))
          
          # Prüfe, ob die Spalte Spectral_Band vorhanden ist
          if ("Spectral_Band" %in% names(verified_correlations)) {
            # Zeige einige Beispiele für Bandnamen
            cat("Beispiele für Bandnamen:\n")
            print(head(unique(verified_correlations$Spectral_Band)))
            
            # Lösche die CSV-Datei, wenn alle wavelength_nm NA sind
            if (all(is.na(verified_correlations$wavelength_nm))) {
              cat("Alle Wellenlängen sind NA! Lösche veraltete CSV-Datei und erstelle neu...\n")
              file.remove(verified_file)
              # Neu erstellen
              verified_correlations <- process_verified_correlations(
                algo_dir, algo_name, correlation_summary, significance_threshold
              )
            }
          } else {
            cat("FEHLER: Die Spalte 'Spectral_Band' fehlt in den verifizierten Korrelationen!\n")
          }
        }
      }
    }

    if (nrow(verified_correlations) == 0) {
      warning("Keine verifizierten signifikanten Korrelationen für", algo_name)
      next
    }

    # Erstelle Visualisierungen für verifizierte Korrelationen
    create_correlation_visualizations(
      verified_correlations, algo_name, boxplot_dir, significance_threshold
    )
  }

  cat("\nVisualisierung aller Algorithmen abgeschlossen!\n")
}

# Funktion zur Erstellung verifizierter Korrelationen
process_verified_correlations <- function(
  algo_dir, algo_name, correlation_summary, significance_threshold
) {
  # Suche nach allen *_correlation_table_corrected.csv Dateien
  corrected_files <- list.files(
    algo_dir,
    pattern = "_correlation_table_corrected\\.csv$",
    full.names = TRUE
  )
  cat(
    "Gefundene _correlation_table_corrected.csv Dateien:",
    length(corrected_files), "\n"
  )
  if (length(corrected_files) == 0) {
    warning(
      "Keine _correlation_table_corrected.csv Dateien gefunden in: ",
      algo_dir
    )
    data.frame()
  }
  all_verified_correlations <- data.frame()
  for (file_path in corrected_files) {
    param <- sub(
      "_correlation_table_corrected\\.csv$", "",
      basename(file_path)
    )
    cat(sprintf(
      "Verarbeite Parameter: %s (Datei: %s)\n",
      param, basename(file_path)
    ))
    correlation_table <- tryCatch({
      read_csv(file_path, show_col_types = FALSE)
    }, error = function(e) {
      warning(
        "Fehler beim Laden der Datei ", file_path, ": ",
        e$message
      )
      NULL
    })
    if (is.null(correlation_table)) next
    
    # Prüfe ob Band-Spalte existiert
    if (!"Band" %in% names(correlation_table)) {
      warning("Spalte 'Band' nicht gefunden in ", file_path)
      next
    }
    
    # Nur Zeilen mit nicht-NA Correlation_corrected verwenden
    correlation_df <- correlation_table %>%
      filter(!is.na(Correlation_corrected)) %>%
      mutate(
        Parameter = param,
        Spectral_Band = as.character(Band),  # Explizit zu character konvertieren
        Correlation_Value = Correlation_corrected
      )
    
    # Füge Wellenlänge hinzu (Bandname ggf. anpassen)
    # Stelle sicher, dass beide Join-Spalten character sind
    wavelength_reference_char <- wavelength_reference %>%
      mutate(band_name = as.character(band_name))
    
    correlation_df <- correlation_df %>%
      mutate(
        Spectral_Band_corr = as.character(ifelse(
          grepl("_corr$", Spectral_Band),
          Spectral_Band,
          paste0(Spectral_Band, "_corr")
        ))
      ) %>%
      left_join(
        wavelength_reference_char,
        by = c("Spectral_Band_corr" = "band_name")
      ) %>%
      select(Parameter, Spectral_Band, wavelength_nm, Correlation_Value)
    all_verified_correlations <- bind_rows(
      all_verified_correlations,
      correlation_df
    )
    cat(sprintf(
      "  - %d verifizierte Korrelationen gefunden\n",
      nrow(correlation_df)
    ))
  }
  # Speichere die verifizierten Korrelationen
  if (nrow(all_verified_correlations) > 0) {
    verified_file <- file.path(
      algo_dir,
      "verified_significant_correlations.csv"
    )
    write_csv(all_verified_correlations, verified_file)
    cat(
      "Verifizierte signifikante Korrelationen gespeichert in:",
      verified_file, "\n"
    )
  }
  all_verified_correlations
}

# Funktion zur Erstellung von Visualisierungen
create_correlation_visualizations <- function(
  verified_correlations, algo_name, boxplot_dir, significance_threshold
) {
  # Prüfe initial auf NA-Werte und filtere sie heraus
  n_total <- nrow(verified_correlations)
  na_rows <- sum(is.na(verified_correlations$wavelength_nm))
  if (na_rows > 0) {
    cat(sprintf("WARNUNG: %d von %d Zeilen enthalten NA-Werte in wavelength_nm und werden entfernt\n", 
                na_rows, n_total))
    verified_correlations <- verified_correlations %>% filter(!is.na(wavelength_nm))
    cat("Nach Filterung verbleiben", nrow(verified_correlations), "Zeilen\n")
  }
  if (nrow(verified_correlations) == 0) {
    warning("Nach Filterung der NA-Werte sind keine Daten mehr übrig!")
    return()
  }

  # Klassifiziere Korrelationsstärke NUR für verifizierte Werte (nur Correlation_Value, da nur diese verwendet werden)
  verified_correlations <- verified_correlations %>%
    mutate(
      Correlation_Strength = case_when(
        abs(Correlation_Value) >= 0.8 ~ "Stark",
        abs(Correlation_Value) >= 0.5 ~ "Moderat",
        abs(Correlation_Value) >= 0.3 ~ "Schwach",
        TRUE ~ "Inexistent"
      ),
      Correlation_Strength = factor(
        Correlation_Strength,
        levels = c("Stark", "Moderat", "Schwach", "Inexistent")
      )
    )

  # 1. Erstelle einen Gesamt-Boxplot für alle Parameter
  all_boxplot <- ggplot(
    verified_correlations,
    aes(x = Parameter, y = Correlation_Value, fill = Correlation_Strength)
  ) +
    geom_boxplot() +
    scale_fill_manual(values = c(
      "Stark" = "#43a2ca",
      "Moderat" = "#7bccc4",
      "Schwach" = "#bae4bc",
      "Inexistent" = "grey80"
    )) +
    labs(
      title = paste0(
        algo_name, ": Verifizierte Korrelationen (nur aus _corrected, |r| > ",
        significance_threshold, ")"
      ),
      x = "Bodenparameter",
      y = "Korrelationskoeffizient",
      fill = "Korrelationsstärke"
    ) +
    theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold"),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  all_boxplot_file <- file.path(boxplot_dir, "all_parameters_boxplot.png")
  ggsave(filename = all_boxplot_file, plot = all_boxplot, width = 12, height = 8, dpi = 300)
  cat("Gesamtboxplot gespeichert in:", all_boxplot_file, "\n")

  # 2. Erstelle für jeden signifikanten Parameter eine Visualisierung
  for (param in unique(verified_correlations$Parameter)) {
    param_data <- verified_correlations %>% filter(Parameter == param)
    if (nrow(param_data) == 0) next
    cat(sprintf("\n--- Debugging für Parameter: %s ---\n", param))
    cat("Anzahl Zeilen:", nrow(param_data), "\n")
    cat("Enthält NA-Werte in wavelength_nm:", sum(is.na(param_data$wavelength_nm)), "\n")
    if (sum(is.na(param_data$wavelength_nm)) > 0) {
      cat("WARNUNG: Parameter", param, "enthält NA-Werte in wavelength_nm\n")
      param_data <- param_data %>% filter(!is.na(wavelength_nm))
      cat("Nach Filterung verbleiben", nrow(param_data), "Zeilen für diesen Parameter\n")
      if (nrow(param_data) == 0) {
        warning("Keine gültigen Wellenlängen für Parameter ", param, ", überspringe...")
        next
      }
    }
    cat(sprintf("Wellenlängen (nm): Min=%s, Max=%s, Avg=%s\n",
                min(param_data$wavelength_nm, na.rm=TRUE),
                max(param_data$wavelength_nm, na.rm=TRUE),
                mean(param_data$wavelength_nm, na.rm=TRUE)))
    param_data <- param_data %>% arrange(wavelength_nm)

    # a) Boxplot für diesen Parameter
    param_boxplot <- ggplot(
      param_data,
      aes(x = "Korrelationen", y = Correlation_Value, fill = Correlation_Strength)
    ) +
      geom_boxplot(width = 0.6) +
      scale_fill_manual(values = c(
        "Stark" = "#43a2ca",
        "Moderat" = "#7bccc4",
        "Schwach" = "#bae4bc",
        "Inexistent" = "grey80"
      )) +
      labs(
        title = paste0(algo_name, ": ", param, " - Korrelationen"),
        x = NULL,
        y = "Korrelationskoeffizient",
        fill = "Korrelationsstärke"
      ) +
      theme_minimal() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold")
      )
    param_boxplot_file <- file.path(boxplot_dir, paste0(param, "_boxplot.png"))
    ggsave(
      filename = param_boxplot_file,
      plot = param_boxplot,
      width = 6, height = 8, dpi = 300
    )
    cat("Boxplot für", param, "gespeichert in:", param_boxplot_file, "\n")

    # b) Balkendiagramm der Korrelationen gegen Wellenlänge
    if (all(is.na(param_data$wavelength_nm))) {
      warning("Alle wavelength_nm-Werte sind NA für Parameter ", param, ", überspringe Balkendiagramm")
    } else {
      param_bars <- ggplot(
        param_data,
        aes(x = wavelength_nm, y = Correlation_Value, fill = Correlation_Strength)
      ) +
        geom_col(width = 5) +
        scale_fill_manual(values = c(
          "Stark" = "#43a2ca",
          "Moderat" = "#7bccc4",
          "Schwach" = "#bae4bc",
          "Inexistent" = "grey80"
        )) +
        scale_x_continuous(
          name = "Wellenlänge (nm)",
          breaks = pretty_breaks(n = 10)
        ) +
        labs(
          title = paste0(algo_name, ": ", param, " - Korrelationen nach Wellenlänge"),
          y = "Korrelationskoeffizient",
          fill = "Korrelationsstärke"
        ) +
        theme_minimal() +
        theme(
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold")
        )
      param_bars_file <- file.path(boxplot_dir, paste0(param, "_wavelength_bars.png"))
      ggsave(
        filename = param_bars_file,
        plot = param_bars,
        width = 12, height = 6, dpi = 300
      )
      cat("Wellenlängen-Diagramm für", param, "gespeichert in:", param_bars_file, "\n")
    }

    # c) Liniendiagramm der Korrelationsintensität gegen Wellenlänge
    if (all(is.na(param_data$wavelength_nm))) {
      warning("Alle wavelength_nm-Werte sind NA für Parameter ", param, ", überspringe Liniendiagramm")
    } else {
      param_line <- ggplot(
        param_data,
        aes(x = wavelength_nm, y = abs(Correlation_Value))
      ) +
        geom_line(linewidth = 1.2, color = "#0868ac") +
        geom_point(aes(color = Correlation_Strength), size = 3) +
        scale_color_manual(values = c(
          "Stark" = "#43a2ca",
          "Moderat" = "#7bccc4",
          "Schwach" = "#bae4bc",
          "Inexistent" = "grey80"
        )) +
        scale_x_continuous(
          name = "Wellenlänge (nm)",
          breaks = pretty_breaks(n = 10)
        ) +
        labs(
          title = paste0(algo_name, ": ", param, " - Korrelationsintensität"),
          y = "Absolute Korrelation",
          color = "Korrelationsstärke"
        ) +
        theme_minimal() +
        theme(
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 16, face = "bold")
        )
      param_line_file <- file.path(boxplot_dir, paste0(param, "_intensity_line.png"))
      ggsave(
        filename = param_line_file,
        plot = param_line,
        width = 12, height = 6, dpi = 300
      )
      cat("Intensitäts-Diagramm für", param, "gespeichert in:", param_line_file, "\n")
    }
  }

  # 3. Vergleich aller Parameter über das Spektrum
  if (all(is.na(verified_correlations$wavelength_nm))) {
    warning("Alle wavelength_nm-Werte sind NA, überspringe Spektrum-Plot")
  } else if (length(unique(verified_correlations$Parameter)) > 1) {
    spectrum_plot <- ggplot(
      verified_correlations,
      aes(x = wavelength_nm, y = Correlation_Value, color = Parameter)
    ) +
      geom_line(linewidth = 1) +
      geom_point(aes(shape = Correlation_Strength), size = 3) +
      scale_x_continuous(
        name = "Wellenlänge (nm)",
        breaks = pretty_breaks(n = 10)
      ) +
      labs(
        title = paste0(algo_name, ": Korrelationen über das Spektrum"),
        y = "Korrelationskoeffizient"
      ) +
      theme_minimal() +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)
      )
    spectrum_plot_file <- file.path(boxplot_dir, "all_parameters_spectrum.png")
    ggsave(
      filename = spectrum_plot_file,
      plot = spectrum_plot,
      width = 14, height = 8, dpi = 300
    )
    cat("Spektrum-Plot gespeichert in:", spectrum_plot_file, "\n")
  }
}

# HINWEIS: Auto-Execution auskommentiert - Script wird nur über master_pipeline.R gesteuert
# # Ausführung der Hauptfunktion
# visualize_correlations(
#   results_base_dir = "results",
#   significance_threshold = 0.3
# )
