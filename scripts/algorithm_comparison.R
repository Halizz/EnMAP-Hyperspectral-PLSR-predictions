# Lade benötigte Pakete
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("gridExtra")) install.packages("gridExtra")
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Script zur Identifikation starker Korrelationen in mehreren Algorithmen
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(scales)

# Hauptfunktion zum Vergleich der Algorithmen
compare_algorithm_correlations <- function(
  results_base_dir = "results",
  correlation_threshold = 0.5,
  min_algorithms = 2
) {
  # Liste der zu analysierenden Algorithmen
  algorithms <- c("BarlowTwin", "Pearson", "Spearman")

  cat("Vergleiche Korrelationen der Algorithmen:", paste(algorithms, collapse = ", "), "\n")

  # Sammle alle verifizierten Korrelationen für jeden Algorithmus
  all_verified_correlations <- list()

  for (algo in algorithms) {
    file_path <- file.path(results_base_dir, algo, "verified_significant_correlations.csv")

    if (file.exists(file_path)) {
      verified_data <- read_csv(file_path, show_col_types = FALSE)

      # Nur starke Korrelationen mit absolutem Wert über dem Schwellenwert behalten
      strong_correlations <- verified_data %>%
        filter(abs(Correlation_Value) >= correlation_threshold) %>%
        # Erstelle einzigartige ID für jede Parameter-Band-Kombination
        mutate(combo_id = paste(Parameter, Spectral_Band, sep = "_"))

      all_verified_correlations[[algo]] <- strong_correlations
      cat(sprintf("  - Algorithmus %s: %d starke Korrelationen (|r| >= %.1f) geladen\n",
                  algo, nrow(strong_correlations), correlation_threshold))
    } else {
      cat(sprintf("  - Warnung: Keine Korrelationsdaten für %s gefunden\n", algo))
      all_verified_correlations[[algo]] <- data.frame()
    }
  }

  # Erstelle Liste aller einzigartigen Parameter-Band-Kombinationen
  all_combos <- unique(c(
    all_verified_correlations$BarlowTwin$combo_id,
    all_verified_correlations$Pearson$combo_id,
    all_verified_correlations$Spearman$combo_id
  ))

  cat(sprintf("\nInsgesamt %d einzigartige Parameter-Band-Kombinationen gefunden\n", length(all_combos)))

  # Erstelle Vergleichsdaten für jede Kombination - behalte die ursprünglichen Werte
  comparison_data <- data.frame(
    combo_id = all_combos,
    stringsAsFactors = FALSE
  )

  # Extrahiere Parameter und Spectral_Band korrekt
  # Suche nach dem Muster "corrected_spectrum_Band" und teile davor und danach
  comparison_data <- comparison_data %>%
    mutate(
      Parameter = sub("_corrected_spectrum_Band[0-9]+$", "", combo_id),
      Spectral_Band = sub("^.*?_", "", sub("^.*?_", "",
                          sub(paste0("^", Parameter, "_"), "", combo_id)))
    )

  # Füge Korrelationswerte für jeden Algorithmus hinzu
  for (algo in algorithms) {
    algo_data <- all_verified_correlations[[algo]]
    if (nrow(algo_data) > 0) {
      # Extrahiere die Korrelationswerte
      algo_values <- algo_data %>%
        select(combo_id, Correlation_Value) %>%
        rename(!!paste0(algo, "_Correlation") := Correlation_Value)

      # Füge durch linken Join hinzu
      comparison_data <- comparison_data %>%
        left_join(algo_values, by = "combo_id")
    } else {
      # Wenn keine Daten vorhanden sind, füge NA-Spalte hinzu
      comparison_data[[paste0(algo, "_Correlation")]] <- NA
    }
  }

  # Ersetze NAs durch 0 für die Zählung
  comparison_data_for_count <- comparison_data %>%
    mutate(
      BarlowTwin_Strong = ifelse(!is.na(BarlowTwin_Correlation), 1, 0),
      Pearson_Strong = ifelse(!is.na(Pearson_Correlation), 1, 0),
      Spearman_Strong = ifelse(!is.na(Spearman_Correlation), 1, 0)
    )

  # Zähle, wie viele Algorithmen eine starke Korrelation zeigen
  comparison_data_for_count <- comparison_data_for_count %>%
    mutate(
      AlgorithmCount = BarlowTwin_Strong + Pearson_Strong + Spearman_Strong
    )

  # Filtere nach Kombinationen, die in mindestens min_algorithms Algorithmen stark korreliert sind
  strong_in_multiple <- comparison_data_for_count %>%
    filter(AlgorithmCount >= min_algorithms) %>%
    select(combo_id, Parameter, Spectral_Band,
           BarlowTwin_Correlation, Pearson_Correlation, Spearman_Correlation,
           AlgorithmCount)

  # Gruppiere nach Parameter und füge Wellenlängeninformationen hinzu
  if (file.exists("data/wavelength_reference.csv")) {
    wavelength_ref <- read_csv("data/wavelength_reference.csv", show_col_types = FALSE)

    # Extrahiere die Band-Nummer aus der Spectral_Band-Spalte
    strong_in_multiple <- strong_in_multiple %>%
      mutate(band_number = as.numeric(sub("corrected_spectrum_Band", "", Spectral_Band))) %>%
      left_join(wavelength_ref, by = c("band_number" = "band_number")) %>%
      arrange(Parameter, wavelength_nm)
  }

  # Zusammenfassung ausgeben
  cat(sprintf("\n%d Parameter-Band-Kombinationen sind in mindestens %d Algorithmen stark korreliert (|r| >= %.1f)\n",
              nrow(strong_in_multiple), min_algorithms, correlation_threshold))

  # Liste der Parameter mit mindestens einer starken Korrelation in mehreren Algorithmen
  unique_params <- unique(strong_in_multiple$Parameter)
  cat("Identifizierte Parameter:", paste(unique_params, collapse = ", "), "\n")

  # Ergebnisse speichern
  output_dir <- file.path(results_base_dir, "algorithm_comparison")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
    cat("Vergleichsverzeichnis erstellt:", output_dir, "\n")
  }

  output_file <- file.path(output_dir, sprintf("strong_correlations_in_%d_of_%d_algorithms.csv",
                                              min_algorithms, length(algorithms)))
  write_csv(strong_in_multiple, output_file)
  cat("Ergebnisse gespeichert in:", output_file, "\n")

  # Erstelle eine Zusammenfassung pro Parameter
  parameter_summary <- strong_in_multiple %>%
    group_by(Parameter) %>%
    summarise(
      TotalBands = n(),
      BarlowTwin_AvgCorrelation = mean(abs(BarlowTwin_Correlation), na.rm = TRUE),
      Pearson_AvgCorrelation = mean(abs(Pearson_Correlation), na.rm = TRUE),
      Spearman_AvgCorrelation = mean(abs(Spearman_Correlation), na.rm = TRUE),
      In_All_3 = sum(AlgorithmCount == 3),
      In_2_Only = sum(AlgorithmCount == 2),
      .groups = "drop"
    ) %>%
    arrange(desc(TotalBands))

  summary_file <- file.path(output_dir, "parameter_correlation_summary.csv")
  write_csv(parameter_summary, summary_file)
  cat("Parameter-Zusammenfassung gespeichert in:", summary_file, "\n")
  
  # Vergleiche die PLSR-Modellleistungen zwischen den Algorithmen
  compare_plsr_results(results_base_dir, output_dir, algorithms)

  # Erstelle ein Plots-Verzeichnis für die Visualisierungen
  plots_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
    cat("Plot-Verzeichnis erstellt:", plots_dir, "\n")
  }

  # Erstelle für jeden Parameter einen Plot der Korrelationsintensität gegen die Wellenlänge
  for (param in unique_params) {
    # Filtere die Daten für diesen Parameter
    param_data <- strong_in_multiple %>% 
      filter(Parameter == param)
    
    # Prüfe, ob der Parameter Daten enthält
    if (nrow(param_data) == 0) {
      cat("Warnung: Keine Daten für Parameter", param, "gefunden\n")
      next
    }
    
    # Verwende Band-Nummer als x-Achse, falls keine Wellenlängendaten verfügbar sind
    has_wavelength <- "wavelength_nm" %in% colnames(param_data) && !all(is.na(param_data$wavelength_nm))
    has_band_number <- "band_number" %in% colnames(param_data) && !all(is.na(param_data$band_number))
    
    if (!has_wavelength && !has_band_number) {
      # Wenn weder Wellenlänge noch Band-Nummer vorhanden sind, dann Band-Index verwenden
      param_data$x_value <- 1:nrow(param_data)
      x_axis_name <- "Spektral-Band Index"
    } else if (!has_wavelength && has_band_number) {
      # Wenn keine Wellenlängendaten verfügbar sind, aber Band-Nummern, dann diese verwenden
      param_data$x_value <- param_data$band_number
      x_axis_name <- "Band-Nummer"
    } else {
      # Wenn Wellenlängendaten verfügbar sind, diese verwenden
      param_data$x_value <- param_data$wavelength_nm
      x_axis_name <- "Wellenlänge (nm)"
    }
    
    # Sortiere nach x-Wert
    param_data <- param_data %>% arrange(x_value)
    
    # Konvertiere zu langem Format für ggplot
    plot_data <- param_data %>%
      select(Parameter, Spectral_Band, x_value, 
             BarlowTwin_Correlation, Pearson_Correlation, Spearman_Correlation) %>%
      pivot_longer(cols = c(BarlowTwin_Correlation, Pearson_Correlation, Spearman_Correlation),
                   names_to = "Algorithm", 
                   values_to = "Correlation") %>%
      mutate(Algorithm = sub("_Correlation$", "", Algorithm))
    
    # Erstelle den Plot für die absolute Korrelation (nur mit vorhandenen Daten)
    plot_data_abs <- plot_data %>% 
      filter(!is.na(Correlation)) %>%
      mutate(Correlation_Abs = abs(Correlation))
    
    if (nrow(plot_data_abs) > 0) {
      p <- ggplot(plot_data_abs, aes(x = x_value, y = Correlation_Abs, color = Algorithm)) +
        geom_point(size = 2, alpha = 0.7) +
        geom_line(aes(group = Algorithm), linewidth = 1) +
        geom_hline(yintercept = correlation_threshold, linetype = "dashed", color = "darkgray") +
        scale_color_manual(values = c("BarlowTwin" = "#0072B2", 
                                      "Pearson" = "#D55E00", 
                                      "Spearman" = "#009E73")) +
        labs(
          title = paste("Korrelationsintensität für Parameter:", param),
          subtitle = paste("Schwellenwert:", correlation_threshold, "- Min. Algorithmen:", min_algorithms),
          x = x_axis_name,
          y = "Korrelationsintensität (|r|)"
        ) +
        scale_y_continuous(limits = c(0, 1)) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 12),
          legend.position = "bottom",
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "lightgray", fill = NA)
        )
      
      # Speichere den Plot
      plot_file <- file.path(plots_dir, paste0(param, "_correlation_intensity.png"))
      ggsave(plot_file, plot = p, width = 10, height = 6, dpi = 300)
      cat("Plot für Parameter", param, "gespeichert in:", plot_file, "\n")
    } else {
      cat("Warnung: Keine gültigen Daten für Intensitätsplot von Parameter", param, "\n")
    }
    
    # Erstelle den Plot für die Korrelation mit Vorzeichen
    if (nrow(plot_data_abs) > 0) {
      p_signed <- ggplot(plot_data_abs, aes(x = x_value, y = plot_data$Correlation[!is.na(plot_data$Correlation)], color = Algorithm)) +
        geom_point(size = 2, alpha = 0.7) +
        geom_line(aes(group = Algorithm), linewidth = 1) +
        geom_hline(yintercept = c(correlation_threshold, -correlation_threshold), 
                   linetype = "dashed", color = "darkgray") +
        scale_color_manual(values = c("BarlowTwin" = "#0072B2", 
                                      "Pearson" = "#D55E00", 
                                      "Spearman" = "#009E73")) +
        labs(
          title = paste("Korrelationskoeffizienten für Parameter:", param),
          subtitle = paste("Schwellenwert: ±", correlation_threshold, "- Min. Algorithmen:", min_algorithms),
          x = x_axis_name,
          y = "Korrelationskoeffizient (r)"
        ) +
        scale_y_continuous(limits = c(-1, 1)) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 12),
          legend.position = "bottom",
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "lightgray", fill = NA)
        )
      
      # Speichere den zweiten Plot
      plot_signed_file <- file.path(plots_dir, paste0(param, "_correlation_signed.png"))
      ggsave(plot_signed_file, plot = p_signed, width = 10, height = 6, dpi = 300)
      cat("Vorzeichen-Plot für Parameter", param, "gespeichert in:", plot_signed_file, "\n")
    } else {
      cat("Warnung: Keine gültigen Daten für Vorzeichen-Plot von Parameter", param, "\n")
    }
  }
  
  # Erstelle einen übergreifenden Plot für alle Parameter
  if (length(unique_params) > 1) {
    # Erstelle Datengrundlage für den Plot über alle Parameter
    all_param_data <- strong_in_multiple
    
    # Verwende Band-Nummer als x-Achse, falls keine Wellenlängendaten verfügbar sind
    has_wavelength <- "wavelength_nm" %in% colnames(all_param_data) && !all(is.na(all_param_data$wavelength_nm))
    
    if (!has_wavelength) {
      # Wenn keine Wellenlängendaten verfügbar sind, können wir keinen kombinierten Plot erstellen
      cat("Warnung: Keine Wellenlängendaten für kombinierten Plot verfügbar. Plot wird übersprungen.\n")
    } else {
      # Filtere nur Daten mit gültigen Wellenlängen
      all_param_data <- all_param_data %>% 
        filter(!is.na(wavelength_nm))
      
      # Konvertiere zu langem Format für ggplot
      all_plot_data <- all_param_data %>%
        select(Parameter, wavelength_nm, 
               BarlowTwin_Correlation, Pearson_Correlation, Spearman_Correlation) %>%
        pivot_longer(cols = c(BarlowTwin_Correlation, Pearson_Correlation, Spearman_Correlation),
                     names_to = "Algorithm", 
                     values_to = "Correlation") %>%
        filter(!is.na(Correlation)) %>%
        mutate(
          Algorithm = sub("_Correlation$", "", Algorithm),
          Correlation_Abs = abs(Correlation)
        )
      
      if (nrow(all_plot_data) > 0) {
        # Erstelle den Plot
        p_all <- ggplot(all_plot_data, aes(x = wavelength_nm, y = Correlation_Abs, 
                                           color = Parameter, linetype = Algorithm)) +
          geom_line(linewidth = 0.8) +
          geom_point(size = 1.2, alpha = 0.7) +
          geom_hline(yintercept = correlation_threshold, linetype = "dashed", color = "darkgray") +
          scale_x_continuous(name = "Wellenlänge (nm)", breaks = pretty_breaks(10)) +
          scale_y_continuous(name = "Korrelationsintensität (|r|)", limits = c(0, 1)) +
          labs(
            title = "Korrelationsintensität aller Parameter",
            subtitle = paste("Schwellenwert:", correlation_threshold, "- Min. Algorithmen:", min_algorithms)
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_text(size = 12),
            legend.position = "right",
            panel.grid.minor = element_blank(),
            panel.border = element_rect(color = "lightgray", fill = NA)
          )
        
        # Speichere den Plot
        all_plot_file <- file.path(plots_dir, "all_parameters_correlation_intensity.png")
        ggsave(all_plot_file, plot = p_all, width = 12, height = 8, dpi = 300)
        cat("Gesamtplot für alle Parameter gespeichert in:", all_plot_file, "\n")
      } else {
        cat("Warnung: Keine gültigen Daten für den Gesamtplot.\n")
      }
    }
  }

  return(strong_in_multiple)
}

# Neue Funktion zum Vergleichen der PLSR-Modellleistungen
compare_plsr_results <- function(results_base_dir, output_dir, algorithms) {
  # Erstelle ein Unterverzeichnis für PLSR-Vergleiche
  plsr_dir <- file.path(output_dir, "plsr_comparison")
  if (!dir.exists(plsr_dir)) {
    dir.create(plsr_dir)
    cat("PLSR-Vergleichsverzeichnis erstellt:", plsr_dir, "\n")
  }
  
  # Sammle PLSR-Ergebnisse für jeden Algorithmus
  all_plsr_results <- list()
  plsr_parameters <- c()
  
  for (algo in algorithms) {
    plsr_summary_file <- file.path(results_base_dir, algo, "cars_plsr", "plsr_summary.csv")
    
    if (file.exists(plsr_summary_file)) {
      plsr_data <- read_csv(plsr_summary_file, show_col_types = FALSE) %>%
        mutate(Algorithm = algo)
      
      all_plsr_results[[algo]] <- plsr_data
      plsr_parameters <- unique(c(plsr_parameters, plsr_data$Parameter))
      cat(sprintf("  - PLSR-Ergebnisse für %s geladen: %d Parameter\n", 
                  algo, nrow(plsr_data)))
    } else {
      cat(sprintf("  - Warnung: Keine PLSR-Ergebnisse für %s gefunden\n", algo))
    }
  }
  
  if (length(all_plsr_results) == 0) {
    cat("Keine PLSR-Ergebnisse gefunden. PLSR-Vergleich wird übersprungen.\n")
    return(invisible(NULL))
  }
  
  # Kombiniere alle Ergebnisse
  combined_plsr <- bind_rows(all_plsr_results)
  
  # Speichere die kombinierten Ergebnisse
  combined_file <- file.path(plsr_dir, "combined_plsr_results.csv")
  write_csv(combined_plsr, combined_file)
  cat("Kombinierte PLSR-Ergebnisse gespeichert in:", combined_file, "\n")
  
  # Erstelle eine breite Version für einfachen Vergleich
  wide_plsr <- combined_plsr %>%
    select(Parameter, Algorithm, R_Squared, RMSE, RPD) %>%
    pivot_wider(
      names_from = Algorithm,
      values_from = c(R_Squared, RMSE, RPD),
      names_sep = "_"
    )
  
  # Speichere die breite Version
  wide_file <- file.path(plsr_dir, "plsr_wide_comparison.csv")
  write_csv(wide_plsr, wide_file)
  cat("Vergleichstabelle für PLSR gespeichert in:", wide_file, "\n")
  
  # Erstelle Vergleichs-Visualisierungen für PLSR-Metriken
  metrics_plot_dir <- file.path(plsr_dir, "metric_plots")
  if (!dir.exists(metrics_plot_dir)) {
    dir.create(metrics_plot_dir)
  }
  
  # R² Vergleich
  if (nrow(combined_plsr) > 0) {
    # Plot für R²
    p_rsq <- ggplot(combined_plsr, aes(x = reorder(Parameter, R_Squared), y = R_Squared, fill = Algorithm)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("BarlowTwin" = "#0072B2", 
                                  "Pearson" = "#D55E00", 
                                  "Spearman" = "#009E73")) +
      labs(
        title = "Vergleich der PLSR R²-Werte nach Algorithmus",
        x = "Parameter",
        y = "R²"
      ) +
      coord_flip() +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom"
      )
    
    # Speichere den Plot
    rsq_plot_file <- file.path(metrics_plot_dir, "rsquared_comparison.png")
    ggsave(rsq_plot_file, plot = p_rsq, width = 10, height = 6 + 0.25 * length(unique(combined_plsr$Parameter)), dpi = 300)
    cat("R²-Vergleichsplot gespeichert in:", rsq_plot_file, "\n")
    
    # Plot für RMSE
    p_rmse <- ggplot(combined_plsr, aes(x = reorder(Parameter, -RMSE), y = RMSE, fill = Algorithm)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("BarlowTwin" = "#0072B2", 
                                   "Pearson" = "#D55E00", 
                                   "Spearman" = "#009E73")) +
      labs(
        title = "Vergleich der PLSR RMSE-Werte nach Algorithmus",
        x = "Parameter",
        y = "RMSE"
      ) +
      coord_flip() +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom"
      )
    
    # Speichere den Plot
    rmse_plot_file <- file.path(metrics_plot_dir, "rmse_comparison.png")
    ggsave(rmse_plot_file, plot = p_rmse, width = 10, height = 6 + 0.25 * length(unique(combined_plsr$Parameter)), dpi = 300)
    cat("RMSE-Vergleichsplot gespeichert in:", rmse_plot_file, "\n")
    
    # Plot für RPD
    p_rpd <- ggplot(combined_plsr, aes(x = reorder(Parameter, RPD), y = RPD, fill = Algorithm)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = c("BarlowTwin" = "#0072B2", 
                                   "Pearson" = "#D55E00", 
                                   "Spearman" = "#009E73")) +
      labs(
        title = "Vergleich der PLSR RPD-Werte nach Algorithmus",
        x = "Parameter",
        y = "RPD"
      ) +
      geom_hline(yintercept = c(1.4, 2.0, 3.0), linetype = "dashed", color = "darkgray") +
      annotate("text", x = 1, y = c(1.4, 2.0, 3.0) + 0.1, label = c("Schwach (1.4)", "Akzeptabel (2.0)", "Gut (3.0)"), hjust = 0) +
      coord_flip() +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom"
      )
    
    # Speichere den Plot
    rpd_plot_file <- file.path(metrics_plot_dir, "rpd_comparison.png")
    ggsave(rpd_plot_file, plot = p_rpd, width = 10, height = 6 + 0.25 * length(unique(combined_plsr$Parameter)), dpi = 300)
    cat("RPD-Vergleichsplot gespeichert in:", rpd_plot_file, "\n")
    
    # Ergebnistabelle: Welcher Algorithmus ist für jeden Parameter am besten?
    best_algo <- combined_plsr %>%
      group_by(Parameter) %>%
      summarize(
        Best_R2_Algo = Algorithm[which.max(R_Squared)],
        Best_R2_Value = max(R_Squared),
        Best_RMSE_Algo = Algorithm[which.min(RMSE)],
        Best_RMSE_Value = min(RMSE),
        Best_RPD_Algo = Algorithm[which.max(RPD)],
        Best_RPD_Value = max(RPD),
        .groups = "drop"
      ) %>%
      arrange(desc(Best_R2_Value))
    
    # Speichere die Best-Algorithmus-Tabelle
    best_algo_file <- file.path(plsr_dir, "best_algorithm_per_parameter.csv")
    write_csv(best_algo, best_algo_file)
    cat("Tabelle der besten Algorithmen pro Parameter gespeichert in:", best_algo_file, "\n")
  } else {
    cat("Warnung: Keine PLSR-Ergebnisse zum Plotten gefunden.\n")
  }
}

# Funktion zum Vergleichen der Algorithmen anhand der PLSR-Ergebnisse
compare_algorithms <- function(plsr_results_file, output_dir) {
  # Erstelle Ausgabeverzeichnis, falls es nicht existiert
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Lade die PLSR-Ergebnisse
  if (!file.exists(plsr_results_file)) {
    stop("PLSR-Ergebnisse nicht gefunden: ", plsr_results_file)
  }
  
  results <- read.csv(plsr_results_file)
  
  # Überprüfe, ob die Daten im erwarteten Format sind
  required_cols <- c("Parameter", "Algorithm", "R_squared", "RMSE", "MAE", "RPD")
  if (!all(required_cols %in% colnames(results))) {
    stop("PLSR-Ergebnisdatei hat nicht das erwartete Format.")
  }
  
  # Überprüfe, ob es mehrere Algorithmen gibt
  if (length(unique(results$Algorithm)) < 2) {
    stop("Mindestens zwei verschiedene Algorithmen werden für den Vergleich benötigt.")
  }
  
  # Formatiere einige Spalten (falls nötig)
  results <- results %>%
    mutate(R_squared = as.numeric(R_squared),
           RMSE = as.numeric(RMSE),
           MAE = as.numeric(MAE),
           RPD = as.numeric(RPD))
  
  # Erstelle eine zusammenfassende Tabelle mit den besten Algorithmen für jeden Parameter
  best_algorithms <- results %>%
    group_by(Parameter) %>%
    summarise(
      Best_R_squared_Algorithm = Algorithm[which.max(R_squared)],
      Best_R_squared = max(R_squared),
      Best_RMSE_Algorithm = Algorithm[which.min(RMSE)],
      Best_RMSE = min(RMSE),
      Best_MAE_Algorithm = Algorithm[which.min(MAE)],
      Best_MAE = min(MAE),
      Best_RPD_Algorithm = Algorithm[which.max(RPD)],
      Best_RPD = max(RPD)
    )
  
  # Speichere die Zusammenfassung
  write.csv(best_algorithms, file.path(output_dir, "best_algorithms.csv"), row.names = FALSE)
  
  # Erstelle einen Barplot für R² (alle Parameter und Algorithmen)
  p1 <- ggplot(results, aes(x = Parameter, y = R_squared, fill = Algorithm)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Vergleich der Algorithmen - R²", 
         x = "Bodenparameter", y = "R²") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, "r_squared_comparison.png"), p1, width = 12, height = 8)
  
  # Erstelle einen Barplot für RMSE
  p2 <- ggplot(results, aes(x = Parameter, y = RMSE, fill = Algorithm)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Vergleich der Algorithmen - RMSE", 
         x = "Bodenparameter", y = "RMSE") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, "rmse_comparison.png"), p2, width = 12, height = 8)
  
  # Erstelle einen Barplot für RPD
  p3 <- ggplot(results, aes(x = Parameter, y = RPD, fill = Algorithm)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Vergleich der Algorithmen - RPD", 
         x = "Bodenparameter", y = "RPD") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, "rpd_comparison.png"), p3, width = 12, height = 8)
  
  # Erstelle eine Heatmap der R²-Werte
  heatmap_data <- results %>%
    select(Parameter, Algorithm, R_squared) %>%
    spread(Algorithm, R_squared)
  
  # Konvertiere zurück zu long format für ggplot
  heatmap_long <- heatmap_data %>%
    gather(key = "Algorithm", value = "R_squared", -Parameter)
  
  p4 <- ggplot(heatmap_long, aes(x = Algorithm, y = Parameter, fill = R_squared)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue") +
    labs(title = "Heatmap der R²-Werte", x = "Algorithmus", y = "Parameter") +
    theme_minimal() +
    geom_text(aes(label = round(R_squared, 3)), color = "black", size = 3)
  
  ggsave(file.path(output_dir, "r_squared_heatmap.png"), p4, width = 10, height = 8)
  
  # Erstelle eine Zusammenfassung für jeden Algorithmus
  algorithm_summary <- results %>%
    group_by(Algorithm) %>%
    summarise(
      Avg_R_squared = mean(R_squared),
      Avg_RMSE = mean(RMSE),
      Avg_MAE = mean(MAE),
      Avg_RPD = mean(RPD),
      Best_Parameters = sum(Algorithm == best_algorithms$Best_R_squared_Algorithm)
    )
  
  write.csv(algorithm_summary, file.path(output_dir, "algorithm_summary.csv"), row.names = FALSE)
  
  # Erstelle Plots für die Zusammenfassung
  p5 <- ggplot(algorithm_summary, aes(x = Algorithm, y = Avg_R_squared)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "Durchschnittliches R² pro Algorithmus", 
         x = "Algorithmus", y = "Durchschnittliches R²") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "avg_r_squared.png"), p5, width = 8, height = 6)
  
  p6 <- ggplot(algorithm_summary, aes(x = Algorithm, y = Best_Parameters)) +
    geom_bar(stat = "identity", fill = "darkgreen") +
    labs(title = "Anzahl der Parameter, bei denen der Algorithmus am besten abschneidet", 
         x = "Algorithmus", y = "Anzahl Parameter") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "best_algorithm_counts.png"), p6, width = 8, height = 6)
  
  # Erstelle einen Boxplot für jede Metrik
  metrics_long <- results %>%
    select(Parameter, Algorithm, R_squared, RMSE, MAE, RPD) %>%
    gather(key = "Metric", value = "Value", -Parameter, -Algorithm)
  
  p7 <- ggplot(metrics_long, aes(x = Algorithm, y = Value)) +
    geom_boxplot(aes(fill = Algorithm)) +
    facet_wrap(~Metric, scales = "free_y") +
    labs(title = "Verteilung der Metriken pro Algorithmus") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, "metrics_boxplot.png"), p7, width = 10, height = 8)
  
  # Erstelle eine Tabelle mit allen Metriken für jeden Parameter und Algorithmus
  all_metrics <- results %>%
    select(Parameter, Algorithm, R_squared, RMSE, MAE, RPD) %>%
    arrange(Parameter, desc(R_squared))
  
  write.csv(all_metrics, file.path(output_dir, "all_metrics.csv"), row.names = FALSE)
  
  cat("Algorithmus-Vergleich abgeschlossen. Ergebnisse gespeichert in", output_dir, "\n")
  
  return(list(best_algorithms = best_algorithms, 
              algorithm_summary = algorithm_summary,
              all_metrics = all_metrics))
}

# Hauptfunktion zum Ausführen des Vergleichs
run_algorithm_comparison <- function(plsr_results_dir, output_dir) {
  # Lade die PLSR-Ergebnisse
  plsr_results_file <- file.path(plsr_results_dir, "all_plsr_results.csv")
  
  # Führe den Vergleich aus
  results <- compare_algorithms(plsr_results_file, output_dir)
  
  return(results)
}

# HINWEIS: Auto-Execution auskommentiert - Script wird nur über master_pipeline.R gesteuert
# # Ausführen der Vergleichsfunktion
# result <- compare_algorithm_correlations(correlation_threshold = 0.6, min_algorithms = 2)
# 
# # Ausgabe der ersten Zeilen des Ergebnisses für schnellen Überblick
# if (nrow(result) > 0) {
#   cat("\nBeispiel-Ergebnisse (erste 10 Zeilen):\n")
#   print(head(result, 10))
# } else {
#   cat("\nKeine Ergebnisse gefunden, die den Kriterien entsprechen.\n")
# }
