# Script für Filter und PLSR über mehrere Algorithmen

library(caret)
library(pls)
library(dplyr)
library(readr)
library(glmnet)  # Für Ridge/Lasso Regression

# --- Hilfsfunktionen für Strata-Qualität ---
compute_balance_score <- function(strata, train_idx, test_idx) {
  train_tab <- prop.table(table(strata[train_idx]))
  test_tab  <- prop.table(table(strata[test_idx]))
  all_levels <- union(names(train_tab), names(test_tab))
  train_vec <- train_tab[all_levels]; train_vec[is.na(train_vec)] <- 0
  test_vec  <- test_tab[all_levels];  test_vec[is.na(test_vec)]  <- 0
  1 - sum(abs(train_vec - test_vec))/2
}

compute_representation_score <- function(strata, train_idx, test_idx) {
  train_levels <- unique(strata[train_idx])
  test_levels  <- unique(strata[test_idx])
  all_levels   <- unique(strata)
  mean(all_levels %in% train_levels & all_levels %in% test_levels)
}

# --- Verbesserte stratified_random_split Funktion ---
stratified_random_split <- function(data, response_var, split_ratio = 1) {
  set.seed(696)
  n <- nrow(data)
  balance_score <- NA
  representation_score <- NA
  stratification_method <- "unknown"
  n_strata <- NA
  
  y <- data[[response_var]]
  
  # Spezielle Behandlung für BS Parameter
  if(response_var == "BS") {
    # Feste Grenzen für BS: 0-33%, 34-66%, 67-100%
    y_clean <- y[!is.na(y)]
    if(length(y_clean) >= 10) {
      bs_breaks <- c(0, 33, 66, 100)
      
      # Prüfe, ob die Daten die Breaks sinnvoll abdecken
      y_range <- range(y_clean)
      if(y_range[1] >= bs_breaks[2] || y_range[2] <= bs_breaks[2]) {
        cat("WARNUNG: BS-Werte liegen außerhalb des erwarteten Bereichs für feste Stratifizierung.\n")
        cat("  BS-Bereich:", y_range[1], "bis", y_range[2], "\n")
        # Fallback: keine BS-spezifische Behandlung, nutze normale Tertile weiter unten
      } else {
        y_strata <- tryCatch({
          cut(y, breaks = bs_breaks, 
              labels = c("niedrig", "mittel", "hoch"),
              include.lowest = TRUE, right = TRUE)
        }, error = function(e) {
          cat("FEHLER bei cut() für BS:", e$message, "\n")
          NULL
        })
        
        # Prüfe ob cut() erfolgreich war und alle Strata besetzt sind
        if(!is.null(y_strata) && !all(is.na(y_strata))) {
          strata_counts <- table(y_strata, useNA = "no")
          if(length(strata_counts) < 2) {
            cat("WARNUNG: BS-Stratifizierung ergibt nur", length(strata_counts), "Stratum/Strata.\n")
            cat("  Fallback auf Tertil-Stratifizierung.\n")
            y_strata <- NULL  # Trigger Fallback
          }
        } else {
          cat("WARNUNG: cut() für BS fehlgeschlagen, Fallback auf Tertil-Stratifizierung.\n")
          y_strata <- NULL
        }
        
        if(!is.null(y_strata)) {
          # NEUE ÄNDERUNG: Manuelle Stratifizierung auch für BS
          strata_table <- table(y_strata)
          train_idx <- c()
          
          # Für jedes Stratum sicherstellen, dass Samples in beiden Splits sind
          for (stratum in names(strata_table)) {
            stratum_indices <- which(y_strata == stratum)
            stratum_size <- length(stratum_indices)
            
            # Mindestens 1 Sample für Validierung reservieren
            train_size <- min(round(stratum_size * split_ratio), stratum_size - 1)
            if (train_size < 1) train_size <- 1  # Mindestens 1 für Training
            
            # Zufällige Auswahl für Training
            train_from_stratum <- sample(stratum_indices, size = train_size)
            train_idx <- c(train_idx, train_from_stratum)
          }
          
          test_idx <- setdiff(seq_len(n), train_idx)
          stratification_method <- "bs_fixed_balanced"
          n_strata <- 3
          
          balance_score <- compute_balance_score(y_strata, train_idx, test_idx)
          representation_score <- compute_representation_score(y_strata, train_idx, test_idx)
          
          cat("BS-Spezifische Stratifizierung (0-33%, 34-66%, 67-100%):\n")
          cat("  Training Strata:"); print(table(y_strata[train_idx], useNA = "ifany"))
          cat("  Test Strata:"); print(table(y_strata[test_idx], useNA = "ifany"))
          cat(sprintf("  Balance Score: %.3f, Representation Score: %.3f\n", 
                     balance_score, representation_score))
          return(list(
            train_indices = train_idx, 
            test_indices = test_idx,
            balance_score = balance_score,
            representation_score = representation_score,
            stratification_method = stratification_method,
            n_strata = n_strata,
            strata = y_strata
          ))
        }
      }
    }
  }
  
  if (is.numeric(y) && length(unique(y[!is.na(y)])) > 6) {
    # Für numerische Variablen mit ausreichend Variabilität: 
    y_clean <- y[!is.na(y)]
    if (length(y_clean) >= 10) {
      
      # ÄNDERUNG: Zurück zu Tertilen für bessere Abdeckung
      tertiles <- quantile(y_clean, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
      
      # Prüfe ob Tertile eindeutig sind
      if(length(unique(tertiles)) < 4) {
        cat("WARNUNG: Tertile nicht eindeutig für", response_var, "- zu wenig Variabilität in den Daten.\n")
        cat("  Tertile:", paste(tertiles, collapse=", "), "\n")
        # Fallback weiter unten
      } else {
        y_strata <- tryCatch({
          cut(y, breaks = tertiles, 
              labels = c("niedrig", "mittel", "hoch"), 
              include.lowest = TRUE, right = TRUE)
        }, error = function(e) {
          cat("FEHLER bei cut() für Tertile:", e$message, "\n")
          NULL
        })
        
        # Prüfe ob cut() erfolgreich war
        if(!is.null(y_strata) && !all(is.na(y_strata))) {
          # ÄNDERUNG: Manuelle Stratifizierung, die Repräsentation in beiden Teilen sicherstellt
          strata_table <- table(y_strata)
          train_idx <- c()
          
          # Für jedes Stratum sicherstellen, dass Samples in beiden Splits sind
          for (stratum in names(strata_table)) {
            stratum_indices <- which(y_strata == stratum)
            stratum_size <- length(stratum_indices)
            
            # Mindestens 1 Sample für Validierung reservieren
            train_size <- min(round(stratum_size * split_ratio), stratum_size - 1)
            if (train_size < 1) train_size <- 1  # Mindestens 1 für Training
            
            # Zufällige Auswahl für Training
            train_from_stratum <- sample(stratum_indices, size = train_size)
            train_idx <- c(train_idx, train_from_stratum)
          }
          
          test_idx <- setdiff(seq_len(n), train_idx)
          stratification_method <- "tertile_balanced"
          n_strata <- 3
          
          # Rest des Codes bleibt gleich
          balance_score <- compute_balance_score(y_strata, train_idx, test_idx)
          representation_score <- compute_representation_score(y_strata, train_idx, test_idx)
          cat("TERTIL-Stratifizierung von", response_var, ":\n")
          cat("  Training Strata:"); print(table(y_strata[train_idx], useNA = "ifany"))
          cat("  Test Strata:"); print(table(y_strata[test_idx], useNA = "ifany"))
          cat(sprintf("  Balance Score: %.3f, Representation Score: %.3f\n", balance_score, representation_score))
          return(list(
            train_indices = train_idx, 
            test_indices = test_idx,
            balance_score = balance_score,
            representation_score = representation_score,
            stratification_method = stratification_method,
            n_strata = n_strata,
            strata = y_strata
          ))
        }
      }
    }
  }
  
  # Fallback: Stratifizierung nach biodiv (falls vorhanden)
  if ("biodiv" %in% colnames(data) && length(unique(data$biodiv[!is.na(data$biodiv)])) > 1) {
    idx <- caret::createDataPartition(data$biodiv, p = split_ratio, list = FALSE)
    train_idx <- as.vector(idx)
    test_idx <- setdiff(seq_len(n), train_idx)
    stratification_method <- "biodiv"
    n_strata <- length(unique(data$biodiv))
    balance_score <- compute_balance_score(data$biodiv, train_idx, test_idx)
    representation_score <- compute_representation_score(data$biodiv, train_idx, test_idx)
    cat("FALLBACK: Stratifizierung nach biodiv (Landnutzung):\n")
    cat("  Training - pf:", sum(data$biodiv[train_idx] == "pf", na.rm = TRUE), 
        ", df:", sum(data$biodiv[train_idx] == "df", na.rm = TRUE), "\n")
    cat("  Test - pf:", sum(data$biodiv[test_idx] == "pf", na.rm = TRUE), 
        ", df:", sum(data$biodiv[test_idx] == "df", na.rm = TRUE), "\n")
    cat(sprintf("  Balance Score: %.3f, Representation Score: %.3f\n", balance_score, representation_score))
    return(list(
      train_indices = train_idx, 
      test_indices = test_idx,
      balance_score = balance_score,
      representation_score = representation_score,
      stratification_method = stratification_method,
      n_strata = n_strata,
      strata = data$biodiv
    ))
  }
  
  # Letzter Fallback: Einfache Zufallsaufteilung
  train_idx <- sample(seq_len(n), size = round(n * split_ratio))
  test_idx <- setdiff(seq_len(n), train_idx)
  stratification_method <- "random"
  n_strata <- 1
  strata_random <- rep("all", n)
  balance_score <- compute_balance_score(strata_random, train_idx, test_idx)
  representation_score <- compute_representation_score(strata_random, train_idx, test_idx)
  cat("WARNUNG: Fallback auf einfachen zufälligen Split (keine Stratifizierung möglich)\n")
  cat(sprintf("  Balance Score: %.3f, Representation Score: %.3f\n", balance_score, representation_score))
  return(list(
    train_indices = train_idx, 
    test_indices = test_idx,
    balance_score = balance_score,
    representation_score = representation_score,
    stratification_method = stratification_method,
    n_strata = n_strata,
    strata = strata_random
  ))
}

# Function to create prediction scatter plots
create_prediction_plots <- function(actual_train, pred_train, actual_valid, pred_valid, 
                                   param, param_dir, algo_name, 
                                   rsquared_train, rmse_train, rpd_train, ncomp_cv,
                                   rsquared_valid, rmse_valid, rpd_valid) {
  
  # Create PNG file
  png(file.path(param_dir, paste0(param, "_prediction_plot.png")), 
     width = 800, height = 600, res = 100)
  
  # Setup layout based on validation data availability
  if(!is.na(rsquared_valid) && length(actual_valid) > 1) {
    par(mfrow=c(1,2))
  }
  
  # Training plot
  plot(actual_train, pred_train, 
      main = paste0(algo_name, ": ", param, " - Training"), 
      xlab = paste("Gemessener Wert:", param), 
      ylab = "Vorhergesagter Wert",
      pch = 19, col = "blue")
  abline(0, 1, col = "red", lwd = 2)
  
  # Add metrics text
  usr <- par("usr")
  x_pos <- usr[1] + 0.05 * (usr[2] - usr[1])
  y_pos <- usr[4] - 0.05 * (usr[4] - usr[3])
  metrics_text <- sprintf("R² = %.3f\nRMSE = %.3f\nRPD = %.3f\nKomp. = %d", 
                         rsquared_train, rmse_train, rpd_train, ncomp_cv)
  text(x = x_pos, y = y_pos, labels = metrics_text, adj = c(0, 1), cex = 1.2, font = 2)
  
  # Validation plot if available
  if(!is.na(rsquared_valid) && length(actual_valid) > 1) {
    plot(actual_valid, pred_valid, 
        main = paste0(algo_name, ": ", param, " - Validierung"), 
        xlab = paste("Gemessener Wert:", param), 
        ylab = "Vorhergesagter Wert",
        pch = 19, col = "green4")
    abline(0, 1, col = "red", lwd = 2)
    
    # Add metrics text
    usr <- par("usr")
    x_pos <- usr[1] + 0.05 * (usr[2] - usr[1])
    y_pos <- usr[4] - 0.05 * (usr[4] - usr[3])
    metrics_text <- sprintf("R² = %.3f\nRMSE = %.3f\nRPD = %.3f", 
                           rsquared_valid, rmse_valid, rpd_valid)
    text(x = x_pos, y = y_pos, labels = metrics_text, adj = c(0, 1), cex = 1.2, font = 2)
  }
  
  dev.off()
  cat("Prediction plot saved to", file.path(param_dir, paste0(param, "_prediction_plot.png")), "\n")
}

# Verbesserte Funktion zur Visualisierung der PLSR-Komponentengewichte
visualize_model_components <- function(pls_model, ncomp_cv, param, param_dir, wavelength_ref) {
  # Extrahiere die Gewichte der Komponenten
  weights_mat <- pls_model$loading.weights[, 1:ncomp_cv, drop = FALSE]
  
  # Speichere Loadings und Gewichte für Referenz
  loadings_mat <- pls_model$loadings[, 1:ncomp_cv, drop = FALSE]
  write.csv(as.data.frame(loadings_mat), 
           file.path(param_dir, paste0(param, "_plsr_loadings_debug.csv")), 
           row.names = TRUE)
  write.csv(as.data.frame(weights_mat), 
           file.path(param_dir, paste0(param, "_plsr_loading_weights_debug.csv")), 
           row.names = TRUE)
  
  # Erstelle Plots für jede Komponente
  for (comp in 1:ncomp_cv) {
    weights_df <- data.frame(
      Band = rownames(weights_mat),
      Weight = as.numeric(weights_mat[, comp])
    )
    
    # Verbinde mit Wellenlängen-Referenz
    if (!is.null(wavelength_ref)) {
      # Extrahiere die Band-Nummer aus dem Band-Namen
      if(grepl("_corr$", weights_df$Band[1])) {
        # Für Band-Namen wie "Band1_corr"
        weights_df$BandNum <- as.numeric(gsub("Band([0-9]+)_corr", "\\1", weights_df$Band))
      } else {
        # Für andere Band-Namen wie "corrected_spectrum_Band1"
        weights_df$BandNum <- as.numeric(gsub(".*Band([0-9]+)$", "\\1", weights_df$Band))
      }
      
      # Verbinde mit Wellenlängen-Referenz - FIXED join syntax
      weights_df <- dplyr::left_join(
        weights_df,
        wavelength_ref %>% mutate(BandNum = band_number),
        by = "BandNum"
      )
      
      # Filtere fehlende Werte und sortiere nach Wellenlänge
      weights_df <- weights_df[!is.na(weights_df$wavelength_nm), ]
      # Hier zusätzliche Filterung hinzufügen, um die gleichen Wellenlängen wie in der Korrelationstabelle auszuschließen
      weights_df <- weights_df %>%
        dplyr::filter(wavelength_nm < 1800 | wavelength_nm > 2010)
      weights_df <- weights_df[order(weights_df$wavelength_nm), ]
      
      xvar <- weights_df$wavelength_nm
      xlab <- "Wellenlänge (nm)"
    } else {
      # Fallback: Verwende BandNum als X-Achse
      weights_df$BandNum <- as.numeric(gsub(".*Band([0-9]+).*", "\\1", weights_df$Band))
      weights_df <- weights_df[!is.na(weights_df$BandNum), ]
      weights_df <- weights_df[order(weights_df$BandNum), ]
      xvar <- weights_df$BandNum
      xlab <- "Band Nummer"
    }
    
    if (nrow(weights_df) == 0) {
      cat("WARNUNG: Keine gültigen Wellenlängendaten für Komponente", comp, "\n")
      next
    }
    
    # Erstelle Plot
    png_file <- file.path(param_dir, paste0(param, "_plsr_component_", comp, ".png"))
    png(png_file, width = 1200, height = 600, res = 120)
    par(mar = c(5, 5, 4, 2))
    plot(xvar, weights_df$Weight, type = "h", lwd = 2, col = "darkblue",
         main = paste("PLSR-Gewichte für", param, "\nKomponente", comp),
         xlab = xlab, ylab = "PLSR-Gewicht", cex.main = 1.3, cex.lab = 1.2)
    abline(h = 0, col = "grey60", lty = 2)
    grid(col = "grey90", lty = 3)
    
    # Füge wichtige Informationen hinzu
    max_weight_idx <- which.max(abs(weights_df$Weight))
    if (length(max_weight_idx) > 0) {
      max_band <- if (!is.null(wavelength_ref)) {
        paste0(round(xvar[max_weight_idx]), " nm")
      } else {
        paste0("Band ", xvar[max_weight_idx])
      }
      
      legend("topright", 
             legend = c(paste("Max:", round(max(abs(weights_df$Weight)), 4)),
                       paste("Bei:", max_band)),
             cex = 0.8, bg = "white")
    }
    dev.off()
    cat("Komponenten-Plot", comp, "gespeichert als", png_file, "\n")
  }
  
  # Erstelle kombinierten Plot für alle Komponenten
  if (ncomp_cv > 1) {
    png_file <- file.path(param_dir, paste0(param, "_plsr_all_components.png"))
    png(png_file, width = 1200, height = 600, res = 120)
    
    # Setze Farben für mehrere Komponenten
    colors <- c("darkblue", "darkred", "darkgreen", "darkorange", "purple")[1:ncomp_cv]
    
    # Bereite Daten für den ersten Plot vor
    weights_df <- data.frame(
      Band = rownames(weights_mat),
      Weight = as.numeric(weights_mat[, 1])
    )
    
    # Verbinde mit Wellenlängen
    if (!is.null(wavelength_ref)) {
      # Extrahiere die Band-Nummer
      if(grepl("_corr$", weights_df$Band[1])) {
        weights_df$BandNum <- as.numeric(gsub("Band([0-9]+)_corr", "\\1", weights_df$Band))
      } else {
        weights_df$BandNum <- as.numeric(gsub(".*Band([0-9]+)$", "\\1", weights_df$Band))
      }
      
      weights_df <- dplyr::left_join(
        weights_df,
        wavelength_ref %>% mutate(BandNum = band_number),
        by = "BandNum"
      )
      weights_df <- weights_df[!is.na(weights_df$wavelength_nm), ]
      # Hier auch für den kombinierten Plot die Filterung hinzufügen
      weights_df <- weights_df %>%
        dplyr::filter(wavelength_nm < 1800 | wavelength_nm > 2010)
      weights_df <- weights_df[order(weights_df$wavelength_nm), ]
      
      xvar <- weights_df$wavelength_nm
      xlab <- "Wellenlänge (nm)"
    } else {
      weights_df$BandNum <- as.numeric(gsub(".*Band([0-9]+).*", "\\1", weights_df$Band))
      weights_df <- weights_df[!is.na(weights_df$BandNum), ]
      weights_df <- weights_df[order(weights_df$BandNum), ]
      xvar <- weights_df$BandNum
      xlab <- "Band Nummer"
    }
    
    # Initialisiere Plot
    plot(xvar, weights_df$Weight, type = "n", 
         main = paste("PLSR-Gewichte für", param, "- Alle Komponenten"),
         xlab = xlab, ylab = "PLSR-Gewicht", cex.main = 1.3, cex.lab = 1.2,
         ylim = c(-max(abs(weights_mat)), max(abs(weights_mat))))
    abline(h = 0, col = "grey60", lty = 2)
    grid(col = "grey90", lty = 3)
    
    # Plotte jede Komponente
    legend_text <- character(ncomp_cv)
    for (comp in 1:ncomp_cv) {
      weights_df <- data.frame(
        Band = rownames(weights_mat),
        Weight = as.numeric(weights_mat[, comp])
      )
      
      # Verbinde mit Wellenlängen
      if (!is.null(wavelength_ref)) {
        if(grepl("_corr$", weights_df$Band[1])) {
          weights_df$BandNum <- as.numeric(gsub("Band([0-9]+)_corr", "\\1", weights_df$Band))
        } else {
          weights_df$BandNum <- as.numeric(gsub(".*Band([0-9]+)$", "\\1", weights_df$Band))
        }
        
        weights_df <- dplyr::left_join(
          weights_df,
          wavelength_ref %>% mutate(BandNum = band_number),
          by = "BandNum"
        )
        weights_df <- weights_df[!is.na(weights_df$wavelength_nm), ]
        # Hier auch für jede Komponentenlinie die Filterung anwenden
        weights_df <- weights_df %>%
          dplyr::filter(wavelength_nm < 1800 | wavelength_nm > 2010)
        weights_df <- weights_df[order(weights_df$wavelength_nm), ]
        
        xvar <- weights_df$wavelength_nm
      } else {
        weights_df$BandNum <- as.numeric(gsub(".*Band([0-9]+).*", "\\1", weights_df$Band))
        weights_df <- weights_df[!is.na(weights_df$BandNum), ]
        weights_df <- weights_df[order(weights_df$BandNum), ]
        xvar <- weights_df$BandNum
      }
      
      lines(xvar, weights_df$Weight, type = "l", lwd = 2, col = colors[comp])
      legend_text[comp] <- paste("Komponente", comp)
    }
    
    # Füge Legende hinzu
    legend("topright", legend = legend_text, col = colors, lwd = 2, cex = 0.8, bg = "white")
    dev.off()
    cat("Plot für alle Komponenten gespeichert als", png_file, "\n")
  }
  
  # Exportiere Gewichte mit Wellenlängen
  weights_export <- data.frame(
    Band = rownames(weights_mat),
    stringsAsFactors = FALSE
  )
  
  for (comp in 1:ncomp_cv) {
    weights_export[[paste0("Component_", comp)]] <- weights_mat[, comp]
  }
  
  if (!is.null(wavelength_ref)) {
    # Extrahiere die Band-Nummer für den Join
    if(grepl("_corr$", weights_export$Band[1])) {
      weights_export$BandNum <- as.numeric(gsub("Band([0-9]+)_corr", "\\1", weights_export$Band))
    } else {
      weights_export$BandNum <- as.numeric(gsub(".*Band([0-9]+)$", "\\1", weights_export$Band))
    }
    
    weights_export <- dplyr::left_join(
      weights_export,
      wavelength_ref %>% mutate(BandNum = band_number),
      by = "BandNum"
    )
    
    # Auch hier die Filterung für die exportierten Daten anwenden
    weights_export <- weights_export %>%
      dplyr::filter(is.na(wavelength_nm) | wavelength_nm < 1800 | wavelength_nm > 2010)
    
    col_order <- c("Band", "BandNum", "wavelength_nm", paste0("Component_", 1:ncomp_cv))
    weights_export <- weights_export[, col_order[col_order %in% colnames(weights_export)]]
  }
  
  write.csv(weights_export, file.path(param_dir, paste0(param, "_plsr_weights.csv")), row.names = FALSE)
  cat("PLSR-Gewichte für alle", ncomp_cv, "Komponenten gespeichert\n")
}

# Hauptfunktion für Filter und PLSR pro Algorithmus
perform_filter_plsr_by_algorithm <- function(
  data_file = "data/final_before_correlation.csv",
  results_base_dir = "results",
  min_bands = 4,
  min_correlation = 0.3,
  use_validation_split = TRUE,    # Option zur externen Validierung
  validation_ratio = 0.3,         # 70% Training, 30% Validierung
  max_components = 5              # Maximale Anzahl an PLSR-Komponenten
) {
  # Daten laden
  cat("Lade Datensatz...\n")
  data <- read_csv(data_file, show_col_types = FALSE)
  
  # Lade Wellenlängen-Referenztabelle
  wavelength_file <- "data/wavelength_reference.csv"
  if (!file.exists(wavelength_file)) {
    wavelength_file <- file.path(dirname(results_base_dir), "data/wavelength_reference.csv")
  }
  if (file.exists(wavelength_file)) {
    wavelength_ref <- read_csv(wavelength_file, show_col_types = FALSE)
    cat("Wellenlängen-Referenz geladen:", nrow(wavelength_ref), "Bänder\n")
  } else {
    cat("WARNUNG: Wellenlängen-Referenzdatei nicht gefunden, fahre ohne Wellenlängeninformationen fort.\n")
    wavelength_ref <- NULL
  }
  
  # Liste aller Algorithmus-Verzeichnisse
  algo_dirs <- list.dirs(results_base_dir, recursive = FALSE)
  algo_dirs <- algo_dirs[basename(algo_dirs) %in% c("BarlowTwin", "Pearson", "Spearman")]
  
  if(length(algo_dirs) == 0) {
    stop("Keine Algorithmus-Verzeichnisse (BarlowTwin, Pearson, Spearman) gefunden.")
  }
  
  cat("Verarbeite folgende Algorithmen:", paste(basename(algo_dirs), collapse=", "), "\n")
  
  # Verarbeite jeden Algorithmus separat
  for(algo_dir in algo_dirs) {
    algo_name <- basename(algo_dir)
    cat("\n\n===== VERARBEITE ALGORITHMUS:", algo_name, "=====\n")
    
    # Lade verifizierte Korrelationen für diesen Algorithmus
    verified_file <- file.path(algo_dir, "verified_significant_correlations.csv")
    
    # Parameter-Liste, die verarbeitet werden soll
    parameters <- c()
    
    if(file.exists(verified_file)) {
      # Wenn verified_file existiert, nutze es
      verified_data <- read_csv(verified_file, show_col_types = FALSE)
      cat("Korrelationen geladen:", nrow(verified_data), "Einträge\n")
      parameters <- unique(verified_data$Parameter)
    } else {
      # Alternative: Suche nach Dateien die auf _correlation_table_corrected.csv enden
      cat("Keine verified_significant_correlations.csv gefunden. Suche nach *_correlation_table_corrected.csv Dateien...\n")
      corr_files <- list.files(algo_dir, pattern = "_correlation_table_corrected\\.csv$", full.names = FALSE)
      
      if(length(corr_files) > 0) {
        # Extrahiere Parameter aus Dateinamen (Format: PARAM_correlation_table_corrected.csv)
        parameters <- sub("_correlation_table_corrected\\.csv$", "", corr_files)
        cat("Parameter aus Korrelationstabellen gefunden:", paste(parameters, collapse=", "), "\n")
      } else {
        cat("WARNUNG: Keine Korrelationstabellen für", algo_name, "gefunden, überspringe.\n")
        next
      }
    }
    
    if(length(parameters) == 0) {
      cat("WARNUNG: Keine Parameter für", algo_name, "gefunden, überspringe.\n")
      next
    }
    
    # Erstelle filter+plsr-Verzeichnis für diesen Algorithmus
    filter_dir <- file.path(algo_dir, "filter+plsr")
    if(!dir.exists(filter_dir)) {
      dir.create(filter_dir, recursive = TRUE)
      cat("Filter-PLSR-Verzeichnis erstellt:", filter_dir, "\n")
    }
    
    cat("Gefundene Parameter:", paste(parameters, collapse=", "), "\n")
    
    # Ergebnisliste für diesen Algorithmus - WICHTIG: Hier initialisieren
    algo_results <- list()
    
    # Verarbeite jeden Parameter einzeln
    for(param in parameters) {
      cat("\n------ Verarbeite Parameter:", param, "------\n")
      
      # Parameter-spezifisches Verzeichnis
      param_dir <- file.path(filter_dir, param)
      if(!dir.exists(param_dir)) {
        dir.create(param_dir, recursive = TRUE)
      }
      
      # Lade die parameterspezifische Korrelationstabelle
      corr_table_file <- file.path(algo_dir, paste0(param, "_correlation_table_corrected.csv"))
      
      if(!file.exists(corr_table_file)) {
        cat("WARNUNG: Keine Korrelationstabelle für", param, "gefunden, überspringe.\n")
        next
      }
      
      cat("Lade Korrelationstabelle von:", corr_table_file, "\n")
      corr_table <- read_csv(corr_table_file, show_col_types = FALSE)
      
      # Extrahiere signifikante Bänder nach bestimmten Kriterien
      # Zuerst: Falls Wellenlängenreferenz vorhanden, filtere problematische Wellenlängen
      if (!is.null(wavelength_ref)) {
        # Join mit wavelength_ref um Wellenlängeninformationen zu erhalten und Filtere problematische Regionen
        # FIXED: changed => to = for join
        corr_table <- corr_table %>%
          dplyr::left_join(wavelength_ref, by = c("Band" = "band_name")) %>%
          dplyr::filter(is.na(wavelength_nm) | wavelength_nm < 1800 | wavelength_nm > 2010)
        cat("Wellenlängenreferenz zugeordnet (Bänder 1800-2010nm ausgeschlossen).\n")
      }
      
      # Dann: Filtere nach Korrelation und Signifikanz in beiden Splits
      significant_bands <- corr_table %>% 
        dplyr::filter(
          abs(Correlation_Split1) >= min_correlation,
          abs(Correlation_Split2) >= min_correlation,
          Split1_p_value < 0.05,
          Split2_p_value < 0.05
        ) %>%
        dplyr::arrange(desc(abs(Correlation)))
      
      if(nrow(significant_bands) < min_bands) {
        cat("WARNUNG: Zu wenige signifikante Bänder für", param, 
            ":", nrow(significant_bands), "<", min_bands, ", überspringe.\n")
        next
      }
      
      # Liste der relevanten Bänder
      relevant_bands <- significant_bands %>% 
        dplyr::pull(Band) %>% 
        unique()
      
      cat("Anzahl signifikanter Bänder vor NA-Filter:", length(relevant_bands), "\n")
      
      # **NEUE FILTERUNG: Entferne Bänder mit NA-Werten**
      # Prüfe, ob alle Bänder in den Daten vorhanden sind
      missing_bands <- relevant_bands[!relevant_bands %in% colnames(data)]
      if(length(missing_bands) > 0) {
        cat("WARNUNG: Einige Bänder fehlen in den Daten:", 
            paste(missing_bands, collapse=", "), "\n")
        relevant_bands <- relevant_bands[relevant_bands %in% colnames(data)]
      }
      
      # Erstelle temporären Datensatz nur mit dem Parameter und relevanten Bändern
      temp_data <- data %>%
        dplyr::select(all_of(c(param, relevant_bands))) %>%
        dplyr::filter(!is.na(!!sym(param)))
      
      cat("Datenpunkte mit gültigem", param, ":", nrow(temp_data), "\n")
      
      # Prüfe jedes Band auf NA-Werte und entferne Bänder mit NAs
      bands_with_na <- character(0)
      for(band in relevant_bands) {
        na_count <- sum(is.na(temp_data[[band]]))
        if(na_count > 0) {
          bands_with_na <- c(bands_with_na, band)
          cat("Band", band, "hat", na_count, "NA-Werte von", nrow(temp_data), "-> wird entfernt\n")
        }
      }
      
      # Entferne Bänder mit NA-Werten
      if(length(bands_with_na) > 0) {
        relevant_bands <- setdiff(relevant_bands, bands_with_na)
        cat("Entfernte Bänder mit NA-Werten:", length(bands_with_na), "\n")
        cat("Verbleibende Bänder:", length(relevant_bands), "\n")
        
        # Aktualisiere significant_bands DataFrame
        significant_bands <- significant_bands %>%
          dplyr::filter(Band %in% relevant_bands)
      }
      
      # Prüfe ob noch genügend Bänder übrig sind
      if(length(relevant_bands) < min_bands) {
        cat("WARNUNG: Nach NA-Filter zu wenige Bänder für", param, 
            ":", length(relevant_bands), "<", min_bands, ", überspringe.\n")
        next
      }
      
      cat("Finale Anzahl verwendbarer Bänder:", length(relevant_bands), "\n")
      
      # Speichere die signifikanten Bänder mit ihren Korrelationswerten
      if (!is.null(wavelength_ref)) {
        # Mit Wellenlängeninformationen
        bands_df <- significant_bands %>%
          dplyr::select(Band, Correlation) %>%
          dplyr::left_join(wavelength_ref, by = c("Band" = "band_name")) %>%
          dplyr::select(Band, wavelength_nm, Correlation) %>%
          dplyr::arrange(desc(abs(Correlation)))
        
        cat("Bänder mit Wellenlängeninformationen gespeichert.\n")
      } else {
        # Ohne Wellenlängeninformationen
        bands_df <- significant_bands %>%
          dplyr::select(Band, Correlation) %>%
          dplyr::arrange(desc(abs(Correlation)))
        
        cat("Bänder ohne Wellenlängeninformationen gespeichert.\n")
      }
      
      write.csv(bands_df, file.path(param_dir, "significant_bands.csv"), row.names = FALSE)
      
      # Prüfe, ob der Parameter in den Daten existiert
      if(!param %in% colnames(data)) {
        cat("WARNUNG: Parameter", param, "nicht in Daten gefunden, überspringe.\n")
        next
      }
      
      # Prüfe, ob alle Bänder in den Daten vorhanden sind
      missing_bands <- relevant_bands[!relevant_bands %in% colnames(data)]
      if(length(missing_bands) > 0) {
        cat("WARNUNG: Einige Bänder fehlen in den Daten:", 
            paste(missing_bands, collapse=", "), "\n")
        relevant_bands <- relevant_bands[relevant_bands %in% colnames(data)]
        
        if(length(relevant_bands) < min_bands) {
          cat("Zu wenige gültige Bänder übrig, überspringe Parameter", param, "\n")
          next
        }
      }
      
      # Erstelle den Datensatz für das Modell
      model_data <- data %>%
        dplyr::select(all_of(c(param, relevant_bands))) %>%
        dplyr::filter(!is.na(!!sym(param)))
      
      # Debugging: Prüfe nochmal auf NA-Werte
      na_check <- sapply(relevant_bands, function(band) sum(is.na(model_data[[band]])))
      if(any(na_check > 0)) {
        cat("WARNUNG: Immer noch NA-Werte gefunden in Bändern:", names(na_check[na_check > 0]), "\n")
        # Entferne diese Bänder auch noch
        clean_bands <- names(na_check[na_check == 0])
        relevant_bands <- clean_bands
        model_data <- model_data %>% dplyr::select(all_of(c(param, relevant_bands)))
      }
      
      cat("Finale Modell-Daten: Parameter =", param, ", Bänder =", length(relevant_bands), ", Datenpunkte =", nrow(model_data), "\n")
      
      # Erstelle zusätzlich einen erweiterten Datensatz mit allen relevanten Spalten für die Dokumentation
      extended_data <- data %>%
        dplyr::select(all_of(c("Name", "biodiv", "X", "Y", param, relevant_bands))) %>%
        dplyr::filter(!is.na(!!sym(param)))
      
      # Auch hier NA-Check für extended_data
      for(band in relevant_bands) {
        if(any(is.na(extended_data[[band]]))) {
          cat("WARNUNG: Extended_data hat noch NA-Werte in Band", band, "\n")
        }
      }
      
      if(nrow(model_data) == 0) {
        cat("WARNUNG: Keine gültigen Daten für", param, ", überspringe.\n")
        next
      }
      
      # **ZUSÄTZLICHE VALIDIERUNG**: Prüfe dass alle Datenpunkte komplett sind
      complete_rows <- complete.cases(model_data)
      if(!all(complete_rows)) {
        cat("WARNUNG: Es gibt immer noch unvollständige Zeilen. Anzahl unvollständig:", sum(!complete_rows), "\n")
        model_data <- model_data[complete_rows, ]
        extended_data <- extended_data[complete_rows, ]
        cat("Nach Entfernung unvollständiger Zeilen:", nrow(model_data), "Datenpunkte\n")
      }
      
      # Stratifizierte Aufteilung in Training und Validierung
      if (use_validation_split) {
        cat("Führe stratifizierte Zufallsaufteilung der Daten durch...\n")
        split <- stratified_random_split(model_data, param, split_ratio = 1 - validation_ratio)
        train_data <- model_data[split$train_indices, ]
        valid_data <- model_data[split$test_indices, ]
        cat("Datenaufteilung: Training =", nrow(train_data), ", Validierung =", nrow(valid_data), "\n")
        
        # Speichere Split-Informationen
        rss_info_dir <- file.path(param_dir, "Split_Info")
        if (!dir.exists(rss_info_dir)) dir.create(rss_info_dir, recursive = TRUE)
        
        split_info <- data.frame(
          Parameter = param,
          N_Train = nrow(train_data),
          N_Valid = nrow(valid_data),
          Train_Mean = mean(train_data[[param]], na.rm = TRUE),
          Valid_Mean = mean(valid_data[[param]], na.rm = TRUE),
          Train_SD = sd(train_data[[param]], na.rm = TRUE),
          Valid_SD = sd(valid_data[[param]], na.rm = TRUE),
          Balance_Score = ifelse(is.null(split$balance_score), NA, split$balance_score),
          Representation_Score = ifelse(is.null(split$representation_score), NA, split$representation_score),
          Stratification_Method = ifelse(is.null(split$stratification_method), "unknown", split$stratification_method),
          N_Strata = ifelse(is.null(split$n_strata), NA, split$n_strata)
        )
        
        split_info_file <- file.path(rss_info_dir, paste0(param, "_split_info.csv"))
        write.csv(split_info, split_info_file, row.names = FALSE)
        cat("Split-Informationen gespeichert in:", split_info_file, "\n")
        
        # Speichere detaillierte Strata-Dokumentation
        detailed_strata <- data.frame(
          Parameter = param,
          Name = if("Name" %in% colnames(extended_data)) extended_data$Name else NA,
          biodiv = if("biodiv" %in% colnames(extended_data)) extended_data$biodiv else NA,
          X = extended_data$X,
          Y = extended_data$Y,
          Parameter_Value = extended_data[[param]],
          Split = ifelse(seq_len(nrow(extended_data)) %in% split$train_indices, "Training", "Validierung"),
          stringsAsFactors = FALSE
        )
        
        detailed_strata$Stratum_Type <- ifelse(is.null(split$stratification_method), "unknown", split$stratification_method)
        
        # Stratum_Value entsprechend der Stratifizierungsmethode
        if (!is.null(split$stratification_method)) {
          if (split$stratification_method == "bs_fixed") {
            # Verwende die festen Grenzen für BS
            bs_breaks <- c(0, 33, 66, 100)
            param_strata <- cut(extended_data[[param]], breaks = bs_breaks,
                              labels = c("niedrig", "mittel", "hoch"),
                              include.lowest = TRUE)
            detailed_strata$Stratum_Value <- as.character(param_strata)
          } else if (split$stratification_method == "tertile") {
            y_clean <- extended_data[[param]][!is.na(extended_data[[param]])]
            if (length(y_clean) >= 10) {
              tertiles <- quantile(y_clean, probs = c(0, 0.33, 0.67, 1))
              param_strata <- cut(extended_data[[param]], breaks = tertiles, 
                                 labels = c("niedrig", "mittel", "hoch"), 
                                 include.lowest = TRUE)
              detailed_strata$Stratum_Value <- as.character(param_strata)
            } else {
              detailed_strata$Stratum_Value <- "insufficient_data"
            }
          } else if (split$stratification_method == "biodiv") {
            detailed_strata$Stratum_Value <- extended_data$biodiv
          } else {
            detailed_strata$Stratum_Value <- "no_stratification"
          }
        } else {
          detailed_strata$Stratum_Value <- "unknown_method"
        }
        
        detailed_strata_file <- file.path(rss_info_dir, paste0(param, "_detailed_strata.csv"))
        write.csv(detailed_strata, detailed_strata_file, row.names = FALSE)
        cat("Detaillierte Strata-Dokumentation gespeichert in:", detailed_strata_file, "\n")
      } else {
        # Wenn keine Validierung gewünscht wird, verwenden wir alle Daten für das Training
        train_data <- model_data
        valid_data <- NULL
      }
      
      # PLSR-Modellierung mit Kreuzvalidierung
      cat("Führe PLSR für", param, "aus...\n")
      set.seed(263) # Für reproduzierbare Ergebnisse
      
      tryCatch({
        # Prüfe Verhältnis von Stichprobengröße zu Prädiktoren
        n_samples <- nrow(train_data)
        n_features <- length(relevant_bands)
        
        cat("DEBUG: n_samples =", n_samples, ", n_features =", n_features, "\n")
        
        # Setze max_comp basierend auf Datenverfügbarkeit
        if (n_features >= n_samples * 0.8) {
          max_comp <- min(2, floor(n_samples * 0.2))  # Sehr konservativ bei vielen Prädiktoren
          cat("WARNUNG: Viele Prädiktoren im Verhältnis zu Stichprobengröße. Komponenten auf", max_comp, "begrenzt.\n")
        } else if (n_features > n_samples * 0.5) {
          max_comp <- min(3, floor(n_samples * 0.3))  # Konservativ
          cat("INFO: Moderate Prädiktor-zu-Stichprobe Ratio. Komponenten auf", max_comp, "begrenzt.\n")
        } else {
          max_comp <- min(max_components, n_samples - 2, n_features - 1)  # Standard
        }
        
        # Mindestens 1 Komponente beibehalten
        max_comp <- max(1, max_comp)
        
        # PLSR-Formel erstellen (Parameter ~ Band1 + Band2 + ...)
        formula_str <- paste(param, "~", paste(relevant_bands, collapse = " + "))
        plsr_formula <- as.formula(formula_str)
        
        # PLSR mit Kreuzvalidierung
        pls_model <- plsr(plsr_formula, data = train_data, 
                          validation = "CV", segments = 5,  # 5-fold kreuzvalidierung
                          ncomp = max_comp)
        
        # Bestimme optimale Anzahl an Komponenten
        rmsep_vals <- RMSEP(pls_model)$val[1,,]
        ncomp_cv <- which.min(rmsep_vals[1:max_comp])[1]
        
        cat("Optimale Anzahl an Komponenten:", ncomp_cv, "\n")
        
        # Vorhersagen auf Trainingsdaten
        pred_train <- predict(pls_model, train_data, ncomp = ncomp_cv)[, 1, 1]
        actual_train <- train_data[[param]]
        
        # Nur Zeilen mit vollständigen Werten verwenden
        valid_idx_train <- complete.cases(pred_train, actual_train)
        pred_train_valid <- pred_train[valid_idx_train]
        actual_train_valid <- actual_train[valid_idx_train]
        
        # Metriken für Trainingsdaten
        rsquared_train <- cor(pred_train_valid, actual_train_valid)^2
        rmse_train <- sqrt(mean((pred_train_valid - actual_train_valid)^2))
        mae_train <- mean(abs(pred_train_valid - actual_train_valid))
        rpd_train <- sd(actual_train_valid) / rmse_train
        
        cat("Trainingsmetriken:", 
            "R² =", round(rsquared_train, 3),
            "RMSE =", round(rmse_train, 3),
            "RPD =", round(rpd_train, 3), "\n")
        
        # Initialisierung der Ergebnisliste
        param_results <- list(
          parameter = param,
          method = "PLSR",
          n_bands = length(relevant_bands),
          ncomp = ncomp_cv,
          r_squared_train = rsquared_train,
          rmse_train = rmse_train,
          mae_train = mae_train,
          rpd_train = rpd_train
        )
        
        # Validierung auf externem Datensatz, wenn aktiviert
        if (use_validation_split && !is.null(valid_data)) {
          # Vorhersagen auf Validierungsdaten
          pred_valid <- predict(pls_model, valid_data, ncomp = ncomp_cv)[, 1, 1]
          actual_valid <- valid_data[[param]]
          
          # Nur vollständige Werte verwenden
          valid_idx_val <- complete.cases(pred_valid, actual_valid)
          pred_valid_valid <- pred_valid[valid_idx_val]
          actual_valid_valid <- actual_valid[valid_idx_val]
          
          if (length(actual_valid_valid) > 1) {
            # Metriken für Validierungsdaten
            rsquared_valid <- cor(pred_valid_valid, actual_valid_valid)^2
            rmse_valid <- sqrt(mean((pred_valid_valid - actual_valid_valid)^2))
            mae_valid <- mean(abs(pred_valid_valid - actual_valid_valid))
            rpd_valid <- sd(actual_valid_valid) / rmse_valid
            
            cat("Validierungsmetriken:", 
                "R² =", round(rsquared_valid, 3),
                "RMSE =", round(rmse_valid, 3),
                "RPD =", round(rpd_valid, 3), "\n")
            
            # Füge Validierungsmetriken zur Ergebnisliste hinzu
            param_results$r_squared_valid <- rsquared_valid
            param_results$rmse_valid <- rmse_valid
            param_results$mae_valid <- mae_valid
            param_results$rpd_valid <- rpd_valid
            
            # CSV mit Training- und Validierungsmetriken
            metrics_df <- data.frame(
              Parameter = param,
              Method = "PLSR",
              N_Bands = length(relevant_bands),
              N_Components = ncomp_cv,
              R_Squared_Train = param_results$r_squared_train,
              RMSE_Train = param_results$rmse_train,
              MAE_Train = mae_train,
              RPD_Train = param_results$rpd_train,
              R_Squared_Valid = param_results$r_squared_valid,
              RMSE_Valid = param_results$rmse_valid,
              MAE_Valid = mae_valid,
              RPD_Valid = param_results$rpd_valid
            )
          } else {
            cat("WARNUNG: Nicht genügend gültige Validierungsdaten für", param, "\n")
            # CSV nur mit Trainingsmetriken
            metrics_df <- data.frame(
              Parameter = param,
              Method = "PLSR",
              N_Bands = length(relevant_bands),
              N_Components = ncomp_cv,
              R_Squared = rsquared_train,
              RMSE = rmse_train,
              MAE = mae_train,
              RPD = rpd_train
            )
          }
        } else {
          # Ohne Validierung - nur Trainingsmetriken
          metrics_df <- data.frame(
            Parameter = param,
            Method = "PLSR",
            N_Bands = length(relevant_bands),
            N_Components = ncomp_cv,
            R_Squared = rsquared_train,
            RMSE = rmse_train,
            MAE = mae_train,
            RPD = rpd_train
          )
        }
        
        # Speichere Metriken und Modell
        write.csv(metrics_df, file.path(param_dir, "plsr_metrics.csv"), row.names = FALSE)
        saveRDS(pls_model, file.path(param_dir, paste0(param, "_plsr_model.rds")))
        
        # Speichere Ergebnisse in der Liste
        algo_results[[param]] <- param_results
        
        # Erstelle Vorhersageplots
        create_prediction_plots(
          actual_train_valid, pred_train_valid,
          if (use_validation_split && exists("pred_valid_valid") && length(pred_valid_valid) > 1) actual_valid_valid else NULL,
          if (use_validation_split && exists("pred_valid_valid") && length(pred_valid_valid) > 1) pred_valid_valid else NULL,
          param, param_dir, algo_name,
          rsquared_train, rmse_train, rpd_train, ncomp_cv,
          if (use_validation_split && exists("rsquared_valid")) rsquared_valid else NA,
          if (use_validation_split && exists("rmse_valid")) rmse_valid else NA,
          if (use_validation_split && exists("rpd_valid")) rpd_valid else NA
        )
        
        # Visualisiere die PLSR-Komponenten mit Wellenlängeninformationen
        visualize_model_components(pls_model, ncomp_cv, param, param_dir, wavelength_ref)
        
        # Erfolgsmeldung
        if (use_validation_split && exists("rsquared_valid")) {
          cat("PLSR für", param, "abgeschlossen: Training R²=", round(rsquared_train, 3), 
              ", Validierung R²=", round(rsquared_valid, 3), "\n")
        } else {
          cat("PLSR für", param, "abgeschlossen: R²=", round(rsquared_train, 3), 
              ", RMSE=", round(rmse_train, 3), ", RPD=", round(rpd_train, 3), "\n")
        }
        
      }, error = function(e) {
        cat("FEHLER bei PLSR für", param, ":", e$message, "\n")
      })
    } # End of parameter loop
    
    # Zusammenfassung für diesen Algorithmus
    if (length(algo_results) > 0) {
      # Unterschiedliche Verarbeitung je nach Validierungsoption
      if (use_validation_split && any(sapply(algo_results, function(x) "r_squared_valid" %in% names(x)))) {
        summary_df <- dplyr::bind_rows(lapply(algo_results, function(x) {
          if ("r_squared_valid" %in% names(x)) {
            data.frame(
              Parameter = x$parameter,
              Method = x$method,
              N_Bands = x$n_bands,
              N_Components = x$ncomp,
              R_Squared_Train = x$r_squared_train,
              RMSE_Train = x$rmse_train,
              MAE_Train = x$mae_train, 
              RPD_Train = x$rpd_train,
              R_Squared_Valid = x$r_squared_valid,
              RMSE_Valid = x$rmse_valid,
              MAE_Valid = x$mae_valid,
              RPD_Valid = x$rpd_valid
            )
          } else {
            data.frame(
              Parameter = x$parameter,
              Method = x$method,
              N_Bands = x$n_bands,
              N_Components = x$ncomp,
              R_Squared_Train = x$r_squared_train,
              RMSE_Train = x$rmse_train,
              MAE_Train = x$mae_train,
              RPD_Train = x$rpd_train,
              R_Squared_Valid = NA,
              RMSE_Valid = NA, 
              MAE_Valid = NA,
              RPD_Valid = NA
            )
          }
        }))
        # Sortiere nach RPD der Validierung oder Training
        valid_indices <- !is.na(summary_df$RPD_Valid)
        if (any(valid_indices)) {
          summary_df <- summary_df[order(-summary_df$RPD_Valid), ]
        } else {
          summary_df <- summary_df[order(-summary_df$RPD_Train), ]
        }
      } else {
        # Ohne Validierung - nur Trainingsmetriken
        summary_df <- dplyr::bind_rows(lapply(algo_results, function(x) {
          data.frame(
            Parameter = x$parameter,
            Method = x$method,
            N_Bands = x$n_bands,
            N_Components = x$ncomp,
            R_Squared = x$r_squared_train,
            RMSE = x$rmse_train,
            MAE = x$mae_train,
            RPD = x$rpd_train
          )
        }))
        # Sortiere nach RPD
        summary_df <- summary_df[order(-summary_df$RPD), ]
      }
      
      # Speichere Zusammenfassung
      write.csv(summary_df, file.path(filter_dir, "plsr_summary.csv"), row.names = FALSE)
      
      cat("\n===== ZUSAMMENFASSUNG für", algo_name, "=====\n")
      print(summary_df)
    } else {
      cat("Keine erfolgreichen PLSR-Modelle für", algo_name, "erstellt.\n")
    }
  } # End of algo_dir loop
  
  cat("\nVerarbeitung aller Algorithmen abgeschlossen!\n")
} # End of function definition

# HINWEIS: Auto-Execution auskommentiert - Script wird nur über master_pipeline.R gesteuert
# # Ausführung der Hauptfunktion mit den gewünschten Parametern
# perform_filter_plsr_by_algorithm(
#   data_file = "data/final_before_correlation.csv",   # Passe an deine Datei an
#   use_validation_split = TRUE,                   # Aktiviere externe Validierung
#   validation_ratio = 0.3,                        # 70% Training, 30% Validierung
#   min_correlation = 0.3                          # Mindestkorrelation für Bandauswahl
# )

