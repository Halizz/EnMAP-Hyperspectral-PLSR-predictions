# ==================================================================================
# VERBESSERTE STRATIFIZIERUNG FÜR KONTINUIERLICHE BODENVARIABLEN
# ==================================================================================

#' Adaptive Stratifizierung für kontinuierliche Bodenvariablen
#' 
#' @param data Dataset
#' @param target_var Zielvariale (kontinuierlich)
#' @param min_strata_size Minimale Größe pro Stratum
#' @param max_strata Maximale Anzahl Strata
#' @param method Stratifizierungsmethode ("adaptive", "quantile", "kmeans", "range")
adaptive_stratification <- function(data, target_var, 
                                  min_strata_size = 5, 
                                  max_strata = 5,
                                  method = "adaptive") {
  
  cat("=== ADAPTIVE STRATIFICATION FOR", target_var, "===\n")
  
  y <- data[[target_var]]
  valid_data <- data[!is.na(y), ]
  y_valid <- y[!is.na(y)]
  n_valid <- length(y_valid)
  
  cat(sprintf("Valid samples: %d\n", n_valid))
  
  # 1. Bestimme optimale Anzahl Strata basierend auf Stichprobengröße
  max_possible_strata <- floor(n_valid / min_strata_size)
  optimal_strata <- min(max_strata, max_possible_strata, 
                       max(2, floor(sqrt(n_valid/10))))  # Heuristik
  
  cat(sprintf("Optimal strata determined: %d (max possible: %d)\n", 
              optimal_strata, max_possible_strata))
  
  if (optimal_strata < 2) {
    cat("⚠️ Insufficient data for stratification, using random split\n")
    return(list(
      strata = rep("all", n_valid),
      method = "random",
      n_strata = 1,
      balance_quality = NA
    ))
  }
  
  # 2. Wähle Stratifizierungsmethode basierend auf Datenverteilung
  strata_result <- switch(method,
    "adaptive" = adaptive_strata_selection(y_valid, optimal_strata),
    "quantile" = quantile_stratification(y_valid, optimal_strata),
    "kmeans" = kmeans_stratification(y_valid, optimal_strata),
    "range" = range_stratification(y_valid, optimal_strata),
    stop("Unknown stratification method:", method)
  )
  
  # 3. Qualitätskontrolle der Stratifizierung
  strata_table <- table(strata_result$strata)
  min_size <- min(strata_table)
  balance_quality <- min_size / max(strata_table)  # Ausgewogenheits-Score
  
  cat("Strata sizes:", paste(as.numeric(strata_table), collapse = ", "), "\n")
  cat(sprintf("Balance quality: %.3f (1.0 = perfect balance)\n", balance_quality))
  
  # 4. Fallback bei unausgewogenen Strata
  if (min_size < min_strata_size || balance_quality < 0.3) {
    cat("⚠️ Poor stratification quality, trying simplified approach\n")
    # Fallback: Einfache Binär-Stratifizierung
    median_val <- median(y_valid)
    strata_result$strata <- ifelse(y_valid <= median_val, "low", "high")
    strata_result$method <- "binary_median"
    strata_result$n_strata <- 2
  }
  
  return(strata_result)
}

#' Adaptive Strata-Selektion basierend auf Datenverteilung
adaptive_strata_selection <- function(y, n_strata) {
  
  # Analysiere Datenverteilung
  skewness <- e1071::skewness(y)
  kurtosis <- e1071::kurtosis(y)
  
  cat(sprintf("Distribution analysis - Skewness: %.3f, Kurtosis: %.3f\n", 
              skewness, kurtosis))
  
  if (abs(skewness) < 0.5) {
    # Normalverteilung: Quantil-basiert
    cat("Using quantile-based stratification (normal distribution)\n")
    return(quantile_stratification(y, n_strata))
    
  } else if (skewness > 1.5) {
    # Stark rechtsschiefe Verteilung: Log-transformation + Quantile
    cat("Using log-transformed quantile stratification (right-skewed)\n")
    y_log <- log(y - min(y) + 1)  # +1 für positive Werte
    return(quantile_stratification(y_log, n_strata, original_values = y))
    
  } else {
    # Moderate Schiefe: K-means clustering
    cat("Using k-means clustering (moderate skewness)\n")
    return(kmeans_stratification(y, n_strata))
  }
}

#' Quantil-basierte Stratifizierung
quantile_stratification <- function(y, n_strata, original_values = NULL) {
  if (is.null(original_values)) original_values <- y
  
  # Berechne gleichmäßige Quantile
  probs <- seq(0, 1, length.out = n_strata + 1)
  breaks <- quantile(y, probs = probs, na.rm = TRUE)
  
  # Handle duplicate breaks (kann bei wenigen unique values passieren)
  if (length(unique(breaks)) < length(breaks)) {
    cat("⚠️ Duplicate quantile breaks detected, adjusting...\n")
    breaks <- unique(breaks)
    n_strata <- length(breaks) - 1
  }
  
  strata <- cut(original_values, breaks = breaks, 
                labels = paste0("Q", 1:n_strata),
                include.lowest = TRUE)
  
  return(list(
    strata = as.character(strata),
    method = "quantile",
    n_strata = n_strata,
    breaks = breaks
  ))
}

#' K-means-basierte Stratifizierung
kmeans_stratification <- function(y, n_strata) {
  
  # K-means clustering auf 1D Daten
  set.seed(42)  # Für Reproduzierbarkeit
  kmeans_result <- kmeans(y, centers = n_strata, nstart = 20)
  
  # Sortiere Cluster nach Zentren
  center_order <- order(kmeans_result$centers)
  cluster_mapping <- setNames(paste0("K", 1:n_strata), center_order)
  
  strata <- cluster_mapping[as.character(kmeans_result$cluster)]
  
  return(list(
    strata = strata,
    method = "kmeans",
    n_strata = n_strata,
    centers = sort(kmeans_result$centers)
  ))
}

#' Range-basierte Stratifizierung (für extreme Wertebereiche)
range_stratification <- function(y, n_strata) {
  
  y_range <- range(y, na.rm = TRUE)
  interval_width <- (y_range[2] - y_range[1]) / n_strata
  
  breaks <- seq(y_range[1], y_range[2], by = interval_width)
  if (length(breaks) != (n_strata + 1)) {
    breaks[length(breaks)] <- y_range[2]  # Ensure exact upper bound
  }
  
  strata <- cut(y, breaks = breaks, 
                labels = paste0("R", 1:n_strata),
                include.lowest = TRUE)
  
  return(list(
    strata = as.character(strata),
    method = "range",
    n_strata = n_strata,
    breaks = breaks
  ))
}

#' Verbesserte stratifizierte Aufteilung mit Qualitätskontrolle
improved_stratified_split <- function(data, target_var, 
                                     train_ratio = 0.7,
                                     stratification_method = "adaptive",
                                     min_strata_size = 5) {
  
  cat("=== IMPROVED STRATIFIED TRAIN/TEST SPLIT ===\n")
  
  # 1. Adaptive Stratifizierung
  strata_info <- adaptive_stratification(
    data, target_var, 
    min_strata_size = min_strata_size,
    method = stratification_method
  )
  
  if (strata_info$method == "random") {
    # Fallback: Standard random split
    set.seed(42)
    train_indices <- sample(nrow(data), size = round(nrow(data) * train_ratio))
    test_indices <- setdiff(1:nrow(data), train_indices)
    
    return(list(
      train_indices = train_indices,
      test_indices = test_indices,
      strata_info = strata_info,
      split_quality = list(balance_score = NA, representation_score = NA)
    ))
  }
  
  # 2. Stratifizierte Aufteilung innerhalb jeder Schicht
  valid_data <- data[!is.na(data[[target_var]]), ]
  strata <- strata_info$strata
  
  train_indices <- c()
  test_indices <- c()
  
  for (stratum in unique(strata)) {
    stratum_indices <- which(strata == stratum)
    n_stratum <- length(stratum_indices)
    
    if (n_stratum < 2) {
      # Zu kleine Schicht: alle zu Training
      train_indices <- c(train_indices, stratum_indices)
      next
    }
    
    # Stratifizierte Aufteilung innerhalb der Schicht
    n_train <- max(1, round(n_stratum * train_ratio))
    set.seed(42)
    stratum_train <- sample(stratum_indices, size = n_train)
    stratum_test <- setdiff(stratum_indices, stratum_train)
    
    train_indices <- c(train_indices, stratum_train)
    test_indices <- c(test_indices, stratum_test)
  }
  
  # 3. Qualitätsbewertung der Aufteilung
  split_quality <- evaluate_split_quality(
    data, target_var, train_indices, test_indices
  )
  
  cat(sprintf("Split quality - Balance: %.3f, Representation: %.3f\n",
              split_quality$balance_score, split_quality$representation_score))
  
  return(list(
    train_indices = train_indices,
    test_indices = test_indices,
    strata_info = strata_info,
    split_quality = split_quality
  ))
}

#' Bewertung der Aufteilungsqualität
evaluate_split_quality <- function(data, target_var, train_indices, test_indices) {
  
  train_values <- data[[target_var]][train_indices]
  test_values <- data[[target_var]][test_indices]
  
  # Balance Score: Ähnlichkeit der Mittelwerte
  train_mean <- mean(train_values, na.rm = TRUE)
  test_mean <- mean(test_values, na.rm = TRUE)
  overall_mean <- mean(data[[target_var]], na.rm = TRUE)
  
  balance_score <- 1 - abs(train_mean - test_mean) / abs(overall_mean)
  
  # Representation Score: Ähnlichkeit der Verteilungen (Kolmogorov-Smirnov)
  ks_test <- ks.test(train_values, test_values)
  representation_score <- 1 - ks_test$statistic  # Higher = more similar
  
  return(list(
    balance_score = max(0, balance_score),
    representation_score = representation_score,
    train_mean = train_mean,
    test_mean = test_mean,
    ks_statistic = ks_test$statistic,
    ks_p_value = ks_test$p.value
  ))
}

# ==================================================================================
# BEISPIEL-VERWENDUNG UND INTEGRATION
# ==================================================================================

#' Integration in bestehende Workflows
demo_improved_stratification <- function() {
  
  # Beispiel mit verschiedenen Bodenparametern
  data_file <- "data/final_before_correlation.csv"
  if (file.exists(data_file)) {
    data <- read.csv(data_file)
    
    # Teste verschiedene Parameter
    test_params <- c("BS", "pH", "CEC", "SWC")
    
    for (param in test_params) {
      if (param %in% colnames(data)) {
        cat("\n", rep("=", 60), "\n")
        cat("TESTING PARAMETER:", param, "\n")
        cat(rep("=", 60), "\n")
        
        # Teste verschiedene Methoden
        methods <- c("adaptive", "quantile", "kmeans")
        
        for (method in methods) {
          cat("\n--- Method:", method, "---\n")
          
          split_result <- improved_stratified_split(
            data, param, 
            stratification_method = method,
            min_strata_size = 3  # Reduziert für kleine Datensätze
          )
          
          cat("Final split sizes: Train =", length(split_result$train_indices),
              ", Test =", length(split_result$test_indices), "\n")
        }
      }
    }
  }
}

# Lade erforderliche Packages
if (!require(e1071)) install.packages("e1071")
library(e1071)

cat("Improved stratification functions loaded.\n")
cat("Usage: improved_stratified_split(data, 'target_variable')\n")
cat("Demo: demo_improved_stratification()\n")
