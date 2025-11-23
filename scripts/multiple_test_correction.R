# Multiple Testkorrektur für Korrelationstabellen (Benjamini-Hochberg FDR) für Pearson und Spearman, nur R2 > 0.4
library(dplyr)
library(readr)
library(stringr)

algorithms <- c("Pearson", "Spearman")
results_dir <- "results"

for (alg in algorithms) {
  alg_dir <- file.path(results_dir, alg)
  files <- list.files(alg_dir, pattern = "_correlation_table.csv$", full.names = TRUE)
  for (f in files) {
    cor_table <- read_csv(f)
    # Prüfe auf Split-Spalten
    has_split <- all(c("Correlation_Split1", "Correlation_Split2") %in% names(cor_table))
    # Probenanzahl automatisch bestimmen, falls vorhanden, sonst Standardwert
    n <- if ("n" %in% names(cor_table)) cor_table$n[1] else 30

    if (has_split) {
      # Split1
      t1 <- cor_table$Correlation_Split1 * sqrt((n - 2) / (1 - cor_table$Correlation_Split1^2))
      p1 <- 2 * pt(-abs(t1), df = n - 2)
      adj_p1 <- p.adjust(p1, method = "BH")
      sig1 <- adj_p1 < 0.05
      # Split2
      t2 <- cor_table$Correlation_Split2 * sqrt((n - 2) / (1 - cor_table$Correlation_Split2^2))
      p2 <- 2 * pt(-abs(t2), df = n - 2)
      adj_p2 <- p.adjust(p2, method = "BH")
      sig2 <- adj_p2 < 0.05
      # Kombiniere
      cor_table <- cor_table %>%
        mutate(
          Split1_p_value = p1,
          Split1_adj_p_value = adj_p1,
          Split1_significant = sig1,
          Split2_p_value = p2,
          Split2_adj_p_value = adj_p2,
          Split2_significant = sig2,
          Both_significant = Split1_significant & Split2_significant,
          Correlation_corrected = ifelse(Both_significant, (Correlation_Split1 + Correlation_Split2)/2, NA_real_)
        )
    } else if ("Correlation" %in% names(cor_table)) {
      # Kein Split: wie bisher
      t_stat <- cor_table$Correlation * sqrt((n - 2) / (1 - cor_table$Correlation^2))
      p_value <- 2 * pt(-abs(t_stat), df = n - 2)
      adj_p_value <- p.adjust(p_value, method = "BH")
      cor_table <- cor_table %>%
        mutate(
          t_stat = t_stat,
          p_value = p_value,
          adj_p_value = adj_p_value,
          significant = adj_p_value < 0.05,
          Correlation_corrected = ifelse(significant, Correlation, NA_real_)
        )
    } else {
      next
    }

    out_file <- str_replace(f, "_correlation_table.csv$", "_correlation_table_corrected.csv")
    write_csv(cor_table, out_file)
    cat(sprintf("[%s] Korrektur für Datei %s abgeschlossen (%d Zeilen).\n", alg, basename(f), nrow(cor_table)))
  }
}
