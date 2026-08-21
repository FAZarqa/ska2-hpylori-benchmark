# Compare SKA2 distances against a benchmark distance matrix (DNAdiff SNP
# distance or ML patristic distance) and plot as a scatter with Pearson r,
# R-squared and n annotated on the panel.
#
# This is meant to be sourced, then called with whichever files/labels
# you need rather than edited per comparison.

library(ggplot2)
library(dplyr)

# strip everything after the accession/version so sample names match
# across SKA2 output, dnadiff output and tree tip labels
shorten_name <- function(x) sub("^(GCA_[0-9]+\\.[0-9]+).*", "\\1", x)

# read an SKA2 distance file and add columns to identify each pair
load_ska2 <- function(filepath) {
  ska <- read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  ska$sample1_id <- shorten_name(ska$Sample1)
  ska$sample2_id <- shorten_name(ska$Sample2)
  ska$genome1 <- pmin(ska$sample1_id, ska$sample2_id)
  ska$genome2 <- pmax(ska$sample1_id, ska$sample2_id)
  ska$pair_id <- paste(ska$genome1, ska$genome2, sep = "_vs_")
  ska
}

# scatter of SKA2 distance vs a benchmark distance, same panel style
# used throughout (black points, red lm line, dashed 1:1 line)
plot_ska2_vs_benchmark <- function(ska_file, benchmark_file, benchmark_col,
                                    k_label, benchmark_label, out_dir) {

  ska <- load_ska2(ska_file)
  valid_genomes <- unique(c(ska$sample1_id, ska$sample2_id))

  benchmark <- read.table(benchmark_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  benchmark <- benchmark %>% filter(genome1 %in% valid_genomes & genome2 %in% valid_genomes)
  benchmark <- benchmark[!duplicated(benchmark$pair_id), ]

  combined <- merge(
    benchmark[, c("pair_id", benchmark_col)],
    ska[, c("pair_id", "Distance")],
    by = "pair_id"
  )
  names(combined)[2:3] <- c("benchmark_distance", "ska_distance")
  combined <- na.omit(combined)

  cor_pearson <- cor(combined$benchmark_distance, combined$ska_distance)
  r_squared <- summary(lm(ska_distance ~ benchmark_distance, data = combined))$r.squared

  cat("\n===", k_label, "vs", benchmark_label, "===\n")
  cat("n pairs:", nrow(combined), "\n")
  cat("Pearson r:", round(cor_pearson, 4), "\n")
  cat("R-squared:", round(r_squared, 4), "\n")

  axis_max <- max(c(combined$benchmark_distance, combined$ska_distance)) * 1.05

  p <- ggplot(combined, aes(benchmark_distance, ska_distance)) +
    geom_point(alpha = 0.4, size = 1.5, color = "black") +
    geom_smooth(method = "lm", color = "#d62728", se = TRUE, linewidth = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
    scale_x_continuous(limits = c(0, axis_max), expand = c(0.02, 0)) +
    scale_y_continuous(limits = c(0, axis_max), expand = c(0.02, 0)) +
    theme_classic(base_size = 14) +
    theme(
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_line(color = "gray95", linewidth = 0.1),
      aspect.ratio = 1,
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    ) +
    labs(
      title = paste0("SKA2 (", k_label, ") vs ", benchmark_label),
      x = paste0(benchmark_label, " distance"),
      y = paste0("SKA2 distance (", k_label, ")")
    ) +
    annotate("text", x = axis_max * 0.95, y = axis_max * 0.05,
             hjust = 1, vjust = 0, size = 5, fontface = "bold",
             label = paste0("r = ", round(cor_pearson, 3),
                            "\nR\u00B2 = ", round(r_squared, 3),
                            "\nn = ", format(nrow(combined), big.mark = ",")))

  out_file <- file.path(out_dir, paste0("ska_", tolower(k_label), "_vs_",
                                         tolower(gsub(" ", "_", benchmark_label)), "_scatter.png"))
  ggsave(out_file, p, width = 8, height = 8, dpi = 300)
  cat("Saved:", out_file, "\n")

  list(plot = p, data = combined, r = cor_pearson, r_squared = r_squared)
}

# example calls, edit paths and rerun as needed
# result_k31 <- plot_ska2_vs_benchmark(
#   ska_file = "SKA distances/K31/hpylori_distance_148_k31.txt",
#   benchmark_file = "SKA distances/dnadiff_cleaned.txt",
#   benchmark_col = "SNP_distance",
#   k_label = "k31",
#   benchmark_label = "DNAdiff",
#   out_dir = "complete genomes"
# )
