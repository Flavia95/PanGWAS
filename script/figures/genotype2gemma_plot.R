##Use:  Rscript genotype2gemma_plot.R -i DBmm10.5.nodup.rpn -n n_rpn -l n_rpnl -o comparison_1M.png -r 1000000
#!/usr/bin/env Rscript
library(ggplot2)
library(data.table)
library(gridExtra)
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", help="RPN input matrix"),
  make_option(c("-n", "--norm"), type="character", help="Simple euclidean normalized RPN"),
  make_option(c("-l", "--length"), type="character", help="Length + euclidean normalized RPN"),
  make_option(c("-o", "--output"), type="character", default="combined_plots.png", help="Output plot file [default=%default]"),
  make_option(c("-r", "--rows"), type="integer", default=-1, help="Number of rows to read (-1 for all) [default=%default]")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

create_plots <- function(input_file, norm_file, normleng_file, output_plot, nrows = -1) {
  read_matrix_values <- function(file) {
    if(nrows == -1) {
      dt <- fread(file, header = TRUE)
    } else {
      dt <- fread(file, header = TRUE, nrows = nrows)
    }
    as.numeric(unlist(dt[, -1]))
  }
  
  rpn_values <- read_matrix_values(input_file)
  euclidean_norm <- read_matrix_values(norm_file)
  length_euclidean_norm <- read_matrix_values(normleng_file)
  
  # Original with linear scale
  p1 <- ggplot(data.frame(value = rpn_values), aes(x = value)) +
    geom_histogram(bins = 50, alpha = 0.7, fill = "orange", color = "black") +
    scale_y_log10(labels = scales::comma, limits = c(1, NA)) +
    scale_x_continuous(labels = scales::comma) +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    ) +
    labs(x = "Value", y = "Frequency", title = "RPN")
  
  # Original with log scale
  p2 <- ggplot(data.frame(value = rpn_values), aes(x = value)) +
    geom_histogram(bins = 50, alpha = 0.7, fill = "orange", color = "black") +
    scale_y_log10(labels = scales::comma, limits = c(1, NA)) +
    scale_x_log10(labels = scales::comma) +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    ) +
    labs(x = "Value (log)", y = "Frequency", title = "RPN (log)")
  
  # Euclidean norm
  p3 <- ggplot(data.frame(value = euclidean_norm), aes(x = value)) +
    geom_histogram(bins = 50, alpha = 0.7, fill = "blue", color = "black") +
    scale_y_log10(labels = scales::comma, limits = c(1, NA)) +
    scale_x_continuous(labels = scales::comma) +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    ) +
    labs(x = "Value", y = "Frequency", title = "n_RPN\n(Euclidean)")
  
  # Length + Euclidean norm
  p4 <- ggplot(data.frame(value = length_euclidean_norm), aes(x = value)) +
    geom_histogram(bins = 50, alpha = 0.7, fill = "green", color = "black") +
    scale_y_log10(labels = scales::comma, limits = c(1, NA)) +
    scale_x_continuous(labels = scales::comma) +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95")
    ) +
    labs(x = "Value", y = "Frequency", title = "n_RPNL\n(Length + Euclidean)")

  # Boxplot comparison norm
  p5 <- ggplot(data.frame(
    value = c(euclidean_norm, length_euclidean_norm),
    type = factor(rep(c("n_rpn", "n_RPNL"),
                     each = length(euclidean_norm)))
  ), aes(x = type, y = value, fill = type)) +
    geom_boxplot(alpha = 0.7) +
    scale_fill_manual(values = c("blue", "green")) +
    scale_y_continuous(labels = scales::comma, limits = c(0, 1)) +
    theme_classic() +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_line(color = "gray90"),
      panel.grid.minor = element_line(color = "gray95"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(x = "", y = "Value", title = "Norm comparison")

  # Arrange plots
  combined <- grid.arrange(p1, p2, p3, p4, p5,
                          layout_matrix = matrix(1:5, nrow = 1),
                          widths = c(1, 1, 1, 1, 1))
  
  ggsave(output_plot, combined, width = 14, height = 2, dpi = 300, bg = "white")
}

create_plots(opt$input, opt$norm, opt$length, opt$output, opt$rows)
