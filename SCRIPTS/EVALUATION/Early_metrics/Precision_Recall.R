# Load required libraries
library(tidyverse)
library(readr)
library(ggplot2)
library(scales)
library(fs)

# Dataset names
datasets <- c("Macrophage_S1", "Macrophage_S2", "K562","E7.5_rep1","E7.5_rep2","E8.5_rep1","E8.5_rep2")

# Method-specific column mappings
column_mappings <- list(
  "CellOracle" = c("source", "target"),
  "FigR" = c("Motif", "DORC"),
  "GRaNIE" = c("TF.name", "gene.name"),
  "LINGER" = c("Source", "Target"),
  "Pando" = c("tf", "target"),
  "SCENIC+" = c("Source", "Target"),
  "TRIPOD" = c("TF", "gene"),
  "DIRECTNET" = c("TF", "Gene")
)

score_mappings <- list(
  "CellOracle" = "coef_abs",
  "FigR" = "Score",
  "GRaNIE" = "TF_gene.r",
  "LINGER" = "Score",
  "Pando" = "estimate",
  "SCENIC+" = "Score",
  "TRIPOD" = "abs_coef"
)

# Color palette
#method_colors <- c(
 # "CELL_ORACLE" = "#1B9E77",
  #"FIGR" = "#D95F02",
  #GRANIE" = "#7570B3",
  #"LINGER" = "#E7298A",
  #"PANDO" = "#66A61E",
  #"SCENIC+" = "#E6AB02",
  #"TRIPOD" = "#A6761D",
  #"DIRECT-NET" = "#666666"
#)
method_colors <- c(
  "CellOracle" = "#1B9E77",  # Teal green
  "FigR"        = "#66A61E",  # Green
  "GRaNIE"      = "#A6761D",  # Brown
  "LINGER"      = "#E7298A",  # Pink/magenta
  "Pando"       = "#7570B3",  # Purple
  "SCENIC+"     = "#D95F02",  # Orange
  "TRIPOD"      = "#E6AB02",  # Mustard yellow
  "DIRECTNET"  = "#666666"   # Grey
)
# Create output folder if it doesn't exist
dir_create("method_wise_metrics_10kEdges")
outdir <- "method_wise_metrics_10kEdges"
# Initialize plot data
plot_data <- tibble(Metric = character(), Score = numeric(), Dataset = character(), Method = character())

# Main loop
for (dataset in datasets) {
  dataset_dir <- dataset
  gt_dir <- file.path(dataset_dir, "GROUND_TRUTH")
  gt_file <- dir_ls(gt_dir, regexp = "RN.*\\.tsv$", recurse = FALSE)[1]

  if (is.na(gt_file)) next

  gt <- tryCatch(read_tsv(gt_file, col_types = cols()), error = function(e) NULL)
  if (is.null(gt)) next

  gt_edges <- unique(tibble(
    source = toupper(gt$Source),
    target = toupper(gt$Target)
  ))
  gt_pairs <- with(gt_edges, paste(source, target, sep = "_"))
  gt_tfs <- unique(gt_edges$source)
  gt_tgs <- unique(gt_edges$target)

  methods <- dir_ls(dataset_dir, type = "directory") %>%
    basename() %>% setdiff("GROUND_TRUTH")

  inferred_networks <- list()
  min_edges <- Inf

  for (method in methods) {
    method_dir <- file.path(dataset_dir, method)
    files <- dir_ls(method_dir, regexp = "\\.(csv|tsv)$", recurse = FALSE)

    if (length(files) == 0 || !(method %in% names(column_mappings))) next

    file <- files[1]
    sep <- if (grepl("\\.csv$", file)) "," else "\t"
    cols <- column_mappings[[method]]
    score_col <- score_mappings[[method]]

    try({
      df <- read_delim(file, delim = sep, col_types = cols(), show_col_types = FALSE)
      df <- df %>% rename(source = all_of(cols[1]), target = all_of(cols[2]))
      if (!is.null(score_col) && score_col %in% names(df)) {
        df[[score_col]] <- abs(df[[score_col]])
        df <- df %>% arrange(desc(.data[[score_col]]))
      }
      inferred_networks[[method]] <- df
      min_edges <- min(min_edges, nrow(df))
    }, silent = TRUE)
  }

  if (is.infinite(min_edges)) next

  summary_list <- list()

  for (method in methods) {
    if (!method %in% names(inferred_networks)) next

    df <- inferred_networks[[method]] %>% head(10000)
    df_edges <- df %>%
      mutate(source = toupper(source), target = toupper(target)) %>%
      transmute(pair = paste(source, target, sep = "_"))

    tp <- intersect(df_edges$pair, gt_pairs)
    fp <- setdiff(df_edges$pair, gt_pairs)
    fn <- setdiff(gt_pairs, df_edges$pair)

    precision <- ifelse((length(tp) + length(fp)) > 0, length(tp) / (length(tp) + length(fp)), 0)
    recall <- ifelse((length(tp) + length(fn)) > 0, length(tp) / (length(tp) + length(fn)), 0)
    f1 <- ifelse((precision + recall) > 0, 2 * precision * recall / (precision + recall), 0)

    plot_data <- bind_rows(plot_data, tibble(
      Metric = c("Precision", "Recall", "F1-score"),
      Score = c(precision, recall, f1),
      Dataset = dataset,
      Method = method
    ))

    inferred_tfs <- unique(toupper(df$source))
    inferred_tgs <- unique(toupper(df$target))

    summary_list[[method]] <- tibble(
      Method = method,
      GT_total = length(gt_pairs),
      Inferred_total = nrow(df),
      Inferred_TFs = length(inferred_tfs),
      Inferred_TGs = length(inferred_tgs),
      TP = length(tp),
      FP = length(fp),
      FN = length(fn),
      TN = NA,
      Precision = round(precision, 4),
      Recall = round(recall, 4),
      F1_score = round(f1, 4),
      Common_TFs = length(intersect(gt_tfs, inferred_tfs)),
      Common_Targets = length(intersect(gt_tgs, inferred_tgs)),
      Inferred_TFs_in_GT_TGs = length(intersect(inferred_tfs, gt_tgs)),
      Inferred_TGs_in_GT_TFs = length(intersect(inferred_tgs, gt_tfs))
    )
  }

  # Save summary
  if (length(summary_list) > 0) {
    df_summary <- bind_rows(summary_list)
    write_tsv(df_summary, file.path(outdir, paste0(dataset, "_method_wise_metrics.tsv")))
  }
}

# Set global method factor levels based on a fixed order (or a chosen metric)
#global_method_order <- plot_data %>%
#  filter(Metric == "F1-score", Dataset == "K562") %>%
#  arrange(desc(Score)) %>%
#  pull(Method) %>%
#  unique()

# Apply consistent factor levels
#plot_data$Method <- factor(plot_data$Method, levels = global_method_order)
# Plot results
library(extrafont)
#font_import(pattern = "Arial", prompt = FALSE)
#loadfonts(device = "all")  # For EPS
loadfonts(device = "pdf")         # For PDF
# List available fonts
#fonts()
#DejaVu Sans Mono
# Clean and set dataset order globally before plotting
plot_data$Dataset <- str_trim(plot_data$Dataset)
plot_data$Dataset <- factor(plot_data$Dataset, levels = datasets)

for (metric in c("Precision", "Recall", "F1-score")) {
  df_metric <- plot_data %>% filter(Metric == metric)
  #plot_data$Dataset <- str_trim(plot_data$Dataset) 
  #plot_data$Dataset <- factor(plot_data$Dataset, levels = datasets)
  # Sort by K562 metric
  method_order <- df_metric %>%
    filter(Dataset == "Macrophage_S1") %>%
    arrange(desc(Score)) %>%
    pull(Method) %>%
    unique()

  df_metric$Method <- factor(df_metric$Method, levels = method_order)

  p <- ggplot(df_metric, aes(x = Dataset, y = Score, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(expand = c(0, 0)) +  # Removes space below bars
  labs(
    title = paste0(metric, " Comparison Across Methods (Sorted by Macrophage_S1 ", metric, ")"),
    x = "Dataset", y = metric
  ) +
  theme_minimal(base_size = 14, base_family = "DejaVu Sans Mono") +
  theme(
    panel.grid = element_blank(),  # Removes all grid lines
    legend.position = "right",
    legend.title = element_text(size = 12),
    axis.line = element_line(color = "black", size = 0.6),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


  #ggsave(paste0(outdir,"/",metric, "_comparison_sorted_by_K562_", metric, ".png"),
   #      plot = p, width = 10, height = 6)
  #ggsave(paste0(outdir,"/",metric, "_comparison_sorted_by_K562_", metric, ".eps"),
   #      plot = p, width = 10, height = 6, device = "eps")
  ggsave(paste0(outdir,"/",metric, "_comparison_sorted_by_Macrophage_S1_", metric, ".pdf"),
         plot = p, width = 10, height = 6, device = "pdf")
}
# Define fixed order of methods
method_order <- c("CellOracle", "FigR", "GRaNIE", "LINGER", "Pando", "SCENIC+", "TRIPOD", "DIRECTNET")

for (metric in c("Precision", "Recall", "F1-score")) {
  df_metric <- plot_data %>% filter(Metric == metric)
  df_metric$Dataset <- factor(df_metric$Dataset, levels = datasets)
  df_metric$Method <- factor(df_metric$Method, levels = method_order)

  # Compute dynamic y-axis limit
  max_score <- max(df_metric$Score, na.rm = TRUE)
  y_limit <- min(1, round(max_score * 1.1, 2))

  p_jitter <- ggplot(df_metric, aes(x = Dataset, y = Score, color = Method)) +
    geom_jitter(
      position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),  # dodging with jitter
      size = 3,
      shape = 16,
      alpha = 0.85
    ) +
    scale_color_manual(values = method_colors) +
    scale_y_continuous(
      limits = c(0, y_limit),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_x_discrete(expand = expansion(add = c(0.9, 0.9))) +
    labs(
      title = paste0(metric, " Jitter Plot Across Methods"),
      x = "Dataset", y = metric
    ) +
    theme_minimal(base_size = 14, base_family = "DejaVu Sans Mono") +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 12),
      axis.line = element_line(color = "black", size = 0.6),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(
    filename = file.path(outdir, paste0(metric, "_jitter_plot.pdf")),
    plot = p_jitter,
    width = 10,
    height = 6,
    device = "pdf"
  )
}

##################lolipop plot#######################
for (metric in c("Precision", "Recall", "F1-score")) { 
  df_metric <- plot_data %>% filter(Metric == metric)
  df_metric$Dataset <- factor(df_metric$Dataset, levels = datasets)

  # Compute dynamic y-axis limit with 10% buffer, capped at 1
  max_score <- max(df_metric$Score, na.rm = TRUE)
  y_limit <- min(1, round(max_score * 1.1, 2))

  p_lollipop <- ggplot(df_metric, aes(x = Dataset, y = Score)) +
    geom_segment(
      aes(x = Dataset, xend = Dataset, y = 0, yend = Score),
      color = "white",  # White lines for lollipop sticks
      position = position_dodge(width = 0.3),
      size = 0.6, alpha = 0.6
    ) +
    geom_point(
      aes(color = Method),
      shape = 16,
      position = position_dodge(width = 0.3),
      size = 3, alpha = 0.9
    ) +
    scale_color_manual(values = method_colors) +
    scale_y_continuous(
      limits = c(0, y_limit),
      expand = expansion(mult = c(0, 0.02))
    ) + scale_x_discrete(expand = expansion(add = c(1.5, 1.5))) +
    labs(
      title = paste0(metric, " Lollipop Plot Across Methods"),
      x = "Dataset", y = metric
    ) +
    theme_minimal(base_size = 14, base_family = "DejaVu Sans Mono") +
    theme(
      panel.grid = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 12),
      axis.line = element_line(color = "black", size = 0.6),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  ggsave(
    filename = file.path(outdir, paste0(metric, "_lollipop_plot.pdf")),
    plot = p_lollipop,
    width = 10,
    height = 6,
    device = "pdf"
  )
}
