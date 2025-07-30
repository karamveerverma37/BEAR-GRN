library(readr)
library(dplyr)
library(pROC)
library(PRROC)
library(ggplot2)
library(stringr)
library(tibble)
library(RColorBrewer)
setwd("/gpfs/Home/kmk7420/Multi_omics_GRN/Inferred_GRNs_all/precision_recall/K562")

# === Load ground truth ===
ground_truth <- read_tsv("filtered_RN117_K562.tsv", col_types = cols()) %>% select(Source, Target)
gt_pairs <- paste(ground_truth$Source, ground_truth$Target, sep = "_")
gt_pairs1 <- unique(gt_pairs)

# === Define method column mappings ===
method_info <- tribble(
  ~method,      ~tf_col,     ~target_col, ~score_col,
  "CellOracle", "source",    "target",    "coef_mean",
  "Scenic+",    "Source",    "Target",    "Score",
  "Pando",      "tf",        "target",    "estimate",
  "LINGER",     "Source",    "Target",    "Score",
  "FigR",       "Motif",     "DORC",      "Score",
  "TRIPOD",     "TF",        "gene",      "abs_coef",
  "GRaNIE",     "TF.name",   "gene.name", "TF_gene.r"
)
#set.seed(42)
input_dir <- "INPUT"
#results <- data.frame(method = character(), AUROC = numeric(), AUPRC = numeric(), AUPRC_random = numeric())
results <- data.frame(
  method = character(),
  AUROC = numeric(),
  AUPRC = numeric(),
  AUPRC_random = numeric(),
  stringsAsFactors = FALSE
)
roc_list <- list()
pr_list <- list()
random_pr_list <- list()

for (i in 1:nrow(method_info)) {
  set.seed(42 + i)
  method_name <- method_info$method[i]
  tf_col <- method_info$tf_col[i]
  target_col <- method_info$target_col[i]
  score_col <- method_info$score_col[i]

  method_path <- file.path(input_dir, method_name)
  file_list <- list.files(method_path, full.names = TRUE)
  file <- file_list[1]

  message("Processing: ", method_name)
  df <- if (str_ends(file, ".csv")) read_csv(file, col_types = cols()) else read_tsv(file, col_types = cols())

    inferred <- df %>%
    select(tf = all_of(tf_col),
           target = all_of(target_col),
           score = all_of(score_col)) %>%
    mutate(score = abs(score),
           tf = str_to_upper(tf),
           target = str_to_upper(target)) %>%
    group_by(tf, target) %>%
    summarise(score = max(score), .groups = "drop") %>%
    mutate(pair = paste(tf, target, sep = "_"),
           label = ifelse(pair %in% gt_pairs1, 1, 0))

  positives <- inferred %>% filter(label == 1)
  negatives <- inferred %>% filter(label == 0)
  if (nrow(positives) == 0) next

  #write.csv(inferred, paste0("labelled_dataframe_", method_name, ".csv"))

  # === AUROC calculation ===
  #set.seed(42)
  neg_sampled <- negatives %>% sample_n(min(nrow(positives), nrow(negatives)))
  auroc_data <- bind_rows(positives, neg_sampled)

  #roc_obj <- roc(auroc_data$label, auroc_data$score)
  roc_obj <- roc(response = auroc_data$label, predictor = auroc_data$score, levels = c(0,1), direction = "<")
  auroc_value <- as.numeric(auc(roc_obj))
  roc_df <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    method = method_name
  )
  roc_df <- roc_df[order(roc_df$fpr), ]
  #print(roc_df)
  roc_list[[method_name]] <- roc_df
  
  # === AUPRC ===
  #set.seed(42)
  neg_sampled_10x <- negatives %>% sample_n(size = min(nrow(positives) * 10, nrow(negatives)))

  pr_obj <- pr.curve(
    scores.class0 = positives$score,
    scores.class1 = neg_sampled_10x$score,
    curve = TRUE
  )
  auprc_value <- pr_obj$auc.integral
  pr_df <- data.frame(
    recall = pr_obj$curve[, 1],
    precision = pr_obj$curve[, 2],
    method = method_name
  )
  pr_list[[method_name]] <- pr_df

  # === Random AUPRC data by score shuffling ===
  #set.seed(42)
  inferred_random <- inferred %>% mutate(score = runif(n()))
  positives_r <- inferred_random %>% filter(label == 1)
  negatives_r <- inferred_random %>% filter(label == 0)
  neg_sampled_10x_r <- negatives_r %>% sample_n(size = min(nrow(positives_r) * 10, nrow(negatives_r)))

  pr_obj_random <- pr.curve(
    scores.class0 = positives_r$score,
    scores.class1 = neg_sampled_10x_r$score,
    curve = TRUE
  )
  auprc_random <- pr_obj_random$auc.integral
  pr_df_random <- data.frame(
    recall = pr_obj_random$curve[, 1],
    precision = pr_obj_random$curve[, 2],
    method = method_name
  )
  random_pr_list[[method_name]] <- pr_df_random

  results <- rbind(results, data.frame(method = method_name,
                                       AUROC = auroc_value,
                                       AUPRC = auprc_value,
                                       AUPRC_random = auprc_random))
}
#results
# === Create output directory ===
dir.create("gap_plots_RN117_K562", showWarnings = FALSE)

# === Custom Theme and Palette ===
custom_theme <- theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )
color_palette <- scale_color_brewer(palette = "Dark2")

# === Bar Plots ===
bar_auroc <- ggplot(results, aes(x = reorder(method, AUROC), y = AUROC, fill = method)) +
  geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Method", y = "AUROC", title = "AUROC Scores") +
  custom_theme

bar_auprc <- ggplot(results, aes(x = reorder(method, AUPRC), y = AUPRC, fill = method)) +
  geom_bar(stat = "identity", width = 0.7, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Method", y = "AUPRC", title = "AUPRC Scores") +
  custom_theme
#print(bar_auprc)
# === Label methods ===
results <- results %>%
  mutate(method_auroc = paste0(method, " (AUROC=", round(AUROC, 3), ")"),
         method_auprc = paste0(method, " (AUPRC=", round(AUPRC, 3), ")"),
         method_random_auprc = paste0(method, " (AUPRC_random=", round(AUPRC_random, 3), ")"))

#results
roc_combined <- bind_rows(roc_list, .id = "method") %>%
  left_join(results, by = "method")
#roc_combined
pr_combined <- bind_rows(pr_list, .id = "method") %>%
  left_join(results, by = "method")

random_pr_combined <- bind_rows(random_pr_list, .id = "method") %>%
  left_join(results, by = "method")
#head(random_pr_combined)
########################################
#print(roc_plot)
# === PR Plot data with Random Baselines ===
# Add curve type and generate combined label using the correct column name `AUROC`
random_pr_combined <- random_pr_combined %>%
  mutate(curve_type = "Randomized",
         legend_label = paste0(method, " ", curve_type, " (AUPRC = ", sprintf("%.3f", AUPRC_random), ")"))

pr_combined <- pr_combined %>%
  mutate(curve_type = "Normal",
         legend_label = paste0(method, " ", curve_type, " (AUPRC = ", sprintf("%.3f", AUPRC), ")"))

# Combine both
combined_pr <- bind_rows(random_pr_combined, pr_combined)

combined_pr <- combined_pr %>%
  mutate(
    method_curve = paste0(method, "_", curve_type),
    legend_label = case_when(
      curve_type == "Normal" ~ paste0(method, " AUPRC=", round(AUPRC, 3)),
      curve_type == "Randomized" ~ paste0(method, " Random AUPRC=", round(AUPRC_random, 3)),
      TRUE ~ method_curve
    )
  )

# Unique methods
unique_methods <- unique(combined_pr$method)
num_methods <- length(unique_methods)

# Assign colors per method
palette_colors <- brewer.pal(min(max(num_methods, 3), 8), "Dark2")
method_colors <- setNames(palette_colors[seq_len(num_methods)], unique_methods)

# Create method_curve list: for each method Normal then Randomized
ordered_method_curves <- unlist(lapply(unique_methods, function(m) {
  c(paste0(m, "_Normal"), paste0(m, "_Randomized"))
}))

# Assign colors and linetypes for these ordered method_curves
method_curve_colors <- method_colors[sub("_(Normal|Randomized)$", "", ordered_method_curves)]
names(method_curve_colors) <- ordered_method_curves

method_curve_linetypes <- ifelse(grepl("_Normal$", ordered_method_curves), "solid", "dotted")
names(method_curve_linetypes) <- ordered_method_curves

# Extract legend labels in this order
unique_labels_df <- combined_pr %>%
  distinct(method_curve, legend_label) %>%
  filter(method_curve %in% ordered_method_curves) %>%
  arrange(factor(method_curve, levels = ordered_method_curves))

unique_legend_labels <- unique_labels_df$legend_label

# Set factor levels on method_curve to fix order in plot and legend
combined_pr$method_curve <- factor(combined_pr$method_curve, levels = ordered_method_curves)
#max(combined_pr$precision)
# Plot
pr_plot <- ggplot(combined_pr, aes(
  x = recall, y = precision,
  color = method_curve,
  linetype = method_curve,
  group = method_curve
)) +
  geom_line(linewidth = 0.8) +
  labs(
    x = "Recall",
    y = "Precision",
    title = "PR Curve",
    color = "Method (AUPRC)",
    linetype = "Method (AUPRC)"
  ) +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_color_manual(values = method_curve_colors, labels = unique_legend_labels) +
  scale_linetype_manual(values = method_curve_linetypes, labels = unique_legend_labels) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 1.5, linetype = method_curve_linetypes)),
    linetype = "none" # hide duplicated linetype legend
  ) +
  custom_theme
#head(combined_pr)
###########################################
library(dplyr)
library(plotrix)
pdf("gap_plots_RN117_K562/PR_plot_RN117_K562_with_gap.pdf", width = 10, height = 8)
#setEPS()
#postscript("plots/PR_plot_RN117_and_RN118_K562_with_gap.ps",
#           width = 10, height = 8, horizontal = FALSE, onefile = FALSE, paper = "special")
#bitmap("plots/PR_plot_bitmap.eps", type = "eps2write",
#       res = 300, width = 10, height = 8, units = "in")
# Save PNG using ragg
#ragg::agg_png("plots/temp_gap.png", width = 10, height = 8, units = "in", res = 300)
#cairo_ps("PR_curve_cairo.eps", width = 7, height = 6, fallback_resolution = 300)
Cairo::CairoPS(file = "gap_plots_RN117_K562/PR_plot_RN117_K562_with_gap.eps",
               width = 10, height = 8)
# Prepare unique method-curve combinations
methods_to_plot <- unique(combined_pr$method_curve)

# Set up color and line type mappings
colors <- method_curve_colors[methods_to_plot]
#linetypes <- method_curve_linetypes[methods_to_plot]
#unique_curves <- unique(combined_pr$method_curve)
#colors
#####################################
# Extract unique labels, colors, and line types per curve
legend_df <- unique(combined_pr[, c("method_curve", "legend_label", "curve_type")])
legend_df <- legend_df[legend_df$method_curve %in% names(colors), ]

# Set line types: dashed for randomized, solid otherwise
legend_df$lty <- ifelse(legend_df$curve_type == "Randomized", 3, 1)
legend_df$col <- colors[legend_df$method_curve]
# Step 1: Get unique base method names
base_methods <- unique(gsub("_(Normal|Randomized)", "", legend_df$method_curve))

# Step 2: Interleave Normal and Randomized for each method
ordered_df <- do.call(rbind, lapply(base_methods, function(method) {
  normal <- legend_df[legend_df$method_curve == paste0(method, "_Normal"), ]
  random <- legend_df[legend_df$method_curve == paste0(method, "_Randomized"), ]
  rbind(normal, random)
}))
par(mar = c(5, 4, 4, 8))  

ordered_df$legend_label
ordered_df$col <-c("#1B9E77","#1B9E77", "#D95F02","#D95F02", "#7570B3","#7570B3", "#E7298A","#E7298A", "#66A61E","#66A61E", "#E6AB02", "#E6AB02", "#A6761D", "#A6761D")

# Set gap region on the y-axis (Precision)
gap_region <- c(0.4, 0.9)
# Set up base plot settings
par(mfrow = c(1, 1))
#par(mar = c(3, 3, 2, 1))
#par(mar = c(5, 4, 4, 2) + 0.1)
par(tck = -0.02)
par(cex.axis = 0.8)
par(mgp = c(2, 1, 0))
# Plot the first curve with gap.plot
first_label <- ordered_df$legend_label[1]
first_color <- ordered_df$col[1]
first_lty <- if (grepl("Random", first_label)) 2 else 1
df_first <- subset(combined_pr, legend_label == first_label)

#df_first <- df_first[is.finite(df_first$precision) & df_first$precision <= 1, ]
print(max(pr_obj$curve[, 2]))
gap.plot(x = df_first$recall, y = df_first$precision,
         gap = gap_region, gap.axis = "y",
         type = "l", col = first_color, lwd = 2, lty = first_lty,
         ylim = c(0, 1), xlim = c(0, 1),
         xlab = "Recall", ylab = "Precision",main = "PR Curve",
         xaxt = "n", yaxt = "n",
         ytics = c(0, 0.1, 0.2, 0.3, 1.02),        # Only desired tick positions
         yticlab = c("0", "0.1", "0.2", "0.3", "1"))  # their labels)
#??gap.plot
# Plot the rest using gap.lines (not available), so instead manually clip and use lines()
for (i in 2:nrow(ordered_df)) {
  df_i <- subset(combined_pr, legend_label == ordered_df$legend_label[i])
  #df_i <- df_i[is.finite(df_i$precision) & df_i$precision <= 1, ]
  this_lty <- if (grepl("Random", ordered_df$legend_label[i])) 2 else 1
  lines(df_i$recall, df_i$precision,
        col = ordered_df$col[i], lwd = 2, lty = this_lty)
}
#print(pr_plot)
# Add axes
#axis(1, at = seq(0, 1, 0.2), labels = TRUE, cex.axis = 0.8)
#axis(2, at = c(0, 0.1, 0.2, 0.3, 1), las = 1, cex.axis = 0.8)

# Add y-axis break marker
abline(h = seq(gap_region[1], gap_region[1] + 0.02, .00001), col = "white")
abline(h = 1, col = "red", lty = 2, lwd = 9)

axis.break(2, 0.405, style = "slash")
#axis.break(2, gap_mid, style = "slash")
# Add legend with correct color and line type
legend("topright",
       legend = ordered_df$legend_label,
       col = ordered_df$col,
       lwd = 2,
       lty = ifelse(grepl("Random", ordered_df$legend_label), 2, 1),
       cex = 1)

dev.off()
# Convert to EPS using ImageMagick
#system("convert plots/temp_gap.png plots/PR_plot_bitmap.eps")

######################
head(combined_pr)


## ROC plot
roc_combined$method <- factor(roc_combined$method, levels = unique_methods)
roc_plot <- ggplot(roc_combined, aes(x = fpr, y = tpr, color = method)) +
  geom_line(linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray") +
  labs(x = "False Positive Rate", y = "True Positive Rate", title = "ROC Curve", color = "Method (AUROC)") +
  coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1),expand = FALSE) +
  scale_color_manual(values = method_colors, labels = results$method_auroc) +
  custom_theme
# Save Plots 
#ggsave("plots/AUROC_bar.pdf", plot = bar_auroc, width = 6, height = 4,device = cairo_pdf)
#ggsave("plots/AUROC_bar.png", plot = bar_auroc, width = 6, height = 4, dpi = 300)

#ggsave("plots/AUPRC_bar.pdf", plot = bar_auprc, width = 6, height = 4,device = cairo_pdf)
#ggsave("plots/AUPRC_bar.png", plot = bar_auprc, width = 6, height = 4, dpi = 300)
library(svglite)

ggsave("gap_plots_RN117_K562/ROC_plot_RN117_K562_with_gap.eps", plot = roc_plot, width = 6, height = 5, device = "eps")
ggsave("gap_plots_RN117_K562/ROC_plot_RN117_K562_with_gap.png", plot = roc_plot, width = 6, height = 5, dpi = 300)

#ggsave("plots/PR_curve.svg", plot = pr_plot, width = 6, height = 5)
#ggsave("plots/PR_curve.png", plot = pr_plot, width = 6, height = 5, dpi = 300)
