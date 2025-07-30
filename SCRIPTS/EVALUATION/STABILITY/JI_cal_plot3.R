# Load necessary libraries
library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(tibble)

# Define directories and methods
sample_dirs <- c("K562", "Macrophage_S1","Macrophage_S2","E7.5_rep1", "E7.5_rep2","E8.5_rep1", "E8.5_rep2")
methods <- c("CellOracle", "FigR", "TRIPOD", "Pando", "LINGER", "Scenic", "GRaNIE","DIRECTNET")
set.seed(123)
# Define sorting column per method
sort_columns <- list(
  CellOracle = "coef_abs",
  FigR = "Score",
  TRIPOD = "abs_coef",
  Pando = "estimate",
  LINGER = "Score",
  Scenic = "Score",
  GRaNIE = "TF_gene.r"
)
# Define source and target columns per method
source_target_columns <- list(
  CellOracle = c("source", "target"),
  FigR = c("Motif", "DORC"),
  TRIPOD = c("TF", "gene"),
  Pando = c("tf", "target"),
  LINGER = c("Source", "Target"),
  Scenic = c("Source", "Target"),
  GRaNIE = c("TF.name", "gene.name"),
  DIRECTNET = c("TF","Gene")
)


# Define file delimiters for each method
delim <- list(
  CellOracle = ",",
  FigR = ",",
  TRIPOD = ",",
  Pando = ",",
  LINGER = ",",
  Scenic = "\t",
  GRaNIE = ",",
  DIRECTNET = ","
)

# Function to read and process a network
read_and_sort_network <- function(file, method) {
  network <- read_delim(file, delim = delim[[method]], show_col_types = FALSE)
  source_col <- source_target_columns[[method]][1]
  target_col <- source_target_columns[[method]][2]
  score_col <- sort_columns[[method]]  # may be NULL for DIRECTNET
  
  if (!(source_col %in% colnames(network)) | !(target_col %in% colnames(network))) {
    stop(paste("Missing expected columns in", file, "for method", method))
  }
  
  # Check if score_col is NULL or not present in columns
  if (is.null(score_col) || !(score_col %in% colnames(network))) {
    message(paste("No score column found or specified for", method, "- skipping sorting"))
    network_long <- network %>%
      select(all_of(source_col), all_of(target_col)) %>%
      rename(Source = all_of(source_col), Target = all_of(target_col))
    network_long$Score <- NA
    return(list(top = network_long, random = network_long %>% sample_n(size = round(0.1 * nrow(network_long)))))
  }
  
  network_long <- network %>%
    select(all_of(source_col), all_of(target_col), all_of(score_col)) %>%
    rename(Source = all_of(source_col), Target = all_of(target_col), Score = all_of(score_col))
  
  # Extract Top 10% and Random 10%
  top_10 <- network_long %>% arrange(desc(Score)) %>% slice_head(n = round(0.1 * nrow(network_long)))
  random_10 <- network_long %>% sample_n(size = round(0.1 * nrow(network_long)), replace = FALSE)
  
  return(list(top = top_10, random = random_10))
}



# Function to compute Jaccard Index
calculate_jaccard <- function(network1, network2) {
  common_edges <- intersect(paste(network1$Source, network1$Target),
                            paste(network2$Source, network2$Target))
  union_edges <- union(paste(network1$Source, network1$Target),
                       paste(network2$Source, network2$Target))
  if (length(union_edges) == 0) return(NA)
  return(length(common_edges) / length(union_edges))
}

# Initialize empty results dataframe
median_ji_results <- data.frame(Sample = character(), Method = character(), Category = character(), Median_JI = numeric(), stringsAsFactors = FALSE)

# Process each sample
for (sample in sample_dirs) {
  for (method in methods) {
    method_path <- file.path(sample, method)
    network_files <- list.files(method_path, pattern = "\\.(csv|tsv|txt)$", full.names = TRUE)

    if (length(network_files) < 2) next

    network_lists <- map(network_files, ~ read_and_sort_network(.x, method))

    # Separate top 10% and random 10% lists
    top_networks <- map(network_lists, "top")
    random_networks <- map(network_lists, "random")

    # Compute Jaccard Index for top 10%
    ji_top <- c()
    for (i in 1:(length(top_networks) - 1)) {
      for (j in (i + 1):length(top_networks)) {
        ji_top <- c(ji_top, calculate_jaccard(top_networks[[i]], top_networks[[j]]))
      }
    }

    # Compute Jaccard Index for random 10%
    ji_random <- c()
    for (i in 1:(length(random_networks) - 1)) {
      for (j in (i + 1):length(random_networks)) {
        ji_random <- c(ji_random, calculate_jaccard(random_networks[[i]], random_networks[[j]]))
      }
    }

    # Append results to the dataframe (with separate sample names for top/random)
    median_ji_results <- rbind(median_ji_results,
                               data.frame(Sample = paste(sample, "(Top 10%)"), Method = method, Category = "Top 10%", Median_JI = median(ji_top, na.rm = TRUE)),
                               data.frame(Sample = paste(sample, "(Random 10%)"), Method = method, Category = "Random 10%", Median_JI = median(ji_random, na.rm = TRUE)))
  }
}
# Explicit factor levels to preserve row/column layout
#desired_order <- c("K562 (Top 10%)", "K562 (Random 10%)", 
 #                  "Macrophage_S1 (Top 10%)", "Macrophage_S1 (Random 10%)", 
  #                 "Macrophage_S2 (Top 10%)", "Macrophage_S2 (Random 10%)", 
   #                "E7.5_rep1 (Top 10%)", "E7.5_rep1 (Random 10%)", 
    #               "E7.5_rep2 (Top 10%)", "E7.5_rep2 (Random 10%)", 
     #              "E8.5_rep1 (Top 10%)", "E8.5_rep1 (Random 10%)", 
      #             "E8.5_rep2 (Top 10%)", "E8.5_rep2 (Random 10%)")

#median_ji_results$Sample <- factor(median_ji_results$Sample, levels = desired_order)

median_ji_results$Sample <- factor(median_ji_results$Sample, levels = rev(unique(median_ji_results$Sample)))
median_ji_results$Method <- factor(median_ji_results$Method, levels = unique(median_ji_results$Method))
library(viridis)

# Create a custom ggplot2 theme
custom_theme <- theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    plot.margin = margin(15, 15, 15, 25, "pt")
  )
# Proper heatmap layout with fixed aspect ratio
ji_heatmap <- ggplot(median_ji_results, aes(x = Method, y = Sample, fill = Median_JI)) +
  geom_tile(color = "grey80", size = 0.2) +
  #geom_text(aes(label = sprintf("%.2f", Median_JI)), size = 4.2, fontface = "bold") +
  geom_text(
  aes(label = sprintf("%.2f", Median_JI), 
      color = ifelse(Median_JI > 0.5, "white", "black")),  # threshold at 0.5
  size = 4.2, fontface = "bold", show.legend = FALSE
  ) + scale_color_identity() +
  scale_fill_viridis(name = "Median JI", option = "D", direction = -1, na.value = "white") +
  coord_fixed(ratio = 0.5) +  # fixes tile to square shape
  labs(
    title = "Median Jaccard Index Heatmap",
    x = "Method",
    y = "Sample"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    panel.grid = element_blank()
  )# Load viridis for a colorblind-friendly palette


# Display the plot
print(ji_heatmap)

# Save plot
ggsave("median_ji_heatmap.png", plot = ji_heatmap, width = 12, height = 7, dpi = 300)

# Save results
write.csv(median_ji_results, "median_ji_results.csv", row.names = FALSE)

