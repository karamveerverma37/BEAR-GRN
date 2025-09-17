# GRN Inference Method Benchmark Dataset Collection

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17087863.svg)](https://doi.org/10.5281/zenodo.17087863)

A comprehensive collection of Gene Regulatory Network (GRN) inference results from multiple state-of-the-art methods across various single-cell multiomics datasets, designed for benchmarking and comparative analysis of new GRN inference methods.

## 📊 Dataset Overview

This collection contains GRN inference results from **8 established methods** across **7 single-cell datasets**, along with corresponding ground truth regulatory networks for rigorous benchmarking. The dataset includes multiple analysis scripts for comprehensive evaluation of GRN methods.

### Included Methods
- **CellOracle**: Dynamic GRN inference from scRNA-seq
- **SCENIC+**: Multi-omics GRN inference (extension of Single-Cell rEgulatory Network Inference and Clustering)
- **Pando**: Lineage-informed GRN reconstruction
- **LINGER**: Lifelong learning for GRN inference
- **FigR**: Functional inference of gene regulation
- **TRIPOD**: Three-step regulatory network inference
- **GRaNIE**: Gene regulatory network inference including enhancers
- **DIRECTNET**: Discover cis-Regulatory Elements and Construct TF regulatory NETwork

### Included Datasets
- **K562**: Human chronic myelogenous leukemia cell line
- **Macrophage_S1**: Mouse macrophage stimulation condition 1
- **Macrophage_S2**: Mouse macrophage stimulation condition 2  
- **mESC_E7.5_rep1**: Mouse embryonic stem cells E7.5 replicate 1
- **mESC_E7.5_rep2**: Mouse embryonic stem cells E7.5 replicate 2
- **mESC_E8.5_rep1**: Mouse embryonic stem cells E8.5 replicate 1
- **mESC_E8.5_rep2**: Mouse embryonic stem cells E8.5 replicate 2

## 📁 Directory Structure

```
dataset_collection/
├── README.md                           # This file
├── LICENSE                             # License information
├── benchmark_new_method.R              # Benchmarking function for new methods
├── reproduce_ROC_PR_plots.R            # ROC and PR curve generation
├── reproduce_early_metrics.R           # Precision/Recall/F1-score analysis
├── analyze_network_stability.R         # Network stability analysis using Jaccard Index
├── example_usage.R                     # Usage examples
├── DATASETS/                           # Main inference results
│   ├── K562/
│   │   ├── CellOracle/
│   │   │   └── inference_results.tsv
│   │   ├── SCENIC+/
│   │   │   └── inference_results.tsv
│   │   ├── Pando/
│   │   │   └── inference_results.tsv
│   │   ├── LINGER/
│   │   │   └── inference_results.tsv
│   │   ├── FigR/
│   │   │   └── inference_results.tsv
│   │   ├── TRIPOD/
│   │   │   └── inference_results.tsv
│   │   ├── GRaNIE/
│   │   │   └── inference_results.tsv
│   │   └── DIRECTNET/
│   │       └── inference_results.tsv
│   ├── Macrophage_S1/
│   │   └── [same method structure]
│   ├── Macrophage_S2/
│   │   └── [same method structure]
│   ├── mESC_E7.5_rep1/
│   │   └── [same method structure]
│   ├── mESC_E7.5_rep2/
│   │   └── [same method structure]
│   ├── mESC_E8.5_rep1/
│   │   └── [same method structure]
│   └── mESC_E8.5_rep2/
│       └── [same method structure]
├── GROUND_TRUTHS/
│   ├── filtered_RN117_K562.tsv
│   ├── filtered_RN204_Buffer1.tsv
│   ├── filtered_RN111_E7.5_rep1.tsv
│   ├── filtered_RN111_E7.5_rep2.tsv
│   ├── filtered_RN111_E8.5_rep1.tsv
│   └── filtered_RN111_E8.5_rep2.tsv
└── STABILITY_GRNS/                     # Additional datasets for stability analysis (optional)
    ├── K562/
    ├── Macrophage_S1/
    └── [other samples with multiple replicates]
```

## 📋 Analysis Functions Overview

This benchmark collection includes four main analysis functions:

### 1. `benchmark_new_method()` - New Method Benchmarking
Evaluates new GRN inference methods against existing methods and ground truth.

### 2. `reproduce_ROC_PR_plots()` - ROC/PR Curve Analysis  
Generates comprehensive ROC and Precision-Recall curves for all methods across datasets.

### 3. `reproduce_early_metrics()` - Top-K Edge Analysis
Computes precision, recall, and F1-scores using filtered networks (top 10K edges).

### 4. `analyze_network_stability()` - Stability Analysis
Analyzes network stability across replicates using Jaccard Index.

## 📋 File Formats and Column Structure

### Method-Specific Column Names

| Method | TF Column | Target Column | Score Column | File Format |
|--------|-----------|---------------|--------------|-------------|
| CellOracle | `source` | `target` | `coef_mean` | CSV |
| SCENIC+ | `Source` | `Target` | `Score` | TSV |
| Pando | `tf` | `target` | `estimate` | CSV |
| LINGER | `Source` | `Target` | `Score` | TSV |
| FigR | `Motif` | `DORC` | `Score` | CSV |
| TRIPOD | `TF` | `gene` | `abs_coef` | CSV |
| GRaNIE | `TF.name` | `gene.name` | `TF_gene.r` | CSV |
| DIRECTNET | `TF` | `Gene` | `NA` | CSV |

### Ground Truth Files
- **Format**: TSV (tab-separated values)
- **Columns**: `Source` (TF), `Target` (target gene)
- **Content**: Experimentally validated regulatory interactions

## 🚀 Quick Start

### 1. Download and Extract

```bash
# Download from Zenodo
wget https://zenodo.org/record/17087863/files/dataset_collection.zip

# Extract
unzip dataset_collection.zip
cd dataset_collection/
```

### 2. Load Required Functions

```r
# Load all analysis functions
source("benchmark_new_method.R")
source("reproduce_ROC_PR_plots.R") 
source("reproduce_early_metrics.R")
source("analyze_network_stability.R")

# Load required libraries (install if needed)
required_packages <- c("readr", "dplyr", "pROC", "PRROC", "ggplot2", 
                      "stringr", "tibble", "RColorBrewer", "plotrix",
                      "tidyr", "purrr", "viridis", "scales")
lapply(required_packages, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
})
```

## 📊 Analysis Examples

### Example 1: Benchmark a New Method

```r
# Benchmark your new method against existing methods
results <- benchmark_new_method(
  dataset_name = "K562",
  new_grn_file = "path/to/your/method_results.tsv",
  new_method_name = "YourNewMethod",
  tf_column = "TF",              # Your TF column name
  target_column = "Gene",        # Your target gene column name  
  score_column = "Confidence",   # Your confidence/score column name
  input_dir = "DATASETS",
  output_dir = "benchmark_results",
  ground_truth_dir = "GROUND_TRUTHS"
)

# View performance summary
print(results$results)
print(results$new_method_performance)
```

### Example 2: Generate ROC and PR Curves for All Methods

```r
# Generate comprehensive ROC and PR curves
roc_pr_results <- reproduce_ROC_PR_plots(
  input_dir = "DATASETS",
  output_dir = "ROC_PR_Results",
  ground_truth_dir = "GROUND_TRUTHS",
  use_filtered_approach = TRUE
)

# Access results for specific dataset
k562_results <- roc_pr_results$results[["K562"]]
print(k562_results$results)  # Performance metrics
```

### Example 3: Analyze Top-K Edge Performance

```r
# Compute precision/recall/F1 for filtered networks (top 10K edges)
metrics_results <- reproduce_early_metrics(
  input_dir = "DATASETS",
  output_dir = "filtered_metrics_results",
  ground_truth_dir = "GROUND_TRUTHS",
  max_edges = 10000
)

# View summary across all datasets
print(metrics_results$combined_summary)

# Access lollipop plots
precision_plot <- metrics_results$plots[["Precision"]]
recall_plot <- metrics_results$plots[["Recall"]]
f1_plot <- metrics_results$plots[["F1-score"]]
```

### Example 4: Network Stability Analysis

```r
# Analyze network stability using Jaccard Index
# Note: This requires multiple replicates per method/dataset
stability_results <- analyze_network_stability(
  base_dir = "STABILITY_GRNS",  # Directory with multiple replicates
  sample_dirs = c("K562", "Macrophage_S1", "Macrophage_S2"),
  methods = c("CellOracle", "FigR", "LINGER", "Pando"),
  top_percent = 0.1,    # Top 10% of edges
  random_percent = 0.1, # Random 10% of edges
  create_plot = TRUE,
  save_plot = TRUE,
  save_results = TRUE
)

# View results
print(stability_results$results)
print(stability_results$plot)  # Jaccard Index heatmap
```

### Example 5: Batch Analysis Across All Datasets

```r
# Run comprehensive analysis across all datasets
datasets <- c("K562", "Macrophage_S1", "Macrophage_S2", 
              "mESC_E7.5_rep1", "mESC_E7.5_rep2", 
              "mESC_E8.5_rep1", "mESC_E8.5_rep2")

# 1. Generate ROC/PR curves for all
all_roc_pr <- reproduce_ROC_PR_plots(
  input_dir = "DATASETS",
  output_dir = "Comprehensive_ROC_PR"
)

# 2. Analyze filtered network metrics
all_metrics <- reproduce_early_metrics(
  input_dir = "DATASETS", 
  output_dir = "Comprehensive_Metrics"
)

# 3. Compare summary statistics
print(all_roc_pr$summary)
print(all_metrics$combined_summary)
```

## 🔧 Advanced Usage

### Custom Method Configuration

```r
# For methods with non-standard column names
custom_method_info <- tribble(
  ~method,        ~tf_col,     ~target_col,   ~score_col,
  "MyMethod",     "regulator", "regulated",   "weight",
  "CustomGRN",    "TF_name",   "target_name", "score"
)

# Use in ROC/PR analysis
results <- reproduce_ROC_PR_plots(
  input_dir = "DATASETS",
  output_dir = "custom_results", 
  method_info = custom_method_info
)
```

### Debug Network Files

```r
# Debug file structure issues
debug_network_files(
  base_dir = "DATASETS",
  sample = "K562", 
  method = "DIRECTNET"
)
```

### Custom Visualization

```r
# Create custom stability analysis with different parameters
stability_custom <- analyze_network_stability(
  base_dir = "STABILITY_GRNS",
  sample_dirs = c("K562"),
  methods = c("CellOracle", "FigR"),
  top_percent = 0.05,     # Top 5% edges
  random_percent = 0.05,  # Random 5% edges
  plot_width = 10,
  plot_height = 8,
  verbose = TRUE          # Enable detailed logging
)
```

## 📈 Output Files

### ROC/PR Analysis Outputs
- **`evaluation_results.csv`**: AUROC, AUPRC metrics for each method
- **`ROC_plot.pdf/png`**: ROC curve comparison plots
- **`PR_gap_plot.pdf`**: Precision-Recall curves with gap visualization
- **`all_datasets_summary.csv`**: Combined results across datasets

### Filtered Metrics Outputs  
- **`[Dataset]_filtered_metrics.tsv`**: Detailed metrics per dataset
- **`all_datasets_filtered_metrics.tsv`**: Combined summary
- **`Precision_filtered_lollipop_plot.pdf/png`**: Precision lollipop plots
- **`Recall_filtered_lollipop_plot.pdf/png`**: Recall lollipop plots  
- **`F1_score_filtered_lollipop_plot.pdf/png`**: F1-score lollipop plots

### Stability Analysis Outputs
- **`median_ji_results.csv`**: Jaccard Index results
- **`median_ji_heatmap.png`**: Stability heatmap visualization

### Benchmark Outputs
- **`benchmark_results.csv`**: Comparative performance metrics
- **`ROC_benchmark.pdf`**: ROC curves including new method
- **`PR_benchmark.pdf`**: PR curves including new method

## 📊 Performance Metrics

The analysis functions calculate:

### ROC/PR Analysis
- **AUROC**: Area Under the ROC Curve
- **AUPRC**: Area Under the Precision-Recall Curve  
- **AUPRC_random**: Random baseline AUPRC
- **Edge statistics**: Total edges, evaluable edges, positives, negatives

### Filtered Network Analysis
- **Precision**: TP / (TP + FP) for top-K edges
- **Recall**: TP / (TP + FN) for top-K edges
- **F1-Score**: Harmonic mean of precision and recall
- **Network statistics**: Original/filtered/final network sizes

### Stability Analysis
- **Jaccard Index**: |A ∩ B| / |A ∪ B| between network pairs
- **Median JI**: Median Jaccard Index across all pairwise comparisons
- **Top vs Random**: Comparison of stability in top-ranked vs random edges

## 🔧 Requirements

### R Version
- R ≥ 4.0.0

### Required R Packages
```r
install.packages(c(
  "readr",        # File reading
  "dplyr",        # Data manipulation  
  "pROC",         # ROC analysis
  "PRROC",        # PR analysis
  "ggplot2",      # Plotting
  "stringr",      # String operations
  "tibble",       # Data frames
  "RColorBrewer", # Color palettes
  "plotrix",      # Gap plots
  "tidyr",        # Data reshaping
  "purrr",        # Functional programming
  "viridis",      # Color scales
  "scales"        # Scale functions
))
```

## 📊 Ground Truth Information

Ground truth regulatory networks are derived from:
- ChIP-seq experiments
- Perturbation studies  
- Literature-curated databases
- Experimental validation studies

Each ground truth file contains high-confidence regulatory interactions specific to the cell type and experimental condition.

## 🎯 Analysis Workflow Recommendations

### For Method Developers
1. **Start with `benchmark_new_method()`** to get initial performance comparison
2. **Use `reproduce_early_metrics()`** to understand top-edge performance  
3. **Apply `reproduce_ROC_PR_plots()`** for comprehensive curve analysis
4. **Employ `analyze_network_stability()`** if you have multiple replicates

### For Comparative Studies
1. **Use `reproduce_ROC_PR_plots()`** for comprehensive method comparison
2. **Apply `reproduce_early_metrics()`** for practical performance assessment
3. **Combine with `analyze_network_stability()`** for robustness evaluation

### For Method Validation
1. **Start with filtered metrics** to understand practical performance
2. **Examine ROC/PR curves** for detailed performance characteristics
3. **Assess stability** across replicates to ensure reproducibility

## 📝 Citation

If you use this dataset collection in your research, please cite:

```bibtex
@dataset{your_dataset_2024,
  title={GRN Inference Method Benchmark Dataset Collection},
  author={Your Name and Collaborators},
  year={2024},
  publisher={Zenodo},
  doi={10.5281/zenodo.XXXXXXX},
  url={https://doi.org/10.5281/zenodo.XXXXXXX}
}
```

## 📄 License

This dataset is released under the [Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

You are free to:
- **Share**: Copy and redistribute the material
- **Adapt**: Remix, transform, and build upon the material

Under the following terms:
- **Attribution**: You must give appropriate credit

## 🤝 Contributing

We welcome contributions to improve this benchmark collection:

1. **New Methods**: Submit results from additional GRN inference methods
2. **New Datasets**: Contribute results on additional single-cell datasets
3. **Analysis Improvements**: Enhance existing analysis functions
4. **Bug Reports**: Report issues with the analysis pipeline
5. **Documentation**: Improve usage examples and documentation

## 📞 Support

For questions, issues, or contributions:
- **GitHub Issues**: [Link to repository issues]
- **Email**: [yuzun@pennstatehealth.psu.edu ; kkaramveer@pennstatehealth.psu.edu]
- **Documentation**: See function documentation and examples for detailed usage

## 🔄 Version History

- **v1.0.0**: Initial release with 8 methods across 7 datasets
- **v1.1.0**: Added comprehensive analysis functions
  - ROC/PR curve generation
  - Filtered network metrics analysis  
  - Network stability assessment
  - Enhanced benchmarking capabilities

## 🙏 Acknowledgments

This benchmark collection builds upon the excellent work of the developers of:
CellOracle, SCENIC+, Pando, LINGER, FigR, TRIPOD, DIRECTNET, and GRaNIE methods.

Special thanks to the single-cell genomics community for providing high-quality datasets and ground truth regulatory networks.

## 🚨 Troubleshooting

### Common Issues

#### Missing Columns
Check that your method's column names match the expected format in the method-specific table above.

#### Empty Results  
Ensure ground truth files exist and contain the expected `Source` and `Target` columns
