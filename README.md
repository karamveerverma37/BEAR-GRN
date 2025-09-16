# GRN Inference Method Benchmark Dataset Collection

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17087863.svg)](https://doi.org/10.5281/zenodo.17087863)

A comprehensive collection of Gene Regulatory Network (GRN) inference results from multiple state-of-the-art methods across various single-cell multiomics datasets, designed for benchmarking and comparative analysis of new GRN inference methods.

## 📊 Dataset Overview

This collection contains GRN inference results from **7 established methods** across **7 single-cell datasets**, along with corresponding ground truth regulatory networks for rigorous benchmarking.

### Included Methods
- **CellOracle**: Dynamic GRN inference from scRNA-seq
- **SCENIC+**: Multi-omics GRN inference 
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
├── benchmark_new_method.R              # Benchmarking function
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
│   │   └── GRaNIE/
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
└── GROUND_TRUTHS/
    ├── filtered_RN117_K562.tsv
    ├── filtered_RN204_Buffer1.tsv
    ├── filtered_RN111_E7.5_rep1.tsv
    ├── filtered_RN111_E7.5_rep2.tsv
    ├── filtered_RN111_E8.5_rep1.tsv
    └── filtered_RN111_E8.5_rep2.tsv
```

## 📋 File Formats and Column Structure

### Method-Specific Column Names

| Method | TF Column | Target Column | Score Column | File Format |
|--------|-----------|---------------|--------------|-------------|
| CellOracle | `source` | `target` | `coef_mean` | TSV |
| SCENIC+ | `Source` | `Target` | `Score` | TSV |
| Pando | `tf` | `target` | `estimate` | TSV |
| LINGER | `Source` | `Target` | `Score` | TSV |
| FigR | `Motif` | `DORC` | `Score` | TSV |
| TRIPOD | `TF` | `gene` | `abs_coef` | TSV |
| GRaNIE | `TF.name` | `gene.name` | `TF_gene.r` | TSV |

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

### 2. Load the Benchmarking Function

```r
# Load the benchmarking function
source("benchmark_new_method.R")

# Load required libraries (install if needed)
required_packages <- c("readr", "dplyr", "pROC", "PRROC", "ggplot2", 
                      "stringr", "tibble", "RColorBrewer", "plotrix")
lapply(required_packages, function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
})
```

### 3. Benchmark Your New Method

```r
# Example: Benchmark your new method against existing methods
results <- benchmark_new_method(
  dataset_name = "K562",
  new_grn_file = "path/to/your/method_results.tsv",
  new_method_name = "YourNewMethod",
  tf_column = "TF",              # Your TF column name
  target_column = "Gene",        # Your target gene column name  
  score_column = "Confidence",   # Your confidence/score column name
  input_dir = "DATASETS",        # Path to extracted datasets
  output_dir = "benchmark_results",
  ground_truth_dir = "GROUND_TRUTHS"
)

# View performance summary
print(results$results)

# Access your method's specific performance
print(results$new_method_performance)
```

### 4. Alternative: R Download and Setup

```r
# Download directly in R
download.file(
  "https://zenodo.org/record/17087863/files/dataset_collection.zip",
  "dataset_collection.zip"
)

# Extract
unzip("dataset_collection.zip")

# Set working directory
setwd("dataset_collection")

# Load function
source("benchmark_new_method.R")
```

## 📊 Usage Examples

### Example 1: Basic Benchmarking (New Method Only)

```r
# Benchmark just your new method against ground truth
results <- benchmark_new_method(
  dataset_name = "mESC_E7.5_rep1",
  new_grn_file = "my_grn_inference.csv",
  new_method_name = "DeepGRN",
  tf_column = "transcription_factor",
  target_column = "target_gene",
  score_column = "edge_weight",
  output_dir = "my_benchmark_results"
)
```

### Example 2: Comprehensive Comparison

```r
# Compare against all existing methods
results <- benchmark_new_method(
  dataset_name = "K562",
  new_grn_file = "novel_method_output.tsv", 
  new_method_name = "NovelGRN",
  tf_column = "source_tf",
  target_column = "target_gene",
  score_column = "regulatory_score",
  input_dir = "DATASETS",
  output_dir = "comprehensive_benchmark",
  ground_truth_dir = "GROUND_TRUTHS"
)

# Results are automatically saved to:
# - comprehensive_benchmark/benchmark_results.csv
# - comprehensive_benchmark/ROC_benchmark.pdf
# - comprehensive_benchmark/PR_benchmark.pdf
```

### Example 3: Batch Processing Multiple Datasets

```r
datasets <- c("K562", "mESC_E7.5_rep1", "mESC_E8.5_rep1")
all_results <- list()

for (dataset in datasets) {
  cat("Processing", dataset, "\n")
  results <- benchmark_new_method(
    dataset_name = dataset,
    new_grn_file = paste0("results_", dataset, ".tsv"),
    new_method_name = "MyMethod",
    tf_column = "TF",
    target_column = "Gene",
    score_column = "Score", 
    input_dir = "DATASETS",
    output_dir = paste0("benchmark_", dataset)
  )
  all_results[[dataset]] <- results
}
```

## 📈 Output Files

The `benchmark_new_method()` function generates:

1. **`benchmark_results.csv`**: Detailed performance metrics (AUROC, AUPRC, etc.)
2. **`ROC_benchmark.pdf`**: ROC curve comparison plot
3. **`PR_benchmark.pdf`**: Precision-Recall curve with gap plot
4. **Console output**: Performance summary and rankings

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
  "plotrix"       # Gap plots
))
```

## 📊 Performance Metrics

The benchmarking function calculates:

- **AUROC**: Area Under the ROC Curve
- **AUPRC**: Area Under the Precision-Recall Curve  
- **AUPRC_random**: Random baseline AUPRC
- **Edge statistics**: Total edges, evaluable edges, positives, negatives

## 🎯 Ground Truth Information

Ground truth regulatory networks are derived from:
- ChIP-seq experiments
- Perturbation studies
- Literature-curated databases
- Experimental validation studies

Each ground truth file contains high-confidence regulatory interactions specific to the cell type and experimental condition.

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
3. **Bug Reports**: Report issues with the benchmarking function
4. **Improvements**: Suggest enhancements to the evaluation pipeline

## 📞 Support

For questions, issues, or contributions:
- **GitHub Issues**: [Link to repository issues]
- **Email**: [yuzun@pennstatehealth.psu.edu ; kkaramveer@pennstatehealth.psu.edu]
- **Documentation**: See `example_usage.R` for additional examples

## 🔄 Version History

- **v1.0.0**: Initial release with 8 methods across 7 datasets
- **v1.1.0**: [Future] Additional methods and improved documentation

## 🙏 Acknowledgments

This benchmark collection builds upon the excellent work of the developers of:
CellOracle, SCENIC+, Pando, LINGER, FigR, TRIPOD, DIRECTNET and GRaNIE methods.

Special thanks to the single-cell genomics community for providing high-quality datasets and ground truth regulatory networks.
