# RNA-seq Analysis Pipeline for Human Liver Dataset

Comprehensive analysis pipeline for human liver RNA-seq gene expression data from ARCHS4.

## Overview

This repository contains a complete analysis pipeline for RNA-seq data analysis, featuring:
- Data quality control and preprocessing
- Dimensionality reduction (PCA)
- Batch effect detection
- Gene expression characterization
- Co-expression network analysis
- Publication-ready visualizations

## Dataset Information

**Source**: ARCHS4 (All RNA-seq and ChIP-seq sample and signature search)
**Tissue**: Human liver
**Expected Size**: 903 samples × 35,238 genes
**Data Format**: Tab-separated values (TSV)
**Kaggle Dataset**: `rana2hin/rna-seq-example-data`

## Requirements

### Python Packages
```bash
pip install kagglehub pandas numpy matplotlib seaborn scikit-learn scipy kaggle
```

### Kaggle API Setup
1. Create a Kaggle account at https://www.kaggle.com
2. Go to Account Settings → API → Create New API Token
3. Place `kaggle.json` in `~/.kaggle/`
4. Set permissions: `chmod 600 ~/.kaggle/kaggle.json`

## Quick Start

### Option 1: Run Complete Analysis Pipeline
```bash
# Run all analyses in sequence
python run_all_analyses.py
```

### Option 2: Run Individual Analysis Steps
```bash
# 1. Download or create sample data
python create_sample_data.py  # For demo
# OR with valid Kaggle credentials:
python download_data.py

# 2. Run analyses
python comprehensive_analysis.py        # Data exploration and QC
python visualization_analysis.py        # Visualizations and PCA
python batch_gene_analysis.py          # Batch effects and gene analysis
python coexpression_report.py          # Co-expression and report
```

## Analysis Components

### 1. Data Exploration (`comprehensive_analysis.py`)
- Basic dataset statistics
- Expression distribution analysis
- Library size quality control
- Gene detection rates
- Sample correlation analysis

**Outputs:**
- `results/qc_metrics.tsv` - Quality metrics for each sample

### 2. Visualizations (`visualization_analysis.py`)
- Library size distributions
- Expression value distributions
- Sample correlation heatmap
- Hierarchical clustering dendrogram
- PCA analysis with scree plots

**Outputs:**
- `results/plots/01_library_sizes.png`
- `results/plots/02_expression_distributions.png`
- `results/plots/03_sample_correlation.png`
- `results/plots/04_hierarchical_clustering.png`
- `results/plots/05_pca_analysis.png`
- `results/pca_coordinates.tsv`

### 3. Batch Effect Detection (`batch_gene_analysis.py`)
- ANOVA tests for batch separation in PC space
- Identification of batch-associated genes
- Library size and detection rate by batch
- PC distribution analysis by batch

**Outputs:**
- `results/plots/06_batch_effects.png`
- `results/plots/07_gene_expression.png`
- `results/gene_statistics.tsv`

### 4. Co-expression Analysis (`coexpression_report.py`)
- Gene-to-gene correlation analysis
- Highly correlated gene pair identification
- Gene clustering by expression patterns
- Comprehensive analysis report generation

**Outputs:**
- `results/plots/08_coexpression.png`
- `results/coexpressed_genes.tsv`
- `results/gene_clusters.tsv`
- `results/ANALYSIS_REPORT.md`
- `results/analysis_summary.tsv`

## Output Structure

```
rna-seq_liver/
├── data/
│   ├── human_liver_expression.tsv
│   └── sample_metadata.tsv
├── results/
│   ├── plots/
│   │   ├── 01_library_sizes.png
│   │   ├── 02_expression_distributions.png
│   │   ├── 03_sample_correlation.png
│   │   ├── 04_hierarchical_clustering.png
│   │   ├── 05_pca_analysis.png
│   │   ├── 06_batch_effects.png
│   │   ├── 07_gene_expression.png
│   │   └── 08_coexpression.png
│   ├── ANALYSIS_REPORT.md
│   ├── qc_metrics.tsv
│   ├── gene_statistics.tsv
│   ├── pca_coordinates.tsv
│   ├── coexpressed_genes.tsv
│   ├── gene_clusters.tsv
│   └── analysis_summary.tsv
├── comprehensive_analysis.py
├── visualization_analysis.py
├── batch_gene_analysis.py
├── coexpression_report.py
└── run_all_analyses.py
```

## Key Findings (Demo Dataset)

### Dataset Characteristics
- **Samples**: 100 (demo) / 903 (full dataset)
- **Genes**: 1,000 (demo) / 35,238 (full dataset)
- **Sparsity**: 0.57% zero values
- **Mean library size**: ~1M counts
- **Sample correlation**: High (mean r = 0.964)

### Quality Control
- 951 genes detected in all samples
- 998 genes detected in >50% of samples
- 1 potential outlier sample identified

### Batch Effects
- **Significant batch effects detected** in PC1 (p < 0.05)
- 9 genes with significant batch association
- Recommended: Apply batch correction (ComBat, limma)

### Gene Expression
- Top expressed genes include liver-specific markers (APOA1, FGA, etc.)
- High variability genes identified for potential biomarker discovery
- 5 distinct expression pattern clusters

## Analysis Methods

### Transformations
- **Log transformation**: log₂(counts + 1) for visualization and PCA
- **Z-score normalization**: For PCA and heatmaps

### Statistical Tests
- **ANOVA**: Batch effect detection
- **Pearson correlation**: Co-expression analysis
- **Spearman correlation**: Sample similarity

### Clustering
- **Hierarchical clustering**: Ward linkage
- **K-means**: Gene cluster identification

## Recommendations for Further Analysis

1. **Differential Expression Analysis**
   - Use DESeq2 or edgeR for robust statistical testing
   - Compare conditions/groups if metadata available

2. **Batch Correction**
   - Apply ComBat (sva package) or limma's removeBatchEffect
   - Re-run PCA after correction

3. **Functional Enrichment**
   - GO/KEGG pathway enrichment on gene clusters
   - Use tools like DAVID, Enrichr, or clusterProfiler

4. **Network Analysis**
   - Construct gene regulatory networks
   - Identify hub genes and key regulators

5. **Biomarker Discovery**
   - Focus on highly variable genes
   - Validate with literature and experimental data

## Troubleshooting

### Kaggle API Issues
- **403 Forbidden**: Check API credentials are valid
- **Dataset not found**: Verify dataset slug is correct
- **Connection timeout**: Check internet connection

### Memory Issues
- For large datasets, process in chunks
- Use sparse matrix representations if needed
- Consider using Dask for out-of-core computing

### Plot Quality
- Increase DPI: `plt.savefig(..., dpi=300)`
- Use vector formats: `.pdf` or `.svg`
- Adjust figure size for publication

## Citation

If you use this pipeline, please cite:

```
ARCHS4: Lachmann, A., Torre, D., Keenan, A.B. et al.
Massive mining of publicly available RNA-seq data from human and mouse.
Nat Commun 9, 1366 (2018). https://doi.org/10.1038/s41597-018-0061-y
```

## License

This analysis pipeline is provided as-is for educational and research purposes.

## Contact

For questions or issues, please open an issue on GitHub.

---

**Note**: The current results are based on a demonstration dataset created due to invalid Kaggle API credentials. To run the analysis on the full dataset (903 samples × 35,238 genes), provide valid Kaggle API credentials and run `download_data.py`.
