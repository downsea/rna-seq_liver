#!/usr/bin/env python3
"""
Comprehensive RNA-seq Analysis for Human Liver Dataset
Performs: data exploration, QC, visualization, batch effect detection,
differential expression, and co-expression analysis
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from scipy.stats import zscore
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
import warnings
warnings.filterwarnings('ignore')

# Set style for better visualizations
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['figure.dpi'] = 100

print("="*80)
print("COMPREHENSIVE RNA-SEQ ANALYSIS - HUMAN LIVER DATASET")
print("="*80)

# ============================================================================
# 1. DATA LOADING
# ============================================================================
print("\n[1] LOADING DATA")
print("-" * 80)

# Load expression data
expr_file = '/home/user/rna-seq_liver/data/human_liver_expression.tsv'
metadata_file = '/home/user/rna-seq_liver/data/sample_metadata.tsv'

print(f"Loading expression data from: {expr_file}")
expr_df = pd.read_csv(expr_file, sep='\t', index_col=0)

print(f"Loading metadata from: {metadata_file}")
metadata = pd.read_csv(metadata_file, sep='\t')

print(f"\n✓ Expression data shape: {expr_df.shape}")
print(f"  - Genes: {expr_df.shape[0]:,}")
print(f"  - Samples: {expr_df.shape[1]:,}")
print(f"\n✓ Metadata shape: {metadata.shape}")
print(f"  - Experiments: {metadata['experiment'].nunique()}")
print(f"  - Batches: {metadata['batch'].nunique()}")

# ============================================================================
# 2. DATA EXPLORATION
# ============================================================================
print("\n[2] DATA EXPLORATION")
print("-" * 80)

# Basic statistics
print("\nExpression Statistics:")
print(f"  - Mean: {expr_df.values.mean():.2f}")
print(f"  - Median: {np.median(expr_df.values):.2f}")
print(f"  - Std Dev: {expr_df.values.std():.2f}")
print(f"  - Min: {expr_df.values.min():.2f}")
print(f"  - Max: {expr_df.values.max():.2f}")

# Sparsity analysis
zero_counts = (expr_df == 0).sum().sum()
total_counts = expr_df.size
print(f"\n  - Zero values: {zero_counts:,} ({100*zero_counts/total_counts:.2f}%)")
print(f"  - Non-zero values: {total_counts - zero_counts:,} ({100*(total_counts-zero_counts)/total_counts:.2f}%)")

# Per-sample statistics
sample_means = expr_df.mean(axis=0)
sample_medians = expr_df.median(axis=0)
print(f"\nPer-Sample Statistics:")
print(f"  - Mean expression range: [{sample_means.min():.2f}, {sample_means.max():.2f}]")
print(f"  - Median expression range: [{sample_medians.min():.2f}, {sample_medians.max():.2f}]")

# Per-gene statistics
gene_means = expr_df.mean(axis=1)
gene_stds = expr_df.std(axis=1)
print(f"\nPer-Gene Statistics:")
print(f"  - Mean expression range: [{gene_means.min():.2f}, {gene_means.max():.2f}]")
print(f"  - Genes with mean > 1000: {(gene_means > 1000).sum()}")
print(f"  - Genes with mean < 100: {(gene_means < 100).sum()}")

# Top expressed genes
print(f"\nTop 10 Most Expressed Genes (by mean):")
top_genes = gene_means.nlargest(10)
for i, (gene, expr) in enumerate(top_genes.items(), 1):
    print(f"  {i:2d}. {gene:15s} - Mean expression: {expr:>10,.2f}")

# Save exploration results
exploration_results = {
    'dataset_shape': expr_df.shape,
    'total_values': total_counts,
    'zero_values': zero_counts,
    'sparsity_percent': 100*zero_counts/total_counts,
    'mean_expression': expr_df.values.mean(),
    'median_expression': np.median(expr_df.values),
    'std_expression': expr_df.values.std(),
    'n_experiments': metadata['experiment'].nunique(),
    'n_batches': metadata['batch'].nunique()
}

# ============================================================================
# 3. QUALITY CONTROL
# ============================================================================
print("\n[3] QUALITY CONTROL")
print("-" * 80)

# Create output directory for plots
import os
os.makedirs('/home/user/rna-seq_liver/results', exist_ok=True)
os.makedirs('/home/user/rna-seq_liver/results/plots', exist_ok=True)

# QC1: Library size distribution
print("\nQC1: Analyzing library sizes...")
library_sizes = expr_df.sum(axis=0)
print(f"  - Mean library size: {library_sizes.mean():,.0f}")
print(f"  - Median library size: {library_sizes.median():,.0f}")
print(f"  - Library size range: [{library_sizes.min():,.0f}, {library_sizes.max():,.0f}]")

# Identify potential outliers (using IQR method)
Q1 = library_sizes.quantile(0.25)
Q3 = library_sizes.quantile(0.75)
IQR = Q3 - Q1
outlier_samples = library_sizes[(library_sizes < Q1 - 1.5*IQR) | (library_sizes > Q3 + 1.5*IQR)]
print(f"  - Potential outlier samples (by library size): {len(outlier_samples)}")
if len(outlier_samples) > 0:
    print(f"    {list(outlier_samples.index)}")

# QC2: Gene detection rate
print("\nQC2: Analyzing gene detection rates...")
detected_genes = (expr_df > 0).sum(axis=1)
detection_rate = 100 * detected_genes / expr_df.shape[1]
print(f"  - Genes detected in all samples: {(detection_rate == 100).sum()}")
print(f"  - Genes detected in >50% samples: {(detection_rate > 50).sum()}")
print(f"  - Genes detected in <10% samples: {(detection_rate < 10).sum()}")

# QC3: Sample correlation
print("\nQC3: Computing sample-to-sample correlations...")
sample_corr = expr_df.corr(method='spearman')
mean_corr = sample_corr.values[np.triu_indices_from(sample_corr.values, k=1)].mean()
print(f"  - Mean pairwise correlation: {mean_corr:.3f}")
print(f"  - Min pairwise correlation: {sample_corr.values[np.triu_indices_from(sample_corr.values, k=1)].min():.3f}")
print(f"  - Max pairwise correlation: {sample_corr.values[np.triu_indices_from(sample_corr.values, k=1)].max():.3f}")

# Save QC results
qc_results = pd.DataFrame({
    'sample_id': expr_df.columns,
    'library_size': library_sizes.values,
    'n_genes_detected': (expr_df > 0).sum(axis=0).values,
    'mean_expression': sample_means.values,
    'median_expression': sample_medians.values
})
qc_results = qc_results.merge(metadata, on='sample_id')
qc_results.to_csv('/home/user/rna-seq_liver/results/qc_metrics.tsv', sep='\t', index=False)
print(f"\n✓ QC metrics saved to: results/qc_metrics.tsv")

print("\n[4] CREATING VISUALIZATIONS")
print("-" * 80)

# Continue to next section...
