#!/usr/bin/env python3
"""
Part 2: Visualization and Advanced Analysis
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
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

print("="*80)
print("VISUALIZATIONS AND ADVANCED ANALYSIS")
print("="*80)

# Load data
expr_df = pd.read_csv('/home/user/rna-seq_liver/data/human_liver_expression.tsv', sep='\t', index_col=0)
metadata = pd.read_csv('/home/user/rna-seq_liver/data/sample_metadata.tsv', sep='\t')

# ============================================================================
# 4. VISUALIZATIONS
# ============================================================================

# VIZ 1: Library size distribution
print("\n[VIZ 1] Creating library size distribution plot...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

library_sizes = expr_df.sum(axis=0)
axes[0].hist(library_sizes, bins=30, edgecolor='black', alpha=0.7)
axes[0].set_xlabel('Library Size (Total Counts)')
axes[0].set_ylabel('Number of Samples')
axes[0].set_title('Distribution of Library Sizes')
axes[0].axvline(library_sizes.mean(), color='red', linestyle='--', label=f'Mean: {library_sizes.mean():,.0f}')
axes[0].axvline(library_sizes.median(), color='green', linestyle='--', label=f'Median: {library_sizes.median():,.0f}')
axes[0].legend()

# Library size by experiment
lib_by_exp = pd.DataFrame({'library_size': library_sizes, 'experiment': metadata['experiment'].values})
axes[1].boxplot([lib_by_exp[lib_by_exp['experiment'] == exp]['library_size'].values
                 for exp in sorted(metadata['experiment'].unique())],
                labels=sorted(metadata['experiment'].unique()))
axes[1].set_xlabel('Experiment')
axes[1].set_ylabel('Library Size')
axes[1].set_title('Library Size Distribution by Experiment')
axes[1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig('/home/user/rna-seq_liver/results/plots/01_library_sizes.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: results/plots/01_library_sizes.png")

# VIZ 2: Expression distribution
print("\n[VIZ 2] Creating expression distribution plots...")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Raw counts distribution
axes[0, 0].hist(expr_df.values.flatten(), bins=100, edgecolor='black', alpha=0.7)
axes[0, 0].set_xlabel('Expression Level')
axes[0, 0].set_ylabel('Frequency')
axes[0, 0].set_title('Raw Expression Distribution (All Values)')
axes[0, 0].set_xlim(0, np.percentile(expr_df.values, 99))

# Log-transformed (non-zero values)
non_zero = expr_df.values[expr_df.values > 0]
axes[0, 1].hist(np.log2(non_zero + 1), bins=100, edgecolor='black', alpha=0.7)
axes[0, 1].set_xlabel('Log2(Expression + 1)')
axes[0, 1].set_ylabel('Frequency')
axes[0, 1].set_title('Log-Transformed Expression Distribution')

# Mean expression per gene
gene_means = expr_df.mean(axis=1)
axes[1, 0].hist(np.log2(gene_means + 1), bins=50, edgecolor='black', alpha=0.7)
axes[1, 0].set_xlabel('Log2(Mean Expression + 1)')
axes[1, 0].set_ylabel('Number of Genes')
axes[1, 0].set_title('Distribution of Mean Gene Expression')

# Variance vs Mean (Dispersion plot)
gene_vars = expr_df.var(axis=1)
axes[1, 1].scatter(np.log2(gene_means + 1), np.log2(gene_vars + 1), alpha=0.3, s=10)
axes[1, 1].set_xlabel('Log2(Mean Expression + 1)')
axes[1, 1].set_ylabel('Log2(Variance + 1)')
axes[1, 1].set_title('Mean-Variance Relationship')
# Add diagonal line
x_range = axes[1, 1].get_xlim()
axes[1, 1].plot(x_range, x_range, 'r--', alpha=0.5, label='y=x')
axes[1, 1].legend()

plt.tight_layout()
plt.savefig('/home/user/rna-seq_liver/results/plots/02_expression_distributions.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: results/plots/02_expression_distributions.png")

# VIZ 3: Sample correlation heatmap
print("\n[VIZ 3] Creating sample correlation heatmap...")
# Use subset for better visualization if many samples
n_samples = min(50, expr_df.shape[1])
sample_subset = expr_df.iloc[:, :n_samples]
sample_corr = sample_subset.corr(method='spearman')

fig, ax = plt.subplots(figsize=(12, 10))
sns.heatmap(sample_corr, cmap='RdYlBu_r', center=0.5,
            square=True, linewidths=0, cbar_kws={"shrink": 0.8},
            xticklabels=False, yticklabels=False)
ax.set_title(f'Sample-to-Sample Correlation Matrix (Spearman, n={n_samples})')
plt.tight_layout()
plt.savefig('/home/user/rna-seq_liver/results/plots/03_sample_correlation.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: results/plots/03_sample_correlation.png")

# VIZ 4: Hierarchical clustering
print("\n[VIZ 4] Creating hierarchical clustering dendrogram...")
# Use log-transformed data for clustering
log_expr = np.log2(expr_df + 1)

# Calculate linkage
linkage_matrix = linkage(log_expr.T, method='ward')

fig, ax = plt.subplots(figsize=(14, 6))
dendrogram(linkage_matrix, labels=expr_df.columns, ax=ax, leaf_font_size=6)
ax.set_xlabel('Samples')
ax.set_ylabel('Ward Distance')
ax.set_title('Hierarchical Clustering of Samples (Ward Linkage)')
plt.xticks([])
plt.tight_layout()
plt.savefig('/home/user/rna-seq_liver/results/plots/04_hierarchical_clustering.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: results/plots/04_hierarchical_clustering.png")

# ============================================================================
# 5. DIMENSIONALITY REDUCTION - PCA
# ============================================================================
print("\n[5] PRINCIPAL COMPONENT ANALYSIS (PCA)")
print("-" * 80)

# Prepare data for PCA (log-transform and standardize)
log_expr = np.log2(expr_df + 1)
scaler = StandardScaler()
scaled_expr = scaler.fit_transform(log_expr.T)

# Perform PCA
pca = PCA()
pca_result = pca.fit_transform(scaled_expr)

# Variance explained
var_explained = pca.explained_variance_ratio_
cum_var_explained = np.cumsum(var_explained)

print(f"  - PC1 variance explained: {100*var_explained[0]:.2f}%")
print(f"  - PC2 variance explained: {100*var_explained[1]:.2f}%")
print(f"  - PC3 variance explained: {100*var_explained[2]:.2f}%")
print(f"  - Cumulative variance (PC1-10): {100*cum_var_explained[9]:.2f}%")

# Create PCA DataFrame
pca_df = pd.DataFrame({
    'PC1': pca_result[:, 0],
    'PC2': pca_result[:, 1],
    'PC3': pca_result[:, 2],
    'sample_id': expr_df.columns,
    'experiment': metadata['experiment'].values,
    'batch': metadata['batch'].values
})

# VIZ 5: Scree plot and PCA
print("\n[VIZ 5] Creating PCA visualizations...")
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

# Scree plot
n_components = min(20, len(var_explained))
axes[0, 0].bar(range(1, n_components + 1), 100*var_explained[:n_components], alpha=0.7)
axes[0, 0].plot(range(1, n_components + 1), 100*cum_var_explained[:n_components],
                'ro-', linewidth=2, markersize=6)
axes[0, 0].set_xlabel('Principal Component')
axes[0, 0].set_ylabel('Variance Explained (%)')
axes[0, 0].set_title('Scree Plot')
axes[0, 0].legend(['Cumulative', 'Individual'])
axes[0, 0].grid(True, alpha=0.3)

# PC1 vs PC2 colored by experiment
for exp in pca_df['experiment'].unique():
    mask = pca_df['experiment'] == exp
    axes[0, 1].scatter(pca_df.loc[mask, 'PC1'], pca_df.loc[mask, 'PC2'],
                      label=exp, alpha=0.7, s=50)
axes[0, 1].set_xlabel(f'PC1 ({100*var_explained[0]:.2f}%)')
axes[0, 1].set_ylabel(f'PC2 ({100*var_explained[1]:.2f}%)')
axes[0, 1].set_title('PCA: PC1 vs PC2 (Colored by Experiment)')
axes[0, 1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
axes[0, 1].grid(True, alpha=0.3)

# PC1 vs PC2 colored by batch
for batch in pca_df['batch'].unique():
    mask = pca_df['batch'] == batch
    axes[1, 0].scatter(pca_df.loc[mask, 'PC1'], pca_df.loc[mask, 'PC2'],
                      label=batch, alpha=0.7, s=50)
axes[1, 0].set_xlabel(f'PC1 ({100*var_explained[0]:.2f}%)')
axes[1, 0].set_ylabel(f'PC2 ({100*var_explained[1]:.2f}%)')
axes[1, 0].set_title('PCA: PC1 vs PC2 (Colored by Batch)')
axes[1, 0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
axes[1, 0].grid(True, alpha=0.3)

# PC2 vs PC3
for batch in pca_df['batch'].unique():
    mask = pca_df['batch'] == batch
    axes[1, 1].scatter(pca_df.loc[mask, 'PC2'], pca_df.loc[mask, 'PC3'],
                      label=batch, alpha=0.7, s=50)
axes[1, 1].set_xlabel(f'PC2 ({100*var_explained[1]:.2f}%)')
axes[1, 1].set_ylabel(f'PC3 ({100*var_explained[2]:.2f}%)')
axes[1, 1].set_title('PCA: PC2 vs PC3 (Colored by Batch)')
axes[1, 1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/user/rna-seq_liver/results/plots/05_pca_analysis.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: results/plots/05_pca_analysis.png")

# Save PCA results
pca_df.to_csv('/home/user/rna-seq_liver/results/pca_coordinates.tsv', sep='\t', index=False)
print("  ✓ PCA coordinates saved to: results/pca_coordinates.tsv")

# Continue to batch effect analysis...
print("\nVisualization and PCA analysis complete!")
