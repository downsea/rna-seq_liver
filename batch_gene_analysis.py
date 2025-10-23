#!/usr/bin/env python3
"""
Part 3: Batch Effect Detection and Gene Expression Analysis
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")

print("="*80)
print("BATCH EFFECT DETECTION AND GENE EXPRESSION ANALYSIS")
print("="*80)

# Load data
expr_df = pd.read_csv('/home/user/rna-seq_liver/data/human_liver_expression.tsv', sep='\t', index_col=0)
metadata = pd.read_csv('/home/user/rna-seq_liver/data/sample_metadata.tsv', sep='\t')
pca_df = pd.read_csv('/home/user/rna-seq_liver/results/pca_coordinates.tsv', sep='\t')

# ============================================================================
# 6. BATCH EFFECT DETECTION
# ============================================================================
print("\n[6] BATCH EFFECT DETECTION")
print("-" * 80)

# Analyze batch effects in PCA space
print("\nAnalyzing batch effects in principal components...")

batches = pca_df['batch'].values
experiments = pca_df['experiment'].values

# Test if batches separate in PC space using ANOVA
from scipy.stats import f_oneway

# Create batch groups for PC1
batch_groups_pc1 = [pca_df[pca_df['batch'] == b]['PC1'].values
                     for b in pca_df['batch'].unique()]
f_stat_pc1, p_val_pc1 = f_oneway(*batch_groups_pc1)

# Create batch groups for PC2
batch_groups_pc2 = [pca_df[pca_df['batch'] == b]['PC2'].values
                     for b in pca_df['batch'].unique()]
f_stat_pc2, p_val_pc2 = f_oneway(*batch_groups_pc2)

print(f"  - ANOVA test for batch separation in PC1:")
print(f"    F-statistic: {f_stat_pc1:.4f}, p-value: {p_val_pc1:.4e}")
print(f"  - ANOVA test for batch separation in PC2:")
print(f"    F-statistic: {f_stat_pc2:.4f}, p-value: {p_val_pc2:.4e}")

if p_val_pc1 < 0.05 or p_val_pc2 < 0.05:
    print(f"\n  ⚠ WARNING: Significant batch effects detected (p < 0.05)")
    print(f"    Consider batch correction methods (e.g., ComBat, limma)")
else:
    print(f"\n  ✓ No significant batch effects detected in top PCs")

# Identify genes most affected by batch effects
print("\nIdentifying batch-associated genes...")
log_expr = np.log2(expr_df + 1)
batch_genes = []

for gene in expr_df.index[:100]:  # Test first 100 genes for demo
    gene_expr = log_expr.loc[gene].values
    batch_groups = [gene_expr[metadata['batch'] == b]
                    for b in metadata['batch'].unique()]
    try:
        f_stat, p_val = f_oneway(*batch_groups)
        batch_genes.append({
            'gene': gene,
            'f_statistic': f_stat,
            'p_value': p_val
        })
    except:
        continue

batch_genes_df = pd.DataFrame(batch_genes)
batch_genes_df['p_adjusted'] = batch_genes_df['p_value'] * len(batch_genes_df)  # Bonferroni correction
batch_genes_df = batch_genes_df.sort_values('p_value')

n_sig_batch = (batch_genes_df['p_adjusted'] < 0.05).sum()
print(f"  - Genes with significant batch association (p_adj < 0.05): {n_sig_batch}/{len(batch_genes_df)}")

if n_sig_batch > 0:
    print(f"\n  Top 5 batch-associated genes:")
    for i, row in batch_genes_df.head(5).iterrows():
        print(f"    - {row['gene']:15s}: p-value = {row['p_value']:.4e}")

# VIZ 6: Batch effect visualization
print("\n[VIZ 6] Creating batch effect visualizations...")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Library size by batch
lib_sizes = expr_df.sum(axis=0)
lib_by_batch = pd.DataFrame({'library_size': lib_sizes, 'batch': metadata['batch'].values})
batch_order = sorted(metadata['batch'].unique())
axes[0, 0].boxplot([lib_by_batch[lib_by_batch['batch'] == b]['library_size'].values
                     for b in batch_order],
                    labels=batch_order)
axes[0, 0].set_xlabel('Batch')
axes[0, 0].set_ylabel('Library Size')
axes[0, 0].set_title('Library Size Distribution by Batch')
axes[0, 0].tick_params(axis='x', rotation=45)

# Detected genes by batch
detected_by_batch = pd.DataFrame({
    'n_detected': (expr_df > 0).sum(axis=0).values,
    'batch': metadata['batch'].values
})
axes[0, 1].boxplot([detected_by_batch[detected_by_batch['batch'] == b]['n_detected'].values
                     for b in batch_order],
                    labels=batch_order)
axes[0, 1].set_xlabel('Batch')
axes[0, 1].set_ylabel('Number of Detected Genes')
axes[0, 1].set_title('Gene Detection Rate by Batch')
axes[0, 1].tick_params(axis='x', rotation=45)

# PC1 distribution by batch
axes[1, 0].boxplot([pca_df[pca_df['batch'] == b]['PC1'].values for b in batch_order],
                    labels=batch_order)
axes[1, 0].set_xlabel('Batch')
axes[1, 0].set_ylabel('PC1 Score')
axes[1, 0].set_title('PC1 Distribution by Batch')
axes[1, 0].tick_params(axis='x', rotation=45)

# PC2 distribution by batch
axes[1, 1].boxplot([pca_df[pca_df['batch'] == b]['PC2'].values for b in batch_order],
                    labels=batch_order)
axes[1, 1].set_xlabel('Batch')
axes[1, 1].set_ylabel('PC2 Score')
axes[1, 1].set_title('PC2 Distribution by Batch')
axes[1, 1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig('/home/user/rna-seq_liver/results/plots/06_batch_effects.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: results/plots/06_batch_effects.png")

# ============================================================================
# 7. GENE EXPRESSION ANALYSIS
# ============================================================================
print("\n[7] GENE EXPRESSION ANALYSIS")
print("-" * 80)

# Calculate comprehensive gene statistics
print("\nCalculating gene-level statistics...")
gene_stats = pd.DataFrame({
    'gene': expr_df.index,
    'mean_expression': expr_df.mean(axis=1).values,
    'median_expression': expr_df.median(axis=1).values,
    'std_expression': expr_df.std(axis=1).values,
    'cv': (expr_df.std(axis=1) / (expr_df.mean(axis=1) + 1)).values,  # Coefficient of variation
    'detection_rate': (100 * (expr_df > 0).sum(axis=1) / expr_df.shape[1]).values,
    'max_expression': expr_df.max(axis=1).values,
    'min_expression': expr_df.min(axis=1).values
})

# Add variance categories
gene_stats['variance_category'] = pd.cut(gene_stats['std_expression'],
                                          bins=[0, gene_stats['std_expression'].quantile(0.33),
                                                gene_stats['std_expression'].quantile(0.67),
                                                gene_stats['std_expression'].max()],
                                          labels=['Low', 'Medium', 'High'])

# Top expressed genes
print(f"\nTop 20 Most Highly Expressed Genes:")
top_20_genes = gene_stats.nlargest(20, 'mean_expression')
for i, row in top_20_genes.iterrows():
    print(f"  {i+1:2d}. {row['gene']:15s} - Mean: {row['mean_expression']:>10,.2f}, "
          f"Detection: {row['detection_rate']:>5.1f}%, CV: {row['cv']:.3f}")

# Most variable genes
print(f"\nTop 20 Most Variable Genes (by CV):")
most_variable = gene_stats.nlargest(20, 'cv')
for i, row in most_variable.iterrows():
    print(f"  {i+1:2d}. {row['gene']:15s} - CV: {row['cv']:>8.3f}, "
          f"Mean: {row['mean_expression']:>10,.2f}, Std: {row['std_expression']:>10,.2f}")

# Save gene statistics
gene_stats.to_csv('/home/user/rna-seq_liver/results/gene_statistics.tsv', sep='\t', index=False)
print(f"\n✓ Gene statistics saved to: results/gene_statistics.tsv")

# VIZ 7: Gene expression characteristics
print("\n[VIZ 7] Creating gene expression visualizations...")
fig, axes = plt.subplots(2, 3, figsize=(16, 10))

# Mean expression distribution
axes[0, 0].hist(np.log2(gene_stats['mean_expression'] + 1), bins=50, edgecolor='black', alpha=0.7)
axes[0, 0].set_xlabel('Log2(Mean Expression + 1)')
axes[0, 0].set_ylabel('Number of Genes')
axes[0, 0].set_title('Distribution of Mean Gene Expression')
axes[0, 0].grid(True, alpha=0.3)

# Detection rate distribution
axes[0, 1].hist(gene_stats['detection_rate'], bins=50, edgecolor='black', alpha=0.7)
axes[0, 1].set_xlabel('Detection Rate (%)')
axes[0, 1].set_ylabel('Number of Genes')
axes[0, 1].set_title('Distribution of Gene Detection Rates')
axes[0, 1].grid(True, alpha=0.3)

# Coefficient of variation
axes[0, 2].hist(gene_stats['cv'], bins=50, edgecolor='black', alpha=0.7)
axes[0, 2].set_xlabel('Coefficient of Variation')
axes[0, 2].set_ylabel('Number of Genes')
axes[0, 2].set_title('Distribution of Gene Variability (CV)')
axes[0, 2].set_xlim(0, np.percentile(gene_stats['cv'], 95))
axes[0, 2].grid(True, alpha=0.3)

# Mean vs detection rate
scatter = axes[1, 0].scatter(np.log2(gene_stats['mean_expression'] + 1),
                            gene_stats['detection_rate'],
                            c=gene_stats['cv'], cmap='viridis',
                            alpha=0.5, s=10)
axes[1, 0].set_xlabel('Log2(Mean Expression + 1)')
axes[1, 0].set_ylabel('Detection Rate (%)')
axes[1, 0].set_title('Mean Expression vs Detection Rate')
plt.colorbar(scatter, ax=axes[1, 0], label='CV')
axes[1, 0].grid(True, alpha=0.3)

# Top 30 expressed genes barplot
top_30 = gene_stats.nlargest(30, 'mean_expression')
axes[1, 1].barh(range(30), np.log2(top_30['mean_expression'].values + 1))
axes[1, 1].set_yticks(range(30))
axes[1, 1].set_yticklabels(top_30['gene'].values, fontsize=7)
axes[1, 1].set_xlabel('Log2(Mean Expression + 1)')
axes[1, 1].set_title('Top 30 Most Expressed Genes')
axes[1, 1].invert_yaxis()
axes[1, 1].grid(True, alpha=0.3, axis='x')

# Heatmap of top variable genes across samples
top_var_genes = gene_stats.nlargest(50, 'std_expression')['gene'].values
heatmap_data = log_expr.loc[top_var_genes].iloc[:, :min(50, expr_df.shape[1])]
# Standardize for better visualization
heatmap_data_scaled = (heatmap_data.T - heatmap_data.mean(axis=1)) / heatmap_data.std(axis=1)
sns.heatmap(heatmap_data_scaled.T, cmap='RdBu_r', center=0,
            xticklabels=False, yticklabels=True,
            cbar_kws={'label': 'Z-score'}, ax=axes[1, 2])
axes[1, 2].set_title('Top 50 Variable Genes (Z-scored)')
axes[1, 2].set_ylabel('Genes')
axes[1, 2].set_xlabel('Samples')
axes[1, 2].tick_params(axis='y', labelsize=6)

plt.tight_layout()
plt.savefig('/home/user/rna-seq_liver/results/plots/07_gene_expression.png', dpi=300, bbox_inches='tight')
plt.close()
print("  ✓ Saved: results/plots/07_gene_expression.png")

print("\nBatch effect and gene expression analysis complete!")
