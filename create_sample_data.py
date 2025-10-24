#!/usr/bin/env python3
"""
Create a realistic sample RNA-seq dataset for demonstration purposes.
This mimics the structure described in README: 903 samples x 35,238 genes
For demo purposes, we'll create a smaller but realistic version.
"""

import numpy as np
import pandas as pd
import os

# Set random seed for reproducibility
np.random.seed(42)

# Dataset parameters (creating a smaller version for demonstration)
n_samples = 100  # Instead of 903 for faster processing
n_genes = 1000   # Instead of 35,238 for faster processing
n_experiments = 10  # Instead of 60

print(f"Creating sample RNA-seq dataset...")
print(f"  - Samples: {n_samples}")
print(f"  - Genes: {n_genes}")
print(f"  - Experiments: {n_experiments}")

# Create realistic gene names (mix of common liver genes and generic names)
common_liver_genes = [
    'ALB', 'AFP', 'APOA1', 'APOB', 'APOC3', 'APOE', 'CYP3A4', 'CYP2E1',
    'G6PC', 'PCK1', 'FBP1', 'HMGCR', 'FASN', 'ACACA', 'SCD', 'FADS2',
    'GCK', 'HNF4A', 'HNF1A', 'CEBPA', 'NR1H4', 'RXRA', 'PPARA', 'PPARG',
    'SLC2A2', 'ALDH2', 'ADH1B', 'ADH1C', 'CPS1', 'OTC', 'ASS1', 'ASL',
    'ABCB11', 'ABCC2', 'SLC10A1', 'SLC22A1', 'SLCO1B1', 'SLCO1B3',
    'TTR', 'TF', 'HPX', 'HP', 'FGA', 'FGB', 'FGG', 'SERPINA1', 'SERPINC1'
]

# Generate gene names
gene_names = common_liver_genes[:min(50, n_genes)]
remaining = n_genes - len(gene_names)
gene_names.extend([f"GENE{i:05d}" for i in range(remaining)])

# Create sample metadata
sample_ids = [f"Sample_{i:03d}" for i in range(n_samples)]
experiment_ids = np.random.choice([f"Exp_{i:02d}" for i in range(n_experiments)],
                                   size=n_samples, replace=True)

# Create batch labels (some correlation with experiments)
batch_ids = []
for exp in experiment_ids:
    exp_num = int(exp.split('_')[1])
    # Experiments tend to cluster in batches
    batch = f"Batch_{(exp_num // 3):02d}"
    batch_ids.append(batch)

# Generate realistic RNA-seq count data
# Using negative binomial distribution which is typical for RNA-seq
expression_data = np.zeros((n_genes, n_samples))

for i in range(n_genes):
    # Mean expression varies by gene (some highly expressed, some lowly expressed)
    mean_expr = np.random.lognormal(mean=5, sigma=2)

    # Add batch effects (some genes affected more than others)
    batch_effect = np.random.random() < 0.3  # 30% of genes have batch effects

    for j in range(n_samples):
        # Base expression
        expr = np.random.negative_binomial(n=10, p=10/(10+mean_expr))

        # Add batch effect
        if batch_effect:
            batch_num = int(batch_ids[j].split('_')[1])
            expr = expr * (1 + 0.3 * np.sin(batch_num))

        # Add some biological variation between experiments
        exp_num = int(experiment_ids[j].split('_')[1])
        expr = expr * (1 + 0.1 * np.random.randn())

        expression_data[i, j] = max(0, expr)

# Create DataFrame
df = pd.DataFrame(expression_data, index=gene_names, columns=sample_ids)

# Save to TSV file
os.makedirs('/home/user/rna-seq_liver/data', exist_ok=True)
output_file = '/home/user/rna-seq_liver/data/human_liver_expression.tsv'
df.to_csv(output_file, sep='\t')

print(f"\n✓ Expression data saved to: {output_file}")
print(f"  Shape: {df.shape}")

# Create metadata file
metadata = pd.DataFrame({
    'sample_id': sample_ids,
    'experiment': experiment_ids,
    'batch': batch_ids
})

metadata_file = '/home/user/rna-seq_liver/data/sample_metadata.tsv'
metadata.to_csv(metadata_file, sep='\t', index=False)

print(f"✓ Sample metadata saved to: {metadata_file}")

# Display summary statistics
print(f"\nDataset Summary:")
print(f"  - Total expression values: {n_samples * n_genes:,}")
print(f"  - Mean expression: {df.values.mean():.2f}")
print(f"  - Median expression: {np.median(df.values):.2f}")
print(f"  - Non-zero values: {(df.values > 0).sum():,} ({100*(df.values > 0).sum()/(n_samples*n_genes):.1f}%)")
print(f"  - Unique experiments: {len(set(experiment_ids))}")
print(f"  - Unique batches: {len(set(batch_ids))}")

print(f"\nNote: This is a demonstration dataset. Replace with actual data using valid Kaggle credentials.")
