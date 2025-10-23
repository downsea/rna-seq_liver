#!/usr/bin/env python3
"""
Master script to run all RNA-seq analyses
Executes the complete analysis pipeline
"""

import subprocess
import time
import sys

scripts = [
    ('comprehensive_analysis.py', 'Basic data exploration and QC'),
    ('visualization_analysis.py', 'Visualizations and PCA'),
    ('batch_gene_analysis.py', 'Batch effects and gene analysis'),
    ('coexpression_report.py', 'Co-expression and report generation')
]

print("="*80)
print("RUNNING COMPLETE RNA-SEQ ANALYSIS PIPELINE")
print("="*80)

total_start = time.time()

for i, (script, description) in enumerate(scripts, 1):
    print(f"\n{'='*80}")
    print(f"Step {i}/{len(scripts)}: {description}")
    print(f"Script: {script}")
    print(f"{'='*80}\n")

    start_time = time.time()

    try:
        result = subprocess.run(['python', script], check=True, capture_output=False)
        elapsed = time.time() - start_time
        print(f"\n✓ {script} completed successfully in {elapsed:.2f} seconds")
    except subprocess.CalledProcessError as e:
        print(f"\n✗ Error running {script}: {e}")
        sys.exit(1)

total_elapsed = time.time() - total_start

print("\n" + "="*80)
print("ALL ANALYSES COMPLETED SUCCESSFULLY!")
print("="*80)
print(f"Total time: {total_elapsed:.2f} seconds ({total_elapsed/60:.2f} minutes)")
print("\nResults are available in: /home/user/rna-seq_liver/results/")
print("  - Review ANALYSIS_REPORT.md for comprehensive findings")
print("  - Check plots/ directory for all visualizations")
print("="*80)
