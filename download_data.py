#!/usr/bin/env python3
"""Download RNA-seq dataset from Kaggle"""

import kagglehub
import os

# Download latest version
print("Downloading dataset from Kaggle...")
path = kagglehub.dataset_download("rana2hin/rna-seq-example-data")

print(f"Path to dataset files: {path}")

# List downloaded files
print("\nDownloaded files:")
for root, dirs, files in os.walk(path):
    for file in files:
        file_path = os.path.join(root, file)
        file_size = os.path.getsize(file_path) / (1024 * 1024)  # Size in MB
        print(f"  - {file} ({file_size:.2f} MB)")

# Save the path for later use
with open('/home/user/rna-seq_liver/dataset_path.txt', 'w') as f:
    f.write(path)

print(f"\nDataset path saved to: /home/user/rna-seq_liver/dataset_path.txt")
