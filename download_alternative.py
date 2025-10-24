#!/usr/bin/env python3
"""Alternative download approach for RNA-seq dataset"""

import requests
import os

# Try downloading using opendatasets which might handle authentication differently
try:
    import opendatasets as od
    print("Using opendatasets to download...")
    od.download("https://www.kaggle.com/datasets/rana2hin/rna-seq-example-data",
                data_dir="/home/user/rna-seq_liver/data")
except ImportError:
    print("opendatasets not installed, installing...")
    os.system("pip install opendatasets -q")
    import opendatasets as od
    od.download("https://www.kaggle.com/datasets/rana2hin/rna-seq-example-data",
                data_dir="/home/user/rna-seq_liver/data")
