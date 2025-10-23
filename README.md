# Human liver RNA-Seq gene expression (903 samples)

ARCHS4 gene expression from human liver tissue

## About Dataset

### Context
RNA sequencing (RNA-seq) is the leading technology for genome-wide transcript quantification. However, publicly available RNA-seq data is currently provided mostly in raw form, a significant barrier for global and integrative retrospective analyses. ARCHS4 is a web resource that makes the majority of published RNA-seq data from human and mouse available at the gene and transcript levels. For developing ARCHS4, available FASTQ files from RNA-seq experiments from the Gene Expression Omnibus (GEO) were aligned using a cloud-based infrastructure. In total 187,946 samples are accessible through ARCHS4 with 103,083 mouse and 84,863 human. Additionally, the ARCHS4 web interface provides an intuitive exploration of the processed data through querying tools, interactive visualization, and gene pages that provide average expression across cell lines and tissues, top co-expressed genes for each gene, and predicted biological functions and protein–protein interactions for each gene based on prior knowledge combined with co-expression.

### Content
This is a subset of the total gene expression contained within ARCHS4. Specifically, this data only contains samples matching human liver samples. The dataset contains 903 unique samples from 60 distinct experiments created by a diverse group of researchers. The data is provided as a simple tab-separated file with the columns representing the samples and the rows are 35238 genes encoded as HUGO gene symbols.

### Inspiration
This is a good example of high dimensional data. It can be used to test visualizations techniques as well as batch effect detection and removal.

### Download

```
import kagglehub

# Download latest version
path = kagglehub.dataset_download("rana2hin/rna-seq-example-data")

print("Path to dataset files:", path)
```
