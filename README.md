# CAGE-seq Data Analysis Guide

This project provides a Python notebook for analyzing **CAGE-seq (Cap Analysis of Gene Expression)** data, focused on transcription start site (TSS) profiling and differential expression.

## Requirements

Before running the notebook, install the following Python packages:

```bash
pip install pandas numpy matplotlib seaborn plotly gseapy
````

For Jupyter Notebook support:

```bash
pip install notebook
```

## Directory Structure

Ensure your working directory contains:

```bash
project_folder/
│
├── CAGE-seq.ipynb         		# The main notebook
├── <your_data_files>.txt.gz            # One or more CAGE-seq data files
├── README.md                           # This file
```


- The script automatically detects and loads **all `.txt.gz` files** in the directory.
- Each file should be a tab-separated text file, optionally compressed with gzip.
- The filename (without extension) is used as a label for each sample or cell line.

### Flexible File Input

The notebook doesn't depends on a file naming format. You can name your CAGE-seq files anything (e.g., `sample1.txt.gz`, `HeLa_rep1.txt.gz`) and they will still be processed, as long as:
- They're compressed (`.txt.gz`)
- They're tab-delimited
- They contain appropriate expression-related columns (`log2FoldChange`, `padj`, etc.)

## Running the Analysis

1. Launch the notebook:

```bash
jupyter notebook CAGE-seq.ipynb
```

2. Step through each cell in order:

   * It will automatically detect and load the CAGE `.txt.gz` files.
   * Visualizations (MA plots, volcano plots, density, etc.) are generated for each cell line.
   * Gene Ontology enrichment is done on significant genes.
   * Interactive plots are rendered for deeper inspection.

## Outputs & Saving

* **Static plots** are shown inline in the notebook.
* **Interactive volcano plots** are saved in:

  ```
  all_volcano_plots.html
  ```
* If needed, save any plot by right-clicking in the notebook or using `plt.savefig()` inside the code.

## Notes

* If your data doesn't follow the naming convention, edit the `cage_files = [...]` line accordingly.
* Ensure gene names, chromosome IDs, and significance columns exist in your dataset.

## Purpose

This tool is ideal for:

* Biologists exploring TSS expression changes
* Bioinformaticians visualizing transcriptional dynamics
* Cross-cell-line comparisons of gene knockdown effects

## Questions?

Open an issue for help with usage or data compatibility.
