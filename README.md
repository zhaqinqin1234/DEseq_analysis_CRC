# 🧬 Gene Expression Analysis in Colorectal Cancer (CRC)

This project demonstrates spatial gene expression analysis in CRC using a simulated dataset, inspired by my NanoString GeoMx DSP poster presentation. The goal is to identify differentially expressed genes (DEGs) between epithelial and stromal regions and explore enriched biological pathways.

---

## 🔍 Project Goals

- Analyze simulated spatial transcriptomics data
- Identify DEGs using `DESeq2`
- Visualize data with volcano plots, heatmaps, and PCA
- Perform pathway enrichment using `clusterProfiler`

---

## 📁 Project Structure

```
CRC_GeneExpression_SpatialDSP/
├── data/                         # input data
├── analysis/                     # R Markdown analysis script
├── figures/                      # Visuals from the analysis
├── README.md                     # Project overview
```

---

## 🔧 Technologies & Tools

- **R**, **DESeq2**, **ggplot2**, **clusterProfiler**
- **R Markdown** for reproducible analysis
- Based on real-world biological use cases

---

## 📊 Sample Visuals

> _Coming soon: volcano plots, heatmaps, PCA, GO enrichment charts._

---

## 🧠 Biological Insights

- This pipeline mimics a real spatial transcriptomics study in CRC
- Enrichment results point to key extracellular matrix remodeling pathways
- Upregulated genes (e.g. `COL1A1`, `S100A9`, `NR4A1`) reflect tumor microenvironment shifts

---

## 🚀 How to Run

1. Clone this repository  
2. Open `CRC_HAnalysis.R` in RStudio  
3. Make sure these R packages are installed:
```r
tidyverse, DESeq2, EnhancedVolcano, pheatmap, clusterProfiler, org.Hs.eg.db
```

---

## 📌 Author

**QinQin Zha**  
Biologist | Data Scientist in Biotech  
Poster presented at: *2023 AACR*  


---

## 📄 License

MIT License (optional)