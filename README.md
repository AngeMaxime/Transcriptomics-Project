
# 🧬 RNA-seq Differential Expression Analysis — Transcriptomics

![R](https://img.shields.io/badge/Language-R-276DC3?style=flat&logo=r&logoColor=white)
![RMarkdown](https://img.shields.io/badge/Report-RMarkdown-blue?style=flat)
![HackBio](https://img.shields.io/badge/HackBio-Internship%202025-orange?style=flat)
![License](https://img.shields.io/badge/License-MIT-green?style=flat)

---

## 📌 Project Overview

This project is part of the **HackBio Internship 2025 — Task Code 2.6: Transcriptomics**.

This analysis explores a **processed RNA-seq dataset** from an experiment comparing:

- 🔴 A **diseased cell line** (control)
- 🟢 A **diseased cell line treated with Compound X** (treatment)

The difference in gene expression between the two conditions is quantified as **Log2 Fold Change (Log2FC)**, and the statistical significance of each change is reported as a **p-value** computed using DESeq2.

The goal is to identify which genes are **differentially expressed** — either activated or suppressed by the compound — and to understand their **biological roles** using GeneCards.

---

## 📂 Dataset

| Resource | Link |
|---|---|
| 📊 RNA-seq Results | [results.txt](https://gist.githubusercontent.com/stephenturner/806e31fce55a8b7175af/raw/1a507c4c3f9f1baaa3a69187223ff3d3050628d4/results.txt) |
| 📖 Original Source | [Stephen Turner — GitHub Gist](https://gist.github.com/stephenturner/806e31fce55a8b7175af) |

### Dataset Structure

The dataset contains **DESeq2 output** with the following columns:

| Column | Description |
|---|---|
| `Gene` | Gene symbol |
| `baseMean` | Average normalised expression across all samples |
| `log2FoldChange` | Log2 ratio of expression: treated vs diseased |
| `lfcSE` | Standard error of the Log2FC estimate |
| `stat` | Wald test statistic |
| `pvalue` | Raw p-value from the Wald test |
| `padj` | Benjamini-Hochberg adjusted p-value |

---

## 🎯 Thresholds Used

| Expression Status | Log2FC Condition | p-value Condition |
|---|---|---|
| 🔴 **Upregulated** | > 1 | < 0.01 |
| 🔵 **Downregulated** | < −1 | < 0.01 |
| ⚪ **Not Significant** | between −1 and 1 | ≥ 0.01 |

> A **Log2FC > 1** means expression at least **doubled** in treated cells.
> A **Log2FC < −1** means expression at least **halved** in treated cells.

---

## ✅ Tasks Completed

### 1. 📦 Package Loading
Packages used: `tidyverse`, `ggplot2`, `ggrepel`, `scales`
Managed via `pacman::p_load()` for automatic installation and loading.

---

### 2. 📥 Data Loading & Inspection
- Dataset loaded directly from the GitHub Gist URL
- Raw data saved locally as `Raw_rnaseq.csv`
- Overview performed: `dim()`, `names()`, `head()`, `summary()`, `colSums(is.na())`

---

### 3. 🧹 Data Cleaning & Preparation
- Rows with `NA` in `pvalue` or `log2FoldChange` removed (genes with insufficient counts)
- New column `neg_log10_pval` computed as `-log10(pvalue)` for volcano plot Y axis
- Each gene classified into `Upregulated`, `Downregulated`, or `Not Significant` using `case_when()`

---

### 4. 🌋 Volcano Plot

A **volcano plot** visualises differential expression simultaneously showing:
- **X axis** → Log2FoldChange (effect size)
- **Y axis** → −log10(p-value) (statistical significance)

Two volcano plots were generated:

| Plot | Description |
|---|---|
| `volcano_plot` | First version with top 5 up and top 5 down genes labelled |
| `Volcano_plot_2` | Enhanced version with quadrant counts annotated |

Both plots include:
- Dashed threshold lines at `|Log2FC| = 1` and `p-value = 0.01`
- Color-coded points: 🔴 red = upregulated, 🔵 blue = downregulated, ⚪ grey = not significant
- `ggrepel` labels for top genes — smart non-overlapping placement

**Saved as:** `Volcano1.png`

---

### 5. 🔼 Upregulated Genes (Log2FC > 1 and p-value < 0.01)

Genes that are **MORE expressed** in the treated condition — pathways **activated** by Compound X.

#### Top 5 Upregulated Genes

| Gene | Full Name | Key Function |
|---|---|---|
| `DTHD1` | Death Domain Containing 1 | Apoptosis signalling via death domain interactions |
| `EMILIN2` | Elastin Microfibril Interfacer 2 | ECM structure, apoptosis activation, tumour suppression |
| `PI16` | Peptidase Inhibitor 16 | Peptidase inhibitor activity, negative regulation of cell growth |
| `C4orf45` | Sperm Microtubule Inner Protein 2 | Structural protein; associated with Hyperekplexia |
| `FAM180B` | Family With Sequence Similarity 180 Member B | Extracellular localisation; associated with Borderline Leprosy |

---

### 6. 🔽 Downregulated Genes (Log2FC < −1 and p-value < 0.01)

Genes that are **LESS expressed** in the treated condition — pathways **suppressed** by Compound X.

#### Top 5 Downregulated Genes

| Gene | Full Name | Key Function |
|---|---|---|
| `TBX5` | T-Box Transcription Factor 5 | Heart development; limb identity; Holt-Oram syndrome |
| `IFITM1` | Interferon Induced Transmembrane Protein 1 | Antiviral defence; restricts influenza A, Ebola, SARS-CoV-2 |
| `TNN` | Tenascin N | ECM remodelling; angiogenesis; cell migration in tumours |
| `COL13A1` | Collagen Type XIII Alpha 1 Chain | Transmembrane collagen; connective tissue function |
| `IFITM3` | Interferon Induced Transmembrane Protein 3 | Antiviral defence; CD225 superfamily member |

---

### 7. 📊 Summary Statistics

| Metric | Count |
|---|---|
| Total genes analysed | computed dynamically |
| Upregulated genes (FC > 1, p < 0.01) | computed dynamically |
| Downregulated genes (FC < −1, p < 0.01) | computed dynamically |
| Not significant | computed dynamically |

---

> **Overall conclusion:**
> Compound X acts through a **dual mechanism** — simultaneously activating
> pro-apoptotic signals while suppressing tumour-promoting ECM remodelling and
> inflammatory responses. This pattern is consistent with a compound that
> **re-sensitises diseased cells to normal regulatory control** — a hallmark
> of effective therapeutic intervention.

---

## 🔧 How to Run

### Prerequisites
Make sure you have **R** (≥ 4.0) and **RStudio** installed.

### Steps

```r
# 1. Clone or download this repository
# 2. Open task2_6_transcriptomics.Rmd in RStudio
# 3. Install pacman if not already installed
install.packages("pacman")

# 4. Click "Knit" to generate the full HTML report
# OR run chunk by chunk interactively
```

All packages are installed automatically via `pacman::p_load()`.

---

## 📦 Dependencies

| Package | Purpose |
|---|---|
| `tidyverse` | Data manipulation (`dplyr`, `tidyr`) + `ggplot2` |
| `ggplot2` | All visualisations including volcano plots |
| `ggrepel` | Smart non-overlapping gene labels on volcano plot |
| `scales` | Axis label formatting |

---

## 🔗 References

- Turner, S. (2015). *DESeq2 RNA-seq results* — [GitHub Gist](https://gist.github.com/stephenturner/806e31fce55a8b7175af)
- [GeneCards Human Gene Database](https://www.genecards.org)
- Love, M.I., Huber, W., Anders, S. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2*. Genome Biology.
- [HackBio Internship 2025](https://github.com/HackBio-Internship)

---

## 👤 Author

**Ange Maxime TCHOUTANG**

**HackBio Internship 2025 — Task Code 2.6: Transcriptomics**

---

## 📜 License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.
