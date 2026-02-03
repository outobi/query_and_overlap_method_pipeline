

# 🧬 Cell Type Inference via Modality-Agnostic Query and Overlap Methods

### Integrating Spatial Transcriptomics/Proteomics with Single-Cell Transcriptomics

------------------------------------------------------------------------

## 📖 Introduction

Recent advancements in human omics technologies—particularly
transcriptomics and proteomics—have enabled high-resolution mapping of
molecular and cellular heterogeneity underlying disease pathogenesis.
Spatial transcriptomics and proteomics allow direct profiling of
histopathological regions from FFPE tissue sections using platforms such
as **Bruker GeoMx**, **CosMx**, and **10X Genomics Visium/Xenium**,
while **laser capture microdissection (LCM)** combined with LC-MS/MS
enables regional proteomic quantification. The greater the alignment
between tissue morphology and molecular readouts, the more biologically
informative the spatial analysis becomes.

This repository presents an R-based workflow to infer **disease-enriched 
or associated cell types** within histopathological regions by **multi-omics
integration**. Conventional deconvolution methods perform well for
transcriptomic integrations (spatial vs single-cell), but they fail
across molecular modalities (e.g., RNA vs protein).

Here, we implement two **modality-agnostic** approaches—**Query** and
**Overlap**—to bridge this gap.

Specifically, I will first combine spatial transcriptomics and single
cell transcriptomics to derive cell type enrichment/depletion in a given
histopathological region from FFPE tissue section in this repository. 
Then I will implement these two methods to combine LCM spatial proteomics 
and single cell transcriptomics in another repository.

These methods were applied to idiopathic pulmonary fibrosis (IPF) datasets
and published in *Proteomes (2025, 13(1):3)*.
👉 DOI:
[10.3390/proteomes13010003](https://doi.org/10.3390/proteomes13010003)

------------------------------------------------------------------------

## 📊 Dataset Overview

-   **Spatial transcriptomics (Bruker GeoMx)**: 1085 genes across 7
    histopathological regions from human IPF and control lungs (Eyres et
    al., *Cell Reports*, 2022; PMID: 35977489).
-   **Single-cell transcriptomics (10x Genomics)**: 33,694 genes across
    89,326 cells from IPF and control donors (Habermann et al., *Science
    Advances*, 2020; PMID: 32832598).
    After filtering, 24,470 genes remained—95.6% overlap with the
    spatial dataset.

Five broad cell classes (epithelial, mesenchymal, myeloid, endothelial,
lymphoid) subdivide into 31 fine-grained cell types, including
IPF-specific **HAS1hi fibroblasts**. Other 30 cell types have both
control donor and IPF patient resources.

------------------------------------------------------------------------

## 📁 Repository Structure



```         
query-and-overlap-method-pipeline/
├── README.md                               # This file
├── IPF multiomics integration Fei Wang 2025.pdf  # Workflow publication
├── Rscripts/
    ├── script 1 region specific gene/
    │   ├── input/
    │   │   ├── Supplementary_Table_1.csv   # sample metadata
    │   │   └── Supplementary_Table_2.xlsx  # expression table
    │   ├── code/
    │   │   └── region_specific_gene.R
    │   └── output/
    │       ├── region_specific_feature_list.rds
    │       ├── region_specific_gene_summary.xlsx
    │       └── region_specfic_gene.RData   # output RData
    │
    ├── script 2 query method enrichment/
    │   ├── input/
    │   │   ├── GSE135893_matrix.mtx            # Raw scRNA-seq expression matrix (large, GSE135893_matrix.mtx.gz PROVIDED IN https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893)
    │   │   ├── GSE135893_genes.tsv             # Gene names
    │   │   ├── GSE135893_barcodes.tsv          # Cell barcodes
    │   │   ├── GSE135893_IPF_metadata.csv      # Cell metadata
    │   │   └── region_specific_feature_list.rds # Output from Step 1
    │   ├── code/
    │   │   └── query_method_enrichment.R       # Main query-based enrichment script
    │   └── output/
    │       ├── Healthy_alveoli_signature.txt
    │       ├── Distal_alveoli_signature.txt
    │       ├── IPF_blood_vessel_signature.txt
    │       ├── Immune_infiltrate_signature.txt
    │       ├── Fibroblast_foci_signature.txt
    │       ├── 24470_gene_annotation.xlsx      # Gene annotation reference after filtering
    │       ├── query_method_enrichment.RData   # Large output; https://zenodo.org/records/17476683)
    │       ├── query_enrichment_score_summary_table.xlsx
    │       ├── query_positive_enrichment_bubble_plot_summary.pdf
    │       └── query_negative_depletion_bubble_plot_summary.pdf
    │
    └── script 3 overlap method enrichment/
        ├── input/
        │   └── query method enrichment.RData                      # Output from Step 2 (large; https://zenodo.org/records/17476683)
        ├── code/
        │   └── overlap method to calculate overlap FDR.R  # Main R script for overlap-based FDR calculation
        └── output/
            ├── overlap enrichment heatmap summary.pdf             # Visualization of overlap enrichment by heatmap
            ├── overlap positive enrichment bubble plot summary.pdf
            ├── overlap negative depletion bubble plot summary.pdf
            ├── summary of d value 30 cell types in control.xlsx   # Control cohort Cohen’s d summary
            ├── summary of d value 31 cell types in IPF.xlsx       # IPF cohort Cohen’s d summary
            └── overlap method to calculate overlap FDR.RData      # Final computed results (large; https://zenodo.org/records/17476683)
```

------------------------------------------------------------------------

## 🚀 Getting Started

### 🧩 Prerequisites

Ensure you have **R ≥ 4.4** installed along with the following R
packages:

``` r
install.packages(c(
  "dplyr", "rstatix", "readxl", "readr", "rio", "ggplot2", 
  "reshape2", "pheatmap", "RColorBrewer", "Matrix", "Seurat"
))
```

------------------------------------------------------------------------

### ⚙️ Installation

Clone this repository:

``` bash
git clone https://github.com/outobi/query-and-overlap-method-pipeline.git
cd query-and-overlap-method-pipeline
```

You can directly run the scripts in RStudio.

------------------------------------------------------------------------

## 🧠 Analysis Workflow

### **Step 1. Region-Specific Gene Extraction**

**Script:**
`Rscripts/script 1 region specific gene/region specific gene.R`

-   Load raw **GeoMx spatial transcriptomics** data and metadata.
-   Perform Wilcoxon rank-sum tests to identify upregulated
    region-specific genes.
-   Export significant region-specific feature lists
    (`region_specific_feature_list.rds`).

------------------------------------------------------------------------

### **Step 2. Query Method: Cell-Type Enrichment**

**Script:**
`Rscripts/script 2 query method enrichment/query enrichment method.R`

-   Import preprocessed **scRNA-seq Seurat object** and
    **`region_specific_feature_list.rds`**.
-   Match region-specific genes with scRNA-seq gene list.
-   Normalize expression across cells and genes.
-   Compute enrichment **z-scores** by summing average expression per
    cell type.
-   Visualize enrichment via **bubble plots**.
-   Output: `query method enrichment.RData`.

------------------------------------------------------------------------

### **Step 3. Overlap Method: Gene-Set Enrichment**

**Script:**
`Rscripts/script 3 overlap method enrichment/overlap method to calculate overlap FDR.R`

-   Load `query method enrichment.RData`.
-   Identify cell type–specific genes via **t-tests** and **Cohen’s d
    values**.
-   Compute **hypergeometric overlaps** between region- and cell-type
    gene sets.
-   Visualize enrichment via **bubble plots**.
-   Output: `overlap method to calculate overlap FDR.RData`.

------------------------------------------------------------------------

## 📈 Interpretation

Unlike classical deconvolution that estimates cell-type proportions,
these methods quantify **relative enrichment**:

-   **Query Method**: Measures summed region gene expression across cell
    types → higher positive *z-score* = stronger enrichment and higher negative *z-score* = stronger depletion. 
-   **Overlap Method**: Evaluates overlap significance between region
    and cell-type gene sets → larger overlap = stronger enrichment and smaller overlap = stronger depletion.

**Caution:**
- These results describe **relative enrichment** rather than absolute cell
type composition derived from common deconvolution methods.
- Compare enrichments **within the same region type**, not across regions.
- The Query method is more quantitative and reliable than the
Overlap method, which is more qualitative.
- Some regions do not have prominent region-specific genes
based on certain standards, like control blood vessel and IPF adjacent
alveoli. In this case, we do not recommend using this method to derive 
cell type enrichment.

**Remarkable biological insights:** 
- **Fibroblast foci** → enriched for
mesenchymal cells (fibroblasts, myofibroblasts, PLIN2+ fibroblasts)
- **Immune infiltrates** → enriched for myeloids and lymphocytes only
- **IPF alveoli** → enriched for epithelial cells compared with control
alveoli from healthy donors, especially transitional AT2 and aberrant
KRT5⁻/KRT17⁺ basaloid cells, implying potential epithelial damage and
dysorganized repair.

------------------------------------------------------------------------

## 📚 References

-   IPF multiomics integration: Wang et al. *Proteomes*, 2025, 13(1):3 — [DOI:
    10.3390/proteomes13010003](https://doi.org/10.3390/proteomes13010003)\
-   GeoMx spatial transcriptomics: Eyres et al., *Cell Reports*, 2022 — [DOI:
    10.1016/j.celrep.2022.111230](https://doi.org/10.1016/j.celrep.2022.111230)\
-   scRNAseq transcriptomics: Habermann et al., *Science Advances*, 2020 — [DOI:
    10.1126/sciadv.aba1972](https://www.science.org/doi/epdf/10.1126/sciadv.aba1972)\
-   Query Method Reference — [BMJ Open Resp. Res. 2023,
    10:e001391](https://bmjopenrespres.bmj.com/content/10/1/e001391)\
-   Overlap Method Reference — [Nat. Biotechnol. 2020,
    38:685–691](https://www.nature.com/articles/s41587-019-0392-8)\
-   General Review — [Nat. Rev. Genet. 2021,
    22:665–681](https://www.nature.com/articles/s41576-021-00370-8)

------------------------------------------------------------------------

## 🤝 Contributing

Contributions are welcome—please submit issues or pull requests on
GitHub.

------------------------------------------------------------------------

## 📧 Contact

For questions or collaborations, please open an issue or contact the
repository maintainer.
