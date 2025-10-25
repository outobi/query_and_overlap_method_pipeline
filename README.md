# Cell type inference by integrating spatial transcriptomics/proteomics and single-cell transcriptomics with modality-agnostic query method

## Introduction 
Recent advancements in human omics technologies, mainly transcriptomics and proteomics, have provided us with powerful tools to map the molecular and cellular landscape of disease pathogenesis. 
Spatial transcriptomics and spatial proteomics directly analyze defined histopathological regions from FFPE tissue sections. Platforms like Bruker GeoMx, CosMx and 10X genomics Visium, Xenium can capture RNA transcripts at regional to subcellular resolution. Similarly, laser capture microdissection (LCM) followed by LC-MS/MS proteomics can quantify protein expression specifically from these regions. The greater the alignment of morphology or histopathological structure with molecular and cellular analyses, the more informative of the study. 

Here, this workflow aims to infer disease-specific or driver cell types in histopathological regions in indications with complex spatial heterogeneity by multi-omics integration.
Conventional deconvolution approaches work well for mini-bulk/spatial transcriptomics and single-cell RNA-seq transcriptomics integration to infer cell type percentage composition in a specific region/spot, but they fail when we cross molecular modalities, such as proteins versus RNA. There has been very limited precedent for directly integrating spatial proteomics with single-cell transcriptomics for cell-type inference. This workflow presents my attempt to tackle that gap.

**Specifically, I focused on two key goals:**
- Developing an appropriate integration framework for LCM spatial proteomics and scRNA-seq data based on an existing spatial-transcriptomics–to–scRNA-seq integration methods
- Comparing the resulting spatial proteomics + scRNA-seq integration against spatial transcriptomics + scRNA-seq integration to assess performance and extract new biological insights






## 📚 References

1. **IPF paper reference:** https://www.mdpi.com/2227-7382/13/1/3

2. **RA paper reference:** https://www.mdpi.com/2227-7382/13/2/17

## 🤝 Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.


## 🙏 Acknowledgments


## 📧 Contact

For questions or collaborations, please open an issue on GitHub.




