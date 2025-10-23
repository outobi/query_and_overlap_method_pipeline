# query-and-overlap-method-pipeline
# 1 Introduction 
This workflow aims to integrate spatial transcriptomics with scRNA-seq transcriptomics to deduce cell type enrichment in histological regions.

Instead of estimating exact cell-type proportions as in conventional spatial deconvolution, our approach highlights relative cell-type enrichment across distinct regions. The advantage is its modality-agnostic potential. Firstly I will apply it in integrating spatial transcriptomics and sc transcriptomics, then more strikingly, in spatial proteomics and sc transcriptomics.
We start with region-specific up-regulated differential genes from a published GeoMx spatial transcriptomics and query their summed expression across cell types defined in the public scRNA-seq dataset. Both two datasets are from human ipf lung samples. For each region, we compute a Z-score to measure how strongly those region-representative genes are expressed within each cell type.  Higher z-score indicated higher enrichment of the given cell type in this region. it evaluates the probabilistic enrichment of cell states: given the transcriptomic program of a histology-defined region, which scRNA clusters show compatible expression signatures? A higher positive score means that cell type is more likely to be active in that region than in others.




