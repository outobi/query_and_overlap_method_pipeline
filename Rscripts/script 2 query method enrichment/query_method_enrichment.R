setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(Matrix)
library(Seurat)
library(rio)
library(readxl)
library(readr)
library(ggplot2)
library(reshape2)




## prepare scRNAseq data--------------------------
#Read in and process single cell data - taken from GEO accession number GSE135893
str(data <- readMM("GSE135893_matrix.mtx"))

#Formal class 'dgTMatrix' [package "Matrix"] with 6 slots
#..@ i       : int [1:338061955] 18 27 45 49 59 61 72 79 86 94 ...
#..@ j       : int [1:338061955] 0 0 0 0 0 0 0 0 0 0 ...
#..@ Dim     : int [1:2] 33694 220213
#..@ Dimnames:List of 2
#.. ..$ : NULL
#.. ..$ : NULL
#..@ x       : num [1:338061955] 1 1 4 2 1 1 6 2 1 3 ...
#..@ factors : list()


metadata <- read.csv("GSE135893_IPF_metadata.csv", stringsAsFactors = F)   # cell meta data
barcodes <- scan("GSE135893_barcodes.tsv", character())  # barcodes
genes <- scan("GSE135893_genes.tsv", character())         # genes
colnames(data) <- barcodes
rownames(data) <- genes

data <- data[,match(metadata$X, colnames(data))]   ## 220213 barcode to 114396 cells
IPF.inds <- grep("IPF", metadata$Diagnosis)  # 57682 cells from 12 IPF patients
Ctrl.inds <- grep("Control", metadata$Diagnosis)  # 31644 cells from 10 healthy control donors

# COMBINED all 57682 + 31644 = 89326 cells from IPF and healthy donors
Combine.inds <- grep(c("IPF|Control"), metadata$Diagnosis)
meta.combine <- metadata[Combine.inds,]
data.combine <- data[,Combine.inds]
rownames(meta.combine) <- meta.combine$X

# create seurat object with 89326 cell column and 33694 gene row
seurat.cts <- CreateSeuratObject(counts = data.combine, project = "Kropski", meta.data = meta.combine)


# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- seurat.cts@assays[["RNA"]]@layers[["counts"]]

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- (counts > 0)

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
sum(keep_genes == TRUE) #24470 genes passing the filter

# filtered gene list 24470 
tmp3 <- data.frame(c(1:33694), keep_genes)
tmp3 <- cbind(tmp3, "Gene.symbol" = genes)
tmp4 <- tmp3[tmp3$keep_genes == T,]
tmp4[[1]] <- c(1:24470) # reorder the index 
export(tmp4,"24470_gene_annotation.xlsx")


# overlap region specific genes with 24470 scRNAseq gene set-------------
region_specific_feature_list <- readRDS("region_specific_feature_list.rds")

#Fibroblast_foci_signature.txt  n = 50 to 50
tmp5 <- region_specific_feature_list$FF.v.all.wilcox$gene.symbol
tmp5 <- intersect(tmp5, tmp4$Gene.symbol)
tmp6 <- tmp4[tmp4$Gene.symbol %in% tmp5, ]
tmp7 <- tmp6[c(1,3)]
colnames(tmp7)[c(1,2)] <- c("index", "gene.symbol")

export(tmp7, "Fibroblast_foci_signature.txt")


#Distal_alveoli_signature.txt n = 41 to 41
tmp5 <- region_specific_feature_list$Dis.Alv.v.all.wilcox$gene.symbol
tmp5 <- intersect(tmp5, tmp4$Gene.symbol)
tmp6 <- tmp4[tmp4$Gene.symbol %in% tmp5, ]
tmp7 <- tmp6[c(1,3)]
colnames(tmp7)[c(1,2)] <- c("index", "gene.symbol")

export(tmp7, "Distal_alveoli_signature.txt")


#Immune_infiltrate_signature.txt n = 111 to 106
tmp5 <- region_specific_feature_list$Imm.v.all.wilcox$gene.symbol
tmp5 <- intersect(tmp5, tmp4$Gene.symbol)
tmp6 <- tmp4[tmp4$Gene.symbol %in% tmp5, ]
tmp7 <- tmp6[c(1,3)]
colnames(tmp7)[c(1,2)] <- c("index", "gene.symbol")

export(tmp7, "Immune_infiltrate_signature.txt")


#IPF_blood_vessel_signature.txt n = 28 to 27
tmp5 <- region_specific_feature_list$IPF.vessel.all.wilcox$gene.symbol
tmp5 <- intersect(tmp5, tmp4$Gene.symbol)
tmp6 <- tmp4[tmp4$Gene.symbol %in% tmp5, ]
tmp7 <- tmp6[c(1,3)]
colnames(tmp7)[c(1,2)] <- c("index", "gene.symbol")

export(tmp7, "IPF_blood_vessel_signature.txt")


#healthy_alveoli_signature.txt n = 32
tmp5 <- region_specific_feature_list$Heal.Alv.v.all.wilcox$gene.symbol
tmp5 <- intersect(tmp5, tmp4$Gene.symbol)
tmp6 <- tmp4[tmp4$Gene.symbol %in% tmp5, ]
tmp7 <- tmp6[c(1,3)]
colnames(tmp7)[c(1,2)] <- c("index", "gene.symbol")

export(tmp7, "Healthy_alveoli_signature.txt")

# Continue to prepare scRNAseq dataset--------------------
# Reassign to filtered Seurat object, 24470 gene expression in 31644 ctrl cells and 57682 IPF cells, 89326 cells in all
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = meta.combine)

# logNormalize in each cell
## normalized data  1. Divide each cell by the total number of molecules measured in the cell
#2. Multiply that number by a scaling factor (i.e. 10000 Add 1, and take a natural e log)

pbmc_combined <- NormalizeData(object = filtered_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
normalized_count <- pbmc_combined@assays[["RNA"]]@layers[["data"]]
remove(counts)
remove(filtered_counts)

# save space
remove(data)
remove(nonzero)
remove(seurat.cts)
remove(pbmc_combined)
remove(data.combine)

# get fibroblast foci region specific gene signature
ff_signature <- read_table("Fibroblast_foci_signature.txt",     
                           col_names = TRUE)

ff_sig_count <- normalized_count[c(ff_signature$index),]  # 50 gene * 89326 cells


# distal alveoli
distal_alv_signature <- read_table("Distal_alveoli_signature.txt", 
                           col_names = TRUE)

distal_alv_sig_count <- normalized_count[c(distal_alv_signature$index),]  ## 41 gene * 89326 cells

# immune infiltrate
imm_infil <- read_table("Immune_infiltrate_signature.txt", 
                                   col_names = TRUE)

imm_infil_sig_count <- normalized_count[c(imm_infil$index),]  ## 106 gene * 89326 cells

# ipf blood vessel
ipf_vessel_signature <- read_table("IPF_blood_vessel_signature.txt", 
                                   col_names = TRUE)

ipf_vessel_sig_count <- normalized_count[c(ipf_vessel_signature$index),]  ## 27 gene * 89326 cells


# healthy alveoli, pay attention, use single cell expression from healthy control donor cell expression profile

healthy_alv_signature <- read_table("Healthy_alveoli_signature.txt", 
                                   col_names = TRUE)

healthy_alv_sig_count <- normalized_count[c(healthy_alv_signature$index),]  ## 32 gene * 89326 cells


# cell type annotation list 31 in all

tmp1 <- c("AT1", "AT2", "B Cells", "Basal", "cDCs", "Ciliated", "Differentiating Ciliated",
          "Endothelial Cells", "Fibroblasts", "HAS1 High Fibroblasts", "KRT5-/KRT17+", 
          "Lymphatic Endothelial Cells", "Macrophages","Mast Cells", "Mesothelial Cells",
          "Monocytes", "MUC5AC+ High", "MUC5B+", "Myofibroblasts", "NK Cells",
          "pDCs" , "Plasma Cells", "PLIN2+ Fibroblasts", "Proliferating Epithelial Cells",
          "Proliferating Macrophages", "Proliferating T Cells", "SCGB3A2+", "SCGB3A2+ SCGB1A1+",
          "Smooth Muscle Cells", "T Cells", "Transitional AT2") ## store each cell type, 31 IN ALL

tmp1_control <- tmp1[-10]  # Control regions do not contain HAS1 high fibroblast cell type


## cell type enrichment score in fibroblast foci---------------------------------
test <- as.matrix(ff_sig_count)    ## dgCmatrix to matrix with HANDLABLE small size

## zi = xi-u/sigma normalization transform
z_test <- test  

test1 <- apply(test, 1, mean)  # calculate each gene mean and sd, like scale by row in heatmap
test2 <- apply(test, 1, sd)

for (i in 1:nrow(z_test)){
  z_test[i,] <- (test[i,] - test1[i])/test2[i]
}

# query in 31 cell types
ff.sig.score <- c()

for (i in 1:31) {
  cell.type.index <- which(meta.combine$celltype == tmp1[i] & meta.combine$Diagnosis == "IPF")
  cell.type.data <- z_test[, cell.type.index]
  ff.sig.score[i] <- sum(cell.type.data)/length(cell.type.index)       # average cell type expression of all 50 genes
}
## get ff.signature score in 31 cell types
ff.signature.score <- as.data.frame(cbind("Number" = c(1:31), "Cell.type" = tmp1, "ff.sig.score" = ff.sig.score))
ff.signature.score$ff.sig.score <- as.numeric(ff.signature.score$ff.sig.score)



## distal alveoli signature score-------------------------------------------------------------

test <- as.matrix(distal_alv_sig_count)    ## dgCmatrix to matrix with HANDLABLE small size, 11 genes 57682 cells

z_test <- test  ## zi = xi-u/sigma transform

test1 <- apply(test, 1, mean)  # calculate each gene mean and sd, like scale by row in heatmap
test2 <- apply(test, 1, sd)

for (i in 1:nrow(z_test)){
  z_test[i,] <- (test[i,] - test1[i])/test2[i]
}

distal.alv.sig.score <- c()
# repeat in 31 cell types
for (i in 1:31) {
  cell.type.index <- which(meta.combine$celltype == tmp1[i] & meta.combine$Diagnosis == "IPF")
  cell.type.data <- z_test[, cell.type.index]
  distal.alv.sig.score[i] <- as.numeric(sum(cell.type.data)/length(cell.type.index))
}
## get distal.alv.signature score in 31 cell types
distal.alv.signature.score <- as.data.frame(cbind("Number" = c(1:31), "Cell.type" = tmp1, "distal.alv.sig.score" = as.numeric(distal.alv.sig.score)))
distal.alv.signature.score$distal.alv.sig.score <- as.numeric(distal.alv.signature.score$distal.alv.sig.score)



## imm infiltrate signature score-------------------------------------------------------------
test <- as.matrix(imm_infil_sig_count)    ## dgCmatrix to matrix with HANDLABLE small size, 11 genes 57682 cells

z_test <- test  ## zi = xi-u/sigma transform

test1 <- apply(test, 1, mean)  # calculate each gene mean and sd, like scale by row in heatmap
test2 <- apply(test, 1, sd)

for (i in 1:nrow(z_test)){
  z_test[i,] <- (test[i,] - test1[i])/test2[i]
}

imm.infil.sig.score <- c()
# repeat in 31 cell types
for (i in 1:31) {
  cell.type.index <- which(meta.combine$celltype == tmp1[i] & meta.combine$Diagnosis == "IPF")
  cell.type.data <- z_test[, cell.type.index]
  imm.infil.sig.score[i] <- as.numeric(sum(cell.type.data)/length(cell.type.index))
}
## get imm infil.signature score in 31 cell types
imm.infil.signature.score <- as.data.frame(cbind("Number" = c(1:31), "Cell.type" = tmp1, "imm.infil.sig.score" = as.numeric(imm.infil.sig.score)))
imm.infil.signature.score$imm.infil.sig.score <- as.numeric(imm.infil.signature.score$imm.infil.sig.score)



## ipf vessel signature score-------------------------------------------------------------
test <- as.matrix(ipf_vessel_sig_count)    ## dgCmatrix to matrix with HANDLABLE small size, 11 genes 57682 cells

z_test <- test  ## zi = xi-u/sigma transform

test1 <- apply(test, 1, mean)  # calculate each gene mean and sd, like scale by row in heatmap
test2 <- apply(test, 1, sd)

for (i in 1:nrow(z_test)){
  z_test[i,] <- (test[i,] - test1[i])/test2[i]
}

ipf.vessel.sig.score <- c()
# repeat in 31 cell types
for (i in 1:31) {
  cell.type.index <- which(meta.combine$celltype == tmp1[i] & meta.combine$Diagnosis == "IPF")
  cell.type.data <- z_test[, cell.type.index]
  ipf.vessel.sig.score[i] <- as.numeric(sum(cell.type.data)/length(cell.type.index))
}
## get ipf vessel signature score in 31 cell types
ipf.vessel.signature.score <- as.data.frame(cbind("Number" = c(1:31), "Cell.type" = tmp1, "ipf.vessel.sig.score" = as.numeric(ipf.vessel.sig.score)))
ipf.vessel.signature.score$ipf.vessel.sig.score <- as.numeric(ipf.vessel.signature.score$ipf.vessel.sig.score)



## healthy alveoli signature score-------------------------------------------------------------

test <- as.matrix(healthy_alv_sig_count)    ## dgCmatrix to matrix with HANDLABLE small size, 5 genes 57682 cells

z_test <- test  ## zi = xi-u/sigma transform

test1 <- apply(test, 1, mean)  # calculate each gene mean and sd, like scale by row in heatmap
test2 <- apply(test, 1, sd)

for (i in 1:nrow(z_test)){
  z_test[i,] <- (test[i,] - test1[i])/test2[i]
}

healthy.alv.sig.score <- c()
# repeat in 30 cell types
for (i in 1:30) {
  cell.type.index <- which(meta.combine$celltype == tmp1_control[i] & meta.combine$Diagnosis == "Control")
  cell.type.data <- z_test[, cell.type.index]
  healthy.alv.sig.score[i] <- as.numeric(sum(cell.type.data)/length(cell.type.index))
}
## get healthy.alv.signature score in 30 cell types
healthy.alv.signature.score <- as.data.frame(cbind("Number" = c(1:30), "Cell.type" = tmp1_control, "healthy.alv.sig.score" = as.numeric(healthy.alv.sig.score)))
healthy.alv.signature.score$healthy.alv.sig.score <- as.numeric(healthy.alv.signature.score$healthy.alv.sig.score)

# export summary signature table
summary_signature_table <- list("ff.signature.score" = ff.signature.score,
                                "distal.alv.signature.score" = distal.alv.signature.score,
                                "imm.infil.signature.score" = imm.infil.signature.score,
                                "ipf.vessel.signature.score" = ipf.vessel.signature.score,
                                "healthy.alv.signature.score" =  healthy.alv.signature.score)


export(summary_signature_table,  "query_enrichment_score_summary_table.xlsx")



# quick plot of cell type enrichment in fibroblast foci---------
ggplot(ff.signature.score, aes(x=Cell.type, y=ff.sig.score)) + 
  geom_bar(stat = "identity")

# prepare bubble plot
query.summary <- cbind(ipf.vessel.signature.score, distal.alv.signature.score$distal.alv.sig.score,
                         ff.signature.score$ff.sig.score, imm.infil.signature.score$imm.infil.sig.score)

query.summary <- query.summary[-10,]
query.summary <- query.summary[,-1]

colnames(query.summary) <- c("cell type", "IPF blood vessel", "IPF distant alveoli", "IPF fibroblast foci",	"IPF immune infiltrate")
query.summary <- cbind(query.summary, "Control alveoli" = healthy.alv.signature.score$healthy.alv.sig.score)
query.summary <- query.summary[c(1,6,2,3,4,5)]

hm.rsg.csg <- query.summary
hm.rsg.csg <- hm.rsg.csg[-1]

hm.rsg.csg <- as.data.frame(hm.rsg.csg)
rownames(hm.rsg.csg) <- tmp1_control

# reorder based on epithelium, mesenchyme, myeloid, endothelium, and lyphoid
new_order <- c(1, 2, 4, 6, 7, 10, 16, 17, 23, 26, 27, 30, 9, 14, 18, 22, 28, 5, 
               12, 13, 15, 20, 24, 8, 11, 3, 19, 21, 25, 29)
hm.rsg.csg <- hm.rsg.csg[new_order, ]

hm.rsg.csg1 <- cbind(hm.rsg.csg, rownames(hm.rsg.csg))
tmp8 <- hm.rsg.csg1[nrow(hm.rsg.csg1):1,]  # reverse the sample order
colnames(tmp8)[6] <- "Sample"


# bubble plot
tmp1 <- hm.rsg.csg

tmp2 <- tmp1[nrow(tmp1):1,]  # reverse the sample order
tmp2 <- cbind(tmp2, "Sample" = rownames(tmp2))




#convert data frame from a "wide" format to a "long" format
pcm = melt(tmp2, id = c("Sample"))
pcm$Sample <- factor(pcm$Sample,levels=unique(pcm$Sample))

# positive signature score, for enrichment
p1 <- ggplot(pcm, aes(y = Sample, x = variable)) + 
  geom_point(aes(size = value, fill = variable), shape = 21) +
  scale_size_continuous(limits = c(0, 160), range = c(1,10), breaks = c(0,50,100,150)) +
  theme(
    legend.key = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(colour = "black", face = "bold", size = 11),
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13)
  ) +
  labs(x = "Regions", y = "Cell Types")
plot(p1)  

# negative signature score, for depletion
p2 <-  ggplot(pcm, aes(y = Sample, x = variable)) + 
  geom_point(aes(size = value, fill = variable),  shape = 21) +
  scale_size_continuous(limits = c(-30, 0), range = c(10,1), breaks = c(-30,-20,-10,0))+
  theme(
    legend.key = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(colour = "black", face = "bold", size = 11),
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13)
  ) +
  labs(x = "Regions", y = "Cell Types")
plot(p2)  


# save bubble plot
ggsave("query_positive_enrichment_bubble_plot_summary.pdf", plot = p1, width = 8.5, height = 11)  
ggsave("query_negative_depletion_bubble_plot_summary.pdf", plot = p2, width = 8.5, height = 11)  


# save data
save.image("query_method_enrichment.RData")




















