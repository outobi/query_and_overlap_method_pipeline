setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Read in data ------------------------------------------------------------

# expression matrix Supplementary Table 2
#Background counts from negative probe, remove this gene, we get 1085 genes in 60 samples. 
# supplementary table 2 sample column align with Supplementary table 1 sample metadata row order 

library(readxl)
library(rio)
library(dplyr)
norm.new <- read_excel("Supplementary Table 2.xlsx")
norm.new <- as.data.frame(norm.new)
rownames(norm.new) <- norm.new[[1]]
norm.new <- norm.new[-1]
colnames(norm.new) <- paste0("S", 1:60)
norm.new <- as.matrix(norm.new)

# sample metadata
meta.new <- read.csv("Supplementary Table 1.csv", stringsAsFactors = F)


## region vs all other regions based on Wilcox rank test
wilcox.comp <- function(group1) {
  wilcox.p <- c()
  wilcox.ff.c <- for(i in 1:nrow(norm.new)) {
    x <- wilcox.test(x = norm.new[i,which(meta.new$Region == group1)], 
                     y = norm.new[i,which(meta.new$Region != group1)])
    p <- x$p.value
    wilcox.p[i] <- p
  }
  wilcox.degs <- data.frame("gene.symbol" = rownames(norm.new), wilcox.p)
  wilcox.degs$p.adj <- p.adjust(p = wilcox.degs$wilcox.p, method = "BH")
  
  wilcox.logFC <- c()
  for(i in 1:nrow(norm.new)) {
    x = norm.new[i,which(meta.new$Region == group1)] %>% mean()
    y = norm.new[i,which(meta.new$Region != group1)] %>% mean()
    p <- log2(x/y)
    wilcox.logFC[i] <- p
  }
  
  wilcox.degs$LogFC <- wilcox.logFC
  return(wilcox.degs)
}

FF.v.all.wilcox <- wilcox.comp("IPF_Fibroblastic_Foci")
Imm.v.all.wilcox <- wilcox.comp("IPF_Immune")
Adj.Alv.v.all.wilcox <- wilcox.comp("IPF_Adjacent_Alveolar")
Dis.Alv.v.all.wilcon <- wilcox.comp(("IPF_Distal_Alveolar"))
IPF.vessel.all.wilcox <- wilcox.comp("IPF_vessel")
Heal.vessel.all.wilcox <- wilcox.comp("Healthy_Vessel")
Heal.Alv.v.all.wilcox <- wilcox.comp("Healthy_Alveolar")

# save output lists
region_specific_feature <- list("FF.v.all.wilcox" = FF.v.all.wilcox,
                                "Imm.v.all.wilcox" = Imm.v.all.wilcox,
                                "Adj.Alv.v.all.wilcox" = Adj.Alv.v.all.wilcox,
                                "Dis.Alv.v.all.wilcon" = Dis.Alv.v.all.wilcon,
                                "IPF.vessel.all.wilcox" = IPF.vessel.all.wilcox,
                                "Heal.vessel.all.wilcox" = Heal.vessel.all.wilcox,
                                "Heal.Alv.v.all.wilcox" = Heal.Alv.v.all.wilcox
                                )
export(region_specific_feature, "region specific gene summary.xlsx")


# Log2FC > 0.25, adjusted p value < 0.05, cutoff to extract up-regualted region-specific gene list
FF.RSG <- FF.v.all.wilcox[(FF.v.all.wilcox$LogFC > 0.25 & FF.v.all.wilcox$p.adj < 0.05), ]   # n = 50
Imm.RSG <- Imm.v.all.wilcox[(Imm.v.all.wilcox$LogFC > 0.25 & Imm.v.all.wilcox$p.adj < 0.05), ]   # n = 111
Adj.Alv.RSG <- Adj.Alv.v.all.wilcox[(Adj.Alv.v.all.wilcox$LogFC > 0.25 & Adj.Alv.v.all.wilcox$p.adj < 0.05), ] # n = 0, no more following analysis
Dis.Alv.RSG <- Dis.Alv.v.all.wilcon[(Dis.Alv.v.all.wilcon$LogFC > 0.25 & Dis.Alv.v.all.wilcon$p.adj < 0.05), ] # n = 41
IPF.vessel.RSG <- IPF.vessel.all.wilcox[(IPF.vessel.all.wilcox$LogFC > 0.25 & IPF.vessel.all.wilcox$p.adj < 0.05), ] # n = 28
Heal.vessel.RSG <- Heal.vessel.all.wilcox[(Heal.vessel.all.wilcox$LogFC > 0.25 & Heal.vessel.all.wilcox$p.adj < 0.05), ] # n = 0
Heal.Alv.RSG <- Heal.Alv.v.all.wilcox[(Heal.Alv.v.all.wilcox$LogFC > 0.25 & Heal.Alv.v.all.wilcox$p.adj < 0.05), ] # n = 32, no more following analysis

# 5 region specific genes
region_specific_feature_list <- list("FF.v.all.wilcox" = FF.RSG,
                                "Imm.v.all.wilcox" = Imm.RSG,
                                "Dis.Alv.v.all.wilcox" = Dis.Alv.RSG,
                                "IPF.vessel.all.wilcox" = IPF.vessel.RSG,
                                "Heal.Alv.v.all.wilcox" = Heal.Alv.RSG
)

# save data
saveRDS(region_specific_feature_list, file = "region_specific_feature_list.rds")

save.image("data summary region vs other region.RData")


