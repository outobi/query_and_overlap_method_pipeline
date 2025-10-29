setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(dplyr)
library(rstatix)
library(readxl)
library(read)
library(rio)

library(ggplot2)
library(reshape2)

library(pheatmap)
library(RColorBrewer)

load("query method enrichment.RData")

## get gene features from 31 cell types-----------------------------------------

summary <- data.frame("index" = 1:24470)    # store p value


for (j in 1:31) {
stat.sig <- c()
d <- c()
cell.type <- tmp1[j]


example1 <- normalized_count[1,]
experiment.index <- which(meta.combine$celltype == cell.type & meta.combine$Diagnosis == "IPF")
reference.index <- which(meta.combine$celltype != cell.type  & meta.combine$Diagnosis == "IPF")

ttest1 <- as.numeric(example1[reference.index])
ttest2 <-as.numeric(example1[experiment.index])

ntest1 <- length(ttest1)
ntest2 <- length(ttest2)
group1 <- rep("control", ntest1)
group2 <- rep("experiment", ntest2)
group <- c(group1, group2)

# differential gene expression in target cell type tmp1[j] versus all non target cells in IPF patients
ttest <- as.numeric(c(ttest1, ttest2))
test.table <- as.data.frame(cbind("expression" = ttest, "group" = group))
test.table$expression <- as.numeric(test.table$expression)

stat.test2 <- test.table %>% 
  t_test(expression ~ group, var.equal = TRUE) %>% ### var.equal = TRUE, student's t test, same variance
  add_significance()
stat.sig[1] <- stat.test2$p

# Cohen's d is a measure of effect size, specifically for t-tests, indicating the standardized difference between two means
d.cohen <- test.table %>%  cohens_d(expression ~ group, var.equal = TRUE)
d[1] <- d.cohen$effsize


#    may take long time to process depending on the dataset size
for (i in 2:24470) {
  example1 <- normalized_count[i,]
  ttest1 <- as.numeric(example1[reference.index])
  ttest2 <-as.numeric(example1[experiment.index])
  test.table$expression <- as.numeric(c(ttest1, ttest2))
  stat.test2 <- test.table %>% 
    t_test(expression ~ group, var.equal = TRUE) %>%
    add_significance()
  if (is.nan(stat.test2$p)) {next}
  if (stat.test2$p < 0.00001) { # only record p < 10e-5, then calculate cohen's d values
    stat.sig[i] <- stat.test2$p
    d.cohen <- test.table %>%  cohens_d(expression ~ group, var.equal = TRUE)
    d[i] <- d.cohen$effsize
  }
}

d.length <- length(d)
name <- paste0("Cell type ", j, " ", tmp1[j])
statis <- as.data.frame(list("index" = c(1:d.length), "d" = d, "stat.sig" = stat.sig))
export(statis, paste0(name, ".xlsx"))

summary <- cbind(summary, statis$d)

}
colnames(summary)[2:32] <- tmp1
export(summary, "summary of d value 31 cell types in IPF.xlsx")
#summary <- read_excel("summary of d value 31 cell types in IPF.xlsx")

# find cell type specific genes with maximum cohen's d value to represent this cell type------------
summary <- summary[-1]
summary[is.na(summary)] <- 0

summary1 <- (-1)*summary   ## d value is inverse in the previous analysis, now it is switch back to experiment - control

summary1$Largest_Column<-colnames(summary1)[apply(summary1,1,which.max)] # cell type with max d value get this cell type specific gene

for (i in 1:24470) {
  if (sd(summary1[i,c(1:31)]) == 0) {summary1$Largest_Column[i] <- NA} # if all cell type d value is 0, sd is 0, then this gene does not belong to any cell type

}

# store cell type specific gene index in each cell type
tmp2 <- list()
tmp4 <- c()
for (i in 1:31){
tmp3 <- summary1[which(summary1$Largest_Column == tmp1[i]),]
tmp2[[i]] <- as.numeric(rownames(tmp3))
tmp4[i] <- nrow(tmp3)
}
sum(tmp4) # 22523
tmp4  # cell type specific gene number in 31 cell types
#[1] 1058  536  385  569  282 5135  300  860  481 1120 1016  594  675  581  452  842  457
#[18]  104  607  708  986  300  576  261  559 1090  211   67  601  727  383



# export protein list from tmp2
names <- c()
for (i in 1:31){
a1 <- as.data.frame(tmp2[[i]])
names <- paste0("IPF cell type specific gene index", i)
export(a1, paste0(names, ".xlsx"))
}


names(tmp2) <- tmp1# cell type names to tmp2 list names

csg <- tmp2
# rsg region specific genes
# csg cell type specific genes

## fibroblast foci rsg overlap with csg---------------------------------------

tmp5 <- ff_signature$index   # ff sig
q <-c()
expected <-c()
direction <- c()
p.value <- c()

for (i in 1:31){
tmp6 <- csg[[i]]
q[i] <- length(intersect(tmp5, tmp6))
if (is.null(q[i])) {q[i] <- 0}
k = length(tmp5)
m = length(tmp6)
n = 24470 - length(tmp6)

expected[i] <- k*m/n
# hypergeometric test
# if q intersect value is larger than expectation, enrichment
# if q intersect value is smaller than expectation, depletion

if (q[i] > expected[i]) {direction[i] <- "enriched"
                        p.value[i] <- phyper(q[i]-1, m, n, k, lower.tail= FALSE)}           

if (q[i] < expected[i]) {direction[i] <- "depleted"
                        p.value[i] <- phyper(q[i], m, n, k, lower.tail= TRUE)}

}

ff.overlap.score <- data.frame("cell.type" = tmp1, "direction" = direction, "p.value" = p.value, "-log10(p.value)" = -1*log10(p.value))

# direction, if enrich, positive; if depletion, negative
tmp7 <- ff.overlap.score$X.log10.p.value.
for (i in 1:31){
if (ff.overlap.score$direction[i] == "depleted") {tmp7[i] <- -1*tmp7[i]}
}
ff.overlap.score <- cbind(ff.overlap.score,"p.value.with.direction" = tmp7)
export(ff.overlap.score,"ff.overlap.score.xlsx")


##  distal alveoli-----------------------
tmp5 <- distal_alv_signature$index   # distal alveoli sig
q <-c()
expected <-c()
direction <- c()
p.value <- c()

for (i in 1:31){
  tmp6 <- csg[[i]]
  q[i] <- length(intersect(tmp5, tmp6))
  if (is.null(q[i])) {q[i] <- 0}
  k = length(tmp5)
  m = length(tmp6)
  n = 24470 - length(tmp6)
  
  expected[i] <- k*m/n
  
  if (q[i] > expected[i]) {direction[i] <- "enriched"
  p.value[i] <- phyper(q[i]-1, m, n, k, lower.tail= FALSE)}           
  
  if (q[i] < expected[i]) {direction[i] <- "depleted"
  p.value[i] <- phyper(q[i], m, n, k, lower.tail= TRUE)}
  
}

distal.alv.overlap.score <- data.frame("cell.type" = tmp1, "direction" = direction, "p.value" = p.value, "-log10(p.value)" = -1*log10(p.value))

tmp7 <- distal.alv.overlap.score$X.log10.p.value.
for (i in 1:31){
  if (distal.alv.overlap.score$direction[i] == "depleted") {tmp7[i] <- -1*tmp7[i]}
}
distal.alv.overlap.score <- cbind(distal.alv.overlap.score,"p.value.with.direction" = tmp7)
export(distal.alv.overlap.score,"distal.alv.overlap.score.xlsx")

##  ipf vessel----------------------------
tmp5 <- ipf_vessel_signature$index   # IPF vessel sig
q <-c()
expected <-c()
direction <- c()
p.value <- c()

for (i in 1:31){
  tmp6 <- csg[[i]]
  q[i] <- length(intersect(tmp5, tmp6))
  if (is.null(q[i])) {q[i] <- 0}
  k = length(tmp5)
  m = length(tmp6)
  n = 24470 - length(tmp6)
  
  expected[i] <- k*m/n
  
  if (q[i] > expected[i]) {direction[i] <- "enriched"
  p.value[i] <- phyper(q[i]-1, m, n, k, lower.tail= FALSE)}           
  
  if (q[i] < expected[i]) {direction[i] <- "depleted"
  p.value[i] <- phyper(q[i], m, n, k, lower.tail= TRUE)}
  
}

ipf.vessel.overlap.score <- data.frame("cell.type" = tmp1, "direction" = direction, "p.value" = p.value, "-log10(p.value)" = -1*log10(p.value))

tmp7 <- ipf.vessel.overlap.score$X.log10.p.value.
for (i in 1:31){
  if (ipf.vessel.overlap.score$direction[i] == "depleted") {tmp7[i] <- -1*tmp7[i]}
}
ipf.vessel.overlap.score <- cbind(ipf.vessel.overlap.score, "p.value.with.direction" = tmp7)
export(ipf.vessel.overlap.score,"ipf.vessel.overlap.score.xlsx")


##  imm infiltrate----------------------------
tmp5 <- imm_infil$index   # imm infiltrate sig
q <-c()
expected <-c()
direction <- c()
p.value <- c()

for (i in 1:31){
  tmp6 <- csg[[i]]
  q[i] <- length(intersect(tmp5, tmp6))
  if (is.null(q[i])) {q[i] <- 0}
  k = length(tmp5)
  m = length(tmp6)
  n = 24470 - length(tmp6)
  
  expected[i] <- k*m/n
  
  if (q[i] > expected[i]) {direction[i] <- "enriched"
  p.value[i] <- phyper(q[i]-1, m, n, k, lower.tail= FALSE)}           
  
  if (q[i] < expected[i]) {direction[i] <- "depleted"
  p.value[i] <- phyper(q[i], m, n, k, lower.tail= TRUE)}
  
}

imm.infil.overlap.score <- data.frame("cell.type" = tmp1, "direction" = direction, "p.value" = p.value, "-log10(p.value)" = -1*log10(p.value))

tmp7 <- imm.infil.overlap.score$X.log10.p.value.
for (i in 1:31){
  if (imm.infil.overlap.score$direction[i] == "depleted") {tmp7[i] <- -1*tmp7[i]}
}
imm.infil.overlap.score <- cbind(imm.infil.overlap.score, "p.value.with.direction" = tmp7)
export(imm.infil.overlap.score,"imm.infil.overlap.score.xlsx")


## control donor's cell type feature extraction FROM 31644 SINGLE CELLS---------
summary2 <- data.frame("index" = 1:24470)    # store p value

for (j in 1:31) {
  if (j ==10) {next}  #no Cell type 10 HAS1 High fibroblast cell type in control 
  stat.sig <- c()
  d <- c()
  cell.type <- tmp1[j]
  
  
  example1 <- normalized_count[1,]
  experiment.index <- which(meta.combine$celltype == cell.type & meta.combine$Diagnosis == "Control")
  reference.index <- which(meta.combine$celltype != cell.type  & meta.combine$Diagnosis == "Control")
  
  ttest1 <- as.numeric(example1[reference.index])
  ttest2 <-as.numeric(example1[experiment.index])
  
  ntest1 <- length(ttest1)
  ntest2 <- length(ttest2)
  group1 <- rep("control", ntest1)
  group2 <- rep("experiment", ntest2)
  group <- c(group1, group2)
  
  
  ttest <- as.numeric(c(ttest1, ttest2))
  test.table <- as.data.frame(cbind("expression" = ttest, "group" = group))
  test.table$expression <- as.numeric(test.table$expression)
  
  stat.test2 <- test.table %>% 
    t_test(expression ~ group, var.equal = TRUE) %>%
    add_significance()
  stat.sig[1] <- stat.test2$p
  

  d.cohen <- test.table %>%  cohens_d(expression ~ group, var.equal = TRUE)
  d[1] <- d.cohen$effsize
  
  
  
  for (i in 2:24470) {
    example1 <- normalized_count[i,]
    ttest1 <- as.numeric(example1[reference.index])
    ttest2 <-as.numeric(example1[experiment.index])
    test.table$expression <- as.numeric(c(ttest1, ttest2))
    stat.test2 <- test.table %>% 
      t_test(expression ~ group, var.equal = TRUE) %>%
      add_significance()
    if (is.nan(stat.test2$p)) {next}
    if (stat.test2$p < 0.00001) { 
      stat.sig[i] <- stat.test2$p
      d.cohen <- test.table %>%  cohens_d(expression ~ group, var.equal = TRUE)
      d[i] <- d.cohen$effsize
    }
  }
  
  d.length <- length(d)
  name <- paste0("Control Cell type ", j, " ", tmp1[j])
  statis <- as.data.frame(list("index" = c(1:d.length), "d" = d, "stat.sig" = stat.sig))
  export(statis, paste0(name, ".xlsx"))
  
  summary2 <- cbind(summary, statis$d)
}

colnames(summary2)[2:31] <- tmp1_control
export(summary2, "summary of d value 30 cell types in control.xlsx")
#summary2 <- read_excel("summary of d value 30 cell types in control.xlsx")

# find cell type specific genes with maximum cohen's d value to represent this cell type------------

summary2 <- summary2[-1]
summary2[is.na(summary2)] <- 0
summary3 <- -1*summary2   ## d value is inverse in the previous analysis, should be experiment - control

summary3$Largest_Column<-colnames(summary3)[apply(summary3,1,which.max)]

for (i in 1:24470) {
  if (sd(summary3[i,c(1:30)]) == 0) {summary3$Largest_Column[i] <- NA}
  
}

tmp2 <- list()
tmp4 <- c()
for (i in 1:30){
  tmp3 <- summary3[which(summary3$Largest_Column == tmp1_control[i]),]
  tmp2[[i]] <- as.numeric(rownames(tmp3))
  tmp4[i] <- nrow(tmp3)
}


sum(tmp4) # 22077
tmp4
#[1] 1073  949  400  497  329 4216  858 1015  687  379  931  902  805  552  850  385  118  506  456  526  624  314
#[23]  617  319 1733  318   95  652  702  269
# store to "csg number index.xlsx

names <- c()
for (i in 1:30){
  a1 <- as.data.frame(tmp2[[i]])
  names <- paste0("Control cell type specific gene index", i)
  export(a1, paste0(names, ".xlsx"))
}


names(tmp2) <- tmp1_control

csg <- tmp2


## control alveoli rsg overlap with control alveoli csg-------------------------------
tmp5 <- healthy_alv_signature$index   # control alveoli signature
q <-c()
expected <-c()
direction <- c()
p.value <- c()

for (i in 1:30){
  tmp6 <- csg[[i]]
  q[i] <- length(intersect(tmp5, tmp6))
  if (is.null(q[i])) {q[i] <- 0}
  k = length(tmp5)
  m = length(tmp6)
  n = 24470 - length(tmp6)
  
  expected[i] <- k*m/n
  
  if (q[i] > expected[i]) {direction[i] <- "enriched"
  p.value[i] <- phyper(q[i]-1, m, n, k, lower.tail= FALSE)}           
  
  if (q[i] < expected[i]) {direction[i] <- "depleted"
  p.value[i] <- phyper(q[i], m, n, k, lower.tail= TRUE)}
  
}

healthy.alv.overlap.score <- data.frame("cell.type" = tmp1_control, "direction" = direction, "p.value" = p.value, "-log10(p.value)" = -1*log10(p.value))

tmp7 <- healthy.alv.overlap.score$X.log10.p.value.
for (i in 1:30){
  if (healthy.alv.overlap.score$direction[i] == "depleted") {tmp7[i] <- -1*tmp7[i]}
}
healthy.alv.overlap.score <- cbind(healthy.alv.overlap.score,"p.value.with.direction" = tmp7)
export(healthy.alv.overlap.score,"healthy.alv.overlap.score.xlsx")



overlap.summary <- cbind(ipf.vessel.overlap.score, distal.alv.overlap.score$p.value.with.direction,
                         ff.overlap.score$p.value.with.direction, imm.infil.overlap.score$p.value.with.direction)
                       
overlap.summary <- overlap.summary[-10,]
overlap.summary <- overlap.summary[,-c(2:4)]

colnames(overlap.summary) <- c("cell type", "IPF blood vessel", "IPF distant alveoli", "IPF fibroblast foci",	"IPF immune infiltrate")
overlap.summary <- cbind(overlap.summary, "Control alveoli" = healthy.alv.overlap.score$p.value.with.direction)
overlap.summary <- overlap.summary[c(1,6,2,3,4,5)]


# make heatmap and bubble plot--------

# heatmap

hm.rsg.csg <- overlap.summary
hm.rsg.csg <- hm.rsg.csg[-1]

hm.rsg.csg <- as.data.frame(hm.rsg.csg)
rownames(hm.rsg.csg) <- tmp1_control

# reorder based on epithelium, mesenchyme, myeloid, endothelium, and lyphoid
new_order <- c(1, 2, 4, 6, 7, 10, 16, 17, 23, 26, 27, 30, 9, 14, 18, 22, 28, 5, 
         12, 13, 15, 20, 24, 8, 11, 3, 19, 21, 25, 29)
hm.rsg.csg <- hm.rsg.csg[new_order, ]



palette_reds <- colorRampPalette(colors = c("red","white"))(250)
scales::show_col(palette_reds[1:250])
palette1 <- palette_reds[1:250]

palette_blues <- colorRampPalette(colors = c("white", "blue"))(250)
scales::show_col(palette_blues[1:90])
palette2 <- palette_blues[1:90]

palette3 <- c(palette1, palette2)



p1 <- pheatmap(
  hm.rsg.csg,
  rev(palette3), 
  #scale = c('row'),
  gaps_col=NULL,fontsize_row=8,
  #border_color = "black",
  #redgreen(75),scale = c('row'),
  show_rownames =T, show_colnames = T,cluster_cols=F, cluster_rows = F,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation", border_color = 'Black',
  #annotation = hm_anno,
  #cutree_rows = 2,
  main = ""
)

ggsave("overlap enrichment heatmap summary.pdf", plot = p1, width = 8.5, height = 11)  

# bubble plot
hm.rsg.csg1 <- cbind(hm.rsg.csg, rownames(hm.rsg.csg))
tmp8 <- hm.rsg.csg1[nrow(hm.rsg.csg1):1,]  # reverse the sample order
colnames(tmp8)[6] <- "Sample"

#convert data frame from a "wide" format to a "long" format
pcm <- melt(tmp8, id = c("Sample"))
pcm$Sample <- factor(pcm$Sample,levels=unique(pcm$Sample))

# enrichment bubble figure
p2 <- ggplot(pcm, aes(y = Sample, x = variable)) + 
  geom_point(aes(size = value, fill = variable), shape = 21) +
  scale_size_continuous(limits = c(0, 25),
                        range = c(1,10),
                        breaks = c(0,5,10,15,20,25)) +
  theme(
    legend.key = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(colour = "black", face = "bold", size = 11),
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13)
  ) +
  labs(x = "Regions", y = "Cell Types")

plot(p2) 



# depletion bubble figure
p3 <- ggplot(pcm, aes(y = Sample, x = variable)) + 
  geom_point(aes(size = value, fill = variable),  shape = 21) +
  scale_size_continuous(limits = c(-10, 0), range = c(10,1), breaks = c(-10,-5, 0))+
  theme(
    legend.key = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(colour = "black", face = "bold", size = 11),
    axis.title.x = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 13)
  ) +
  labs(x = "Regions", y = "Cell Types")


plot(p3)


ggsave("overlap positive enrichment bubble plot summary.pdf", plot = p2, width = 8.5, height = 11)  
ggsave("overlap negative depletion bubble plot summary.pdf", plot = p3, width = 8.5, height = 11)  

# save data
save.image("overlap method to calculate overlap FDR.RData")


















