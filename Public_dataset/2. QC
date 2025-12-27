#### 20251226 质控与标准化操作 ####

library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(patchwork)
library(RColorBrewer)
library(harmony)
library(remotes) 
library(ggplot2)  # 3.5.2
install_version("ggplot2", version = "3.5.2", repos = "https://cran.r-project.org")
install_version("ggplot2", version = "4.0.1", repos = "https://cran.r-project.org")

#### 1.对每个seurat对象的features进行统一------------------------------------------------------------------------------------------------------------------------------
# 0️⃣ 读入单细胞数据
setwd('G:/课题二 成纤维细胞相关/Results/Public data/Result/2.QC')
dir <-  c("GSM7307094","GSM7307095","GSM7307096","GSM7307097","GSM7307098",
  "GSM7307099","GSM7307100","GSM7307101","GSM7307102","GSM7307103","GSM7307104","GSM7307105")
dir <- paste0("G:/课题二 成纤维细胞相关/Results/Public data/GSE231993_data1/", dir)
sample_name <-  c(
  "GSM7307094","GSM7307095","GSM7307096","GSM7307097","GSM7307098",
  "GSM7307099","GSM7307100","GSM7307101","GSM7307102","GSM7307103",
  "GSM7307104","GSM7307105"
)
scRNAlist <- list()

setwd('G:/课题二 成纤维细胞相关/Results/Public data/Result/2.QC')
pdf("QC_violin_all_samples.pdf", width = 7, height = 5)
for(i in 1:length(dir)) {
  counts <- Read10X(data.dir = dir[i], gene.column = 2)  # 读取数据
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=sample_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-") 
  # p <- VlnPlot(scRNAlist[[i]], features = c("nCount_RNA", "nFeature_RNA","percent.mt"), ncol = 3)
  # print(p)
}
dev.off()

# 0️⃣ 指定你要操作的 assay（推荐 SYMBOL）
assay_use <- "RNA"
layer_use <- "counts"

# 1️⃣ 计算所有对象的 feature 交集
common_genes <- Reduce(
  intersect,
  lapply(scRNAlist, function(obj) {
    rownames(GetAssayData(obj, assay = assay_use, layer = layer_use))
  })
)
length(common_genes)  # 看交集大小, 14314
dir <-  c(
  "GSM7307094","GSM7307095","GSM7307096","GSM7307097","GSM7307098",
  "GSM7307099","GSM7307100","GSM7307101","GSM7307102","GSM7307103","GSM7307104","GSM7307105")


# 2️⃣ 核心步骤：构建“新的 Seurat 对象”（不修改原对象）
scRNAlist_common <- lapply(scRNAlist, function(obj) {
  
  # 取 counts 并裁剪到交集基因
  mat <- GetAssayData(obj, assay = assay_use, layer = layer_use)
  mat <- mat[common_genes, , drop = FALSE]
  
  # 用裁剪后的 counts 新建 Seurat 对象
  new_obj <- CreateSeuratObject(
    counts = mat,
    project = dir[i]
  )
  
  # 拷贝 meta.data（按 cell 名对齐）
  new_obj@meta.data <- obj@meta.data[colnames(mat), , drop = FALSE]
  
  # 拷贝细胞身份（如果你之前设置过）
  Idents(new_obj) <- Idents(obj)
  
  new_obj
})

sapply(scRNAlist_common, function(obj) {
  c(
    nFeature = nrow(obj),
    nCell    = ncol(obj)
  )
})

#### 2.QC------------------------------------------------------------------------------------------------------------------------------
sceList<-list()
sceList = lapply(scRNAlist_common, function(x) {
  subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
           nCount_RNA > 400 & nCount_RNA < 40000 & percent.mt < 15)
})

sapply(sceList, function(obj) {
  c(
    nFeature = nrow(obj),
    nCell    = ncol(obj)
  )
})

pdf("after.QC_violin_all_samples.pdf", width = 7, height = 5)
for(i in 1:length(dir)) {
  p <- VlnPlot(sceList[[i]], features = c("nCount_RNA", "nFeature_RNA","percent.mt"), ncol = 3)
  print(p)
}
dev.off()

# 为每个细胞名加后缀
for (i in 1:length(sceList)) {
  suffix <- substr(sample_name[i], nchar(sample_name[i]) - 2, nchar(sample_name[i]))  # substr(x, start, end)：截取字符串子串
  new_cellnames <- paste0(colnames(sceList[[i]]), "_", suffix)  # 修改细胞名：原始细胞名 + "_后缀"
  sceList[[i]] <- RenameCells(sceList[[i]], new.names = new_cellnames)
}

scRNA <- merge(sceList[[1]], y = c(sceList[2:12])) #合并所有样本
scRNA <- JoinLayers(scRNA)
saveRDS(scRNA,"merged.scRNA.RDS")


