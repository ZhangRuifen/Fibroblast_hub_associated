#### 20251225 将human中的ENSEMBL转为SYMBOL ####

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(Matrix)

obj<-scRNAlist[[1]]

#### 1.读入数据------------------------------------------------------------------------------------------------------------------------------

# 1.构建数据的目标路径
dir <-  c(
  "GSM7307094","GSM7307095","GSM7307096","GSM7307097","GSM7307098",
  "GSM7307099","GSM7307100","GSM7307101","GSM7307102","GSM7307103","GSM7307104","GSM7307105")

dir <- paste0("G:/课题二 成纤维细胞相关/Results/Public data/GSE231993_data1/", dir)

sample_name <-  c(
  "GSM7307094","GSM7307095","GSM7307096","GSM7307097","GSM7307098",
  "GSM7307099","GSM7307100","GSM7307101","GSM7307102","GSM7307103",
  "GSM7307104","GSM7307105"
)
scRNAlist <- list()

# 2.读入

for(i in 1:length(dir)) {
  counts <- Read10X(data.dir = dir[i], gene.column = 1)  # 读取数据
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=sample_name[i], min.cells=3, min.features=200)
}
# 检查数据的一致性
sapply(scRNAlist, function(obj) {
  c(
    nFeature = nrow(obj),
    nCell    = ncol(obj)
  )
})

#### 2.构建转换单细胞数据counts信息的函数------------------------------------------------------------------------------------------------------------------------------
convert_ensembl_to_symbol <- function(seurat_obj) {
  
  counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  
  ensembl_ids <- rownames(counts)
  
  # 1.ENSEMBL -> SYMBOL 映射
  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = ensembl_ids,
    keytype = "ENSEMBL",
    columns = c("SYMBOL")
  )
  dim(mapping) # 37119     2
  
  # 2.去掉 NA，得到可用的SYMBOL
  mapping <- mapping[!is.na(mapping$SYMBOL), ]
  dim(mapping)  # 24838     2
  
  length(unique(mapping$ENSEMBL))  # 24320
  length(unique(mapping$SYMBOL))  # 24608
  
  # 每个 ENSEMBL 只保留第一次出现的 SYMBOL
  mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
  
  # 3.保证顺序一致;按照mapping的顺序将counts顺序重排，如果某个rownames不在ENSEMBL中，就会match后就会返回NA值
  #   通过match函数：将mapping的行顺序调整成与counts或valid_ensembl完全一致的顺序
  mapping <- mapping[match(rownames(counts), mapping$ENSEMBL), ]   # 使用ENSEMBL将mapping的顺序按照counts顺序排列为mapping
  dim(mapping)  # 36601     2
  
  # 找出mapping中SYMBOL不为NA的那些行号
  valid_idx <- which(!is.na(mapping$SYMBOL))
  
  counts <- counts[valid_idx, ]
  dim(counts) # 24320  4890
  
  symbols <- mapping$SYMBOL[valid_idx]
  length(symbols) # 24320
    
  # ！！！至此counts跟symbol的顺序是一致的
  length(symbols) # 24320
  length(unique(symbols)) # 24274
  # 这里表示存在多个ENSEMBL对应一个SYMBOL的情况！因此下一步通过稀疏矩阵将构建将重复的SYMBOl进行合并
  
  # 4.合并重复 SYMBOL（关键步骤）：
  # 生成一个：行数 = length(symbols)；列数 = length(unique(symbols))
  # 每一行只有一个 1，其余都是 0：表示“这一行属于哪个 SYMBOL”
  agg_mat <- sparse.model.matrix(~ symbols - 1)  # agg_mat 是 (genes × symbols)
  
  # model.matrix的列名规则是变量名+level名；这里将变量名“symbols”去除
  colnames(agg_mat) <- sub("^symbols", "", colnames(agg_mat))
  
  # 矩阵计算
  # 维度变化：
  # counts 是 (genes × cells)
  # agg_mat 是 (genes × symbols)
  # t(agg_mat) 是 (symbols × genes)
  # 所以乘起来：
  # t(agg_mat) %*% counts 得到 (symbols × cells)
  new_counts <- t(agg_mat) %*% counts  # %*%: 矩阵乘法运算符

  # 写回 Seurat
  seurat_obj[["RNA_symbol"]] <- CreateAssayObject(counts = new_counts)
  DefaultAssay(seurat_obj) <- "RNA_symbol"
  return(seurat_obj)
}

#### 3.批量处理并保存------------------------------------------------------------------------------------------------------------------------------
scRNAlist <- lapply(scRNAlist, convert_ensembl_to_symbol)

# 1️⃣ 准备输出目录
out_dir <- "G:/课题二 成纤维细胞相关/Results/Public data/Result/1.ENSEMBL_to_SYMBOL"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 检查数据的一致性
sapply(scRNAlist, function(obj) {
  c(
    nFeature = nrow(obj),
    nCell    = ncol(obj)
  )
})

# 对不同样本的features取交集
# 取所有样本 RNA_symbol 的基因交集
common_genes <- Reduce(
  intersect,
  lapply(scRNAlist, function(obj) {
    rownames(GetAssayData(obj, assay = "RNA_symbol", layer = "counts"))
  })
)
length(common_genes)  # 13468

# 对每个对象 subset 到共同基因集
scRNAlist <- lapply(scRNAlist, function(obj) {
  # 取 RNA_symbol 的 counts
  mat <- GetAssayData(obj, assay = "RNA_symbol", layer = "counts")
  
  # 只保留 common_genes
  mat2 <- mat[common_genes, , drop = FALSE]
  
  # 重新写回 RNA_symbol（只保留 counts；后面你反正会 NormalizeData）
  obj[["RNA_symbol"]] <- CreateAssayObject(counts = mat2)
  DefaultAssay(obj) <- "RNA_symbol"
  obj
})

# 检查统一后是否一致
sapply(scRNAlist, function(obj) c(nFeature = nrow(obj), nCell = ncol(obj)))


# 2️⃣ 依次保存每个 Seurat 对象
for (i in seq_along(scRNAlist)) {
  saveRDS(
    scRNAlist[[i]],
    file = file.path(
      out_dir,
      paste0(sample_name[i], ".rds")
    )
  )
}


