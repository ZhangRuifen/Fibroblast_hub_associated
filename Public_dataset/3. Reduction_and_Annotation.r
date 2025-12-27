#### 20251226 对质控后的样本进行降维聚类注释

library(clustree)
library(dplyr)
library(Seurat)
library(patchwork)


setwd('G:/课题二 成纤维细胞相关/Results/Public data/Result/3.降维聚类')

#### 1.降维聚类------------------------------------------------------------------------------------------------------------------------------

DefaultAssay(scRNA)<-"RNA"
scRNA[["RNA"]] <- split(scRNA[["RNA"]], f = scRNA$orig.ident)
Layers(scRNA[["RNA"]])

scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA)
scRNA <- ScaleData(scRNA)
scRNA <- RunPCA(scRNA)
ElbowPlot(scRNA, ndims = 50)
### 筛选贡献度最大的PC数
pct <- scRNA[["pca"]]@stdev / sum(scRNA[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
component1 <- which(cumu > 90 & pct < 5)[1] # determine the point where the principal component contributes < 5% of standard deviation and the principal components so far have cumulatively contributed 90% of the standard deviation.
component2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1 # identify where the percent change in variation between consecutive PCs is less than 0.1%
prin_comp <- min(component1, component2) 
prin_comp

colnames(scRNA@meta.data)
scRNA <- IntegrateLayers(
  object = scRNA, 
  method = HarmonyIntegration,
  orig.reduction = "pca", 
  new.reduction = "integrated.harmony",
  group.by.vars = "orig.ident", 
  verbose = FALSE
)

#Dimensional reduction and Cell cluster identification
scRNA[["RNA"]] <- JoinLayers(scRNA[["RNA"]])

scRNA <- FindNeighbors(scRNA,  reduction = "integrated.harmony", dims = 1:prin_comp)
scRNA <- FindClusters(scRNA, resolution = c(0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7))  
scRNA <- RunUMAP(scRNA,  reduction = "integrated.harmony", dims = 1:prin_comp)
c<-clustree(scRNA)
c
ggsave(c,file= 'clustree.tree.pdf',width = 7,height = 8)

a<-DimPlot(scRNA, reduction = "umap", group.by = c("RNA_snn_res.0.05","RNA_snn_res.0.1", "RNA_snn_res.0.15", "RNA_snn_res.0.2", 
                                                  "RNA_snn_res.0.3","RNA_snn_res.0.4","RNA_snn_res.0.5"))
ggsave(a,file= 'umap.pdf',width = 15,height = 12)


#### 2.注释------------------------------------------------------------------------------------------------------------------------------
# 高变基因进行定义
colnames(scRNA@meta.data)

list_genes <- list(
  # Stem_cells = c('LGR5','ASCL2','SMOC2'),
  # TA = c('MKI67','UBE2C','TOP2A','PCNA'),
  # Mature_colonocyte = c('AQP8','FABP2'),
  # Goblets = c('TFF3','MUC2','SPDEF','WFDC2'),
  # EEC  = c('CHGA','CHGB','NEUROD1'), 
  # NK_cells = c('GNLY'),
  # Neutrophils = c('CSF3R','CXCR2','S100A8', 'S100A9', 'LCN2', 'FCGR3B'),
  # Myeloid_cells = c("LYZ","APOE","LAMP3","C1QB",'ITGAM')
  
  Epithelial = c('EPCAM','KRT8','KRT18','KRT19','TFF3'),  
  Fibroblasts = c('COL1A1','COL1A2','ADAMDEC1','COL3A1','LUM','SPARC','DCN','SPP1'),
  Endothelial = c('CDH5', 'CLDN5', 'PECAM1'),

  B_cells = c('CD79A','CD79B','CD19'),
  Plasma = c('MZB1','JCHAIN'),   
  T_cells = c('CD3E','CD3G','CD3D','NKG7'),
  
  Mast_cells = c('CPA3','MS4A2','TPSAB1','TPSB2'),
  Macrophage = c('CD68','CD14','CTSD','LYZ')
  )

all_genes <- unlist(list_genes)
dup_genes <- all_genes[duplicated(all_genes)]

DefaultAssay(scRNA)<-'RNA'
colnames(scRNA@meta.data)
Idents(scRNA)<-scRNA$RNA_snn_res.0.05

p<-DotPlot(scRNA,features = list_genes)+
   theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10)
  )
p
ggsave(plot=p,file="final.0.05.dotplot.pdf",width = 15,height = 3)


### 注释结果
# 1.分辨率0.05：9个cluster
0:T_cells
1:Plasma
2:B_cells
3:Epithelial
4:Fibroblast
5:Macrophage
6:Endothelial
7:Mast_cells
8:Fibroblast_SPP1

### 强制按照顺序绘制 ###
ord <- as.character(0:8)

p <- DotPlot(
  scRNA,
  features = list_genes,
  cols = c("grey", "#9d0208"),
  scale = TRUE,
  cluster.idents = TRUE
)

# 强制 y 轴顺序
p$data$id <- factor(p$data$id, levels = ord)
dotplot <- p +
  scale_y_discrete(limits = ord) +
  RotatedAxis() +
  theme(
    panel.border = element_rect(color="black"),
    panel.spacing = unit(1, "mm"),
    strip.text = element_text(margin=margin(b=3, unit="mm")),
    strip.placement = "outlet",
    axis.line = element_blank()
  ) +
  labs(x="", y="")
ggsave(plot=dotplot,file="final_0.05.dotplot_红色.pdf",width = 13,height = 3.5)


#### 3.定义细胞类型+可视化------------------------------------------------------------------------------------------------------------------------------
# 1.建立cluster→celltype的映射表
cluster2type <- c(
  "0" = "T_cells",
  "1" = "Plasma",
  "2" = "B_cells",
  "3" = "Epithelial",
  "4" = "Fibroblast",
  "5" = "Macrophage",
  "6" = "Endothelial",
  "7" = "Mast_cells",
  "8" = "Fibroblast_SPP1"
)

# 2.写入 Seurat 并设置身份
# 确保 cluster ID 来自当前 Idents（你的情况就是 0~8）
scRNA$cluster_id <- as.character(Idents(scRNA))
colnames(scRNA@meta.data)

# 写入细胞类型:将聚类结果映射为细胞类型注释
scRNA$celltype <- unname(cluster2type[scRNA$cluster_id]) # 从Seurat对象中提取聚类标签，并从cluster2type中查找对应的细胞类型
# 可选：检查是否有 NA（说明映射漏了某个 cluster）
table(scRNA$celltype, useNA = "ifany")
# 把身份切换为 celltype（后续画图、DE 都会按 celltype 分组）
Idents(scRNA) <- factor(scRNA$celltype, levels = unique(cluster2type))
DimPlot(scRNA, group.by = "celltype", label = TRUE)

# 1️⃣ 明确 celltype 的顺序（这是核心）
celltype_levels <- c(
  "Epithelial",
  "T_cells",
  "B_cells",
  "Plasma",
  "Macrophage",
  "Mast_cells",
  "Fibroblast",
  "Fibroblast_SPP1",
  "Endothelial"
)
scRNA$celltype <- factor(scRNA$celltype, levels = celltype_levels)
Idents(scRNA) <- scRNA$celltype

# 2️⃣ 定义颜色（顺序必须一一对应）
celltype_colors <- c(
  "#1F77B4",  # Epithelial - 深蓝
  "#D62728",  # T_cells - 鲜红
  "#17BECF",  # B_cells - 绿色
  "#9467BD",  # Plasma - 紫色
  "#FF7F0E",  # Macrophage - 橙色
  "#E377C2",  # Mast_cells - 洋红
  "#8C564B",  # Fibroblast - 深棕
  "#2CA02C",  # Fibroblast_SPP1 - 青蓝
  "#BCBD22"   # Endothelial - 黄绿
)
names(celltype_colors) <- celltype_levels

# 3️⃣ 作图（关键：不让 ggplot 重新排序）
p <- DimPlot(
  scRNA,
  reduction = "umap",
  group.by = "celltype",
  cols = celltype_colors,
  label = TRUE,
  repel = TRUE
) +
  scale_color_manual(
    values = celltype_colors,
    breaks = celltype_levels,   # 控制 legend 顺序
    limits = celltype_levels    # 防止 ggplot 自动重排
  ) +
  theme_classic(base_size = 14)
p
ggsave(p, file="DimPlot.pdf",width = 6,height = 4.5)
saveRDS(scRNA,"G:/课题二 成纤维细胞相关/Results/Public data/Result/3.降维聚类/scRNA.RDS")
