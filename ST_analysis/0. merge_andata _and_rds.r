#### 20251229 合并python中的counts+coordinate信息 & R中SingleR的注释信息 ####

library(reticulate)
Sys.setenv(RETICULATE_CONDA = path.expand("D:\\computation-software\\Anaconda3\\install\\condabin\\conda.bat"))
use_condaenv("scanpy310", required = TRUE)
library(anndata)
library(Seurat)

days = c("d0", "d1", "d3", "d5")

for (day in days) {
  rds_file_name <- paste0(day, "_singler.bigref.rds")
  adata_file_name <- paste0(day, ".adata_assigned.h5ad")
  adata <- read_h5ad(adata_file_name)
  pred <- readRDS(rds_file_name)
  celltypes <- pred$pruned.labels
  obj <- CreateSeuratObject(
    counts = t(adata$X),
    meta.data = data.frame(x = adata$obs$x, y = adata$obs$y, celltype = celltypes)
  )
  saveRDS(obj, file = paste0(day, ".rds"))
}
