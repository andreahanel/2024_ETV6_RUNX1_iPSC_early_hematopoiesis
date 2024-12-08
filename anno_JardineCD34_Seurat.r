library(Seurat)
library(SeuratObject)
library(dplyr)
#library(SeuratDisk)
#library(patchwork)
#library(ggplot2)
library(anndata)

options(future.globals.maxSize = 1e9)
fig <- function(width, heigth){
options(repr.plot.width = width, repr.plot.height = heigth) }


## reference --------
adata_jardine <- read_h5ad('/research/groups/sysgen/PROJECTS/sysgen_team/andrea_work/VAE-Andrea/output/jardine_processed_withrawcounts.h5ad')

mat <- t(adata_jardine$X)


jardine.ref <- CreateSeuratObject(counts = mat, 
                             meta.data = adata_jardine$obs)

# Identify highly variable genes from adata_jardine
#hvg <- rownames(adata_jardine$var)[adata_jardine$var$highly_variable]

# Set highly variable genes in the Seurat object
#VariableFeatures(jardine.ref) <- hvg

umap_coords <- adata_jardine$obsm$X_umap

# Ensure row names in UMAP coordinates match the Seurat object cell names
rownames(umap_coords) <- rownames(jardine.ref@meta.data)

# Add UMAP coordinates to the Seurat object as a new DimReduc object
jardine.ref[["X_umap"]] <- CreateDimReducObject(embeddings = umap_coords, key = "UMAP_")

# pre-process dataset (without integration)
jardine.ref <- NormalizeData(jardine.ref)
jardine.ref <- FindVariableFeatures(jardine.ref)
jardine.ref <- ScaleData(jardine.ref)
jardine.ref <- RunPCA(jardine.ref)

adata = read_h5ad ("/research/groups/sysgen/PROJECTS/sysgen_team/andrea_work/2023_Fischer_iPSC/revision_final/data/adata2_PCAharmony_leiden_LS.h5ad")
mat <- t(adata$X)
query <- CreateSeuratObject(counts = mat, 
                             meta.data = adata$obs)

query[["RNA"]] <- split(query[["RNA"]], f = query$sample)
query <- NormalizeData(query)

query.anchors <- FindTransferAnchors(reference = jardine.ref, 
                                        query = query, 
                                        dims = 1:30,
                                        reference.reduction = "pca")


predictions <- TransferData(anchorset = query.anchors, 
                            refdata = jardine.ref$cell.labels, 
                            dims = 1:30)
#
query <- AddMetaData(query, metadata = predictions)
query <- JoinLayers(query)
saveRDS(object = query, 
          file = "/research/groups/sysgen/PROJECTS/sysgen_team/andrea_work/2023_Fischer_iPSC/revision_final/data/seurat_mapping_nonscaled_defaultHVG.rds"
)
adata = read_h5ad ("/research/groups/sysgen/PROJECTS/sysgen_team/andrea_work/2023_Fischer_iPSC/revision_final/data/adata2_PCAharmony_leiden_LS.h5ad")

  # Extract the Seurat object metadata (obs)
obs <- query[[]]
ad <- AnnData(X = adata$X,  # original count matrix
              obs = obs,    # updated metadata
              var = adata$var#,  # original variable metadata
             )

ad$write_h5ad("/research/groups/sysgen/PROJECTS/sysgen_team/andrea_work/2023_Fischer_iPSC/revision_final/data/seurat_mapping_nonscaled_defaultHVG.h5ad")




