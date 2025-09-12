# T-cell analysis
# publication: https://doi.org/10.1016/j.immuni.2020.05.010
# data access: GSE139042
library(harmony)
library(Seurat)
library(tidyverse)
library(patchwork)
library(DoubletFinder)
# use version3 assay structure for this analysis
options(Seurat.object.assay.version = "v3")

# read count matrices
m_count_T1_34p <- as.matrix(read.table('data/t_cell_data/GSE139042_RAW/GSM4127993_Thymus1counts_34positive.txt',
                                       header = TRUE, row.names = 1,sep = '\t'))
m_count_T1_34n <- as.matrix(read.table('data/t_cell_data/GSE139042_RAW/GSM4127994_Thymus1counts_34negative.txt',
                                       header = TRUE, row.names = 1,sep = '\t'))
m_count_T2_34p <- as.matrix(read.table('data/t_cell_data/GSE139042_RAW/GSM4505165_Thymus2counts_34positive.txt',
                                       header = TRUE, row.names = 1,sep = '\t'))
m_count_T3_34p <- as.matrix(read.table('data/t_cell_data/GSE139042_RAW/GSM4505166_Thymus3counts_34positive.txt',
                                       header = TRUE, row.names = 1,sep = '\t'))
m_count_T3_34n <- as.matrix(read.table('data/t_cell_data/GSE139042_RAW/GSM4505167_Thymus3counts_34negative.txt',
                                       header = TRUE, row.names = 1,sep = '\t'))
m_count_T7_34n <- as.matrix(read.table('data/t_cell_data/GSE139042_RAW/GSM4505171_Thymus7counts_34negative.txt',
                                       header = TRUE, row.names = 1,sep = '\t'))
m_count_T7_CD4SP <- as.matrix(read.table('data/t_cell_data/GSE139042_RAW/GSM4505172_Thymus7counts_CD4SP.txt',
                                         header = TRUE, row.names = 1,sep = '\t'))
m_count_T7_CD8SP <- as.matrix(read.table('data/t_cell_data/GSE139042_RAW/GSM4505173_Thymus7counts_CD8SP.txt',
                                         header = TRUE, row.names = 1,sep = '\t'))


# generate seurat objects for each gene
# we do some QC and filtering
# thresholds are in the S1 of the paper

count_T1_34p <- CreateSeuratObject(m_count_T1_34p,min.cells = 3)
count_T1_34n <- CreateSeuratObject(m_count_T1_34n,min.cells = 3)
count_T2_34p <- CreateSeuratObject(m_count_T2_34p,min.cells = 3)
count_T3_34p <- CreateSeuratObject(m_count_T3_34p,min.cells = 3)
count_T3_34n <- CreateSeuratObject(m_count_T3_34n,min.cells = 3)
count_T7_34n <- CreateSeuratObject(m_count_T7_34n,min.cells = 3)
count_T7_CD4SP <- CreateSeuratObject(m_count_T7_CD4SP,min.cells = 3)
count_T7_CD8SP <- CreateSeuratObject(m_count_T7_CD8SP,min.cells = 3)

# we will use this list later
dataset_list=list(count_T1_34p, count_T2_34p, count_T3_34p,
                  count_T1_34n, count_T3_34n, count_T7_34n,
                  count_T7_CD4SP, count_T7_CD8SP)

# first, let's start with the percentage of mt genes
# calculate mt gene percentage for each dataset

count_T1_34p[["percent.mt"]] <- PercentageFeatureSet(count_T1_34p, pattern = "^MT[-\\.]")
count_T1_34n[["percent.mt"]] <- PercentageFeatureSet(count_T1_34n, pattern = "^MT[-\\.]")
count_T2_34p[["percent.mt"]] <- PercentageFeatureSet(count_T2_34p, pattern = "^MT[-\\.]")
count_T3_34p[["percent.mt"]] <- PercentageFeatureSet(count_T3_34p, pattern = "^MT[-\\.]")
count_T3_34n[["percent.mt"]] <- PercentageFeatureSet(count_T3_34n, pattern = "^MT[-\\.]")
count_T7_34n[["percent.mt"]] <- PercentageFeatureSet(count_T7_34n, pattern = "^MT[-\\.]")
count_T7_CD4SP[["percent.mt"]] <- PercentageFeatureSet(count_T7_CD4SP, pattern = "^MT[-\\.]")
count_T7_CD8SP[["percent.mt"]] <- PercentageFeatureSet(count_T7_CD8SP, pattern = "^MT[-\\.]")

count_T1_34p <- subset(count_T1_34p,subset = nFeature_RNA > 200 & nCount_RNA < 30000 & percent.mt < 10) #nFeatureRNA > 200 can be skipped
count_T1_34n <- subset(count_T1_34n,subset = nFeature_RNA > 200 & nCount_RNA < 30000 & percent.mt < 10)
count_T2_34p <- subset(count_T2_34p,subset = nFeature_RNA > 200 & nCount_RNA < 40000 & percent.mt < 20)
count_T3_34p <- subset(count_T3_34p,subset = nFeature_RNA > 200 & nCount_RNA < 50000 & percent.mt < 20)
count_T3_34n <- subset(count_T3_34n,subset = nFeature_RNA > 200 & nCount_RNA < 35000 & percent.mt < 20)
count_T7_34n <- subset(count_T7_34n,subset = nFeature_RNA > 200 & nCount_RNA < 40000 & percent.mt < 20)
count_T7_CD4SP <- subset(count_T7_CD4SP,subset = nFeature_RNA > 200 & nCount_RNA < 20000 & percent.mt < 20)
count_T7_CD8SP <- subset(count_T7_CD8SP,subset = nFeature_RNA > 200 & nCount_RNA < 30000 & percent.mt < 20)

# Normalize the data
count_T1_34p <- NormalizeData(count_T1_34p)
count_T1_34n <- NormalizeData(count_T1_34n)
count_T2_34p <- NormalizeData(count_T2_34p)
count_T3_34p <- NormalizeData(count_T3_34p)
count_T3_34n <- NormalizeData(count_T3_34n)
count_T7_34n <- NormalizeData(count_T7_34n)
count_T7_CD4SP <- NormalizeData(count_T7_CD4SP)
count_T7_CD8SP <- NormalizeData(count_T7_CD8SP)


# regress out cell cycle effects. use the gene list provided in the Supp. File 1

gene_list <- readxl::read_excel('data/t_cell_data/NIHMS1599710-supplement-2.xlsx',sheet = 2)
s_genes <- unique(na.omit(trimws(gene_list$`cell cycle (S phase genes)`)))
g2m_genes <- unique(na.omit(gene_list$`cell cycle (G2/M phase genes)`))



count_T1_34p <- CellCycleScoring(count_T1_34p, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
count_T1_34n <- CellCycleScoring(count_T1_34n, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
count_T2_34p <- CellCycleScoring(count_T2_34p, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
count_T3_34p <- CellCycleScoring(count_T3_34p, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
count_T3_34n <- CellCycleScoring(count_T3_34n, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
count_T7_34n <- CellCycleScoring(count_T7_34n, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
count_T7_CD4SP <- CellCycleScoring(count_T7_CD4SP, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
count_T7_CD8SP <- CellCycleScoring(count_T7_CD8SP, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)

# feature selection (find highly variable genes / hvgs)

count_T1_34p <- FindVariableFeatures(count_T1_34p)
count_T1_34n <- FindVariableFeatures(count_T1_34n)
count_T2_34p <- FindVariableFeatures(count_T2_34p)
count_T3_34p <- FindVariableFeatures(count_T3_34p)
count_T3_34n <- FindVariableFeatures(count_T3_34n)
count_T7_34n <- FindVariableFeatures(count_T7_34n)
count_T7_CD4SP <- FindVariableFeatures(count_T7_CD4SP)
count_T7_CD8SP <- FindVariableFeatures(count_T7_CD8SP)


# scaling
count_T1_34p <- ScaleData(count_T1_34p, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","percent.mt"))
count_T1_34n <- ScaleData(count_T1_34n, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","percent.mt"))
count_T2_34p <- ScaleData(count_T2_34p, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","percent.mt"))
count_T3_34p <- ScaleData(count_T3_34p, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","percent.mt"))
count_T3_34n <- ScaleData(count_T3_34n, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","percent.mt"))
count_T7_34n <- ScaleData(count_T7_34n, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","percent.mt"))
count_T7_CD4SP <- ScaleData(count_T7_CD4SP, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","percent.mt"))
count_T7_CD8SP <- ScaleData(count_T7_CD8SP, vars.to.regress = c("S.Score","G2M.Score","nCount_RNA","percent.mt"))

# PCA / dimension reduction
count_T1_34p <- RunPCA(count_T1_34p, npcs = 50)
count_T1_34n <- RunPCA(count_T1_34n, npcs = 50)
count_T2_34p <- RunPCA(count_T2_34p, npcs = 50)
count_T3_34p <- RunPCA(count_T3_34p, npcs = 50)
count_T3_34n <- RunPCA(count_T3_34n, npcs = 50)
count_T7_34n <- RunPCA(count_T7_34n, npcs = 50)
count_T7_CD4SP <- RunPCA(count_T7_CD4SP, npcs = 50)
count_T7_CD8SP <- RunPCA(count_T7_CD8SP, npcs = 50)

# we can have a look at the elbow plots to see the explained variation per PC 
ElbowPlot(count_T1_34p, ndims = ncol(Embeddings(count_T1_34p, "pca")))

# clustering

count_T1_34p <- FindNeighbors(count_T1_34p, dims = 1:30)
count_T1_34p <- FindClusters(count_T1_34p, resolution = 1)
count_T1_34p <- RunUMAP(count_T1_34p, dims = 1:30)
count_T1_34p <- RunTSNE(count_T1_34p, dims = 1:30)

count_T1_34n <- FindNeighbors(count_T1_34n, dims = 1:30)
count_T1_34n <- FindClusters(count_T1_34n, resolution = 1)
count_T1_34n <- RunUMAP(count_T1_34n, dims = 1:30)
count_T1_34n <- RunTSNE(count_T1_34n, dims = 1:30)

count_T2_34p <- FindNeighbors(count_T2_34p, dims = 1:30)
count_T2_34p <- FindClusters(count_T2_34p, resolution = 1)
count_T2_34p <- RunUMAP(count_T2_34p, dims = 1:30)
count_T2_34p <- RunTSNE(count_T2_34p, dims = 1:30)

count_T3_34p <- FindNeighbors(count_T3_34p, dims = 1:30)
count_T3_34p <- FindClusters(count_T3_34p, resolution = 1)
count_T3_34p <- RunUMAP(count_T3_34p, dims = 1:30)
count_T3_34p <- RunTSNE(count_T3_34p, dims = 1:30)

count_T3_34n <- FindNeighbors(count_T3_34n, dims = 1:30)
count_T3_34n <- FindClusters(count_T3_34n, resolution = 1)
count_T3_34n <- RunUMAP(count_T3_34n, dims = 1:30)
count_T3_34n <- RunTSNE(count_T3_34n, dims = 1:30)

count_T7_34n <- FindNeighbors(count_T7_34n, dims = 1:30)
count_T7_34n <- FindClusters(count_T7_34n, resolution = 1)
count_T7_34n <- RunUMAP(count_T7_34n, dims = 1:30)
count_T7_34n <- RunTSNE(count_T7_34n, dims = 1:30)

count_T7_CD4SP <- FindNeighbors(count_T7_CD4SP, dims = 1:30)
count_T7_CD4SP <- FindClusters(count_T7_CD4SP, resolution = 1)
count_T7_CD4SP <- RunUMAP(count_T7_CD4SP, dims = 1:30)
count_T7_CD4SP <- RunTSNE(count_T7_CD4SP, dims = 1:30)

count_T7_CD8SP <- FindNeighbors(count_T7_CD8SP, dims = 1:30)
count_T7_CD8SP <- FindClusters(count_T7_CD8SP, resolution = 1)
count_T7_CD8SP <- RunUMAP(count_T7_CD8SP, dims = 1:30)
count_T7_CD8SP <- RunTSNE(count_T7_CD8SP, dims = 1:30)

p1 <- DimPlot(count_T1_34p,reduction = 'umap', label = TRUE) + ggtitle("T1_34p")
p2 <- DimPlot(count_T1_34n,reduction = 'umap', label = TRUE) + ggtitle("T1_34n")
p3 <- DimPlot(count_T2_34p,reduction = 'umap', label = TRUE) + ggtitle("T2_34p")
p4 <- DimPlot(count_T3_34p,reduction = 'umap', label = TRUE) + ggtitle("T3_34p")
p5 <- DimPlot(count_T3_34n,reduction = 'umap', label = TRUE) + ggtitle("T3_34n")
p6 <- DimPlot(count_T7_34n,reduction = 'umap', label = TRUE) + ggtitle("T7_34n")
p7 <- DimPlot(count_T7_CD4SP,reduction = 'umap', label = TRUE) + ggtitle("T7_CD4SP")
p8 <- DimPlot(count_T7_CD8SP,reduction = 'umap', label = TRUE) + ggtitle("T7_CD8SP")
(p1|p2) / (p3|p4) / (p5|p6) / (p7|p8)

# most of the steps above are just applying the same function with the same parameters to different datasets
# instead of doing this many times for each dataset, we can also write a function including all steps from normalization to clustering
# then we can apply this to all our objects either in a loop or with using lapply.
# below is a function to do that

# process_seurat_object <- function(seurat_obj, s_genes, g2m_genes) {
#   seurat_obj <- NormalizeData(seurat_obj)
#   seurat_obj <- CellCycleScoring(seurat_obj,s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)#   
#   seurat_obj <- FindVariableFeatures(seurat_obj)
#   seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt"))
#   seurat_obj <- RunPCA(seurat_obj, npcs = 50)
#   seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
#   seurat_obj <- FindClusters(seurat_obj, resolution = 1)
#   seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
#   seurat_obj <- RunTSNE(seurat_obj, dims = 1:30)
#   return(seurat_obj)
# }

# alternatively, we can use pipe "%>%" and apply all commands to a single dataset and repeat this for all datasets.
# this would be more useful if we need to tweak parameters for different datasets


# based on clustering, authors applied several more cut-offs and removed some cluster
# first, a cluster in Thymus 2 CD34+ sample had increased expression mitochondrial genes. Let's find out which.
# we already have a percent.mt column that gives us the percentage of mt genes within each cell.
# if we get the mean for each cluster, we will see how clusters differ from each other with respect to the expression of mt genes

mt_by_cluster <- tapply(count_T2_34p$percent.mt, Idents(count_T2_34p), mean)
sort(mt_by_cluster)
# all clusters are around 4 except for one cluster which has mean of 9.97. 
# this is most likely the one authors mentioned. you can see it with a boxplot as well. #boxplot(mt_by_cluster) 
# let's remove it.

drop_mt_cluster <- names(mt_by_cluster)[mt_by_cluster>9]
count_T2_34p <- subset(count_T2_34p, subset = seurat_clusters!=drop_mt_cluster)


# next, in all CD34- cells, there is a cluster with high expression of PAX5, CD19, and CD22
# this cluster likely represents thymic B cells and will be excluded

bcell_markers <- c("PAX5","CD19","CD22")

avg_b <- AverageExpression(count_T1_34n, features = bcell_markers)$RNA # gives average per cluster for each gene in the list
b_scores <- colMeans(as.matrix(avg_b), na.rm = TRUE) # mean expression of three genes per cluster
drop_B <- gsub("^g", "", names(b_scores)[b_scores>1])
count_T1_34n <- subset(count_T1_34n, subset = seurat_clusters!=drop_B)


# count_T1_34n_updatedB <- subset(count_T1_34n, subset = seurat_clusters!=drop_B)
# #unique(count_T1_34n@meta.data$seurat_clusters) # sanity check
# ptmp <- DimPlot(count_T1_34n_updatedB,reduction = 'umap', label = TRUE) + ggtitle("T1_34n")
# 
# umap_coords <- Embeddings(count_T1_34n, reduction = "umap")
# 
# # Look at the first rows to confirm
# head(umap_coords)
# 
# # Find the cells where UMAP_1 < -10
# cells_subset <- rownames(umap_coords)[umap_coords[,1] < -10]
# clusters_vec <- Idents(count_T1_34n)
# clusters_of_subset <- clusters_vec[cells_subset]
# 
# 
# count_T1_34n_updatedB <- FindNeighbors(count_T1_34n_updatedB, dims = 1:30)
# count_T1_34n_updatedB <- FindClusters(count_T1_34n_updatedB, resolution = 1)
# count_T1_34n_updatedB <- RunUMAP(count_T1_34n_updatedB, dims = 1:30)
# count_T1_34n_updatedB <- RunTSNE(count_T1_34n_updatedB, dims = 1:30)
# ptmp2 <- DimPlot(count_T1_34n_updatedB,reduction = 'umap', label = TRUE) + ggtitle("T1_34n")


# repeat for other CD34- samples

avg_b <- AverageExpression(count_T3_34n, features = bcell_markers)$RNA # gives average per cluster for each gene in the list
b_scores <- colMeans(as.matrix(avg_b), na.rm = TRUE) # mean expression of three genes per cluster
drop_B <- gsub("^g", "", names(b_scores)[b_scores>1])
count_T3_34n <- subset(count_T3_34n, subset = seurat_clusters!=drop_B)


avg_b <- AverageExpression(count_T7_34n, features = bcell_markers)$RNA # gives average per cluster for each gene in the list
b_scores <- colMeans(as.matrix(avg_b), na.rm = TRUE) # mean expression of three genes per cluster
drop_B <- gsub("^g", "", names(b_scores)[b_scores>1])
count_T7_34n <- subset(count_T7_34n, subset = seurat_clusters!=drop_B)


# now we cleaned the data, we re-run the clustering for the samples that we changed

count_T1_34n <- FindNeighbors(count_T1_34n, dims = 1:30)
count_T1_34n <- FindClusters(count_T1_34n, resolution = 1)
count_T1_34n <- RunUMAP(count_T1_34n, dims = 1:30)
count_T1_34n <- RunTSNE(count_T1_34n, dims = 1:30)

count_T2_34p <- FindNeighbors(count_T2_34p, dims = 1:30)
count_T2_34p <- FindClusters(count_T2_34p, resolution = 1)
count_T2_34p <- RunUMAP(count_T2_34p, dims = 1:30)
count_T2_34p <- RunTSNE(count_T2_34p, dims = 1:30)

count_T3_34n <- FindNeighbors(count_T3_34n, dims = 1:30)
count_T3_34n <- FindClusters(count_T3_34n, resolution = 1)
count_T3_34n <- RunUMAP(count_T3_34n, dims = 1:30)
count_T3_34n <- RunTSNE(count_T3_34n, dims = 1:30)

count_T7_34n <- FindNeighbors(count_T7_34n, dims = 1:30)
count_T7_34n <- FindClusters(count_T7_34n, resolution = 1)
count_T7_34n <- RunUMAP(count_T7_34n, dims = 1:30)
count_T7_34n <- RunTSNE(count_T7_34n, dims = 1:30)



# let's have a look at the data again

p1 <- DimPlot(count_T1_34p,reduction = 'umap', label = TRUE) + ggtitle("T1_34p")
p2 <- DimPlot(count_T1_34n,reduction = 'umap', label = TRUE) + ggtitle("T1_34n")
p3 <- DimPlot(count_T2_34p,reduction = 'umap', label = TRUE) + ggtitle("T2_34p")
p4 <- DimPlot(count_T3_34p,reduction = 'umap', label = TRUE) + ggtitle("T3_34p")
p5 <- DimPlot(count_T3_34n,reduction = 'umap', label = TRUE) + ggtitle("T3_34n")
p6 <- DimPlot(count_T7_34n,reduction = 'umap', label = TRUE) + ggtitle("T7_34n")
p7 <- DimPlot(count_T7_CD4SP,reduction = 'umap', label = TRUE) + ggtitle("T7_CD4SP")
p8 <- DimPlot(count_T7_CD8SP,reduction = 'umap', label = TRUE) + ggtitle("T7_CD8SP")
(p1|p2) / (p3|p4) / (p5|p6) / (p7|p8)



# now we are done with the cleaning, we will merge datasets
# unlike the paper, we will merge all datasets together

all_merged <- merge(count_T1_34p,
                    y = list(count_T2_34p, count_T3_34p,
                             count_T1_34n,count_T3_34n,count_T7_34n,
                             count_T7_CD4SP,count_T7_CD8SP),
                    add.cell.ids = c("T1_34p","T2_34p","T3_34p",
                                     "T1_34n","T3_34n","T7_34n",
                                     "CD4SP", "CD8SP"),
                    project = "all_merged")

all_merged <- CellCycleScoring(all_merged, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
all_merged <- ScaleData(all_merged,vars.to.regress = c("S.Score","G2M.Score"))

# this was a heavy step and it is better to save the object now just in case
saveRDS(all_merged, file = "data/t_cell_data/all_merged_scaled.rds")

# for the variable feature set, we follow the paper's approach (including T7 as well)
# first, find the top 1000 HVGs for each set and find the union of them.

hvgs_list <- list(
  T1_34p = head(VariableFeatures(count_T1_34p), 1000),
  T2_34p = head(VariableFeatures(count_T2_34p), 1000),
  T3_34p = head(VariableFeatures(count_T3_34p), 1000),
  T1_34n = head(VariableFeatures(count_T1_34n), 1000),
  T3_34n = head(VariableFeatures(count_T3_34n), 1000),
  T7_34n = head(VariableFeatures(count_T7_34n), 1000),
  CD4SP  = head(VariableFeatures(count_T7_CD4SP), 1000),
  CD8SP  = head(VariableFeatures(count_T7_CD8SP), 1000)
)

hvgs_union <- (unique(unlist(hvgs_list)))

# now find genes common in all datasets
genes_common <- rownames(dataset_list[[1]])
for (i in 2:length(dataset_list)){genes_common <- intersect(genes_common,rownames(dataset_list[[i]]))}

# find the intersect of these two datasets, giving us the HVGs present in all samples
hvgs_final <- intersect(hvgs_union, genes_common)
VariableFeatures(all_merged) <- hvgs_final

# run PCA on common HVGs
all_merged <- RunPCA(all_merged, features = hvgs_final, npcs = 50, verbose = FALSE)


# Now we run harmony to control for the batch effect introduced by different thymus samples

cn <- colnames(all_merged)

all_merged$thymus <- dplyr::case_when(
  grepl("^T1_", cn) ~ "T1",
  grepl("^T2_", cn) ~ "T2",
  grepl("^T3_", cn) ~ "T3",
  grepl("^T7_", cn) ~ "T7",
  grepl("^CD4SP_", cn) ~ "T7",
  grepl("^CD8SP_", cn) ~ "T7",
  TRUE ~ NA_character_
)


all_merged <- RunHarmony(
  all_merged,
  group.by.vars = "thymus",
  reduction.use = "pca",
  dims.use = 1:20,
  theta= 2
)


all_merged <- FindNeighbors(all_merged, reduction = "harmony", dims = 1:20, k.param = 30)
all_merged <- FindClusters(all_merged, resolution = 0.5)  # tune if needed
all_merged <- RunUMAP(all_merged, reduction = "harmony", dims = 1:20)
#all_merged <- RunTSNE(all_merged, reduction = "harmony", dims = 1:20)


# pre/post mixing check
#DimPlot(all_merged, reduction = "pca", group.by = "thymus")      # BEFORE Harmony
# (after running Harmony below)
DimPlot(all_merged, reduction = "umap", group.by = "thymus")     # AFTER Harmony
#DimPlot(all_merged, reduction = "tsne", group.by = "thymus")     # AFTER Harmony



# 1) Define panels
markers <- list(
  DN1_ETP = c("CD34","KIT","IL7R","CCR9","NOTCH1","HES1","TCF7","LYL1","LMO2","MEF2C","HHEX"),
  DN2     = c("BCL11B","GATA3","TCF7","IL7R","CCR9","NOTCH1","HES1","DTX1"),
  DN3     = c("PTCRA","RAG1","RAG2","DNTT","CD3D","CD3E","CD3G","LCK","ZAP70"),
  DN4     = c("MKI67","TOP2A","UBE2C","HMGB2","TYMS"),
  DP      = c("CD4","CD8A","CD8B","TRAC","CD3D","CD3E","CD3G","THEMIS","LCK","ZAP70"),
  CD4SP   = c("CD4","IL7R","CCR7","SELL","SATB1","GATA3","ZBTB7B","BCL11B"),
  CD8SP   = c("CD8A","CD8B","RUNX3","IKZF1","ITGAE")
)

# 2) Keep only genes present in your matrix to avoid warnings
markers <- lapply(markers, function(gs) Seurat::CaseMatch(gs, rownames(all_merged)))

# 3) Add module scores (one column per panel, suffixed with '1' by Seurat)
for (nm in names(markers)) {
  all_merged <- AddModuleScore(all_merged, features = list(markers[[nm]]), name = paste0(nm, "_score"))
}

# 4) Visualize on Harmony UMAP
FeaturePlot(all_merged,reduction = 'umap',
            features = c("DN1_ETP_score1","DN2_score1","DN3_score1","DN4_score1",
                         "DP_score1","CD4SP_score1","CD8SP_score1"),
            min.cutoff = "q10", max.cutoff = "q90", ncol = 3
)


# 5) (Optional) Derive a single stage label by argmax of scores
score_cols <- c("DN1_ETP_score1","DN2_score1","DN3_score1","DN4_score1","DP_score1","CD4SP_score1","CD8SP_score1")
S <- as.matrix(FetchData(all_merged, vars = score_cols))
lab <- gsub("_score1$", "", colnames(S))[max.col(S, ties.method = "first")]
all_merged$Stage_refined <- factor(lab, levels = c("DN1_ETP","DN2","DN3","DN4","DP","CD4SP","CD8SP"))
DimPlot(all_merged, group.by = "Stage_refined", label = TRUE,reduction = 'umap') + ggtitle("DN → DP → SP continuum (refined)")
DimPlot(all_merged, group.by = "Stage_refined", label = TRUE,reduction = 'tsne') + ggtitle("DN → DP → SP continuum (refined)")


FeaturePlot(all_merged, features = c("CD34", "CD4", "CD8A", "CD8B"), 
            reduction = 'umap')
DotPlot(all_merged, features = c("CD34", "CD44", "CD25", "CD4", "CD8A", "CD8B"), 
        group.by = "Stage_refined") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



