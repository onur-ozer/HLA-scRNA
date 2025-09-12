library(Seurat)
library(tidyverse)
library(patchwork)
library(SeuratData)
library(harmony)
#InstallData("pbmc3k")
#InstallData("pbmcsca")
options(timeout = 1000)
#InstallData("lungref")
#InstallData("ifnb")

tmp <- AvailableData()

pbmc3k <- UpdateSeuratObject(pbmc3k)
pbmcsca <- UpdateSeuratObject(pbmcsca)
ifnb <- UpdateSeuratObject(ifnb)
lungref <- UpdateSeuratObject(object = 'lungref')

hla1  <- c("HLA-A","HLA-B","HLA-C")
hla2 <- c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1")
hla_genes <- c(hla1,hla2)


seurat_dataset <- pbmcsca
seurat_dataset <- seurat_dataset %>% 
  subset(subset = nCount_RNA < 20000) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
#  RunHarmony(group.by.vars = c('Experiment',"Method"), dims.use = 1:30) %>%
  #  FindNeighbors(dims = 1:20) %>%
#  FindClusters(resolution = 1, verbose = FALSE) %>%
  RunUMAP(dims = 1:30, reduction = 'pca', verbose = FALSE) #%>%
#  RunUMAP(dims = 1:30, reduction = 'harmony', verbose = FALSE)

DimPlot(seurat_dataset, reduction = "umap", label = TRUE, group.by = "CellType", label.size = 3)

DotPlot(seurat_dataset, features = hla_genes, group.by = "CellType") + RotatedAxis()# + ggtitle("MHC-I (control)")

DimPlot(seurat_dataset, reduction = "umap", label = TRUE, group.by = c('Experiment',"Method"), label.size = 3)




##### PBMC3K DATASET
seurat_dataset <- pbmc3k

seurat_dataset <- seurat_dataset %>% 
  subset(subset = nCount_RNA < 7000) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
#  FindNeighbors(dims = 1:20) %>%
#  FindClusters(resolution = 1, verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)

DimPlot(seurat_dataset, reduction = "umap", label = TRUE, group.by = "seurat_annotations", label.size = 3)

Idents(seurat_dataset)

# we need to tweak this a bit because NA values break the plots
seurat_dataset$seurat_annotations <- as.character(seurat_dataset$seurat_annotations)
seurat_dataset$seurat_annotations[is.na(seurat_dataset$seurat_annotations)] <- "Unknown"
DotPlot(seurat_dataset, features = hla_genes,  group.by = "seurat_annotations") + RotatedAxis()

seurat_dataset$seurat_annotations <- factor(seurat_dataset$seurat_annotations,
                                            levels = rev(c("Unknown", "Platelet","Memory CD4 T","Naive CD4 T","CD8 T","NK","CD14+ Mono","FCGR3A+ Mono","DC","B")))


##### IFN-B STIMULATED PBMC DATASET

seurat_dataset <- ifnb

# let's get rid of erythrocytes and megakaryocytes because very low expression of HLA in erythrocytes skews the results 
seurat_dataset <- subset(seurat_dataset,subset = !seurat_annotations %in% c("Eryth", "Mk"))
seurat_dataset$seurat_annotations <- droplevels(seurat_dataset$seurat_annotations)

seurat_dataset <- seurat_dataset %>% 
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunHarmony(group.by.vars = "stim", dims.use = 1:30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30,reduction.name = 'UMAPharmony') %>%
  RunUMAP(reduction = "pca", dims = 1:30,reduction.name = 'UMAPpca')

p1 <- DimPlot(seurat_dataset, reduction = "UMAPpca", label = TRUE, group.by = "seurat_annotations", label.size = 3)
p2 <- DimPlot(seurat_dataset, reduction = "UMAPharmony", label = TRUE, group.by = "seurat_annotations", label.size = 3)
p3 <- DimPlot(seurat_dataset, reduction = "UMAPpca", label = FALSE, group.by = "stim")
p4 <- DimPlot(seurat_dataset, reduction = "UMAPharmony", label = FALSE, group.by = "stim")

(p1 + p2) / (p3 + p4)
seurat_dataset$seurat_annotations <- factor(seurat_dataset$seurat_annotations,
                                            levels =rev(c("CD4 Naive T","CD4 Memory T","CD8 T","T activated","NK",
                                                          "CD16 Mono","CD14 Mono","DC","pDC", "B", "B Activated")))

DotPlot(seurat_dataset, features = hla_genes, group.by = "seurat_annotations", scale = TRUE) + 
  RotatedAxis()+ggtitle("MHC gene expression by cell type (complete dataset)")

DotPlot(seurat_dataset, features = hla_genes, group.by = "seurat_annotations", split.by = "stim",scale = TRUE) + 
  RotatedAxis()+ggtitle("MHC gene expression by cell type (CTRL vs IFN-B stimulated)")

DotPlot(seurat_dataset, features = c(hla1,hla2), group.by = "seurat_annotations", 
        split.by = "stim",scale = TRUE,cols = 'RdGy') + 
  RotatedAxis()+ggtitle("MHC genes by cell type (CTRL vs IFN-B)")


#seurat_dataset$celltype.stim <- paste(seurat_dataset$seurat_annotations, seurat_dataset$stim, sep = "_")
#Idents(seurat_dataset) <- "celltype.stim"

cell_type <- "B"
mono.de <- FindMarkers(seurat_dataset, 
                       ident.1 = paste(cell_type, "STIM", sep = "_"), 
                       ident.2 = paste(cell_type, "CTRL", sep = "_"), 
                       features = hla_genes,
                       test.use = "MAST",
                       # thresholds below are normally used to exclude rare genes or genes with similar expression
                       # as we do not want to exclude any of our HLA genes, we set them to zero
                       # that way even if there is no change in the expression, it will be included in the final dataset
                       logfc.threshold=0,min.cells.group=0,min.pct=0)


mono.de <- mono.de %>% mutate(gene = rownames(.), cell_type = cell_type)

ggplot(mono.de, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(data = mono.de, size = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  ggrepel::geom_text_repel(data = mono.de,aes(label = gene), max.overlaps = 20) +
  labs(title = paste("Volcano:", cell_type, "STIM vs CTRL"),
       x = "log2FC (STIM/CTRL)", y = "-log10(adj p)")




ct <- "T activated"
b.de <- FindMarkers(seurat_dataset, 
                    ident.1 = paste(ct, "STIM", sep = "_"), 
                    ident.2 = paste(ct, "CTRL", sep = "_"), 
                    features = hla_genes,
                    test.use = "MAST",
                    logfc.threshold=0,min.cells.group=0,min.pct=0)
#b.de <- b.de[rownames(b.de)%in%c(hla_genes),] %>% mutate(gene = rownames(.), ct = ct)
b.de <- b.de %>% mutate(gene = rownames(.), ct = ct)
ggplot(b.de, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(data = subset(b.de, gene %in% c(hla_genes)), size = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  ggrepel::geom_text_repel(data = subset(b.de, gene %in% c(hla_genes)),
                           aes(label = gene), max.overlaps = 20) +
  labs(title = paste0("STIM vs CTRL ",'(',ct,')'),
       x = "log2FC (STIM/CTRL)", y = "-log10(adj p)")



cell_types <- levels(seurat_dataset$seurat_annotations)

de_results_list <- list()

for (cell_type in cell_types) {

  subset_seu <- seurat_dataset
  Idents(subset_seu) <- subset_seu$seurat_annotations
  subset_seu <- subset(subset_seu, idents = cell_type)
  Idents(subset_seu) <- 'stim'
  
  de_genes <- FindMarkers(subset_seu,
                          ident.1 = "STIM",
                          ident.2 = "CTRL",
                          features = hla_genes,
                          test.use = "MAST",,
                          logfc.threshold=0,min.cells.group=0,min.pct=0)
  
  tmp_de <- de_genes %>% mutate(gene = rownames(.), cell_type = cell_type)
  vol_plot <- ggplot(tmp_de, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    geom_point(data = subset(tmp_de, gene %in% c(hla_genes)), size = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    ggrepel::geom_text_repel(data = subset(tmp_de, gene %in% c(hla_genes)),
                             aes(label = gene), max.overlaps = 20) +
    labs(title = paste0("STIM vs CTRL ",'(',cell_type,')'),
         x = "log2FC (STIM/CTRL)", y = "-log10(adj p)")
  print(vol_plot)
  
  de_genes$cell_type <- cell_type
  de_genes$gene <- rownames(de_genes)
  de_results_list[[cell_type]] <- de_genes
  }

combined_de_results <- do.call(rbind, de_results_list)
rownames(combined_de_results) <- NULL

ggplot(combined_de_results, aes(x = gene, y = cell_type, fill = avg_log2FC)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "lightgray", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "HLA Gene", y = "Cell Type", fill = "Log2FC\n(STIM vs CTRL)",
       title = "Differential Expression of HLA Genes after IFN-B Stimulation")














de_one_ct <- function(obj, ct, feats = NULL, assay = DefaultAssay(obj)) {
  x <- subset(obj, seurat_annotations == ct)
  Idents(x) <- "stim"
  FindMarkers(x, ident.1 = "STIM", ident.2 = "CTRL",
              features = feats, logfc.threshold = 0, min.pct = 0,
              test.use = if (assay == "RNA") "MAST" else "wilcox") %>%
    mutate(gene = rownames(.), ct = ct)
}
ct <- "CD14 Mono"    # pick a cell type to highlight
de <- de_one_ct(seurat_dataset, ct, feats = hla2)

p2 <- ggplot(de, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
  geom_point(alpha = 0.6) +
  geom_point(data = subset(de, gene %in% c(hla1,hla2)), size = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  ggrepel::geom_text_repel(data = subset(de, gene %in% c(hla1,hla2)),
                           aes(label = gene), max.overlaps = 20) +
  labs(title = paste("Volcano:", ct, "STIM vs CTRL"),
       x = "log2FC (STIM/CTRL)", y = "-log10(adj p)")

p1 + p2





ct <- "B"

# 1) Work on the same object & assay for both analyses
DefaultAssay(seurat_dataset) <- "RNA"   # ensure consistency for this comparison

# 2) Subset to one cell type
x <- subset(seurat_dataset, seurat_annotations == ct)
Idents(x) <- "stim"

# 3) What does the volcano use? Recalculate DE (RNA):
de_mast <- FindMarkers(
  x, ident.1 = "STIM", ident.2 = "CTRL",
  logfc.threshold = 0, min.pct = 0,
  test.use = "MAST"
)
de_mast[rownames(de_mast) == "HLA-DQB1", , drop = FALSE]

# 4) Simple group means that DotPlot is based on (RNA/data slot):
m <- AverageExpression(x, features = "HLA-DQB1", group.by = "stim", assays = "RNA")$RNA
m  # columns CTRL, STIM

# Also inspect per-cell stats directly:
df <- FetchData(x, vars = c("HLA-DQB1","stim"))
aggregate(`HLA-DQB1` ~ stim, df, function(z) c(mean = mean(z), pct = mean(z > 0)))

