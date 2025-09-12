HLA expression on scRNA-seq data
================
Onur Özer
2025-09-08

There is a well-known text book knowledge: HLA-classI genes are
expressed in all cells containing nucleus, while HLA-classII genes are
mainly expressed in antigen-presenting cells (APCs).

Now that we have many available single-cell RNA sequencing data, we can
have a look at this again.

First, let’s load some libraries that will be used.

``` r
library(Seurat)
library(tidyverse)
library(patchwork)
library(SeuratData)
library(harmony)
```

Now we will install some interesting scRNA-seq datasets from the
SeuratData package. These are small datasets but if the installation
fails, it is most likely because of R’s limit of 60 seconds of download.
Try setting `options(timeout = 600)`

``` r
InstallData("pbmc3k")
InstallData("pbmcsca")
InstallData("ifnb")

pbmc3k <- UpdateSeuratObject(pbmc3k)
pbmcsca <- UpdateSeuratObject(pbmcsca)
ifnb <- UpdateSeuratObject(ifnb)
```

Finally, we define HLA genes to be analyzed. If the dataset use a
different naming for genes (for example ENSEMBL IDs), you may either
change this list or -better- change the ENSEMBL IDs in your dataset to
gene symbols.

``` r
hla1  <- c("HLA-A","HLA-B","HLA-C")
hla2 <- c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1")
hla_genes <- c(hla1,hla2)
```

Now we can start our analysis. Let’s begin with the `pbmc3k` and get a
quick peek. This is a commonly used dataset of peripheral blood
mononuclear cells.

``` r
head(pbmc3k@meta.data)
```

    ##                orig.ident nCount_RNA nFeature_RNA seurat_annotations
    ## AAACATACAACCAC     pbmc3k       2419          779       Memory CD4 T
    ## AAACATTGAGCTAC     pbmc3k       4903         1352                  B
    ## AAACATTGATCAGC     pbmc3k       3147         1129       Memory CD4 T
    ## AAACCGTGCTTCCG     pbmc3k       2639          960         CD14+ Mono
    ## AAACCGTGTATGCG     pbmc3k        980          521                 NK
    ## AAACGCACTGGTAC     pbmc3k       2163          781       Memory CD4 T

Above, as rownames we have cell barcodes. Columns give us the dataset,
the number of RNA reads (nCount_RNA) and the number of features
(nFeature_RNA; in our case, genes). We also have already available
annotations for each cell, telling us the cell type.

We start with regular pre-processing and quality control.

``` r
pbmc3k_dataset <- pbmc3k

pbmc3k_dataset <- pbmc3k_dataset %>% 
  subset(subset = nCount_RNA < 7000) %>%
# 7000 here is somehow arbitrary. 
# The idea is to remove cells with very high RNA counts as these are likely doublets
  NormalizeData() %>% # normalize across each cell to account for library size differences
  FindVariableFeatures() %>% # focus only on the variable genes. default is to use top 2000
  ScaleData() %>% # scale and center across each gene to account for normal expression differences between genes
  RunPCA() %>%
  RunUMAP(dims = 1:20) %>%
  RunTSNE(dims = 1:20)


# normally, the aim of scRNA analysis is to find clusters of cells and annotate them. 
# For this we would need functions below to identify clusters and later annotate these clusters.
# However, we already have our annotations in the dataset, so we skip this step.
#  FindNeighbors(dims = 1:20) %>% 
#  FindClusters(resolution = 1, verbose = FALSE)
```

Here is how the data looks with the UMAP and TSNE

``` r
p1 <- DimPlot(pbmc3k_dataset, reduction = "umap", label = TRUE, group.by = "seurat_annotations", label.size = 3)
p2 <- DimPlot(pbmc3k_dataset, reduction = "tsne", label = TRUE, group.by = "seurat_annotations", label.size = 3)
p1+p2
```

<img src="HLA_expression_files/figure-gfm/unnamed-chunk-6-1.png" width="100%" />

We have both APCs and other cell types so we can look at how HLA gene
expression differs between these.

``` r
# we need to tweak annotations a bit because NA values break the plot
pbmc3k_dataset$seurat_annotations <- as.character(pbmc3k_dataset$seurat_annotations)
pbmc3k_dataset$seurat_annotations[is.na(pbmc3k_dataset$seurat_annotations)] <- "Unknown"

# sort cells so that the plot looks neat
pbmc3k_dataset$seurat_annotations <- factor(pbmc3k_dataset$seurat_annotations,
                                            levels = rev(c("Unknown", "Platelet","Memory CD4 T","Naive CD4 T","CD8 T","NK","CD14+ Mono","FCGR3A+ Mono","DC","B")))

DotPlot(pbmc3k_dataset, features = c(hla_genes),  group.by = "seurat_annotations") + RotatedAxis()
```

![](HLA_expression_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Right away, we can see the pattern for HLA-classII genes. APCs robustly
express them but expression in other cell types is sporadic. So we can
confirm our beloved Janeway’s Immunobiology was not lying to us.

There are other differences that catches the eye as well. For example
monocytes appear to have lower expression of HLA-DQ compared to DCs and
B cells. There are examples from the literature confirming low
expression of DQ in monocytes.

Another interesting observation is the slightly higher expression of HLA
classII genes in CD8+ T cells. This may suggest that some these cells
were closer to an activated state.

Nice, let’s try another dataset `pbmcsca`. We do the same processing.

``` r
pbmcsca_dataset <- pbmcsca
pbmcsca_dataset <- pbmcsca_dataset %>% 
  subset(subset = nCount_RNA < 20000) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30, reduction = 'pca', verbose = FALSE)

DimPlot(pbmcsca_dataset, reduction = "umap", label = TRUE, group.by = "CellType", label.size = 3)
```

![](HLA_expression_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

That looks weird. It looks like the same cell types are all over the
UMAP with no apparent clustering. Let’s check the metadata to see what’s
going on.

``` r
unique(pbmcsca@meta.data$Experiment)
```

    ## [1] "pbmc1" "pbmc2"

``` r
unique(pbmcsca@meta.data$Method)
```

    ## [1] "Smart-seq2"          "CEL-Seq2"            "10x Chromium (v2) A"
    ## [4] "10x Chromium (v2) B" "10x Chromium (v3)"   "Drop-seq"           
    ## [7] "Seq-Well"            "inDrops"             "10x Chromium (v2)"

As we can see, `pbmcsca` dataset contains multiple experiments ran on
different scRNAseq technologies. This means, the different clusters we
see on the UMAP plot most likely caused by the batch effects introduced
by these differences. We can quickly check that on UMAP by grouping with
experiment or method.

``` r
DimPlot(pbmcsca_dataset, reduction = "umap", label = TRUE, group.by = c('Experiment',"Method"), label.size = 3)
```

<img src="HLA_expression_files/figure-gfm/unnamed-chunk-10-1.png" width="100%" />

It is clear that the clustering is mainly driven by these two
categories, rather than the cell type as we expect. One straightforward
way to deal with this is to use the `harmony` package. It works on PC
scores and tries removing the effect of batches. Here we have two
sources of batches, the experiment and the method (i.e. technology).

Let’s repeat the steps above by adding harmony.

``` r
pbmcsca_dataset <- pbmcsca
pbmcsca_dataset <- pbmcsca_dataset %>% 
  subset(subset = nCount_RNA < 20000) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunHarmony(group.by.vars = c('Experiment',"Method"), dims.use = 1:30) %>%
  RunUMAP(dims = 1:30, reduction = 'harmony', verbose = FALSE) # important to set reduction to harmony

DimPlot(pbmcsca_dataset, reduction = "umap", label = TRUE, group.by = "CellType", label.size = 3)
```

![](HLA_expression_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

That looks better. Now, for our purpose, which is looking at the
expression of HLA genes, we actually do not need this. Again,
annotations are already provided with the dataset. But if that was a
novel dataset without any annotation information, then we would need to
do a clustering and later identify which clusters are which type of
cells. Without batch correction (either using harmony or other methods),
this would be impossible.

We can now look at the HLA expression.

``` r
# some sorting to make plot look better
pbmcsca_dataset$CellType <- factor(pbmcsca_dataset$CellType,
                                            levels = rev(c("Unassigned","Megakaryocyte", "CD4+ T cell" ,"Cytotoxic T cell",
                                                           "Natural killer cell","CD16+ monocyte","CD14+ monocyte",
                                                           "Plasmacytoid dendritic cell","Dendritic cell" ,"B cell")))

DotPlot(pbmcsca_dataset, features = hla_genes, group.by = "CellType") + RotatedAxis()
```

![](HLA_expression_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Very similar to the plot from the `pbmc3k`, we again see the difference
in HLA-classII expression and low expression of HLA-DQ in some APCs.

Let’s spice it up a bit. We will look at the `ifnb` dataset which is
also PBMCs but one group is treated with interferon beta, while the
other group is control/untreated.

``` r
ifnb_dataset <- ifnb

# let's get rid of erythrocytes and megakaryocytes because very low expression of HLA in erythrocytes skews the results 
ifnb_dataset <- subset(ifnb_dataset,subset = !seurat_annotations %in% c("Eryth", "Mk"))
ifnb_dataset$seurat_annotations <- droplevels(ifnb_dataset$seurat_annotations)

ifnb_dataset <- ifnb_dataset %>% 
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunHarmony(group.by.vars = "stim", dims.use = 1:30) %>%
  RunUMAP(reduction = "harmony", dims = 1:30,reduction.name = 'UMAPharmony') %>%
  RunUMAP(reduction = "pca", dims = 1:30,reduction.name = 'UMAPpca')

p1 <- DimPlot(ifnb_dataset, reduction = "UMAPpca", label = TRUE, group.by = "seurat_annotations", label.size = 3)
p2 <- DimPlot(ifnb_dataset, reduction = "UMAPharmony", label = TRUE, group.by = "seurat_annotations", label.size = 3)
p3 <- DimPlot(ifnb_dataset, reduction = "UMAPpca", label = FALSE, group.by = "stim")
p4 <- DimPlot(ifnb_dataset, reduction = "UMAPharmony", label = FALSE, group.by = "stim")

(p1 + p2) / (p3 + p4)
```

<img src="HLA_expression_files/figure-gfm/unnamed-chunk-13-1.png" width="100%" />

Similar to the `pbmcsca` dataset, we have more than one conditions here.
Does that mean we need to control for batch effects?

In this dataset we actually expect to see differences between two groups
as one of them is treated with IFN-B. This is very nicely reflected in
the UMAP above (left panels).It looks as if there is a mirror running
across UMAPpca_2=0 and we have very similar clusters on each side.
Bottom left figure shows that these clustered are formed by the
treatment rather than cell type.

On the right panels, clustering based on treatment is no longer the case
and we see that datasets are harmonized, clusters are formed by cell
types. But this time we lose the effect of treatment. Well, we did not
lose it as we still have original RNA counts in the dataset but we just
do not see it in the UMAP figure as `harmony` deals with the batch
effects.

So, which one should we use? The rule of thumb is that with the batch
correction, we want to get rid of technical differences. For example
when we have samples from the same tissue of different healthy
individuals, we expect them to be similar. So we should opt in for
correction. But if we have samples from healthy and sick individuals, a
batch correction will obscure the real biological difference. In our
case, treatment with IFN-B leads to a biological difference, so we
should avoid using harmony.

Let’s move one with our HLA question. How does the expression look in
the combined dataset?

``` r
ifnb_dataset$seurat_annotations <- factor(ifnb_dataset$seurat_annotations,
                                            levels =rev(c("CD4 Naive T","CD4 Memory T","CD8 T","T activated","NK",
                                                      "CD16 Mono","CD14 Mono","DC","pDC", "B", "B Activated")))

DotPlot(ifnb_dataset, features = hla_genes, group.by = "seurat_annotations", scale = TRUE) + 
  RotatedAxis()+ggtitle("MHC gene expression by cell type (complete dataset)")
```

<img src="HLA_expression_files/figure-gfm/unnamed-chunk-14-1.png" width="100%" />

Nothing out of place. Constitutive expression of HLA classI in all cells
and expression of HLA classII in APCs. But this is the combined dataset.
We can look at it again, separated based on the IFN-B stimulation.

``` r
DotPlot(ifnb_dataset, features = hla_genes, group.by = "seurat_annotations", split.by = "stim",scale = TRUE) + 
  RotatedAxis()+ggtitle("MHC gene expression by cell type (CTRL vs IFN-B stimulated)")
```

<img src="HLA_expression_files/figure-gfm/unnamed-chunk-15-1.png" width="100%" />

Well, I don’t know what you think but when I saw the above plot for the
first time, I said “WHOAA!”. It looks as if IFN-B stimulation immensely
increases the expression of each gene in every cell, because light gray
is low expression and blue is high expression, right?

No. When `split.by` is provided to the `DotPlot`, it assigns colors to
each group separately. So what we see is actually shades of gray
representing differences for the IFN-B stimulation group and shades of
blue representing differences for the control group. Therefore, in the
figure above, comparison between groups is not possible.

We can make this possible by providing a continuous color scale to the
`DotPlot` function. Pick your favorite divergent scale from
`RColorBrewer::brewer.pal.info`

``` r
DotPlot(ifnb_dataset, features = c(hla1,hla2), group.by = "seurat_annotations", 
              split.by = "stim",scale = TRUE,cols = 'RdGy') + 
  RotatedAxis()+ggtitle("MHC genes by cell type (CTRL vs IFN-B)")
```

<img src="HLA_expression_files/figure-gfm/unnamed-chunk-16-1.png" width="100%" />

Now the color of circles are comparable to each other.

If you look carefully, an interesting pattern appears. For almost every
cell type, IFN-B stimulation increases the expression of HLA-classI
molecules. This is a well-known effect of type-I interferons which
include IFN-B as potent antiviral molecules. High expression of
HLA-classI increases the peptide presentation on infected cells, hence
increasing the chance of cytotoxic T-cell activation.

IFN-B’s effect on HLA-classII molecules appears more subtle. In fact,
some of the APCs appear to have slightly decreased expression after
stimulation. We can look at it in different ways.

For example, we can try a volcano plot. Normally, when running formal
tests of differential expression, this is done for all genes but we will
now limit our plot to HLA genes.

``` r
# first, define the cell type that we want to look at
ifnb_dataset$celltype_stim <- paste(ifnb_dataset$seurat_annotations,ifnb_dataset$stim, sep = "_")
Idents(ifnb_dataset) <- "celltype_stim"

cell_type <- "B"
mono.de <- FindMarkers(ifnb_dataset, 
                       ident.1 = paste0(cell_type, "_STIM"), 
                       ident.2 = paste0(cell_type, "_CTRL"), 
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
```

![](HLA_expression_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

We can confirm the dotplot that in B cells HLA-classI genes are
upregulated after IFN-B treatment and HLA-DQ seems downregulated. Let’s
calculate log2FC values for all cell types and make a heatmap to get a
better overview.

``` r
cell_types <- levels(ifnb_dataset$seurat_annotations)
de_results_list <- list()

for (cell_type in cell_types) {
  # first, we will subset the dataset based on the cell type. So out Idents should be annotations without the treatment information
  subset_seu <- ifnb_dataset
  Idents(subset_seu) <- subset_seu$seurat_annotations
  subset_seu <- subset(subset_seu, idents = cell_type)
  # now we switch back to treatment information and run the DE analysis
  Idents(subset_seu) <- 'stim'

  de_genes <- FindMarkers(subset_seu,
                          ident.1 = "STIM",
                          ident.2 = "CTRL",
                          features = hla_genes,
                          test.use = "MAST",
                          logfc.threshold=0,min.cells.group=0,min.pct=0)
  
  de_genes$cell_type <- cell_type
  de_genes$gene <- rownames(de_genes)
  de_results_list[[cell_type]] <- de_genes
  }

combined_de_results <- bind_rows(de_results_list)
rownames(combined_de_results) <- NULL

ggplot(combined_de_results, aes(x = gene, y = cell_type, fill = avg_log2FC)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "lightgray", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "HLA Gene", y = "Cell Type", fill = "Log2FC\n(STIM vs CTRL)",
       title = "Differential Expression of HLA Genes after IFN-B Stimulation")
```

![](HLA_expression_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

Now the difference between class1 and class2 genes is more clear and the
upregulating effect of IFN-B on class1 across most cell types stands
out, whereas class2 changes are more subtle and cell-type specific.
