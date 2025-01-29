library(Seurat)
library(harmony)
library(viridis)
library(ggplot2)

set.seed(123)
options(future.globals.maxSize = 10 * 1024^3) # Set limit to 10 GB

# starting from 21 individual Seurat objects, merge into complete dataset
object_list <- mget(paste0("Object", 1:21)) # Dynamically retrieve objects
merged_ST <- do.call(merge, c(object_list, list(add.cell.ids = as.character(1:21))))

# remove spots with less than 500 genes detected
filtered_seurat <- subset(x = merged_ST, subset= (nGene >= 500))

# apply SCTransform normalization to samples
all_samples.list <- SplitObject(filtered_seurat, split.by = "ID")
all_samples.list <- lapply(all_samples.list, SCTransform, assay = "Spatial", verbose = FALSE)
all_samples <- do.call(merge, c(all_samples.list, list(add.cell.ids = as.character(1:21))))

# add metadata to Seurat object (e.g. library batch, RIN, age, fluorescence IF intensities, etc.)
meta <- read.csv(file = "ST_metadata.csv") # read in CSV file containing per-spot metadata values, rows are spotIDs, columns are metadata values

# harmonize samples
features  <- SelectIntegrationFeatures(object.list =all_samples.list, nfeatures = 3000)
all_samples <- RunPCA(object =all_samples, features = features, assay = "SCT", npcs = 50)
ElbowPlot(all_samples, ndims = 50, reduction = "pca")
all_samples <- RunHarmony(object =all_samples, assay.use = "SCT", reduction = "pca", dims.use = 1:50, group.by.vars = c("ID", "librarybatch", "RIN", "slide_nr", "cDNAbatch"), reduction.save = "harmony")
ElbowPlot(all_samples, ndims = 20, reduction = "harmony")
all_samples <- FindNeighbors(object =all_samples, assay = "SCT", reduction = "harmony", dims = 1:11)
all_samples <- FindClusters(object =all_samples, resolution = 0.3)
all_samples <- RunUMAP(object =all_samples, assay = "SCT", reduction = "harmony", dims = 1:10)
save(all_samples, file = "Seurat_object.rds")
# load("Seurat_object.rds")

# annotate clusters
levels(all_samples@active.ident)  # view current annotations
new.cluster.ids <- c(
  "L5", "L3-5", "WM", "L1-5 Low UMI", "L2-3", "L6", 
  "L1/Astrocytes", "Meninges", "Blood Vessel", "Interneuron", "Mixed glia")

names(new.cluster.ids) <- levels(all_samples)
all_samples <- RenameIdents(all_samples, new.cluster.ids)

# Set the desired factor levels
all_samples@active.ident <- factor(all_samples@active.ident, 
                                   levels = c("L1/Astrocytes", "L2-3", "L3-5", "L5", 
                                              "L1-5 Low UMI", "L6", "WM", "Meninges", 
                                              "Blood Vessel", "Interneuron", "Mixed Glia"))

# identifty cluster-enriched genes (Table S2)
cluster_names <- levels(all_samples@active.ident)
marker_list <- lapply(cluster_names, function(cluster) {
  FindMarkers(all_samples, ident.1 = cluster, min.pct = 0.25, only.pos = TRUE)
})
names(marker_list) <- cluster_names

# Save all marker lists to CSV
lapply(names(marker_list), function(cluster) {
  write.csv(marker_list[[cluster]], paste0(cluster, "_markers.csv"))
})

# custom color palettes
cluster_colors <- c(
  "L5" = "#0A7B83", "L3-5" = "#7ABAF2", "WM" = "#FFD265", 
  "L1-5 Low UMI" = "#595959", "L2-3" = "#2AA876", "L6" = "#F19C65", 
  "L1/Astrocytes" = "#8B63A6", "Meninges" = "#607D8B", 
  "Blood Vessel" = "#ff5252", "Interneuron" = "#4C1B1B", 
  "Mixed Glia" = "#fb9a99"
)
cluster_colors_AB <- c("#ffffff","#E8C51E","#364eA1")
cluster_colors_RGN_comb <- c("#d7301f","#0570b0")

# plot UMAP (Figure 1B)
DimPlot(all_samples, cols = cluster_colors, reduction = "umap", label = TRUE, label.size = 7, label.box = T, label.color = "white", shuffle=TRUE, repel = T) + NoLegend()

# Figure 1C; replace "section_image" with name of desired section
SpatialDimPlot(all_samples, images = "section_image", label = F, label.size = 10, crop = FALSE) + scale_fill_manual(values = cluster_colors)

# Figure 1D
SpatialFeaturePlot(all_samples, images= "section_image", features = c("SNAP25","MBP","PCP4"), alpha = c(0.2,1))

# plot UMAP by AD diagnosis, nGene, nUMI, or mito% (Supplementary Fig 1D-G)
DimPlot(all_samples,reduction = "umap", group.by="ad_reagan",label = FALSE,  shuffle=TRUE, repel = T)
FeaturePlot(all_samples,features = "nGene",cols=viridis(10000))
FeaturePlot(all_samples,features = "log10_nUMI",cols=viridis(10000))
FeaturePlot(all_samples,features = "mitoRatio", cols=magma(256,begin=0,end=60))

# plot UMAP cluster distribution by case or AD diagnosis (Supplementary Fig 1H)

# by individual ID
pct <- table(Idents(all_samples), all_samples$ID)
pct <- as.data.frame(pct)
pct$Var1 <- as.character(pct$Var1)

ggplot(pct, aes(x = Var2, y = Freq, fill = fct_inorder(Var1))) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") +
  scale_fill_manual(values = cluster_colors) +
  theme(legend.title = element_blank())

# by AD diagnosis
pctAD <- table(Idents(all_samples), all_samples$ad_reagan)
pctAD <- as.data.frame(pctAD)
pctAD$Var1 <- as.character(pctAD$Var1)

ggplot(pctAD, aes(x = Var2, y = Freq, fill = fct_inorder(Var1))) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("AD status") +
  ylab("Proportion") +
  scale_fill_manual(values = cluster_colors) +
  theme(legend.title = element_blank())

# IF image quant
# read in CSV file containing plaque type annotations for 781 spots (column called 'strat' includes annotations)
spot_check_quant <- read.csv(file = "TableS5.csv")

# A-beta, GFAP and IBA1 violin plots (Figure 2D)
spot_check_quant$TYPE <- factor(x = spot_check_quant$TYPE, levels = c('diffuse', 'compact','dense core'))
ggplot(spot_check_quant, aes(x=TYPE, y=log_amyloid, fill = TYPE)) + geom_jitter(shape=16, position=position_jitter(0.01)) + theme_classic() + geom_violin(width=1) + geom_boxplot(outlier.shape = NA, color = "black", width = 0.2) + scale_fill_manual(values = c("#ffffff", "#808184", "#414042"))
ggplot(spot_check_quant, aes(x=TYPE, y=log_GFAP, fill = TYPE)) + geom_jitter(shape=16, position=position_jitter(0.01)) + theme_classic() + geom_violin(width=1) + geom_boxplot(outlier.shape = NA, color = "black", width = 0.2) + scale_fill_manual(values = c("#ffffff", "#808184", "#414042"))
ggplot(spot_check_quant, aes(x=TYPE, y=log_IBA1, fill = TYPE)) + geom_jitter(shape=16, position=position_jitter(0.01)) + theme_classic() + geom_violin(width=1) + geom_boxplot(outlier.shape = NA, color = "black", width = 0.2) + scale_fill_manual(values = c("#ffffff", "#808184", "#414042"))

# Figure 2E proportion bar plot
# subset by low or high A-beta for stacked bar plot
high_quant <- subset(spot_check_quant, strat == "abeta_high")
low_quant <- subset(spot_check_quant, strat == "abeta_low")
# high A-beta
pctAD <- table(high_quant$TYPE, high_quant$strat)
pctAD <- as.data.frame(pctAD)
pctAD$Var1 <- as.character(pctAD$Var1)

ggplot(pctAD, aes(x = Var2, y = Freq, fill = fct_inorder(Var1))) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("abeta strat") +
  ylab("Proportion") +
  scale_fill_manual(values = c("#ffffff", "#808184", "#414042")) +
  theme(legend.title = element_blank())

# low A-beta
pctAD <- table(low_quant$TYPE, low_quant$strat)
pctAD <- as.data.frame(pctAD)
pctAD$Var1 <- as.character(pctAD$Var1)

ggplot(pctAD, aes(x = Var2, y = Freq, fill = fct_inorder(Var1))) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("abeta strat") +
  ylab("Proportion") +
  scale_fill_manual(values = c("#ffffff", "#808184", "#414042")) +
  theme(legend.title = element_blank())

# A-beta stratification violin plot (Figure 2F)
VlnPlot(all_samples,features=c("log_amyloid_avg"), split.by = "abeta_avg_strat2", group.by = "abeta_avg_strat2", pt.size = 0, cols=cluster_colors_AB) + geom_boxplot()

# Spatial plots for Figures 3F, 4C, 5F and 6B
# UMAP clusters overlaid on representative section
SpatialDimPlot(all_samples, images = "section_image", label = F, label.size = 10, crop = FALSE, pt.size.factor = 1.25, stroke = 0) + scale_fill_manual(values = cluster_colors)
# by average log2intensity of A-beta 
SpatialFeaturePlot(all_samples, images= "section_image", features = "log_amyloid_avg", pt.size.factor = 1.25, stroke = 0)
# by A-beta stratification
SpatialPlot(all_samples, images= "section_image", group.by = "abeta_avg_strat2", pt.size.factor = 1.25, stroke = 0) + scale_fill_manual(values = cluster_colors_AB)
# by glia stratification
SpatialPlot(all_samples, images= "section_image", group.by = "RGN_comb2", pt.size.factor = 1.25, stroke = 0) + scale_fill_manual(values = cluster_colors_RGN_comb)

# subset to 4 individuals, 4 samples each using Subset
# Define the sections you want to subset
sections_to_keep <- c("6_A1", "6_B1", "6_C1", "6_D1", 
                      "15_A1", "15_B1", "15_C1", "15_D1", 
                      "8_A1", "8_B1", "8_C1", "8_D1", 
                      "17_A1", "17_B1", "17_C1", "17_D1")

# Get cell/spot names for sections to keep
spots_to_keep <- rownames(all_samples@meta.data)[all_samples@meta.data$section %in% sections_to_keep]
# subset
subset_samples <- all_samples[, spots_to_keep]
# save/load
save(subset_samples, file = "subset_samples.rds")
load("subset_samples.rds")

# re-run normalization, PCA, harmony, and UMAP on subsetted samples (Supplementary Fig. 3)
subset_samples <- SCTransform(subset_samples, assay = "Spatial", verbose = FALSE, new.assay.name = "subset_SCT")
subset_samples <- RunPCA(subset_samples, assay = "subset_SCT", npcs = 50, reduction.name = "subset_pca")
ElbowPlot(subset_samples, ndims = 50, reduction = "subset_pca")
subset_samples <- RunHarmony(subset_samples, assay.use = "subset_SCT", reduction = "subset_pca", dims.use = 1:50, group.by.vars = c("projID", "librarybatch", "RIN", "slide_nr", "cDNAbatch"), reduction.save = "subset_harmony")
ElbowPlot(subset_samples, ndims = 20, reduction = "subset_pca")
subset_samples <- FindNeighbors(subset_samples, assay = "subset_SCT", reduction = "subset_harmony", dims = 1:11)
subset_samples <- FindClusters(subset_samples, resolution = 0.3)
subset_samples <- RunUMAP(subset_samples, assay = "subset_SCT", reduction = "subset_harmony", dims = 1:10, nn.name = , reduction.name = "subset_umap")
save(subset_samples, file = "harmonized_subset_samples.rds")

# custom color palette that maps subsetted Louvain cluster colors similarly to those of original UMAP
subset_colors <- c("#2AA876","#7ABAF2","#0A7B83","#FFD265","#F19C65","#8B63A6","#ff5252","#607D8B","#4C1B1B","#595959")

# plot UMAP - Supplementary Fig. 3E
DimPlot(subset_samples, group.by = "seurat_clusters", cols = subset_colors, reduction = "umap", shuffle=TRUE, repel = T) + NoLegend() 

# custom color palette that maps BayesSpace cluster colors similarly to those of original UMAP
palette_8 <- c(
  "3" = "#0A7B83",
  "2" = "#7ABAF2",
  "1" = "#FFD265",
  "4" = "#2AA876",
  "5" = "#F19C65",
  "8" = "#8B63A6",
  "7" = "#607D8B",
  "6" = "#00BFC4"
)

# plot UMAP - Supplementary Fig. 3F
DimPlot(subset_samples, group.by = "Bayes_cluster", cols = palette_8, reduction = "umap", shuffle=TRUE, repel = T) + NoLegend()
