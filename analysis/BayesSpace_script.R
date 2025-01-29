# install BayeSpace
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BayesSpace")

# load packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(scater)
  library(harmony)
  library(BayesSpace)
  library(SingleCellExperiment)
  library(Matrix)
  library(Seurat)
})

options(matrixStats.useNames.NA = "deprecated")

# custom color palette
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

# read in Seruat metadata
meta <- read.csv(file = "path/to/ST_metadata.csv")
filtered_spot_ids <- meta$cells # Get the filtered spot IDs from metadata; enable direct comparison with Seurat clusters

# Code for importing raw data from SpaceRanger output
Sample8A <- readVisium("path/to/outs/") # repeat for remaining sections

# B17012 (Sample 8, NCI)
Sample8A <- readVisium("path/to/outs/")
Sample8B <- readVisium("path/to/outs/")
Sample8C <- readVisium("path/to/outs/")
Sample8D <- readVisium("path/to/outs/")

# B16002 (Sample 15, AD)
Sample15A <- readVisium("path/to/outs/")
Sample15B <- readVisium("path/to/outs/")
Sample15C <- readVisium("path/to/outs/")
Sample15D <- readVisium("path/to/outs/")

# B17078 (Sample 6, NCI)
Sample6A <- readVisium("path/to/outs/")
Sample6B <- readVisium("path/to/outs/")
Sample6C <- readVisium("path/to/outs/")
Sample6D <- readVisium("path/to/outs/")

# B17062 (Sample 17, AD)
Sample17A <- readVisium("path/to/outs/")
Sample17B <- readVisium("path/to/outs/")
Sample17C <- readVisium("path/to/outs/")
Sample17D <- readVisium("path/to/outs/")

# filter spots so only considering same spots used for Seurat/Louvain clusters
# B17012 (Sample 8)
colnames(Sample8A) <- paste0("8_A1_", colnames(Sample8A)) # rename spot IDs in BayesSpace object to match Seurat
Sample8A <- Sample8A[, colnames(Sample8A) %in% filtered_spot_ids]
colnames(Sample8B) <- paste0("8_B1_", colnames(Sample8B)) # rename spot IDs in BayesSpace object to match
Sample8B <- Sample8B[, colnames(Sample8B) %in% filtered_spot_ids]
colnames(Sample8C) <- paste0("8_C1_", colnames(Sample8C)) # rename spot IDs in BayesSpace object to match
Sample8C <- Sample8C[, colnames(Sample8C) %in% filtered_spot_ids]
colnames(Sample8D) <- paste0("8_D1_", colnames(Sample8D)) # rename spot IDs in BayesSpace object to match
Sample8D <- Sample8D[, colnames(Sample8D) %in% filtered_spot_ids]
# B16002 (Sample 15)
colnames(Sample15A) <- paste0("15_A1_", colnames(Sample15A)) # rename spot IDs in BayesSpace object to match Seurat
Sample15A <- Sample15A[, colnames(Sample15A) %in% filtered_spot_ids]
colnames(Sample15B) <- paste0("15_B1_", colnames(Sample15B)) # rename spot IDs in BayesSpace object to match
Sample15B <- Sample15B[, colnames(Sample15B) %in% filtered_spot_ids]
colnames(Sample15C) <- paste0("15_C1_", colnames(Sample15C)) # rename spot IDs in BayesSpace object to match
Sample15C <- Sample15C[, colnames(Sample15C) %in% filtered_spot_ids]
colnames(Sample15D) <- paste0("15_D1_", colnames(Sample15D)) # rename spot IDs in BayesSpace object to match
Sample15D <- Sample15D[, colnames(Sample15D) %in% filtered_spot_ids]
# B17078 (Sample 6)
colnames(Sample6A) <- paste0("6_A1_", colnames(Sample6A)) # rename spot IDs in BayesSpace object to match Seurat
Sample6A <- Sample6A[, colnames(Sample6A) %in% filtered_spot_ids]
colnames(Sample6B) <- paste0("6_B1_", colnames(Sample6B)) # rename spot IDs in BayesSpace object to match
Sample6B <- Sample6B[, colnames(Sample6B) %in% filtered_spot_ids]
colnames(Sample6C) <- paste0("6_C1_", colnames(Sample6C)) # rename spot IDs in BayesSpace object to match
Sample6C <- Sample6C[, colnames(Sample6C) %in% filtered_spot_ids]
colnames(Sample6D) <- paste0("6_D1_", colnames(Sample6D)) # rename spot IDs in BayesSpace object to match
Sample6D <- Sample6D[, colnames(Sample6D) %in% filtered_spot_ids]
# B17062 (Sample 17)
colnames(Sample17A) <- paste0("17_A1_", colnames(Sample17A)) # rename spot IDs in BayesSpace object to match Seurat
Sample17A <- Sample17A[, colnames(Sample17A) %in% filtered_spot_ids]
colnames(Sample17B) <- paste0("17_B1_", colnames(Sample17B)) # rename spot IDs in BayesSpace object to match
Sample17B <- Sample17B[, colnames(Sample17B) %in% filtered_spot_ids]
colnames(Sample17C) <- paste0("17_C1_", colnames(Sample17C)) # rename spot IDs in BayesSpace object to match
Sample17C <- Sample17C[, colnames(Sample17C) %in% filtered_spot_ids]
colnames(Sample17D) <- paste0("17_D1_", colnames(Sample17D)) # rename spot IDs in BayesSpace object to match
Sample17D <- Sample17D[, colnames(Sample17D) %in% filtered_spot_ids]

# add $sample.name metadata before merging
Sample8A$sample_name <- "Sample8A"
Sample8B$sample_name <- "Sample8B"
Sample8C$sample_name <- "Sample8C"
Sample8D$sample_name <- "Sample8D"
Sample15A$sample_name <- "Sample15A"
Sample15B$sample_name <- "Sample15B"
Sample15C$sample_name <- "Sample15C"
Sample15D$sample_name <- "Sample15D"
Sample6A$sample_name <- "Sample6A"
Sample6B$sample_name <- "Sample6B"
Sample6C$sample_name <- "Sample6C"
Sample6D$sample_name <- "Sample6D"
Sample17A$sample_name <- "Sample17A"
Sample17B$sample_name <- "Sample17B"
Sample17C$sample_name <- "Sample17C"
Sample17D$sample_name <- "Sample17D"

#Combine into 1 SCE and preprocess (log-normalize and run PCA)
sce.combined = cbind(Sample8A, Sample8B, Sample8C, Sample8D, Sample15A, Sample15B, Sample15C, Sample15D, Sample6A, Sample6B, Sample6C, Sample6D, Sample17A, Sample17B, Sample17C, Sample17D, deparse.level = 1)
sce.combined = spatialPreprocess(sce.combined, n.PCs = 50) #lognormalize, PCA

# determine optimal cluster number (q)
sce.combined <- qTune(sce.combined, qs=seq(2, 10), platform="Visium", d=7)
qPlot(sce.combined) # Supplementary Fig. 3C

# Batch Correction and Clustering
sce.combined = runUMAP(sce.combined, dimred = "PCA")
colnames(reducedDim(sce.combined, "UMAP")) = c("UMAP1", "UMAP2")

ggplot(data.frame(reducedDim(sce.combined, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(sce.combined$sample_name))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()

# Add Sample metadata, want to run Harmony by sample rather than section
n8A <- ncol(Sample8A)
n8B <- ncol(Sample8B)
n8C <- ncol(Sample8C)
n8D <- ncol(Sample8D)
n15A <- ncol(Sample15A)
n15B <- ncol(Sample15B)
n15C <- ncol(Sample15C)
n15D <- ncol(Sample15D)
n6A <- ncol(Sample6A)
n6B <- ncol(Sample6B)
n6C <- ncol(Sample6C)
n6D <- ncol(Sample6D)
n17A <- ncol(Sample17A)
n17B <- ncol(Sample17B)
n17C <- ncol(Sample17C)
n17D <- ncol(Sample17D)

# Get the number of columns (cells) in each original object
n8 <- n8A + n8B + n8C + n8D
n15 <- n15A + n15B + n15C + n15D
n6 <- n6A + n6B + n6C + n6D
n17 <- n17A + n17B + n17C + n17D

# Assign the 'sample' metadata column based on column indices
sce.combined$sample <- rep(
  c("8", "15", "6", "17"),
  times = c(n8, n15, n6, n17)
)

# Verify the result
table(sce.combined$sample)

# Harmony batch correction by sample
sce.combined = RunHarmony(sce.combined, "sample", verbose = F)
sce.combined = runUMAP(sce.combined, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(sce.combined, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")

ggplot(data.frame(reducedDim(sce.combined, "UMAP.HARMONY")), 
       aes(x = UMAP1, y = UMAP2, color = factor(sce.combined$sample_name))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()

# Harmony batch correction by section
sce.combined = RunHarmony(sce.combined, "sample_name", verbose = F, reduction.save = "HARMONY.section")
sce.combined = runUMAP(sce.combined, dimred = "HARMONY.section", name = "UMAP.HARMONY.section")
colnames(reducedDim(sce.combined, "UMAP.HARMONY.section")) = c("UMAP1", "UMAP2")

ggplot(data.frame(reducedDim(sce.combined, "UMAP.HARMONY.section")), 
       aes(x = UMAP1, y = UMAP2, color = factor(sce.combined$sample_name))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()

# Offset sections to enable joint clustering/visualization
sce.combined$row[sce.combined$sample_name == "Sample8B"] = 
  100 + sce.combined$row[sce.combined$sample_name == "Sample8B"]
sce.combined$col[sce.combined$sample_name == "Sample8C"] = 
  150 + sce.combined$col[sce.combined$sample_name == "Sample8C"]
sce.combined$row[sce.combined$sample_name == "Sample8D"] = 
  100 + sce.combined$row[sce.combined$sample_name == "Sample8D"]
sce.combined$col[sce.combined$sample_name == "Sample8D"] = 
  150 + sce.combined$col[sce.combined$sample_name == "Sample8D"]

sce.combined$row[sce.combined$sample_name == "Sample15A"] = 
  200 + sce.combined$row[sce.combined$sample_name == "Sample15A"]
sce.combined$row[sce.combined$sample_name == "Sample15B"] = 
  300 + sce.combined$row[sce.combined$sample_name == "Sample15B"]
sce.combined$row[sce.combined$sample_name == "Sample15C"] = 
  200 + sce.combined$row[sce.combined$sample_name == "Sample15C"]
sce.combined$col[sce.combined$sample_name == "Sample15C"] = 
  150 + sce.combined$col[sce.combined$sample_name == "Sample15C"]
sce.combined$row[sce.combined$sample_name == "Sample15D"] = 
  300 + sce.combined$row[sce.combined$sample_name == "Sample15D"]
sce.combined$col[sce.combined$sample_name == "Sample15D"] = 
  150 + sce.combined$col[sce.combined$sample_name == "Sample15D"]

sce.combined$col[sce.combined$sample_name == "Sample6A"] = 
  300 + sce.combined$col[sce.combined$sample_name == "Sample6A"]
sce.combined$row[sce.combined$sample_name == "Sample6B"] = 
  100 + sce.combined$row[sce.combined$sample_name == "Sample6B"]
sce.combined$col[sce.combined$sample_name == "Sample6B"] = 
  300 + sce.combined$col[sce.combined$sample_name == "Sample6B"]
sce.combined$col[sce.combined$sample_name == "Sample6C"] = 
  450 + sce.combined$col[sce.combined$sample_name == "Sample6C"]
sce.combined$row[sce.combined$sample_name == "Sample6D"] = 
  100 + sce.combined$row[sce.combined$sample_name == "Sample6D"]
sce.combined$col[sce.combined$sample_name == "Sample6D"] = 
  450 + sce.combined$col[sce.combined$sample_name == "Sample6D"]

sce.combined$row[sce.combined$sample_name == "Sample17A"] = 
  200 + sce.combined$row[sce.combined$sample_name == "Sample17A"]
sce.combined$col[sce.combined$sample_name == "Sample17A"] = 
  300 + sce.combined$col[sce.combined$sample_name == "Sample17A"]
sce.combined$row[sce.combined$sample_name == "Sample17B"] = 
  300 + sce.combined$row[sce.combined$sample_name == "Sample17B"]
sce.combined$col[sce.combined$sample_name == "Sample17B"] = 
  300 + sce.combined$col[sce.combined$sample_name == "Sample17B"]
sce.combined$row[sce.combined$sample_name == "Sample17C"] = 
  200 + sce.combined$row[sce.combined$sample_name == "Sample17C"]
sce.combined$col[sce.combined$sample_name == "Sample17C"] = 
  450 + sce.combined$col[sce.combined$sample_name == "Sample17C"]
sce.combined$row[sce.combined$sample_name == "Sample17D"] = 
  300 + sce.combined$row[sce.combined$sample_name == "Sample17D"]
sce.combined$col[sce.combined$sample_name == "Sample17D"] = 
  450 + sce.combined$col[sce.combined$sample_name == "Sample17D"]

clusterPlot(sce.combined, "sample_name", color = NA) + #make sure no overlap between samples
  labs(fill = "Sample", title = "Offset check")

# Cluster harmonized offset sections, test q=8 or q=10
sce.combined_q8 = spatialCluster(sce.combined, use.dimred = "HARMONY", q = 8, nrep = 10000) #use HARMONY; started @ ~12:00; takes several hours for joint 16 sect
clusterPlot(sce.combined_q8, color = NA, palette = palette_8) + #plot clusters
  labs(title = "BayesSpace joint clustering q8")

# read joint BayesSpace object
sce.combined <- readRDS(file = "Joint_16sect_q8_BayesSpace.rds")
palette_10 <- c("#F8766D", "#00BFC4", "#A3A500", "#619CFF", "#E49B36", "#00BF7D", "#C77CFF", "#A05195", "#FF9DA7", "#5A8F29")

# check ?spatialCluster should be able to store q8 and q10 in the same object under different names
sce.combined_q10 = spatialCluster(sce.combined, use.dimred = "HARMONY", q = 10, nrep = 10000) #use HARMONY

clusterPlot(sce.combined, color = NA, palette = palette_10) + #plot clusters
  labs(title = "BayesSpace joint clustering")

sce.combined$spatialCluster_q8 <- sce.combined$cluster.init
sce.combined$spatialCluster_q10 <- sce.combined_q10$cluster.init

# save/load joint BayesSpace object
saveRDS(sce.combined, file = "Joint_16sect_q8_q10_BayesSpace.rds")
sce.combined <- readRDS(file = "Joint_16sect_q8_q10_BayesSpace.rds")

# Ensure matching spot IDs
spot_ids <- sce.combined@colData@rownames
seurat_df <- subset_samples[subset_samples$cells %in% spot_ids, ] # subsetted Seurat object

# Spot IDs in sce.combined but not in seurat_df
missing_spots <- setdiff(sce.combined@colData@rownames, seurat_df$cells)
length(missing_spots)  # Should match the number of NA values
print(missing_spots)  # Inspect the IDs

# Match Seurat clusters to spot IDs
seurat_clusters <- seurat_df$annotated_clusters[match(spot_ids, seurat_df$cells)]

# Assign "Unknown" to unmatched spots (NA values)
seurat_clusters[is.na(seurat_clusters)] <- "Unknown"

# Add the clusters to the BayesSpace object
sce.combined$seurat_clusters <- factor(
  x = seurat_clusters,
  levels = c(unique(seurat_clusters))  # Automatically include "Unknown"
)

# Re-order clusters as factor
sce.combined$seurat_clusters <- factor(x = sce.combined$seurat_clusters, levels = c("L1/Astrocytes","L2-3","L3-5","L5","L1-5 Low UMI","L6","WM","meninges","Blood Vessel","Interneuron","Mixed Glia","Unknown"))

# convert back to character
sce.combined$seurat_clusters <- as.character(sce.combined$seurat_clusters)

# Ensure spot IDs are present in sce.combined@colData@rownames
spot_ids <- rownames(sce.combined@colData)

# Add rownames to seurat_clusters
names(sce.combined$seurat_clusters) <- spot_ids

#custom colors
cluster_colors <- c(
  "L5" = "#0A7B83",
  "L3-5" = "#7ABAF2",
  "WM" = "#FFD265",
  "L1-5 Low UMI" = "#595959",
  "L2-3" = "#2AA876",
  "L6" = "#F19C65",
  "L1/Astrocytes" = "#8B63A6",
  "meninges" = "#607D8B",
  "Blood Vessel" = "#ff5252",
  "Interneuron" = "#4C1B1B",
  "Mixed Glia" = "#fb9a99",
  "Unknown" = "gray"  # Add a color for "Unknown"
)

# to plot, temporarily replace $spatial.cluster, which is pulled by default for clusterPlot, don't save new BayesSpace object!
# for Seurat clusters; Supplementary Fig. 3A
sce.combined$spatial.cluster <- sce.combined$seurat_clusters

clusterPlot(sce.combined, color = NA, palette = cluster_colors) + #plot clusters
  labs(title = "Seurat Louvain clusters")
# change back
sce.combined$spatial.cluster <- sce.combined$spatialCluster_q8

# for q8 clusters; Supplementary Fig. 3D
clusterPlot(sce.combined, color = NA, palette = palette_8) + #plot clusters
  labs(title = "BayesSpace q8 clusters")


# plot subsetted seurat clusters on BayesSpace clusterplot; Supplementary Fig. 3B
# Ensure matching spot IDs
spot_ids <- sce.combined@colData@rownames

spots_to_keep <- rownames(subset_samples@meta.data)[subset_samples@meta.data$cells %in% spot_ids]
seurat_df_2 <- subset_samples[, spots_to_keep]

# Spot IDs in sce.combined but not in seurat_df
missing_spots <- setdiff(sce.combined@colData@rownames, seurat_df_2$cells)
length(missing_spots)  # Should match the number of NA values
print(missing_spots)  # Inspect the IDs


# Match Seurat clusters to spot IDs
seurat_clusters <- seurat_df_2$seurat_clusters[match(spot_ids, seurat_df_2$cells)]

# Assign "Unknown" to unmatched spots (NA values)
seurat_clusters[is.na(seurat_clusters)] <- "Unknown"

# Add the clusters to the BayesSpace object
sce.combined$subset_seurat_clusters <- factor(
  x = seurat_clusters,
  levels = c(unique(seurat_clusters))  # Automatically include "Unknown"
)

# Re-order clusters as factor
sce.combined$subset_seurat_clusters <- factor(x = sce.combined$subset_seurat_clusters, levels = c("0","1","2","3","4","5","6","7","8","9","NA"))

# convert back to character
sce.combined$subset_seurat_clusters <- as.character(sce.combined$subset_seurat_clusters)

# temporarily replace $spatial.cluster, which is pulled by default for clusterPlot, don't save new BayesSpace object!
sce.combined$spatial.cluster <- sce.combined$seurat_clusters

# custom colors (roughly matching all-samples cluster colors based on manual inspection of cluster marker genes)
subset_colors <- c("#2AA876","#7ABAF2","#0A7B83","#FFD265","#F19C65","#8B63A6","#ff5252","#607D8B","#4C1B1B","#595959","lightgray")

clusterPlot(sce.combined, color = NA, palette = subset_colors) + #plot clusters
  labs(title = "Seurat Louvain clusters (downsampled spots)")
