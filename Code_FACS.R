


library(devtools)
source_url("https://github.com/DimitriAnderfuhren/hdFunctions/blob/master/flowFunctions.R?raw=TRUE")


# Adapted workflow from Nuñez et al Nat Imm 2023

# Set your working directory with raw fcs files, metadata and panel file
setwd("/Users/nico/Desktop/TPFACS")

# Check if you are in the correct folder
getwd()
list.files()

# Read metadata (xlsx or csv or any other table format)
metadata_filename = "matadata.xlsx"
library(readxl)

md = read_excel(metadata_filename)
View(md)

# Make condition variables to factor variables
md$Condition = factor(md$Condition, levels = c("SINO_SINO", "MOD_MOD"))


# Define colors for conditions (The number of colors should correspond to n conditions)
color_conditions = c("Red","Black") 
names(color_conditions) = levels(md$Condition)

# Read your raw fcs files
## try http:// if https:// URLs are not supported
#BiocManager::install("flowCore")
library(flowCore)
fcs_raw = read.flowSet(md$file_name, transformation = FALSE, truncate_max_range = FALSE)

# Rename and Replace problematic characters
channels <- as.data.frame(pData(parameters(fcs_raw[[1]])))
channels <- channels$desc
colnames(fcs_raw) <-channels

# Obtén los nombres de las columnas actuales
current_colnames <- colnames(fcs_raw)

# Reemplaza los caracteres problemáticos
new_colnames <- gsub("-", "_", current_colnames)
new_colnames <- gsub(" ", "_", new_colnames)

# Asigna los nuevos nombres de columna al objeto fcs_raw
colnames(fcs_raw) <- new_colnames

# Extract expression matrix
fcs = fsApply(fcs_raw, exprs)


### Constants
minCofactor = 0
maxCofactor = 10000

color1 = "cadetblue"
color2 = "chocolate"
head(fcs)

### Select the columns to use
new_colnames1<- colnames(fcs_raw[, c(7:16, 18:42)])

cofactorValues = createCofactorVector(fcs, colsToUse = new_colnames1)


library(manipulate)
### interactive plot: lets add 4660
manipulate(interactiveScatterPlot(dataFrame = fcs,colsToUse = new_colnames1,
                                  markerX = X,
                                  markerY = Y,
                                  asinh_bool = asinh,
                                  histogram_bool = histogram,
                                  cofactor_int = cofactor,
                                  sampleSize = 5000,
                                  saveValue_bool = saveButton,
                                  transformAll_bool = transformAllButton,
                                  clearAll_bool = clearAllButton),
           Y = picker(as.list(new_colnames1)),
           X = picker(as.list(new_colnames1)),
           histogram = checkbox(initial = F,label = "Histogram"),
           cofactor = slider(min = minCofactor,max = maxCofactor,step = 10,label = "Cofactor"),
           asinh = checkbox(initial = F,label = "Asinh"),
           saveButton = button(label = "Save value"),
           transformAllButton = button(label = "Transform all"),
           clearAllButton = button(label = "Clear all"))
print(cofactorValues)
write.csv(cofactorValues, file = "saved_cofactors.csv")

### In case that we would like to modify the cofactors manually, 
###they can be read and changed in the file "saved_cofactors.csv"

#cf = read.csv("saved_cofactors.csv")
### interactive plot
#cofactorValues = cofactorValues[1,]
#cofactorValues
#names(cofactorValues) = cf[,1]

# Save cofactor values
### Actual Transformation and normalization of the dataframe
fcs1 = asinhTransform(fcs,cofactorValues)
### Ensure that these are not large values similar to the flowJO fluorescence values
head(fcs1)

# Adjust data to start from (approximately) 0
q.vector <- apply(fcs1[,new_colnames1], 2, function(x) quantile(x, 0.01, names = F))

data.shift <- fcs1[,new_colnames1]
data.shift[,new_colnames1] <- sweep(data.shift[,new_colnames1], 2, q.vector)
head(data.shift)


# Normalize to have everything from 0 to 1 
per.vector <- apply(data.shift[,new_colnames1], 2, function(x) quantile(x, 0.9995, names = F))

data.shift[,new_colnames1] <- t(t(data.shift[,new_colnames1]) / as.numeric(per.vector))
head(data.shift)

### Check successful transformation for all markers via biaxial plot, Example "CD4" vs "CD8" -> name axis
library(ggplot2)
ggplot(data = data.frame(data.shift[1:10000,]), aes(x =CD3, y = CD45)) + geom_point(alpha=1/2) + ylim(-0.05,1) + xlim(-0.05,1)


# Generate sample IDs corresponding to each cell per sample in the 'expr' matrix
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))

#delete duplicates (events coming from FlowJo export) from data.shift, the normalized 
dups <- which(!duplicated(data.shift[, new_colnames1]))

## Data subsampling: create indices by sample
inds <- split(1:length(sample_ids), sample_ids)

## How many cells to downsample per-sample to run the histogram (plot1 to check all markers)
tsne_ncells <- pmin(table(sample_ids), 10000)

## Get subsampled indices
set.seed(1234)
tsne_inds <- lapply(names(inds), function(i){
  s <- sample(inds[[i]], tsne_ncells[i], replace = FALSE)
  intersect(s, dups)
})

tsne_inds <- unlist(tsne_inds)

#function before the histogram: to track the single cells in each samples (P01, P02, ...)

sample_ids <- data.frame(rep(md$sample_id, fsApply(fcs_raw, nrow)))

ggdf <- data.frame(sample_id = sample_ids[tsne_inds,], data.shift[tsne_inds, new_colnames1])
library(reshape2)
ggdf <- melt(ggdf, id.var = "sample_id",
             value.name = "expression", variable.name = "antigen")
#add the variable sample ID

mm <- match(ggdf$sample_id, md$sample_id)

ggdf$condition <- md$Condition [mm]


#### Plot the histogram

plot1 <- ggplot(ggdf, aes(x = expression, color = condition, group = sample_id)) +
  geom_density() +
  facet_wrap(~ antigen, nrow = 4, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
  guides(color = guide_legend(ncol = 1)) +
  scale_color_manual(values = color_conditions) + xlim(-0.2,1.2)
plot1


## Data subsampling: create indices by sample
expr <- as.data.frame(data.shift)
sample_ids <- rep(md$sample_id, fsApply(fcs_raw, nrow))

expr$cell.id <- as.factor(1:nrow(expr))

tsne_expr <- expr[tsne_inds, new_colnames1]

library(umap)
data_umap <-  umap(tsne_expr, random_state=123, verbose =T)

Umap <- as.data.frame(data_umap$layout)
colnames(Umap) <- c("Umap1", "Umap2")


## Plot t-SNE colored by single expression expression instead new_colnames1, 
## Select the markers you wish to include in the UMAP plot.
dr_umap <- data.frame(Umap, expr[tsne_inds, new_colnames1])

## UMAP color palette

actual2 <- c('#2F2C62', '#42399B', '#4A52A7', '#A7DA64',
             '#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131')


# prepare the expression data


data.ix.df <- data.frame(expr[tsne_inds,c(new_colnames1, "cell.id")])

data.melt <- melt(data.ix.df, variable.name = "antigen", value.name = "expression")

dr_umap$cell.id <-expr[tsne_inds,c("cell.id")]

# bind with last number of colunms

dr_umap2 <-dr_umap[,c("Umap1", "Umap2", "cell.id")]

joined.expr_umap <- merge(data.melt, dr_umap2, by = "cell.id")


#or


p7_umap <- ggplot(joined.expr_umap, aes(x = Umap1, y = Umap2, color = expression)) +
  geom_point(size = 0.04) +
  scale_color_gradientn(colours = actual2,
                        limits = c(-0, 1)) +  
  facet_wrap(~ antigen, ncol = 6, scales = "free") 

p7_umap


## FlowSOM
BiocManager::install("FlowSOM")
library(FlowSOM)

fsom <- ReadInput(flowFrame(data.shift, desc = list(FIL = 1)), transform = FALSE, scale = FALSE)

set.seed(1234)
new_colnames1
lineage_markers = c("CD1C","CD16","CD3","CD8","CD19","TCRGD", "CD56", "CD4", "CD14", "CD123")

som <- BuildSOM(fsom, colsToUse = lineage_markers)

## Metaclustering into 40 clusters with ConsensusClusterPlus

library(ConsensusClusterPlus)
codes <- som$map$codes
plot_outdir <- "consensus_plots"


nmc <- 30
mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100,
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png", clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average", distance = "euclidean", seed = 1234)


## Get cluster ids for each cell
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[som$map$mapping[,1]]


color_clusters <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00",
                    "#8B3800", "#0a0d05", "#4E8500", "#2F2C62", "#D70131",
                    "#E64B35B2","#4DBBD5B2", "#00A087B2","#3C5488B2", "#F39B7FB2", 
                    "#4E8500", "#4C00FFFF", "#000DFFFF", "#0068FFFF", "#00C1FFFF", 
                    "#00FF24FF", "#42FF00FF", "#A8FF00FF", "#FFF90CFF", "#FFE247FF", 
                    "#FFDB83FF", "#4a5b67","#91D1C2B2","#DC0000B2","#7E6148B2",
                    "#B09C85B2","#15b2d3","#236e96","#ffd700","#f3872f","#ff598f",
                    "#4E8500", "#4C00FFFF", "#000DFFFF", "#0068FFFF", "#00C1FFFF", 
                    "#00FF24FF", "#42FF00FF", "#A8FF00FF", "#FFF90CFF", "#FFE247FF", 
                    "#FFDB83FF", "#4a5b67","#91D1C2B2","#DC0000B2","#7E6148B2",
                    "#B09C85B2","#15b2d3","#236e96","#ffd700","#f3872f","#ff598f")

library(dplyr)
library(RColorBrewer)
library(pheatmap)

plot_clustering_heatmap_wrapper <- function(fcs1, data.shift, 
                                            cell_clustering, color_clusters, cluster_merging = NULL) {
  
  # Calculate the median expression
  expr_median <- data.frame(fcs1, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_each(funs(median))
  data.shift_median <- data.frame(data.shift, cell_clustering = cell_clustering) %>%
    group_by(cell_clustering) %>%
    summarize_each(funs(median))
  
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering))
  
  # This clustering is based on the markers that were used for the main clustering
  d <- dist(expr_median[, colnames(fcs1)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(data.shift_median[, colnames(data.shift)])
  rownames(expr_heat) <- data.shift_median$cell_clustering
  
  labels_row <- paste0(rownames(expr_heat), " (",
                       round(clustering_table / sum(clustering_table) * 100, 2), "%)")
  labels_col <- colnames(expr_heat)
  
  # Row annotation for the heatmap
  annotation_row <- data.frame(cluster = factor(data.shift_median$cell_clustering)) 
  rownames(annotation_row) <- rownames(expr_heat)
  
  color_clusters <- color_clusters[1:nlevels(annotation_row$cluster)]
  names(color_clusters) <- levels(annotation_row$cluster)
  annotation_colors <- list(cluster = color_clusters)
  annotation_legend <- FALSE
  
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$cluster_merging <- cluster_merging$new_cluster
    color_clusters <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters) <- levels(cluster_merging$new_cluster)
    annotation_colors$cluster_merging <- color_clusters
    annotation_legend <- TRUE
  }  
  
  # Colors for the heatmap
  color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  pheatmap(expr_heat, color = color,
           cluster_cols = FALSE, cluster_rows = cluster_rows,
           labels_col = labels_col, labels_row = labels_row,
           display_numbers = TRUE, number_color = "black",
           fontsize = 8, fontsize_number = 0.001,
           annotation_row = annotation_row, annotation_colors = annotation_colors,
           annotation_legend = annotation_legend)
}

plot_clustering_heatmap_wrapper(fcs1 = fcs1[, lineage_markers],
                                data.shift = data.shift[, lineage_markers],
                                cell_clustering = cell_clustering1, color_clusters = color_clusters)


## Plot FlowSOM 

dr_umap$sample_id <- sample_ids[tsne_inds]

#or
mm <- match(dr_umap$sample_id, md$sample_id)
dr_umap$condition <- md$Condition[mm]
dr_umap$cell_clustering1 <- factor(cell_clustering1[tsne_inds], levels = 1:nmc)

ggp_umap <-ggplot(dr_umap, aes(x = Umap1, y = Umap2,color = cell_clustering1)) +
  geom_density_2d(data = dr_umap[,c(1,2)], aes(x = Umap1, y = Umap2), colour = "lightgrey", size =0.3, bins=10) +
  geom_point(size = 0.02) + coord_fixed(ratio = 1)+  
  ylim(-12,13) + xlim(-12,13) +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))
ggp_umap + facet_wrap(~ condition)


###################################################################################


## Manual metaclustering to annotate all the clusters manually to show the general overview of the main populations
#Load the content of the cluster_mergings.xlsx file into R (see Nowicka et al for description of "cluster_mergings")

cluster_merging_filename <- "cluster_mergings.xlsx" 
cluster_merging1 <- read_excel(cluster_merging_filename) 
data.frame(cluster_merging1)
mm <- match(cell_clustering1, cluster_merging1$original_cluster)
cell_clustering1m <- cluster_merging1$new_cluster[mm]

# You update the t-SNE plot with the new annotated cell populations.
dr_umap$cell_clustering1m <- factor(cell_clustering1m[tsne_inds])


############### 


plot5 <-ggplot(dr_umap, aes(x = Umap1, y = Umap2,color = cell_clustering1m)) +
  geom_density_2d(data = dr_umap[,c(1,2)], aes(x = Umap1, y = Umap2), colour = "lightgrey", size =0.2, bins=10) +
  geom_point(size = 0.1) + coord_fixed(ratio = 1)+  
  ylim(-12,15) + xlim(-12,12) +
  scale_color_manual(values = color_clusters) +
  guides(color = guide_legend(override.aes = list(size = 4)))
plot5

# Facet per condition
plot5 + facet_wrap(~ condition)

#towards the final heatmap (black and white) ---> better visualization????????
dr_sub <-dr_umap[,c(lineage_markers, "cell_clustering1m")]
head(dr_sub)

median_subset <- dr_sub %>%
  group_by(cell_clustering1m) %>%
  summarize_all(funs(median(., na.rm=TRUE)))

class(median_subset)
median_subset <- as.data.frame(median_subset)
median_subset_2 <- as.matrix(sapply(median_subset[, -1], as.numeric))
rownames(median_subset_2) <- median_subset$cell_clustering1m
class(median_subset_2)

# Colors for the heatmap
color <- colorRampPalette(brewer.pal(n = 9, name = "Blues"))(10)

pheatmap(median_subset_2, color = color,
         number_color = "black", fontsize_number = 5, clustering_method = "average",
         cluster_cols = FALSE,
         border_color = "black",fontsize = 10,
         cellwidth = 20, cellheight = 15,
         display_numbers = matrix(ifelse(median_subset_2 > 5, "*", ""), nrow(median_subset_2)))

#dev.off()

# towards plotting the frequencies of each population per sample
counts_table <- table(cell_clustering1m, sample_ids)
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

ggdf <- melt(data.frame(cluster = rownames(props), props),
             id.vars = "cluster", value.name = "Frequency", variable.name = "sample_id")

# Add condition info
mm <- match(ggdf$sample_id, md$sample_id)
ggdf$condition <- factor(md$Condition[mm])

ggplot(ggdf, aes(x = sample_id, y = Frequency, fill = cluster)) + geom_bar(stat = "identity") +
  facet_wrap(~ condition, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_manual(values = color_clusters)

ggdf$sample_id <- factor(md$sample_id[mm])

library(ggpubr)

# Boxplot
p <- ggboxplot(ggdf, x = "condition", y = "Frequency", color = "cluster",
               add = "jitter")+
  scale_color_manual(values = color_clusters)+
  labs(title = "Sino vs MOD")+
  ylim(c(0,75))+
  facet_wrap(~cluster)
p 
