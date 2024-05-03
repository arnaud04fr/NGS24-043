library(tidyverse)
library(pheatmap)

Data_Norm_DEG <- read.table("./results/DESeq_Table_norm_count_DEG.tsv", header = TRUE, row.names = 1)
Data_Norm_DEG_all <-  read.table("./results/DESeq_Table_norm_count_all.tsv", header = TRUE, row.names = 1)
metadata <- read.csv2("Metadata_NGS24-043.csv")

nom <- c("888_CM", "888_ctrl", "889_CM", "889_ctrl","894_CM", "894_ctrl","897_CM", "897_ctrl","913_CM", "913_ctrl")
colnames(Data_Norm_DEG) <- nom
colnames(Data_Norm_DEG_all) <- nom
row.names(metadata) <- nom
LIF_MYF <- read.csv("LIF_MYF_signature.csv")
LIF <- LIF_MYF[LIF_MYF$signature == "LIF", ]
MYF <- LIF_MYF[LIF_MYF$signature == "MYF", ]

Data_Norm_LIF <- subset(Data_Norm_DEG, row.names(Data_Norm_DEG) %in% LIF$gene_name)
Data_Norm_MYF <- subset(Data_Norm_DEG, row.names(Data_Norm_DEG) %in% MYF$gene_name)
Data_Norm_sign <- subset(Data_Norm_DEG, row.names(Data_Norm_DEG) %in% LIF_MYF$gene_name)

Data_Norm_LIF_all <- subset(Data_Norm_DEG, row.names(Data_Norm_DEG) %in% LIF$gene_name)
Data_Norm_MYF_all <- subset(Data_Norm_DEG, row.names(Data_Norm_DEG) %in% MYF$gene_name)
Data_Norm_sign_all <- subset(Data_Norm_DEG_all, row.names(Data_Norm_DEG_all) %in% LIF_MYF$gene_name)

# Get the order of row names in LIF_MYF table
order <- match(LIF_MYF$gene_name, rownames(Data_Norm_sign_all))

# Reorder row names of Data_Norm_sign_all ans remove NA rows
table1_reordered <- Data_Norm_sign_all[order, , drop = FALSE]
table1_reordered <- na.omit(table1_reordered)

#Clustering and heatmap

#annotations des colonnes pour graph par traitement
my_sample_col <- data.frame(metadata$treatment)
row.names(my_sample_col) <- row.names(metadata)
colnames(my_sample_col) <- "Treatment"

my_sample_row <- data.frame(LIF_MYF$signature)
row.names(my_sample_row) <- LIF_MYF$gene_name
colnames(my_sample_row) <- "signature"


# Specify colors for annotations
my_colour = list(
  Treatment = c("cont" = "#5977ff", "stim" = "#f74747"),
  signature = c("LIF" = "Green", MYF = "Orange" )
)

method <- c("euclidean", "maximum", "manhattan", "canberra", "binary","minkowski")
n <- length(method)


pheatmap(table1_reordered,
         #kmeans_k = 4,
         annotation_col = my_sample_col,
         annotation_row = my_sample_row,
         annotation_colors = my_colour,
         #clustering_distance_rows = method[i],
         #clustering_distance_cols = "correlation",
         cluster_rows = FALSE,  # Cluster rows
         cluster_cols = TRUE,  # Cluster columns,
         cutree_rows = 4,
         cutree_cols =2,
         treeheight_col = 20,
         treeheight_row = 0, # hauteur de l'arbre pour ligne
         scale = "row",        # Scale rows (-> Z-score)
         main = "Heatmap of Differentially Expressed Genes",
         fontsize = 8,         # Adjust the font size for row and column names
         show_rownames = TRUE, # Show row names
         show_colnames = TRUE,  # Show column names
)


for (i in 1:n) {
   pheatmap(Data_Norm_sign_all,
           #kmeans_k = 4,
           annotation_col = my_sample_col,
           annotation_row = my_sample_row,
           annotation_colors = my_colour,
           #clustering_distance_rows = method[i],
           #clustering_distance_cols = "correlation",
           cluster_rows = FALSE,  # Cluster rows
           cluster_cols = ,  # Cluster columns,
           cutree_rows = 4,
           cutree_cols =2,
           treeheight_col = 20,
           treeheight_row = 0, # hauteur de l'arbre pour ligne
           scale = "row",        # Scale rows (-> Z-score)
           main = paste("Heatmap of Differentially Expressed Genes - ",method[i]),
           fontsize = 8,         # Adjust the font size for row and column names
           show_rownames = TRUE, # Show row names
           show_colnames = TRUE,  # Show column names
            )

}

