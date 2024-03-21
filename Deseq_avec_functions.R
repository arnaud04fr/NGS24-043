# Chargement des librairies
library(DESeq2)
library(tidyverse)
library(FactoMineR)
library(pheatmap)


######### functions ###################

get_geneID <- function(results) {
  # load biomart
  library(biomaRt)

  #create empty vectors 
  geneID <- NULL
  list_names <- NULL
  
  # Creation of the list of referencegene names
  list_names <- row.names(results)
  
  # Define the Ensembl database and create a biomart object
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Specify the attributes to retrieve
  attributes <- c("ensembl_gene_id", "ensembl_gene_id_version", "external_gene_name")
  
  # Retrieve information from Biomart
  geneID <- getBM(attributes = attributes, filters = "ensembl_gene_id_version", values = list_names, mart = ensembl)
  
  # Get the number of genes
  n <- nrow(geneID)
  
  # Replace empty cells with the stable ID
  for (i in 1:n) {
    if (geneID[i, 3] == "") {
      geneID[i, 3] <- geneID[i, 1]
    }
  }
  #change the name of the columns to ease table union
  colnames(geneID)[2] <- "gene"
  colnames(geneID)[3] <- "gene_name"
  
  return(geneID)
}

#cette fonction pour insérer en premier une colonne avec les transcript IDs avant GeneID
column_geneName <- function(results,p) {
  geneName <- NULL  
  library(dplyr)
    geneName <- row.names(results) #extraction de la liste des gènes
    results <- as.data.frame(cbind(results,geneName))
    results <- results %>% relocate(geneName, .before = p)
    return(results)
}

########### DESeq2 analysis #########

#chargement de la table de comptage (dans wd)
countdata <- read.table("count_NHLF_CM.txt", sep = "\t", header = T, row.names = 1, check.names=F)

#metadata:
coldata <- read.csv(file = "Metadata_NGS24-043.csv", header = T, sep =";")
str(coldata)
coldata$CM <- factor(coldata$CM) #cree les différents niveaux pour les différents CM
coldata$treatment <- factor(coldata$treatment) #cree les différents niveaux pour treatment (cont or stim)
str(coldata)

#Paired analysis on patient and time 
dds <- DESeqDataSetFromMatrix(countData = countdata, 
                                            colData = coldata, 
                                            design = ~ CM + treatment)
keep <- rowSums(counts(dds)) > 0 # added 8/03 ->  genes with counts
dds <- dds[keep,] # added 8/03 -> to remove genes with no counts
dds$treatment <- relevel(dds$treatment, "cont") # added 03/21 -> confirm control conditions

# Run the DESeq analysis and obtain the results table
dds <- DESeq(dds)
results_table <- results(dds)

# Extract the results for all genes
log2FoldChange <- results_table$log2FoldChange
pvalue <- results_table$pvalue
padj <- results_table$padj

# Create a new data frame with the results, including row names
results <- data.frame(
  log2FoldChange = log2FoldChange,
  pvalue = pvalue,
  padj = padj
)
rownames(results) <- rownames(results_table)

##### rajout des noms officiels ########
geneID <- get_geneID(results)
results$gene <- row.names(results)
results <- inner_join(results,geneID, by = "gene")
row.names(results) <- results$ensembl_gene_id
results <- results[,-c(4:5)] # remove the ensembl_gene_id
results <- as.data.frame(results)
results <- results %>% relocate(gene_name, .before = 1) #changer l'ordre des colonnes
results <- results[!duplicated(results$gene_name), ]
head(results)


######## Calcul et rajout du score GSEA ##########
n <- nrow(results)
score_GSEA <- matrix(0,n,1)
for (i in 1:n){
  score_GSEA[i] <- sign(results[i,2])*(2^abs(results[i,2]))*-log10(results[i,4])
}
results <- cbind(results,score_GSEA)
head(results)

#préparation tableau avec GSEA
table_GSEA <- as.data.frame(matrix(0,n,2))
table_GSEA[,1] <- results[,1]
table_GSEA[,2] <- as.numeric(results[,5])
colnames(table_GSEA) <- c("Gene_name","Rank")
table_GSEA <- table_GSEA %>% drop_na() #remove lines with NA value
head(table_GSEA)

#ecriture des tables au format tsv dans le dossier results
write.table(results,"./results/DESeq-table_gene_all.tsv",sep='\t', row.names = F, col.names=T, quote = F)
write.table(table_GSEA,"./results/DESeq-gene_GSEAonly_all.tsv",sep='\t', row.names = F, col.names=T, quote = F)


###### selection des gènes avec LogFC >=2 et padj <= 0.05
result_diff <- subset(results, padj <= 0.05 & abs(log2FoldChange) >= 1)
Down_expressed_genes <- subset(result_diff, log2FoldChange <= -1)
UP_expressed_genes <- subset(result_diff, log2FoldChange >= 1)
nr <- nrow(result_diff)
print(paste("Total DEG" , nrow(result_diff), "/ UP ", nrow(UP_expressed_genes), " / DOWN ", nrow(Down_expressed_genes)))
write.table(result_diff,"./results/DESeq-table_gene_rank_DEG.tsv",sep='\t', row.names = T, col.names=T)


############ extraction tables comptages normalisées ###########
ddsNorm <- estimateSizeFactors(dds)
Data_Norm <- (counts(ddsNorm, normalized=TRUE))
Data_Norm_all <- as.data.frame(Data_Norm)
head(Data_Norm_all)
rm(geneID)
geneID <- get_geneID(Data_Norm_all) #rajout nom officiel
Data_Norm_all$gene <- row.names(Data_Norm_all) 
Data_Norm_all <- inner_join(Data_Norm_all,geneID, by = "gene")
Data_Norm_all <- Data_Norm_all[!duplicated(Data_Norm_all$gene_name), ]
row.names(Data_Norm_all) <- Data_Norm_all$gene_name
Data_Norm_all <- Data_Norm_all[,-c(11:13)] # remove the ensembl_gene_id
colnames(Data_Norm_all) <- sub("_.*", "", colnames(Data_Norm_all))
write.table(Data_Norm_all,"./results/DESeq_Table_norm_count_all.tsv",sep='\t', row.names = T, col.names=T) #sauvegarde table comptages normalisées

# table will only DEG:
Data_Norm_DEG <- as.data.frame(Data_Norm_all[result_diff[,1], ])
write.table(Data_Norm_DEG,"./results/DESeq_Table_norm_count_DEG.tsv",sep='\t', row.names = T, col.names=T) #sauvegarde table comptages normalisées


###################" graph generation ###################################

#MA-plots 
plotMA(results_table, ylim=c(-6,6), alpha = 0.1)

#Volcano plot generation

# Extract the relevant information from Results
results <- as.data.frame(results)
results_df <- data.frame(
  gene = results$gene,
  log2FoldChange = results$log2FoldChange,
  neg_log10_pvalue = -log10(results$pvalue),
  padj = results$padj
)

# Create the volcano plot using ggplot2
ggplot(results_df, aes(x = log2FoldChange, y = neg_log10_pvalue)) +
  geom_point(aes(color = ifelse(neg_log10_pvalue > -log10(0.05) & abs(log2FoldChange) > 1, "Significant", "Non-significant")), show.legend =F) +
  scale_color_manual(values = c("Significant" = "cyan3", "Non-significant" = "grey")) +
  xlim(c(-10, 15)) + ylim(c(0, 100)) +
  labs(x = "Log2 Fold Change", y = "-log10(p-value)",
       title = "Volcano Plot of Differentially Expressed Genes") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  theme_minimal() +
  geom_text(aes(label = ifelse(neg_log10_pvalue > 10 & abs(log2FoldChange) > 3, gene, "")), 
            vjust = -0.9, hjust = 0.5, size = 2, angle = 20) +
    geom_text(aes(label = ifelse(neg_log10_pvalue > 40 & log2FoldChange < 2, gene, "")), 
            vjust = -0.9, hjust = 0.5, size = 2, angle = 20) 
  
#Clustering and heatmap

pheatmap(Data_Norm_DEG ,
         #kmeans_k = 4,
         clustering_distance_rows = "correlation",
         cluster_rows = TRUE,  # Cluster rows
         cluster_cols = TRUE,  # Cluster columns,
         scale = "row",        # Scale rows (you can change this if needed)
         main = "Heatmap of Differentially Expressed Genes",
         fontsize = 8,         # Adjust the font size for row and column names
         show_rownames = FALSE, # Show row names
         show_colnames = TRUE  # Show column names
)

######### test PCA analysis from DESeq2 norm data

# Step 1: Prepare ans transpose the data 
data_for_pca2 <- t(Data_Norm) #transpose le tableau de données
#row.names(data_for_pca2) <- c("ASF859_D0","ASF859_D14","ASF876_D0","ASF876_D14","ASF879-2_D0","ASF879-2_D14","ASF880_D0","ASF880_D14")


# Step 2: Perform PCA
pca_result <- PCA(data_for_pca2, graph = FALSE)

# Step 3: Extract PCA results for plotting
pc_df <- as.data.frame(pca_result$ind$coord[, 1:2]) # Extract the first two principal components
pc_df$Sample <- rownames(pc_df)                     # Add a "Sample" column for labels

# Step 4: Create the ggplot visualization with colors
ggplot(pc_df, aes(x = Dim.1, y = Dim.2, label = Sample)) +
  geom_point(aes(color = grepl("MC", Sample)), size = 2, shape = 16) +  # Grey dots for columns with "D0"
  geom_point(aes(color = grepl("ctrl", Sample)), size = 2, shape = 16) + # Black dots for columns with "D14"
  geom_text(aes(y = Dim.2), 
            position = position_nudge(y = 10), # Adjust the labels to be above the points
            #vjust = -0.1, # Center the labels vertically
            #hjust = -0.6,
            size = 2, angle = -45) +  # Add labels for each point with adjusted size and 45° angle
  labs(x = "PC1 (Principal Component 1)", y = "PC2 (Principal Component 2)") +
  ggtitle("PCA Analysis - First Two Principal Components") +
  theme_minimal() +
  scale_color_manual(name = "Treatment", values = c("grey", "black"), labels = c("ctrl", "stim")) +  # Custom color scale with legend labels
  coord_cartesian(xlim = c(-150, 150))

########### PCA version 2 #########

#load data for PCA using normalized data dss
#df <- read.table("./results/DESeq_Table_norm_count_all.tsv", sep = "\t", header = T, row.names = 1, check.names=F)
df <- as.data.frame(t(Data_Norm_all))


# Assuming df.pca is the output of prcomp
df.pca <- prcomp(df,scale = TRUE)
percentage_variance <- (df.pca$sdev^2 / sum(df.pca$sdev^2)) * 100

## Create a data frame with PC and Percentage Variance
variance_table <- data.frame(
  PC = paste0("PC", 1:length(df.pca$sdev)),
  Percentage_Variance = paste0(round(percentage_variance, 2), "%")
)

PC1 <- df.pca$x[,1]
PC2 <- df.pca$x[,2]

# Define custom colors for treatment + metadata
treatment <- rep(c("stim","cont"),5)
treatment_color <- c("cont" = "#619CFF", "stim" = "#F8766D")

ggplot(df, 
       aes(x = PC1, 
           y = PC2, 
           color = treatment, 
           label = rownames(df))) +
  geom_point(size = 3, shape = 16) +
  geom_text(vjust = -2, size = 3, angle = -20) +
  stat_ellipse(type = "t", level = 0.9, segments = 51, linetype = "dashed") +
  labs(title = "PCA Analysis",
       x = paste0("PC1 (", variance_table$Percentage_Variance[1], ")"),
       y = paste0("PC2 (", variance_table$Percentage_Variance[2], ")")) +
  ggtitle("PCA Analysis - First Two Principal Components") +
  coord_cartesian(xlim = c(-210, 210)) +
  theme_minimal() +
  scale_color_manual(values = treatment_color)




