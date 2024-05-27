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

########### DESeq2 analysis #########

#chargement de la table de comptage (dans wd)
countdata <- read.table("count_NHLF_CM.txt", sep = "\t", header = T, row.names = 1, check.names=F)

#remove CM from ASF894 -> to test.
#colnames(countdata) #ASF398 -> columns 5 and 6
#countdata <- countdata[,-(5:6)]

#metadata:
coldata <- read.csv(file = "Metadata_NGS24-043.csv", header = T, sep =";")
#remove metadata of ASF894.
#coldata[,1] #ASF398 -> rows 5 and 6
#coldata <- coldata[-(5:6),]
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
dds$treatment <- relevel(dds$treatment, "cont") # added 03/21 -> define clearly the control condition

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
row.names(results) <- results$gene
results <- results[,-c(4:5)] # remove the ensembl_gene_id
results <- as.data.frame(results)
results <- results %>% relocate(gene_name, .before = 1) #changer l'ordre des colonnes
#results <- results[!duplicated(results$gene_name), ]
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
table_GSEA  <- arrange(table_GSEA, desc(Rank)) #order genes by rank
head(table_GSEA)

#ecriture des tables au format tsv dans le dossier results
results$gene <- row.names(results)
results <- results %>% relocate(gene, .before = 1)
#write.table(results,"./results/DESeq-table_gene_all.tsv",sep='\t', row.names = F, col.names=T, quote = F)
#write.table(table_GSEA,"./results/DESeq-gene_GSEAonly_all.tsv",sep='\t', row.names = F, col.names=T, quote = F)


###### selection des gènes avec LogFC >=2 et padj <= 0.05
result_diff <- subset(results, padj <= 0.05 & abs(log2FoldChange) >= 1)
Down_expressed_genes <- subset(result_diff, log2FoldChange <= -1)
UP_expressed_genes <- subset(result_diff, log2FoldChange >= 1)
nr <- nrow(result_diff)
print(paste("Total DEG" , nrow(result_diff), "/ UP ", nrow(UP_expressed_genes), " / DOWN ", nrow(Down_expressed_genes)))
#write.table(result_diff,"./results/DESeq-table_gene_rank_DEG.tsv",sep='\t', row.names = T, col.names=T)


############ extraction tables comptages normalisées ###########
ddsNorm <- estimateSizeFactors(dds)
Data_Norm <- (counts(ddsNorm, normalized=TRUE))
Data_Norm_all <- as.data.frame(Data_Norm)
head(Data_Norm_all)
rm(geneID)
geneID <- get_geneID(Data_Norm_all) #rajout nom officiel
Data_Norm_all$gene <- row.names(Data_Norm_all) 
Data_Norm_all <- inner_join(Data_Norm_all,geneID, by = "gene")
Data_Norm_all <- Data_Norm_all[!duplicated(Data_Norm_all$gene_name), ] #remove duplicated symbols
row.names(Data_Norm_all) <- Data_Norm_all$gene_name
Data_Norm_all <- Data_Norm_all[,1:10] # remove the ensembl_gene_id -> check the number of columns before
colnames(Data_Norm_all) <- sub("_.*", "", colnames(Data_Norm_all))
head(Data_Norm_all)
#write.table(Data_Norm_all,"./results/DESeq_Table_norm_count_all.tsv",sep='\t', row.names = T, col.names=T) #sauvegarde table comptages normalisées

# table will only DEG:
Data_Norm_DEG <- as.data.frame(Data_Norm_all[result_diff$gene_name, ])
#write.table(Data_Norm_DEG,"./results/DESeq_Table_norm_count_DEG.tsv",sep='\t', row.names = T, col.names=T) #sauvegarde table comptages normalisées



###################" graph generation ###################################
#Data_Norm_DEG <- read.table("./results/DESeq_Table_norm_count_DEG.tsv",sep='\t', header =T)
#Data_Norm_all <- read.table("./results/DESeq_Table_norm_count_all.tsv",sep='\t', header =T)

## MA-plots ## 
plotMA(results_table, ylim=c(-6,6), alpha = 0.05)


## Scatter plot ##

norm_counts <- Data_Norm_all
colnames(norm_counts) <- c("MC1","CT1","MC2",'CT2',"MC3",'CT3',"MC4",'CT4',"MC5",'CT5')
norm_counts <- as.data.frame(norm_counts[,order(colnames(norm_counts))])

n <- length(rownames(norm_counts))
counts_mean <- as.data.frame(matrix(0 ,nrow = n, ncol = 2 ))
colnames(counts_mean) <- c("cont", "stim")
row.names(counts_mean) <- row.names(norm_counts)

for (i in 1:n) {
  counts_mean[i,1] <- log10(mean(as.numeric(norm_counts[i,1:5])) +1)
  counts_mean[i,2] <- log10(mean(as.numeric(norm_counts[i,6:10])) +1)
}

up_genes <- subset(results, padj < 0.05 & log2FoldChange > 1)
up_genes <- up_genes[!duplicated(up_genes$gene_name), ] #remove duplicated symbols
down_genes <- subset(results, padj < 0.05 & log2FoldChange < -1)
down_genes <- down_genes[!duplicated(up_genes$gene_name), ] #remove duplicated symbols


# Add a new column to 'counts_mean' indicating the regulation status of each gene
counts_mean$regulation <- ifelse(rownames(counts_mean) %in% up_genes$gene_name, "Up",
                                 ifelse(rownames(counts_mean) %in% down_genes$gene_name, "Down", "NS"))

counts_mean$genes <- row.names(counts_mean)
gene_list <- data.frame(
  genes = results$gene_name,
  log2FoldChange = results$log2FoldChange,
  neg_log10_pvalue = -log10(results$pvalue),
  padj = results$padj,
  GSEA = score_GSEA
)
counts_mean <- inner_join(counts_mean,gene_list, by= "genes")
counts_mean <- counts_mean[!duplicated(counts_mean$genes), ] #remove duplicated symbols
row.names(counts_mean) <- counts_mean$genes
 
#Create the scatter plot
ggplot(counts_mean, aes(x = cont, y = stim, color = regulation)) +
  geom_point() +
  scale_color_manual(values = c("blue", "black", "red")) + # Set the colors for up, down, and not changed genes
  xlab("Log10 Normalized Counts (Control)") +
  ylab("Log10 Normalized Counts (Stim)") +
  ggtitle("Scatter plot of Control vs Treated conditions") +
  geom_text(aes(label = ifelse(GSEA > 550 | GSEA < -225, genes, "")), 
          vjust = -0.9, hjust = 0.5, size = 2, angle = 20) 


## Volcano plot generation ##

# Extract the relevant information from Results
results <- as.data.frame(results)
results_df <- data.frame(
  gene = results$gene_name,
  log2FoldChange = results$log2FoldChange,
  neg_log10_pvalue = -log10(results$pvalue),
  padj = results$padj
)

# Create the volcano plot using ggplot2
ggplot(results_df, aes(x = log2FoldChange, y = neg_log10_pvalue)) +
  geom_point(aes(color = ifelse(neg_log10_pvalue > -log10(0.05) & abs(log2FoldChange) > 1, "Significant", "Non-significant")), show.legend =F) +
  scale_color_manual(values = c("Significant" = "cyan3", "Non-significant" = "grey")) +
  xlim(c(-10, 15)) + ylim(c(0, 150)) +
  labs(x = "Log2 Fold Change", y = "-log10(p-value)",
       title = "Volcano Plot of Differentially Expressed Genes") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  theme_minimal() +
  geom_text(aes(label = ifelse(neg_log10_pvalue > 10 & abs(log2FoldChange) > 3, gene, "")), 
            vjust = -0.9, hjust = 0.5, size = 2, angle = 20) +
    geom_text(aes(label = ifelse(neg_log10_pvalue > 40 & log2FoldChange < 2, gene, "")), 
            vjust = -0.9, hjust = 0.5, size = 2, angle = 20) +
  geom_text(aes(label = ifelse(neg_log10_pvalue > 40 & log2FoldChange > 1, gene, "")), 
            vjust = -0.9, hjust = 0.5, size = 2, angle = 20)

  
#######Clustering and heatmap##########


#annotations des colonnes pour graph par traitement
my_sample_col <- data.frame(coldata$treatment)
row.names(my_sample_col) <- colnames(Data_Norm_DEG)
colnames(my_sample_col) <- "Treatment"

# Specify colors for annotations
my_colour = list(
  Treatment = c("cont" = "#5977ff", "stim" = "#f74747"))

pheatmap(Data_Norm_DEG ,
         #kmeans_k = 4,
         annotation_col = my_sample_col,
         clustering_distance_rows = "correlation",
         cluster_rows = TRUE,  # Cluster rows
         cluster_cols = TRUE,  # Cluster columns,
         scale = "row",        # Scale rows (-> zscore)
         cutree_rows = 2,
         cutree_cols =2,
         main = "Heatmap of Differentially Expressed Genes",
         fontsize = 8,         # Adjust the font size for row and column names
         show_rownames = FALSE, # Show row names
         show_colnames = TRUE  # Show column names
)

######### test PCA analysis from DESeq2 norm data

# Step 1: Prepare ans transpose the data 
data_for_pca2 <- t(Data_Norm) #transpose le tableau de données


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
str(df)

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
treatment <- NULL
treatment <- rep(c("stim","cont"),5) #modify according to the number of samples 
treatment_color <- c("cont" = "#619CFF", "stim" = "#F8766D")

ggplot(df, 
       aes(x = PC1, 
           y = PC2, 
           color = treatment, 
           label = rownames(df))) +
  geom_point(size = 3, shape = 16) +
  geom_text(vjust = -2, size = 3, angle = -20) +
  stat_ellipse(type = "t", level = 0.9, segments = 20, linetype = "dashed") +
  labs(title = "PCA Analysis",
       x = paste0("PC1 (", variance_table$Percentage_Variance[1], ")"),
       y = paste0("PC2 (", variance_table$Percentage_Variance[2], ")")) +
  ggtitle("PCA Analysis - First Two Principal Components") +
  coord_cartesian(xlim = c(-210, 210)) +
  theme_minimal() +
  scale_color_manual(values = treatment_color)



######################## Pathway analysis ##################

library(pathview)
library(GO.db)
library(enrichplot)
library(clusterProfiler)
library(biomaRt)


# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

#reading data
df <- read.table("./results/DESeq-table_gene_all.tsv", header = T, sep = "\t", quote = "")

##### rajout des noms officiels ########

#creation de la liste:
list_names <- df$gene

# Définissez la base de données Ensembl et créez un objet biomart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Précisez les attributs que vous souhaitez récupérer
attributes <- c("ensembl_gene_id", "ensembl_gene_id_version")

# Récupérez les informations à partir de Biomart
gene_info <- getBM(attributes = attributes, filters = "ensembl_gene_id_version", values = list_names, mart = ensembl)


#fusionner les tables avec Tidyverse:
df <- inner_join(df, gene_info, by = c("gene" = "ensembl_gene_id_version"))


######################################################
############### GO ANALYSIS ##########################
#####################################################

# we want the log2FoldChange score
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$ensembl_gene_id
head(original_gene_list)

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results (padj < 0.05)
sig_genes_df = subset(df, padj < 0.05)

# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$log2FoldChange

# Name the vector
names(genes) <- sig_genes_df$ensembl_gene_id

# omit NA values
genes <- na.omit(genes)

# filter on min log2fold change (log2FoldChange > 2)
genes <- names(genes)[abs(genes) > 2]

onto  <- c("MF", "BP", "CC") 
n <- length(onto)
go_enrich_list <- NULL 


for (i in 1:n) {
  print(onto[i])
  go_enrich <- enrichGO(gene = genes,
                        universe = names(gene_list),
                        OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
                        keyType = 'ENSEMBL',
                        readable = T,
                        ont = onto[i],
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)
  #graph generation
  titre <- c(paste("GO",onto[i]))
 print(upsetplot(go_enrich),
       title = titre )
 
 print(barplot(go_enrich, 
          drop = TRUE, 
          showCategory = 10, 
          title = titre,
          font.size = 8))
 
  print(dotplot(go_enrich),
        title = titre)
  
  # Assign the enrichGO result to the list with a dynamic name
  go_enrich_list[[paste("go_enrich_", onto[i], sep = "")]] <- go_enrich
  
 }


######################################################
############### GSEA ANALYSIS #######################
#####################################################

# we want the GSEA score rank from df
original_gene_list2 <- df$score_GSEA

# name the vector
names(original_gene_list2) <- df$ensembl_gene_id

# omit any NA values 
gene_list<-na.omit(original_gene_list2)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
             pAdjustMethod = "none")


#Dotplot
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

ridgeplot(gse) + labs(x = "enrichment distribution")

# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)





####################################################
############ KEGG Pathway Enrichment    ############
####################################################


# Convert gene IDs for enrichKEGG function
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db") # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 <- inner_join(df, dedup_ids, by = c("ensembl_gene_id" = "ENSEMBL"))

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$ENTREZID

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

# Exctract significant results from df2
kegg_sig_genes_df = subset(df2, padj < 0.05)

# From significant results, we want to filter on log2fold change
kegg_genes <- kegg_sig_genes_df$log2FoldChange

# Name the vector with the CONVERTED ID!
names(kegg_genes) <- kegg_sig_genes_df$ENTREZID

# omit NA values
kegg_genes <- na.omit(kegg_genes)
head(kegg_genes)

# filter on log2fold change (PARAMETER)
kegg_genes <- (kegg_genes)[abs(kegg_genes) > 1]
head(kegg_genes)

##### Create enrichKEGG object
kegg_organism = "hsa"
kk <- enrichKEGG(gene=names(kegg_genes), universe=names(kegg_gene_list),organism=kegg_organism, pvalueCutoff = 0.05, keyType = "ncbi-geneid")

# omit any NA values 
gene_list<-na.omit(original_gene_list)

barplot(kk, 
        showCategory = 40, 
        title = "Enriched Pathways",
        font.size = 8)

dotplot(kk, 
        showCategory = 40, 
        title = "Enriched Pathways",
        font.size = 8)

####### UPSET PLOT
upsetplot(kk)

####### CNET PLOT

# categorySize can be either 'pvalue' or 'geneNum'
cnetplot(kk, categorySize="pvalue", foldChange=gene_list)

####### PATHVIEW

#change dme according to gseKEGG output table
KEGGf <- as.data.frame(kk)
ID <- KEGGf$ID
n <- length(ID)

#PNG
for (i in 1:n) {
  dme <- pathview(gene.data=kegg_genes, 
                  pathway.id=ID[i], 
                  species = kegg_organism,
                  kegg.native = TRUE,
                  same.layer = F) 
}

#PDF
for (i in 1:n) {
  dme <- pathview(gene.data=kegg_genes, 
                  pathway.id=ID[i], 
                  species = kegg_organism,
                  kegg.native = F,
                  same.layer = T) 
}

############## generating heatmap using GoI#######
print(KEGGf$Description)
# Assuming KEGGf$geneID is a character vector
geneID_list <- strsplit(KEGGf$geneID, "/")
names(geneID_list) <- KEGGf$ID
# Convert the split strings to numeric vectors
geneID_list <- lapply(geneID_list, function(x) as.numeric(x))

#liste cytokine
cytokine <- geneID_list$hsa04060
cytokine <- bitr(cytokine, fromType = "ENTREZID", toType = "SYMBOL", OrgDb="org.Hs.eg.db") #rajoute le gène symbol
#list PPARg
PPARg <- geneID_list$hsa03320
PPARg <- bitr(PPARg, fromType = "ENTREZID", toType = "SYMBOL", OrgDb="org.Hs.eg.db") #rajoute le gène symbol

#filter DEG:
Data_Norm_DEG_cytokine <- Data_Norm_DEG[row.names(Data_Norm_DEG) %in% cytokine$SYMBOL,]
Data_Norm_DEG_PPARg <- Data_Norm_DEG[row.names(Data_Norm_DEG) %in% PPARg$SYMBOL,]
colnames(Data_Norm_DEG_PPARg)
nom <- c("MC-888", "ctrl-888","MC-889", "ctrl-889","MC-894", "ctrl-894","MC-897", "ctrl-897","MC-913", "ctrl-913")
colnames(Data_Norm_DEG_PPARg) <- nom
Data_Norm_DEG_PPARg <- sort.DataFrame(Data_Norm_DEG_PPARg)
Data_Norm_DEG_PPARg <- Data_Norm_DEG_PPARg[, order(names(Data_Norm_DEG_PPARg))]
my_sample_col2 <- as.data.frame(colnames(Data_Norm_DEG_PPARg))
colnames(my_sample_col2) <- "Treatment"
row.names(my_sample_col2) <- my_sample_col2$Treatment
my_sample_col2[1:5,] <- "cont"
my_sample_col2[6:10,] <- "stim"

#generate heatmap


pheatmap(Data_Norm_DEG_cytokine ,
         #kmeans_k = 4,
         annotation_col = my_sample_col,
         clustering_distance_rows = "correlation",
         cluster_rows = TRUE,  # Cluster rows
         cluster_cols = TRUE,  # Cluster columns,
         scale = "row",        # Scale rows (-> zscore)
         cutree_rows = 2,
         cutree_cols =2,
         main = "hsa04060 - Cytokine-cytokine receptor interaction",
         fontsize = 8,         # Adjust the font size for row and column names
         show_rownames = TRUE, # Show row names
         show_colnames = FALSE  # Show column names
)

pheatmap(Data_Norm_DEG_PPARg ,
         #kmeans_k = 4,
         annotation_col = my_sample_col2,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "canberra",
         cluster_rows = TRUE,  # Cluster rows
         cluster_cols = FALSE,  # Cluster columns,
         scale = "row",        # Scale rows (-> zscore)
         cutree_rows = 2,
         cutree_cols =2,
         main = "hsa03320 - PPAR signaling pathway",
         fontsize = 8,         # Adjust the font size for row and column names
         show_rownames = TRUE, # Show row names
         show_colnames = TRUE  # Show column names
)



