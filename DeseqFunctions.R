#################
#### gene_ID ####
#################

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
  #end of function get_geneID
}

############################################
#### Graph_Matrice - scatterplot Matrix ####
############################################

Graph_Matrice <- function(dds_subset,mot) {

# Load necessary libraries
library(DESeq2)
library(ggplot2)

# Normalize the counts (get variance stabilized or rlog-transformed counts) for the subset
design(dds_subset) <- ~1  # Simplifying the design for the subset
vsd_subset <- vst(dds_subset, blind = TRUE) #Use rlog() instead of vst() if preferred
  
# Extract the transformed counts for visualization
normalized_counts_subset <- assay(vsd_subset)

pairs(normalized_counts_subset, 
      upper.panel = function(x, y) {
        # Correlation text
        usr <- par("usr")
        on.exit(par(usr))
        par(usr = c(0, 1, 0, 1))
        text(0.5, 0.5, round(cor(x, y), 2), cex = 1.5)
      },
      lower.panel = function(x, y) {
        # Scatter plot with regression line
        points(x, y, pch = 20, col = scales::alpha("black", 0.3))
        abline(lm(y ~ x), col = "red")
      },
      diag.panel = function(x) {
        # Histogram
        par(new = TRUE)
        hist(x, breaks = 20, col = "cyan", border = "white", main = "", xlab = "", ylab = "", axes = FALSE)
      },
      main = paste("Scatterplot matrix of", mot, "samples"))

#end of function Graph_Matrice
}  

##########################################################
###### bar plot -  before / after Normalization ##########
##########################################################

Graph_norm <- function(dds) {
  
# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(reshape2)  # For reshaping data

# Assuming you have a DESeqDataSet object called `dds`
# Get the raw counts from the DESeqDataSet
raw_counts <- counts(dds)

# Log2 transform the raw counts (adding 1 to avoid log(0))
log2_raw_counts <- log2(raw_counts + 1)

# Normalize the counts (using variance stabilizing transformation or rlog)
vsd <- vst(dds, blind = FALSE)  # Use rlog(dds, blind = FALSE) if preferred

# Extract normalized counts
normalized_counts <- assay(vsd)

# Log2 transform the normalized counts (if using vst, already transformed)
# Log2 is typically already accounted for in vst, so this step may be redundant.
# Uncomment the line below if needed for your transformation method
#log2_normalized_counts <- log2(normalized_counts + 1)

# Prepare data for plotting
# Melt the raw and normalized counts for ggplot
raw_counts_df <- melt(log2_raw_counts)
normalized_counts_df <- melt(normalized_counts)  # If log2 already applied

# Add a new column to differentiate between raw and normalized counts
raw_counts_df$Normalization <- "Before Normalization"
normalized_counts_df$Normalization <- "After Normalization"

# Combine the two data frames
combined_df <- rbind(raw_counts_df, normalized_counts_df)

# Correct the column labels (Sample and Gene are inverted)
colnames(combined_df) <- c("Gene", "Sample", "Log2_Counts", "Normalization")

# Separate the combined data into two data frames for individual plots
before_norm_df <- combined_df[combined_df$Normalization == "Before Normalization", ]
after_norm_df <- combined_df[combined_df$Normalization == "After Normalization", ]

# Plot for raw counts (before normalization)
figure_before <- ggplot(before_norm_df, aes(x = Sample, y = Log2_Counts, fill = Sample)) +
  geom_boxplot(outlier.shape = TRUE) +
  labs(title = "Log2(Counts) Before Normalization",
       x = "Sample",
       y = "Log2(Counts)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot for normalized counts (after normalization)
figure_after <- ggplot(after_norm_df, aes(x = Sample, y = Log2_Counts, fill = Sample)) +
  geom_boxplot(outlier.shape = TRUE) +
  labs(title = "Log2(Counts) After Normalization",
       x = "Sample",
       y = "Log2(Counts)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combined plot for comparison
figure_combined <- ggplot(combined_df, aes(x = Sample, y = Log2_Counts, fill = Normalization)) +
  geom_boxplot(outlier.shape = TRUE) +
  scale_fill_manual(values = c("Before Normalization" = "lightblue", "After Normalization" = "lightgreen")) +
  labs(title = "Log2(Counts) Before and After Normalization",
       x = "Sample",
       y = "Log2(Counts)",
       fill = "Normalization") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plots
print(figure_before)
print(figure_after)
print(figure_combined)

#end of function Graph_norm 
}  

##########################################################
######            Scatter plot count            ##########
##########################################################

Graph_scatter<- function(norm_counts) {
  
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
  
  #end of function Graph_scatter 
}  





