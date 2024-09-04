library(ggVennDiagram)
library(ggvenn)
library(VennDiagram)

#Best expressed gene (counts)

Data_Norm <- (counts(estimateSizeFactors(dds), normalized=TRUE))
colnames(Data_Norm)
n <- nrow(Data_Norm)
group <- c("ctrl", "MC") #indicate sample groups to analyze.

for (i in 1:length(group)) {
  mot <- group[i]
  data_subset <- as.data.frame(Data_Norm[, colnames(dds) %in% grep(mot, colnames(dds), value = TRUE)])
  data_subset$mean_count <- rowMeans(data_subset[, 1:5])
  data_subset <- arrange(data_subset,desc(mean_count))
  Top30 <- data_subset[1:31,]
  Gene_ID <- get_geneID(Top30)
  Top30$gene <- row.names(Top30)
  Top30 <- inner_join(Top30,Gene_ID, "gene")
  Top30 <- Top30 %>% relocate(gene_name, .before = 1)
  Top30 <- Top30 %>% relocate(mean_count, .before = 2)
  Top30 <- Top30[,1:2]
  Top30$gene_name <- factor(Top30$gene_name, levels = Top30$gene_name)
  assign(paste("Top30_", sep ="", mot), Top30$gene_name)
  
  
  ggplot(Top30, aes(x = mean_count, y = fct_rev(factor(gene_name)))) +
    geom_bar(stat = "identity", fill = "lightgray", color = "black", alpha = 0.5) +
    labs(title = paste("Best expressed gene for the", mot, "conditions"),
         x = "mean count", y = "") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10, face = "bold.italic"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.ticks.x = element_line(color = "black"),
      axis.line.x = element_line(color = "black"),  # Add line for x-axis
      panel.grid.major = element_blank(),           # Remove major grid lines
      panel.grid.minor = element_blank()            # Remove minor grid lines
    )  
}

# Create a list containing your gene lists
gene_lists <- list(Basal = Top30_ctrl, CM = Top30_MC)


#Draw Venn diagram
ggVennDiagram(
  gene_lists, label_alpha = 0,
  category.names = c("Basal","CM")
) +
  ggplot2::scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")

# Create a Venn diagram with embedded data:
v.table <- venn(gene_lists,intersection=TRUE)
isect <- attr(v.table, "intersection")
isect$CM
isect$Basal
