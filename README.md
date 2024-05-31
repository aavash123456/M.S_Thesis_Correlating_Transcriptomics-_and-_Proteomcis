# M.S_Thesis_Correlating_Transcriptomics-_and-_Proteomcis
part of my thesis where I calculated the correlation between transcriptomics and proteomics datasets.


protein_data <- "/Users/aavashadhikari/desktop/data/mergedResultsmaster.xlsx"
protein_df<-read_excel(protein_data, sheet = "mergedResultsMaster")
transcriptomics <- read.csv("Users/aavashadhikari/desktop/data/res0.01.csv")
colnames(transcriptomics)[1] <- "gene" 

#dim(protein_df)==many

pro_df <- protein_df[1:3729, ]

#dim(pro_df)
#colnames(pro_df)

protein_df_1<- pro_df[,c(37:60, 121)]

#dim(protein_df_1)

col_names <- names(protein_df_1)

new_col_order <- c(col_names[25], col_names[-25])

protein_df_1 <- protein_df_1[, new_col_order]

names(protein_df_1)[1] <- "gene"

protein_df_1$gene <- gsub(";", "", protein_df_1$gene)

names(protein_df_1) <- gsub("\\.y$", "", names(protein_df_1))

# Remove specified columns
protein_df_1 <- protein_df_1[, !names(protein_df_1) %in% c("FA4", "SA5", "FB5", "SB5", "FC4", "SC2")]

dim(protein_df_1)

# List of column prefixes for each group
groups <- c("FA", "FB", "FC", "SA", "SB", "SC")

# Loop through each group
for (group in groups) {
  # Extract columns for the current group
  group_cols <- grep(paste0("^", group), names(protein_df_1), value = TRUE)
  
  # Convert columns to numeric
  protein_df_1[, group_cols] <- sapply(protein_df_1[, group_cols], as.numeric)
  
  # Calculate row means for the current group
  mean_values <- rowMeans(protein_df_1[, group_cols], na.rm = TRUE)
  
  # Attach mean values as new columns to the original dataframe
  protein_df_1[[paste0("mean_", group)]] <- mean_values
}

# Calculate mean of mean values for groups FA, FB, and FC
mean_F <- rowMeans(protein_df_1[, c("mean_FA", "mean_FB", "mean_FC")], na.rm = TRUE)

# Add mean_F as a new column "F" to the dataframe
protein_df_1$F <- mean_F

# Calculate mean of mean values for groups SA, SB, and SC
mean_S <- rowMeans(protein_df_1[, c("mean_SA", "mean_SB", "mean_SC")], na.rm = TRUE)

# Add mean_S as a new column "S" to the dataframe
protein_df_1$S <- mean_S



# Calculate fold change
protein_df_1$FOLD_CHANGE <- protein_df_1$S / protein_df_1$F

# Calculate log2 fold change
protein_df_1$log2FC <- log2(protein_df_1$S / protein_df_1$F)


# Scatter plot of F and S
plot(protein_df_1$F, protein_df_1$S, 
     xlab = "F", ylab = "S",
     main = "Scatter Plot of F and S")


# Calculate correlation between F and S
correlation <- cor(protein_df_1$F, protein_df_1$S, use= "pairwise.complete.obs")
#corr = 0.97 



# Merge datasets based on the common 'gene' column
common_df <- merge(transcriptomics, protein_df_1, by = "gene")

dim(common_df)
3647 35



#checking if data are normally distributes

# Histogram for log2FC
ggplot(common_df, aes(x = log2FC)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black") +
  labs(title = "Histogram of log2FC")


# Histogram for log2FoldChange
ggplot(common_df, aes(x = log2FoldChange)) +
  geom_histogram(binwidth = 0.5, fill = "blue", color = "black") +
  labs(title = "Histogram of log2FoldChange")


#computing the correlation
correlation <- cor(common_df$log2FoldChange, common_df$log2FC, method = "pearson")
#0.2262
spearman_correlation <- cor(common_df$log2FoldChange, common_df$log2FC, method = "spearman")
#0.212137


# Load ggplot2 for visualization
library(ggplot2)

# Create scatter plot
ggplot(common_df, aes(x = log2FoldChange, y = log2FC)) +
  geom_point() +
  geom_smooth(method = "lm", col = "blue") +
  labs(title = "Scatter Plot of log2FoldChange vs log2FC",
       x = "log2FoldChange (Transcriptomics)",
       y = "log2FC (Proteomics)") +
  theme_minimal()


#both F AND S column have equal no of elements(row) and both do not have missing values. 

# Initialize a vector to store p-values
p_values <- numeric(nrow(protein_df_1))

# Loop through each row and perform a paired t-test
for (i in 1:nrow(protein_df_1)) {
  # Extract mean FA, FB, FC values for F and mean SA, SB, SC values for S
  F_values <- as.numeric(protein_df_1[i, c("mean_FA", "mean_FB", "mean_FC")])
  S_values <- as.numeric(protein_df_1[i, c("mean_SA", "mean_SB", "mean_SC")])
  
  # Perform paired t-test
  if (!any(is.na(F_values)) && !any(is.na(S_values))) {
    t_test_result <- t.test(F_values, S_values, paired = TRUE)
    p_values[i] <- t_test_result$p.value
  } else {
    p_values[i] <- NA  # Assign NA if there are missing values
  }
}

# Add p-values as a new column to the dataset
protein_df_1$p_value <- p_values

# Adjust p-values for multiple testing using the Benjamini-Hochberg method
protein_df_1$p_adj_value <- p.adjust(protein_df_1$p_value, method = "BH")


# Check the distribution of raw p-values
summary(protein_df_1$p_value)

# Visualize the distribution of raw p-values
hist(protein_df_1$p_value, breaks = 50, main = "Distribution of Raw P-Values", xlab = "P-Value")

# Print proteins with raw p-values less than 0.05
signi_protein_raw <- protein_df_1 %>%
  filter(p_value < 0.05)

# Print the filtered dataset with raw p-values less than 0.05
print(signi_protein_raw)


protein_comp1 <- c("FBgn0020513", 
                   "FBgn0036837",
                   "FBgn0033907",
                   "FBgn0261862",
                   "FBgn0019982")

comm_protein_comp1 <- intersect(protein_comp1, signi_protein_raw$gene)
# all protein selected from comp 1 were present in significant protein datsets.


protein_comp_2 <- c("FBgn0011205",
"FBgn0038870",
"FBgn0038319",
"FBgn0024285",
"FBgn0027567")


comm_protein_comp2 <- intersect(protein_comp_2, signi_protein_raw$gene)

#not matching anything

Protein_comp_3 <- c("FBgn0027866", 	"FBgn0001218", "FBgn0023529", 	"FBgn0263260",
                    "FBgn0052069")

comm_protein_comp3 <- intersect(Protein_comp_3, signi_protein_raw$gene)

#not matching any proteins

#dataframe containing only protein from comp 1

comp1_protein_dataframe <- signi_protein_raw %>%
  filter(gene %in% protein_comp1)


genes_comp1 <- c("FBgn0032078",  
"FBgn0028424",  
"FBgn0052549",               
"FBgn0051100",  
"FBgn0261283",  
"FBgn0052758",  
"FBgn0260945",  
"FBgn0028988",  
"FBgn0267668",  
"FBgn0001248")


comm_genes_comp1 <- intersect(genes_comp1, signi_protein_raw$gene)
#no match

# do a PCA for protein_df_1
colnames(protein_df_1)

pca_dataframe <- protein_df_1[, 1:19]

pca_dataframe2 <- t(pca_dataframe)

new_colnames <- as.character(pca_dataframe2[1, ])

# Remove the first row from the dataframe
pca_dataframe2 <- pca_dataframe2[-1, ]

pca_dataframe2 <- as.data.frame(pca_dataframe2)

# Assign the new column names to the dataframe
colnames(pca_dataframe2) <- new_colnames

pca_raw_prot <- pca(pca_dataframe2)

plotIndiv(pca_raw_prot)


#shows a clear batch effect. Might need to address

Genes_comp2 <- c("FBgn0261041", 
                 "FBgn0027836",
                 "FBgn0039584", 
                 "FBgn0053679",  
                 "FBgn0263405", 
                 "FBgn0267734", 
                 "FBgn0259241", 
                 "FBgn0267124", 
                 "FBgn0259222", 
                 "FBgn0039178")

comm_genes_comp2 <- intersect(Genes_comp2, signi_protein_raw$gene)

#FBgn0039178

Genes_comp_3 <- c("FBgn0027590", 
                  "FBgn0034290",
                  "FBgn0038250",
                  "FBgn0031315",
                  "FBgn0037063",  
                  "FBgn0031633",
                  "FBgn0032900",
                  "FBgn0033179",
                  "FBgn0034936",
                  "FBgn0031639")



comm_genes_comp3 <- intersect(Genes_comp_3, signi_protein_raw$gene)

#FBgn0031639"













