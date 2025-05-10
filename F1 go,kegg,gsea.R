library(ggplot2)
library(limma)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(dplyr)
library(pathview)

packageVersion("clusterProfiler")
packageVersion("org.Sc.sgd.db")

#read data
data_80 <- read.csv("G:/项目/data_80.csv")
#select the "Ca" group and "Na" group
data_82 <- data_80[, c(2:7, 16:19)]
colnames(data_82)

#Filtering of data removes expression in Ca or Na groups with absolute values less than 0.8
filtered_data <- data_82[!(apply(abs(data_82[, 3:6]), 1, function(x) all(x < 0.8)) &
                             apply(abs(data_82[, 7:10]), 1, function(x) all(x < 0.8))), ]
cat("保留的行数:", nrow(filtered_data), "\n")

#save the file
write.csv(filtered_data, "G:/项目/filter_data.csv")

# Extract YORF columns as row names
rownames(filtered_data) <- filtered_data$YORF
filtered_data1 <- filtered_data[, c(3:10)]

#Transposing data in preparation for PCA analysis
data_for_pca <- t(filtered_data1)

#pCA analysis
pca_result <- prcomp(data_for_pca, scale. = TRUE, center = FALSE)
pca_df <- as.data.frame(pca_result$x)
#Extract row names as Sample columns
pca_df$Sample <- rownames(pca_df)

#Characters with “Na” are extracted as Na group and characters with “Ca” are extracted as Ca group.
pca_df$Group <- ifelse(grepl("^Ca", pca_df$Sample), "Ca", "Na")

#Plotting the results of PCA analysis
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
geom_point(size = 4) +
geom_text(vjust = -1, size = 3) +
labs(
 title = "PCA of Yeast under Ca/Na Treatments",
  x = paste0("PC1 (", round(summary(pca_result)$importance[2,1] * 100, 1), "%)"),
  y = paste0("PC2 (", round(summary(pca_result)$importance[2,2] * 100, 1), "%)")
) +
theme_minimal() +
scale_color_manual(values = c("Ca" = "blue", "Na" = "orange"))


# Select the expression matrices of the four time points processed by Ca
log_matrix <- filtered_data1[, c("Ca5.", "Ca15.", "Ca30.", "Ca60.")]

# Define the time variable vector
time1 <- c(5, 15, 30, 60)

# Construct a design matrix for the linear model analysis of the expression trend changing over time
design <- model.matrix(~ time1)
print(design)

# Fit the linear model with limma
fit <- lmFit(log_matrix, design)

# Robust estimation of model parameters using empirical Bayesian method
fit <- eBayes(fit)

#Extract the results of the linear model and retrieve all genes according to the time1 coefficient
results <- topTable(fit, coef = "time1", number = Inf)
View(results)

# Significant genes were screened based on the adj. p.al and logFC thresholds
Ca_genes4pa <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 0.05, ]
view(Ca_genes4pa)
Ca_genes4p <- results[results$P.Value < 0.01 & abs(results1$logFC) > 0.05, ]
view(Ca_genes4p)

count(Ca_genes4pa)
count(Ca_genes4p)

# Select the expression matrices of the four time points processed by Na
log_matrix1 <- filtered_data1[, c("Na15.", "Na30.", "Na45.", "Na60.")]

# Define the vector of time variables
time1 <- c(15, 30, 45, 60)

# Build a design matrix for linear model analysis of the expression trend that changes over time
design <- model.matrix(~ time1)
print(design)

# Fit the linear model with limma
fit <- lmFit(log_matrix1, design)

# Robust estimation of model parameters using empirical Bayesian method
fit <- eBayes(fit)

#Extract the results of the linear model and retrieve all genes according to the time1 coefficient
results1 <- topTable(fit, coef = "time1", number = Inf)
View(results1)

# Significant genes were screened based on the adj. p.al and logFC thresholds
Na_genes4pa <- results1[results1$adj.P.Val < 0.05 & abs(results1$logFC) > 0.05, ]
Na_genes4p <- results1[results1$P.Value < 0.01 & abs(results1$logFC) > 0.05, ]

count(Na_genes4p)

#save the file
write.csv(Na_genes4p, "G:/项目/时间趋势分析和log fold change/时间趋势分析/Na_genes4p.csv")
write.csv(Ca_genes4pa, "G:/项目/时间趋势分析和log fold change/时间趋势分析/Ca_genes4pa.csv")
write.csv(Ca_genes4p, "G:/项目/时间趋势分析和log fold change/时间趋势分析/Ca_genes4p.csv")
write.csv(de_genes1, "G:/项目/时间趋势分析和log fold change/log fold change/Ca60-Ca5.csv")
write.csv(de_genes, "G:/项目/时间趋势分析和log fold change/log fold change/Na60-Na15.csv")


#log FC Na60-Na15
expr1 <- data_82[, "Na15."]
expr2 <- data_82[, "Na60."]

# 合并结果为数据框
gseaNa <- data.frame(
  Gene = data_82$YORF,
  sample1 = expr1,
  sample2 = expr2,
  log2FC = expr1 -expr2
)



#log FC Ca60-Ca5
expr1 <- data_82[, "Ca5."]
expr2 <- data_82[, "Ca60."]

gseaCa <- data.frame(
  Gene = data_82$YORF,
  sample1 = expr1,
  sample2 = expr2,
  log2FC = expr1 -expr2
)



#go,kegg lfc
# Extract the expression data at the two time points of Na60 and Na15
expr1 <- filtered_data[, "Na15."]
expr2 <- filtered_data[, "Na60."]

# Construct the result data frame, including the gene name, expression values at two time points, and log2FC (Na15-Na60)
fc_resultNa <- data.frame(
  Gene = filtered_data$YORF,
  sample1 = expr1,
  sample2 = expr2,
  log2FC = expr1 -expr2
)

# Set the threshold of log2 fold change
fc_threshold <- 1  

# Initialize the expression direction of the difference as "NotDE"
fc_resultNa$Direction <- "NotDE"

# Mark as "Up" or "Down" based on the log2FC direction and threshold
fc_resultNa$Direction[fc_resultNa$log2FC >  fc_threshold] <- "Up"
fc_resultNa$Direction[fc_resultNa$log2FC < -fc_threshold] <- "Down"

de_genesNa <- fc_resultNa[fc_resultNa$Direction == "Up" | fc_resultNa$Direction == "Down", ]
count(de_genesNa)

# Extract the expression data at the two time points of Ca60 and Ca5
expr1 <- filtered_data[, "Ca5."]
expr2 <- filtered_data[, "Ca60."]

# Construct the result data frame, including the gene name, expression values at two time points, and log2FC (Ca5-Ca60)
fc_resultCa <- data.frame(
  Gene = filtered_data$YORF,
  sample1 = expr1,
  sample2 = expr2,
  log2FC = expr1 -expr2
)
view(fc_resultCa)

# Set the fold change threshold
fc_threshold <- 1  

# Initialize the expression direction of the difference as "NotDE"
fc_resultCa$Direction <- "NotDE"

# Mark as "Up" or "Down" based on the log2FC direction and threshold
fc_resultCa$Direction[fc_resultCa$log2FC >  fc_threshold] <- "Up"
fc_resultCa$Direction[fc_resultCa$log2FC < -fc_threshold] <- "Down"

de_genesCa <- fc_resultCa[fc_resultCa$Direction == "Up" | fc_resultCa$Direction == "Down", ]
count(de_genesCa)

#Extract the lfc column
de_Ca <- de_genesCa[, c("Gene", "log2FC", "Direction")]
de_Na <- de_genesNa[, c("Gene", "log2FC", "Direction")]

#Add the row name as the Gene column
Ca_genes4pa$Gene <- rownames(Ca_genes4pa)
Ca_genes4p$Gene <- rownames(Ca_genes4p)
Na_genes4pa$Gene <- rownames(Na_genes4pa)
Na_genes4p$Gene <- rownames(Na_genes4p)


#Extract the "Gene" column
Ca_4pa <- Ca_genes4pa$Gene
Ca_4p  <- Ca_genes4p$Gene
Na_4pa <- Na_genes4pa$Gene
Na_4p <- Na_genes4p$Gene

Ca_4p <- as.data.frame(Ca_4p)
colnames(Ca_4p) <- "Gene"
Na_4p <-as.data.frame(Na_4p)
colnames(Na_4p) <- "Gene"

#Extract the intersection and union of the lfc and time trend results
Ca4pa_intersection <- merge(de_Ca, Ca_4pa, by = "Gene")
Ca4p_intersection <- merge(de_Ca, Ca_4p, by = "Gene")
Na4p_intersection <- merge(de_Na, Na_4p, by = "Gene")
Na4p_union <- full_join(de_Na, Na_4p, by = "Gene")
Ca4P_union <- full_join(Ca_4p, de_Ca, by = "Gene")


#save the file
write.csv(Na4p_union, "G:/项目/Na4p_union.csv")
write.csv(Ca4P_union, "G:/项目/Ca4P_union.csv")
write.csv(Ca4p_intersection, "G:/项目/Ca4p_intersection.csv")
write.csv(Na4p_intersection, "G:/项目/Na4p_intersection.csv")

#GSEA
#Extract the gene columns and log2FC columns for preprocessing
fcCa <- gseaCa[, c("Gene", "log2FC")]
fcNa <- gseaNa[, c("Gene", "log2FC")]

# Extract the log2FC vector processed by Ca and set the gene name as the vector name. Sort from high to low
fc_Ca <- fcCa$log2FC
names(fc_Ca ) <- fcCa$Gene
fc_Ca <- sort(fc_Ca, decreasing = TRUE)


fc_Na <- fcNa$log2FC
names(fc_Na) <- fcNa$Gene
fc_Na <- sort(fc_Na, decreasing = TRUE)

#gsea
result_gsea_fcCa <- gseGO(
  geneList = fc_Ca,
  OrgDb = org.Sc.sgd.db, 
  keyType = "ORF",     
  ont = "BP",            
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE,
  eps = 0
)


result_gsea_fcNa <- gseGO(
  geneList = fc_Na,
  OrgDb = org.Sc.sgd.db,  
  keyType = "ORF",    
  ont = "BP",         
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE,
  eps = 0
)

# Convert the GSEA result into a data frame
result_gsea_fcCa <- as.data.frame(result_gsea_fcCa)
result_gsea_fcNa <- as.data.frame(result_gsea_fcNa)

# Save as a CSV file
write.csv(result_gsea_fcCa, file = "G:/项目/result_gsea_fcCa.csv", row.names = TRUE)
write.csv(result_gsea_fcNa, file = "G:/项目/result_gsea_fcNa.csv", row.names = TRUE)


# Draw a specific path
gseaplot2(
  result_gsea_fcCa,
  geneSetID = "GO:0034727",  # 替换为对应的通路 ID
  title = "GSEA Plot for piecemeal microautophagy of the nucleus",
  color = "blue")

gseaplot2(
  result_gsea_fcNa,
  geneSetID = "GO:0034727",  # 替换为对应的通路 ID
  title = "GSEA Plot for piecemeal microautophagy of the nucleus",
  color = "blue")


gseaplot2(
  result_gsea_fcCa,
  geneSetID = "GO:0071038",  # 替换为对应的通路 ID
  title = "GSEA Plot for TRAMP-dependent tRNA surveillance pathway",
  color = "blue")



gseaplot2(
  result_gsea_fcNa,
  geneSetID = "GO:0006450",  # 替换为对应的通路 ID
  title = "GSEA Plot for regulation of translational fidelity",
  color = "blue")

#Extract the results and select the 20 pathways with the highest NES
gsea_table <- result_gsea_fcCa@result
top20 <- gsea_table[order(-abs(gsea_table$NES)), ][1:20, ]

#Draw a bar chart of the 20 pathways with the highest NES
ggplot(top20, aes(x = reorder(Description, NES), y = setSize, fill = NES)) +
  geom_bar(stat = "identity") +  # Draw a bar chart
  coord_flip() +  # Flip the coordinate axes
  scale_fill_gradient(low = "red", high = "blue") +  # Set the color gradient
  labs(
    title = "GSEA Ca Top 20 Enriched Pathways",
    x = "Pathways",
    y = "Set Size",
    fill = "NES"
  ) +
  theme_minimal()

# Calculate the total number of rows in the result data frame
total_rows <- nrow(gsea_table)  
print(paste("总行数:", total_rows))

#The number of rows with positive values in the NES column
positive_nes <- sum(gsea_table$NES > 0) 
#The number of rows with negative values in the NES column
negative_nes <- sum(gsea_table$NES < 0)  

print(paste("NES 列中正值的行数:", positive_nes))
print(paste("NES 列中负值的行数:", negative_nes))

#Read the results and extract the 20 pathways with the highest NES
gsea_table1 <- result_gsea_fcNa@result
top20 <- gsea_table1[order(-abs(gsea_table1$NES)), ][1:20, ]

#Draw a bar chart
ggplot(top20, aes(x = reorder(Description, NES), y = setSize, fill = NES)) +
  geom_bar(stat = "identity") +  
  coord_flip() +  
  scale_fill_gradient(low = "red", high = "blue") +  
  labs(
    title = "GSEA Na Top 20 Enriched Pathways",
    x = "Pathways",
    y = "Set Size",
    fill = "NES"
  ) +
  theme_minimal()

# Calculate the total number of rows in the result data frame
total_rows <- nrow(gsea_table1)  
print(paste("总行数:", total_rows))


#The number of rows with positive values in the NES column
positive_nes <- sum(gsea_table1$NES > 0)
#The number of rows with negative values in the NES column
negative_nes <- sum(gsea_table1$NES < 0)  

print(paste("NES 列中正值的行数:", positive_nes))
print(paste("NES 列中负值的行数:", negative_nes))

#Go,kegg
#Preprocessing, extraction of genes and direction
gsea_Ca4pa <- Ca4pa_intersection[, c("Gene", "Direction")]
gsea_Ca4p <- Ca4p_intersection[, c("Gene", "Direction")]
gsea_Nain <- Na4p_intersection[, c("Gene", "Direction")]
gsea_Naun <- Na4p_union[, c("Gene", "Direction")]


#Separate up and down
up_Ca4pin <- gsea_Ca4p[gsea_Ca4p$Direction == "Up", ]
down_Ca4pin <- gsea_Ca4p[gsea_Ca4p$Direction == "Down", ]

up_Na4pin <- gsea_Nain[gsea_Nain$Direction == "Up", ]
down_Na4pin <- gsea_Nain[gsea_Nain$Direction == "Down", ]


go_upCa4pin <- enrichGO(
  gene          = up_Ca4pin$Gene,      
  OrgDb         = org.Sc.sgd.db,   
  keyType       = "ORF",            # Gene ID type, yeast system name “ORF”
  ont           = "All",            # GO classification with options “BP” (biological process), “MF” (molecular function), “CC” (cellular component)
  pAdjustMethod = "BH",            # p-value correction method 
  pvalueCutoff  = 0.05,           # p-value threshold
  qvalueCutoff  = 0.2,             # q-value threshold
)


go_downCa4pin <- enrichGO(
  gene          = down_Ca4pin$Gene,      
  OrgDb         = org.Sc.sgd.db,   
  keyType       = "ORF",          
  ont           = "All",            
  pAdjustMethod = "BH",          
  pvalueCutoff  = 0.05,            
  qvalueCutoff  = 0.2,          
)

go_upNa4pin <- enrichGO(
  gene          = up_Na4pin$Gene,      
  OrgDb         = org.Sc.sgd.db,   
  keyType       = "ORF",           
  ont           = "All",           
  pAdjustMethod = "BH",           
  pvalueCutoff  = 0.05,          
  qvalueCutoff  = 0.2,        
)

go_downNa4pin <- enrichGO(
  gene          = down_Na4pin$Gene,     
  OrgDb         = org.Sc.sgd.db,   
  keyType       = "ORF",          
  ont           = "All",           
  pAdjustMethod = "BH",          
  pvalueCutoff  = 0.05,          
  qvalueCutoff  = 0.2,           
)

go_upCa4pin_result <- go_upCa4pin@result
go_downCa4pin_result <- go_downCa4pin@result
go_upNa4pin_result <- go_upNa4pin@result
go_downNa4pin_result <- go_downNa4pin@result

write.csv(go_upCa4pin_result, "G:/项目/GO富集/go_upCa4pin_result.csv", row.names = FALSE)
write.csv(go_downCa4pin_result, "G:/项目/GO富集/go_downCa4pin_result.csv", row.names = FALSE)
write.csv(go_upNa4pin_result, "G:/项目/GO富集/go_upNa4pin_result.csv", row.names = FALSE)
write.csv(go_downNa4pin_result, "G:/项目/GO富集/go_downNa4pin_result.csv", row.names = FALSE)

#Plotting the results bar graph
ggplot(go_upCa4pin_result, aes(x = reorder(Description, pvalue), y = Count, fill = ONTOLOGY)) +
  geom_bar(stat = "identity") + 
  coord_flip() +  
  scale_fill_manual(values = c("red", "blue", "yellow")) +  
  labs(
    title = "Go up Ca Enriched Pathways",
    x = "Pathways",
    y = "Gene Count",
    fill = "ONTOLOGY"
  ) +
  theme_bw()

ggplot(go_downCa4pin_result, aes(x = reorder(Description, pvalue), y = Count, fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +  
  coord_flip() +
  scale_fill_manual(values = c("red", "blue", "yellow")) + 
  labs(
    title = "Go down Ca Enriched Pathways",
    x = "Pathways",
    y = "Gene Count",
    fill = "ONTOLOGY"
  ) +
  theme_bw()

ggplot(go_downNa4pin_result, aes(x = reorder(Description, pvalue), y = Count, fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +  
  coord_flip() +  
  scale_fill_manual(values = c("red", "blue", "yellow")) +  
  labs(
    title = "Go down Na Enriched Pathways",
    x = "Pathways",
    y = "Gene Count",
    fill = "ONTOLOGY"
  ) +
  theme_bw()

ggplot(go_upNa4pin_result, aes(x = reorder(Description, pvalue), y = Count, fill = ONTOLOGY)) +
  geom_bar(stat = "identity") + 
  coord_flip() +  
  scale_fill_manual(values = c("red", "blue", "yellow")) +  
  labs(
    title = "Go up Na Enriched Pathways",
    x = "Pathways",
    y = "Gene Count",
    fill = "ONTOLOGY"
  ) +
  theme_bw()

#Kegg
kegg_upCa4pin <- enrichKEGG(
  gene = up_Ca4pin$Gene,        # List of input gene IDs
  organism = 'sce',         # Species code
  pvalueCutoff = 0.05,      #P value threshold
  qvalueCutoff = 0.2        # q-value threshold
)

kegg_downCa4pin <- enrichKEGG(
  gene = down_Ca4pin$Gene,        
  organism = 'sce',        
  pvalueCutoff = 0.05,    
  qvalueCutoff = 0.2      
)

kegg_upNa4pin <- enrichKEGG(
  gene = up_Na4pin$Gene,       
  organism = 'sce',        
  pvalueCutoff = 0.05,     
  qvalueCutoff = 0.2       
)

kegg_downNa4pin <- enrichKEGG(
  gene = down_Na4pin$Gene,        
  organism = 'sce',        
  pvalueCutoff = 0.05,    
  qvalueCutoff = 0.2       
)

#save the file
write.csv(kegg_upCa4pin, "G:/项目/Kegg分析/kegg_upCa4pin.csv", row.names = FALSE)
write.csv(kegg_downCa4pin, "G:/项目/Kegg分析/kegg_downCa4pin.csv", row.names = FALSE)
write.csv(kegg_upNa4pin, "G:/项目/Kegg分析/kegg_upNa4pin.csv", row.names = FALSE)
write.csv(kegg_downNa4pin, "G:/项目/Kegg分析/kegg_downNa4pin.csv", row.names = FALSE)