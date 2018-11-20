
# HD_RNASeq_NormalizationAndDE_CD44.R
# For publication by Osipovitch et al. 2018 on Huntington's disease (HD)-derived glia
# Author and contact: Mikhail Osipovitch (mxo@sund.ku.dk)
#
# The study used two glial progenitor cell types: the bi-potential Glial Progenitor Cells (GPCs)
# identified by the CD140a cell surface marker and the astrocyte-restricted Astrocyte Progenitor
# Cells (APCs) identified by the CD44 cell surface marker. The following script preforms the
# analysis for the CD44+ APCs
#
# The following script performs RUVSeq-based normalization and differential expression analysis
# of RNA-Seq data produced from HD- and control-derived human glial cells. We recommend creating
# a new R project in a separate folder. The output will then be saved to that folder. The procedure
# is accomplished in three steps:
#
# 1 - first-pass differential expression analysis for determination of in silico negative control
#     genes (FDR-corrected P Value > 0.75) that are not affected by the condition of interest
# 2 - calculation of variance normalization factors by RUVg function of RUVSeq
# 3 - second-pass differential expression analysis (5% FDR and and no fold change threshold) for
#     determination of disease-dysregulated genes using the original counts but with adjusting for
#     RUVg-calculated variance factors by multi-factor GLM models implemented in edgeR and DESeq2

### load required libraries:
library(gplots)
library(ape)
library(scatterplot3d)
library(RColorBrewer)
library(DESeq2)
library(edgeR)
library(RUVSeq)
library(VennDiagram)

##########################################################################################################################
### DATASET SPECIFICATION ################################################################################################
##########################################################################################################################

experiment_title = "HD Differentiation Delay - CD44+ APCs"

### read counts and sample data sheet
countData_allSamples  <- read.table("countData_GeneaHD_CD44.txt", header = TRUE, sep = "\t", row.names = 1)
sampleData_allSamples <- read.table("sampleData_GeneaHD_CD44.txt", header = TRUE, sep = "\t", row.names = 1)

### specify colors for cell lines to be used in plots and figures:
col_02 <- rep("#636363", 4)
col_19 <- rep("#cccccc", 4)
col_17 <- rep("#6497b1", 5)
col_18 <- rep("#005b96", 3)
col_20 <- rep("#03396c", 4)


for (which_comparison in 1:4) {
  
  if (which_comparison == 1) {
    
    ### suffix is used in filenames to distinguish output for different comparisons:
    suffix <- " - HD17vCTR_CD44"
    
    ### subset countData:
    countData <- countData_allSamples[,c(1:8, 9:13)]
    
    ### subset sampleData:
    sampleData <- sampleData_allSamples[c(1:8, 9:13),]
    
    ### combined vector of colors to be used in plots and figures:
    cell_line_colors <- c(col_02, col_19, col_17)
    
    ### factors for edgeR differential expression analysis (1 = CTR, 2 = HD):
    x <- factor(c(rep(1, 8), rep(2, 5)))
    
  } else if (which_comparison == 2) {
    
    suffix <- " - HD18vCTR_CD44"
    countData <- countData_allSamples[,c(1:8, 14:16)]
    sampleData <- sampleData_allSamples[c(1:8, 14:16),]
    cell_line_colors <- c(col_02, col_19, col_18)
    x <- factor(c(rep(1, 8), rep(2, 3)))
    
  } else if (which_comparison == 3) {
    
    suffix <- " - HD20vCTR_CD44"
    countData <- countData_allSamples[,c(1:8, 17:20)]
    sampleData <- sampleData_allSamples[c(1:8, 17:20),]
    cell_line_colors <- c(col_02, col_19, col_20)
    x <- factor(c(rep(1, 8), rep(2, 4)))
    
  } else if (which_comparison == 4) {
    
    suffix <- " - HD20vCTR19_CD44"
    countData <- countData_allSamples[,c(5:8, 17:20)]
    sampleData <- sampleData_allSamples[c(5:8, 17:20),]
    cell_line_colors <- c(col_19, col_20)
    x <- factor(c(rep(1, 4), rep(2, 4)))
    
  } # end else if
  
  
  ##########################################################################################################################
  ### START of ANALYSIS: FILTER DATA and MAKE PRELIMINARY PLOTS  ###########################################################
  ##########################################################################################################################
  
  ### make necessary folders to store output:
  
  output_folder <- paste("Output", suffix, "/", sep = "")
  dir.create(output_folder)
  
  output_rle <- paste(output_folder, "Output - RLE Plots/", sep = "")
  output_pca <- paste(output_folder, "Output - PCA Plots/", sep = "")
  output_matrices <- paste(output_folder, "Output - Normalized Matrices/", sep = "")
  dir.create(output_rle)
  dir.create(output_pca)
  dir.create(output_matrices)
  
  ### filter out lowly expressed transcripts:
  countData <- countData[apply(countData, 1, function(x) length(x[x >= 5]) >= 3),]
  
  ### write out filtered count matrix (NN = No Normalization):
  write.table(countData, paste(output_matrices,"countData - NN", suffix, ".txt", sep = ""), sep = "\t", quote = FALSE)
  
  ### UQ Normalization to account for library sizes:
  uq <- as.data.frame(betweenLaneNormalization(as.matrix(countData), which = "upper"))
  
  ### write out UQ count matrix (UQ = Upper Quartile):
  write.table(uq, paste(output_matrices, "countData - UQ", suffix, ".txt", sep = ""), sep = "\t", quote = FALSE)
  
  ### RLE plot (NN = No Normalization):
  png(paste(output_rle, "RLE - NN", suffix, ".png", sep = ""), width = 660, height = 440)
  par(cex = 1.2)
  plotRLE(as.matrix(countData), outline = FALSE, ylim = c(-2, 2.5), col = cell_line_colors, las = 2, ylab = "Relative Log Expression")
  graphics.off()
  
  ### 2D PCA Plot (NN = No Normalization):
  png(paste(output_pca, "PCA 2D - NN", suffix, ".png", sep = ""), width = 660, height = 440)
  par(cex = 1.3)
  plotPCA(as.matrix(countData), col = cell_line_colors, xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
  graphics.off()
  
  ### RLE plot (UQ = Upper Quartile):
  png(paste(output_rle, "RLE - UQ", suffix, ".png", sep = ""), width = 660, height = 440)
  par(cex = 1.2)
  plotRLE(as.matrix(uq), outline = FALSE, ylim = c(-2, 2.5), col = cell_line_colors, las = 2, ylab = "Relative Log Expression")
  graphics.off() 
  
  ### 2D PCA Plot (UQ = Upper Quartile):
  png(paste(output_pca, "PCA 2D - UQ", suffix, ".png", sep = ""), width = 660, height = 440)
  par(cex = 1.3)
  plotPCA(as.matrix(uq), col = cell_line_colors, xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
  graphics.off()
  
  
  ##########################################################################################################################
  ### STEP 1: FIRST PASS DIFFERENTIAL EXPRESSION ANALYSIS ##################################################################
  ##########################################################################################################################
  
  ### make necessary folders to store output:
  dir.create(paste(output_folder, "Output - DE First Pass/", sep  = ""))
  
  output_path_int <- paste(output_folder, "Output - DE First Pass/DEGs_INTs/", sep = "")
  output_path_deseq2 <- paste(output_folder, "Output - DE First Pass/DEGs_DESeq2/", sep = "")
  output_path_edger <- paste(output_folder, "Output - DE First Pass/DEGs_edgeR/", sep = "")
  output_path_allde <- paste(output_folder, "Output - DE First Pass/All DE Results Tables/", sep = "")
  output_path_pvalhist <- paste(output_folder, "Output - P Value Histograms/", sep = "")
  
  dir.create(output_path_int)
  dir.create(output_path_deseq2)
  dir.create(output_path_edger)
  dir.create(output_path_allde)
  dir.create(output_path_pvalhist)
  
  ##### DESeq2 Differential Expresstion:
  
  ### Create a DESeq2 dataset - no normalization
  labels <- unique(sampleData$condition)
  dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ condition)
  dds$condition <- factor(dds$condition, levels = labels)
  
  ### DESEq2 single factor: ~ condition
  dds <- DESeq(dds)
  res <- results(dds)
  res <- res[order(res$padj),]
  res_matrix <- as.matrix(res)
  
  ### write out complete differential expression results from DESeq2:
  write.table(res_matrix, paste(output_path_allde, "All DE Results - DESeq2 - First Pass", suffix, ".txt"), sep = "\t", quote = FALSE)
  
  degs_deseq2_001 <- as.data.frame(res[which(res$padj < 0.01), ])
  degs_deseq2_005 <- as.data.frame(res[which(res$padj < 0.05), ])
  degs_deseq2_010 <- as.data.frame(res[which(res$padj < 0.10), ])
  
  ### P Value histogram for DESeq2 results:
  h1 <- hist(res_matrix[,6], breaks=0:50/50, plot=FALSE)
  png(paste(output_path_pvalhist ,"PVal Hist DESeq2 - NN", suffix, ".png", sep = ""),
      width = 660, height = 440)
  par(cex = 1.3)
  barplot(height = h1$counts, beside = FALSE, ylim = c(0, 5000),
          col = "#d2d4dc", space = 0, ylab="Frequency")
  text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
       adj = c(0.5,1.7), xpd=NA)
  graphics.off()
  
  ##### edgeR Differential Expression:
  
  ### edgeR single factor: ~ condition (~ x)
  design <- model.matrix(~ x, data = sampleData)
  y <- DGEList(counts = uq, group = x)
  y <- calcNormFactors(y)
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit)
  all_tags <- topTags(lrt, n = nrow(uq))
  all_tags <- all_tags$table
  all_tags <- all_tags[order(all_tags$FDR),]
  
  ### write out complete differential expression results from edgeR:
  write.table(all_tags, paste(output_path_allde, "All DE Results - edgeR - First Pass", suffix, ".txt"), sep = "\t", quote = FALSE)
  
  degs_edger_001 <- all_tags[which(all_tags$FDR < 0.01), ]
  degs_edger_005 <- all_tags[which(all_tags$FDR < 0.05), ]
  degs_edger_010 <- all_tags[which(all_tags$FDR < 0.10), ]
  
  ### P Value histogram for edgeR results:
  h1 <- hist(all_tags$FDR, breaks=0:50/50, plot=FALSE)
  png(paste(output_path_pvalhist, "PVal Hist edgeR - NN", suffix, ".png", sep = ""),
      width = 660, height = 440)
  par(cex = 1.3)
  barplot(height = h1$counts, beside = FALSE, ylim = c(0, 5000),
          col = "#d2d4dc", space = 0, ylab="Frequency")
  text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
       adj = c(0.5,1.7), xpd=NA)
  graphics.off()
  
  ### write out individual and intersection tables from DESeq2 and edgeR:  
  write.table(as.data.frame(degs_deseq2_001), paste(output_path_deseq2, "degs_deseq2_001", suffix,".txt", sep = ""), sep = "\t", quote = FALSE)
  write.table(as.data.frame(degs_deseq2_005), paste(output_path_deseq2, "degs_deseq2_005", suffix,".txt", sep = ""), sep = "\t", quote = FALSE)
  write.table(as.data.frame(degs_deseq2_010), paste(output_path_deseq2, "degs_deseq2_010", suffix,".txt", sep = ""), sep = "\t", quote = FALSE)
  
  write.table(degs_edger_001, paste(output_path_edger, "degs_edgeR_001", suffix,".txt", sep = ""), sep = "\t", quote = FALSE)
  write.table(degs_edger_005, paste(output_path_edger, "degs_edgeR_005", suffix,".txt", sep = ""), sep = "\t", quote = FALSE)
  write.table(degs_edger_010, paste(output_path_edger, "degs_edgeR_010", suffix,".txt", sep = ""), sep = "\t", quote = FALSE)
  
  int_001 <- intersect(rownames(degs_deseq2_001), rownames(degs_edger_001))
  int_005 <- intersect(rownames(degs_deseq2_005), rownames(degs_edger_005))
  int_010 <- intersect(rownames(degs_deseq2_010), rownames(degs_edger_010))
  
  write.table(cbind(degs_deseq2_001[int_001,], degs_edger_001[int_001,]), paste(output_path_int, "DEGs_001", suffix, ".txt", sep = ""),
              quote = FALSE, sep = "\t")
  write.table(cbind(degs_deseq2_005[int_005,], degs_edger_005[int_005,]), paste(output_path_int, "DEGs_005", suffix, ".txt", sep = ""),
              quote = FALSE, sep = "\t")
  write.table(cbind(degs_deseq2_010[int_010,], degs_edger_010[int_010,]), paste(output_path_int, "DEGs_010", suffix, ".txt", sep = ""),
              quote = FALSE, sep = "\t")
  
  ##### select Negative Controls:
  
  PVal_T = 0.75
  
  cutoff_deseq2 <- length(which(res_matrix[,6] < PVal_T)) 
  cutoff_edger <- length(which(all_tags$FDR < PVal_T))
  
  ### sorted by increasing p value, take the bottom rows:
  emp_deseq2 <- row.names(res_matrix[c(cutoff_deseq2:nrow(res_matrix)),])
  emp_edger <- row.names(all_tags[c(cutoff_edger:nrow(all_tags)),])
  
  ### final list of control genes, intersection of DESeq2 and edgeR:
  empirical <- intersect(emp_deseq2, emp_edger)
  
  ### write out the empirical control genes with DESeq2 (res_matrix) and edgeR (all_tags) DE data:
  write.table(cbind(res_matrix[empirical,], all_tags[empirical,]),
              paste(output_folder, "Output - DE First Pass/Negative Empiricals", suffix, ".txt", sep = ""),
              quote = FALSE, sep = "\t")
  
  
  ##########################################################################################################################
  ### STEP 2: RUVg NORMALIZATION BY NEGATIVE CONTROLS ######################################################################
  ##########################################################################################################################
  
  ### make necessary folders to store output:
  dir.create(paste(output_folder, "Output - DE Second Pass/", sep = ""))
  
  output_path_int <- paste(output_folder, "Output - DE Second Pass/DEGs_INTs/", sep = "")
  output_path_deseq2 <- paste(output_folder, "Output - DE Second Pass/DEGs_DESeq2/", sep = "")
  output_path_edger <- paste(output_folder, "Output - DE Second Pass/DEGs_edgeR/", sep = "")
  output_path_allde <- paste(output_folder, "Output - DE Second Pass/All DE Results Tables/", sep = "")
  
  dir.create(output_path_int)
  dir.create(output_path_deseq2)
  dir.create(output_path_edger)
  dir.create(output_path_allde)
  
  
  i = 1 # value for parameter k in RUVg function, also used in output file names
  
  ### call to RUVg function for normalization:
  s <- RUVg(as.matrix(uq), empirical, k = i)
  
  ### 2D PCA plot:
  png(paste(output_pca,"PCA 2D - UQ_RUVg_", i, suffix, ".png", sep = ""), width = 660, height = 440)
  par(cex = 1.3)
  plotPCA(s$normalizedCounts, col = cell_line_colors, xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))
  graphics.off()
  
  ### RLE Plot:
  png(paste(output_rle, "RLE - UQ_RUVg_", i, suffix, ".png", sep = ""), width = 660, height = 440)
  par(cex = 1.2)
  plotRLE(s$normalizedCounts, outline = FALSE, ylim = c(-2, 2.5), col = cell_line_colors, las = 2, ylab = "Relative Log Expression")
  graphics.off()
  
  ### write out RUVg-normalized matrix:
  write.table(s$normalizedCounts, paste(output_matrices, "countData - UQ_RUVg_", i, suffix, ".txt", sep = ""), sep = "\t", quote = FALSE)
  
  
  ##########################################################################################################################
  ### STEP 3: SECOND PASS DIFFERENTIAL EXPRESSION ##########################################################################
  ##########################################################################################################################
  
  ### make a copy of sampleData with appended normalization factores from RUVg:
  sampleDataW <- cbind(sampleData, s$W)
  
  ##### DESeq2 Differential Expresstion:
  
  ### create a DESeq2 dataset
  ddsW <- DESeqDataSetFromMatrix(countData = countData, colData = sampleDataW, design = ~ condition)
  ddsW$condition <- factor(ddsW$condition, levels = labels)
  
  ### DESeq2 multi factor:
  design(ddsW) <- formula(~ W_1 + condition)
  ddsW <- DESeq(ddsW)
  res <- results(ddsW)
  res <- res[order(res$padj),]
  res_matrix <- as.matrix(res)
  
  ### write out complete differential expression results from DESeq2:
  write.table(res_matrix, paste(output_path_allde, "All DE Results - DESeq2 - UQ_RUVg_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
  
  degs_deseq2_001 <- as.data.frame(res[which(res$padj < 0.01), ])
  degs_deseq2_005 <- as.data.frame(res[which(res$padj < 0.05), ])
  degs_deseq2_010 <- as.data.frame(res[which(res$padj < 0.10), ])
  
  ### P Value histogram for DESeq2 results:
  h1 <- hist(res_matrix[,6], breaks=0:50/50, plot=FALSE)
  png(paste(output_path_pvalhist, "PVal Hist DESeq2 - UQ_RUVg_", i, ".png", sep = ""),
      width = 660, height = 440)
  par(cex = 1.3)
  barplot(height = h1$counts, beside = FALSE, ylim = c(0, 5000),
          col = "#d2d4dc", space = 0, ylab="Frequency")
  text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
       adj = c(0.5,1.7), xpd=NA)
  graphics.off()
  
  ##### edgeR Differential Expression:
  
  ### edgeR multi factor
  design <- model.matrix(~x + W_1, data = sampleDataW)
  y <- DGEList(counts = uq, group = x)
  y <- calcNormFactors(y)
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef = 2)
  all_tags <- topTags(lrt, n = nrow(uq))
  all_tags <- all_tags$table
  all_tags <- all_tags[order(all_tags$FDR),]
  
  ### write out complete differential expression results from edgeR:
  write.table(all_tags, paste(output_path_allde, "All DE Results - edgeR - UQ_RUVg_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
  
  degs_edger_001 <- all_tags[which(all_tags$FDR < 0.01), ]
  degs_edger_005 <- all_tags[which(all_tags$FDR < 0.05), ]
  degs_edger_010 <- all_tags[which(all_tags$FDR < 0.10), ]
  
  ### p value histogram for edgeR results:
  h1 <- hist(all_tags$FDR, breaks=0:50/50, plot=FALSE)
  png(paste(output_path_pvalhist, "PVal Hist edgeR - UQ_RUVg_", i, ".png", sep = ""),
      width = 660, height = 440)
  par(cex = 1.3)
  barplot(height = h1$counts, beside = FALSE, ylim = c(0, 5000),
          col = "#d2d4dc", space = 0, ylab="Frequency")
  text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
       adj = c(0.5,1.7), xpd=NA)
  graphics.off()
  
  ### write out individual and intersection tables from DESeq2 and edgeR:
  write.table(as.data.frame(degs_deseq2_001), paste(output_path_deseq2, "degs_deseq2_001 - W_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
  write.table(as.data.frame(degs_deseq2_005), paste(output_path_deseq2, "degs_deseq2_005 - W_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
  write.table(as.data.frame(degs_deseq2_010), paste(output_path_deseq2, "degs_deseq2_010 - W_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
  
  write.table(degs_edger_001, paste(output_path_edger, "degs_edgeR_001 - W_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
  write.table(degs_edger_005, paste(output_path_edger, "degs_edgeR_005 - W_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
  write.table(degs_edger_010, paste(output_path_edger, "degs_edgeR_010 - W_", i, ".txt", sep = ""), sep = "\t", quote = FALSE)
  
  ### obtain intersects:
  int_001 <- intersect(rownames(degs_deseq2_001), rownames(degs_edger_001))
  int_005 <- intersect(rownames(degs_deseq2_005), rownames(degs_edger_005))
  int_010 <- intersect(rownames(degs_deseq2_010), rownames(degs_edger_010))
  
  write.table(cbind(degs_deseq2_001[int_001,], degs_edger_001[int_001,]), paste(output_path_int, "DEGs_001 - W_", i, ".txt", sep = ""),
              quote = FALSE, sep = "\t")
  write.table(cbind(degs_deseq2_005[int_005,], degs_edger_005[int_005,]), paste(output_path_int, "DEGs_005 - W_", i, ".txt", sep = ""),
              quote = FALSE, sep = "\t")
  write.table(cbind(degs_deseq2_010[int_010,], degs_edger_010[int_010,]), paste(output_path_int, "DEGs_010 - W_", i, ".txt", sep = ""),
              quote = FALSE, sep = "\t")
  
  
  ### write out the sampleData with appended normalization factors from RUVg:
  write.table(sampleDataW, paste(output_folder, "Output - DE Second Pass/sampleDataW", suffix, ".txt", sep = ""),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
} # end for

##########################################################################################################################
### GENE SET INTERSECTIONS - FDR 5% - NO FC THRESHOLD  ###################################################################
##########################################################################################################################

degs_005_17vsCTR <- read.table("Output - HD17vCTR_CD44/Output - DE Second Pass/DEGs_INTs/DEGs_005 - W_1.txt", sep = "\t")
degs_005_18vsCTR <- read.table("Output - HD18vCTR_CD44/Output - DE Second Pass/DEGs_INTs/DEGs_005 - W_1.txt", sep = "\t")
degs_005_20vsCTR <- read.table("Output - HD20vCTR_CD44/Output - DE Second Pass/DEGs_INTs/DEGs_005 - W_1.txt", sep = "\t")
degs_005_20vs19  <- read.table("Output - HD20vCTR19_CD44/Output - DE Second Pass/DEGs_INTs/DEGs_005 - W_1.txt", sep = "\t")


degs_005_fc0_17vsCTR_up <- degs_005_17vsCTR[which(degs_005_17vsCTR$logFC >  0),]
degs_005_fc0_17vsCTR_dn <- degs_005_17vsCTR[which(degs_005_17vsCTR$logFC < -0),]

degs_005_fc0_18vsCTR_up <- degs_005_18vsCTR[which(degs_005_18vsCTR$logFC >  0),]
degs_005_fc0_18vsCTR_dn <- degs_005_18vsCTR[which(degs_005_18vsCTR$logFC < -0),]

degs_005_fc0_20vsCTR_up <- degs_005_20vsCTR[which(degs_005_20vsCTR$logFC >  0),]
degs_005_fc0_20vsCTR_dn <- degs_005_20vsCTR[which(degs_005_20vsCTR$logFC < -0),]

degs_005_fc0_20vs19_up <- degs_005_20vs19[which(degs_005_20vs19$logFC >  0),]
degs_005_fc0_20vs19_dn <- degs_005_20vs19[which(degs_005_20vs19$logFC < -0),]


output_folder <- "Output - DE Genes Intersections - CD44/"
dir.create(output_folder)

labels <- c("17vsCTR","18vsCTR","20vsCTR","20vs19")


###### Gene set intersection and Venn diagram for up-regulated genes, 5% FDR and FC > 0:

l1 <- rownames(degs_005_fc0_17vsCTR_up)
l2 <- rownames(degs_005_fc0_18vsCTR_up)
l3 <- rownames(degs_005_fc0_20vsCTR_up)
l4 <- rownames(degs_005_fc0_20vs19_up)

png(paste(output_folder, "Venn 005_FC0 - 4 Lists - UP.png", sep = ""), width = 800, height = 700, units = "px")  
vd <- venn.diagram(list(v1 = l1, v2 = l2, v3 = l3, v4 = l4),
                   main = "Up-regulated genes in each of 3 HD hGPC lines\ncompared to control hGPC lines\n5% FDR and FC > 0",
                   main.cex = 2.5, main.fontfamily = "Verdana",
                   category.names = labels, fill = c("#adcbe3","#adcbe3","#adcbe3","#adcbe3"),
                   fontfamily = "Verdana", cat.fontfamily = "Verdana", cex = 3, cat.cex = 2.7,
                   margin = c(0.1,0.1), filename = NULL)
grid.draw(vd)
graphics.off()

### compile the intersection gene list with fold changes and P Values calculated by edgeR:
degs_int_005_fc0_up <- cbind(degs_005_fc0_17vsCTR_up[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)],
                             degs_005_fc0_18vsCTR_up[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)],
                             degs_005_fc0_20vsCTR_up[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)],
                             degs_005_fc0_20vs19_up[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)])

###### Gene set intersection and Venn diagram for down-regulated genes, 5% FDR and FC > 0:

l1 <- rownames(degs_005_fc0_17vsCTR_dn)
l2 <- rownames(degs_005_fc0_18vsCTR_dn)
l3 <- rownames(degs_005_fc0_20vsCTR_dn)
l4 <- rownames(degs_005_fc0_20vs19_dn)

png(paste(output_folder, "Venn 005_FC0 - 4 Lists - DN.png", sep = ""), width = 800, height = 700, units = "px")  
vd <- venn.diagram(list(v1 = l1, v2 = l2, v3 = l3, v4 = l4),
                   main = "Down-regulated genes in each of 3 HD hGPC lines\ncompared to control hGPC lines\n5% FDR and FC > 0",
                   main.cex = 2.5, main.fontfamily = "Verdana",
                   category.names = labels, fill = c("#adcbe3","#adcbe3","#adcbe3","#adcbe3"),
                   fontfamily = "Verdana", cat.fontfamily = "Verdana", cex = 3, cat.cex = 2.7,
                   margin = c(0.1,0.1), filename = NULL)
grid.draw(vd)
graphics.off()

### compile the intersection gene list with fold changes and P Values calculated by edgeR:
degs_int_005_fc0_dn <- cbind(degs_005_fc0_17vsCTR_dn[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)],
                             degs_005_fc0_18vsCTR_dn[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)],
                             degs_005_fc0_20vsCTR_dn[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)],
                             degs_005_fc0_20vs19_dn[Reduce(intersect, list(v1 = l1, v2 = l2, v3 = l3, v4 = l4)),c(7,11)])


degs_int_005_fc0 <- rbind(degs_int_005_fc0_up, degs_int_005_fc0_dn)

colnames(degs_int_005_fc0) <- c("logFC_HD17vsCTR","FDR_HD17vsCTR",
                                "logFC_HD18vsCTR","FDR_HD18vsCTR",
                                "logFC_HD20vsCTR","FDR_HD20vsCTR",
                                "logFC_HD20vs19","FDR_HD20vs19")

### write out the compiled intersection DE gene list:
write.table(degs_int_005_fc0, paste(output_folder, "DEGs 4 Lists - FDR005_FC0.txt", sep = ""), sep = "\t", quote = FALSE)

