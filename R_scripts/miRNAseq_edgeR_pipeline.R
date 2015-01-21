#############################################################
#############################################################
##                                                         ##
##  miRNA-seq data analysis of sense counts (paired data)  ##
##                                                         ##
#############################################################
#############################################################

# Analysis of 10 animals infected over a 15 weeks time course for which 
# miRNA-seq libraries where prepared from serum samples

#############################
# List of required packages #
#############################

# Source the common functions used across this script
source(file = paste("C:/Users/nnalpas/Dropbox/Home_work_sync/Work", 
                     "/Bioinformatics/R/General_function.R", sep = ""))
source(file = paste("C:/Users/Nicolas Nalpas/Dropbox/Home_work_sync/Work", 
                    "/Bioinformatics/R/General_function.R", sep = ""))

# Load the required packages
loadpackage(package = edgeR)
loadpackage(package = MASS)
loadpackage(package = ggplot2)
loadpackage(package = PerformanceAnalytics)
loadpackage(package = VennDiagram)
loadpackage(package = reshape2)

#########################################################################
# Analysis of counts data obtained via Novoalign-featureCounts pipeline #
#########################################################################

################################################
# Read in and concatenate input files within R #
################################################

# Reads and merges a set of files containing counts
Count <- readDGE(files = list.files(
                 path = paste("C:/Users/nnalpas/Dropbox/Home_work_sync/Work/",
                              "TIDA/miRNA-seq/Results/Counts/Novo-feature/", 
                              "mature_miRNA", sep = ""), pattern = "6"),
                 path = paste("C:/Users/nnalpas/Dropbox/Home_work_sync/Work/",
                              "TIDA/miRNA-seq/Results/Counts/Novo-feature/", 
                              "mature_miRNA", sep = ""), columns = c(1,7), 
                 skip = 1, header = TRUE)

# Check the readDGE variable output
names(Count)
Count$samples
head(Count$counts)

############################################################
# Gene annotation using information obtained from GTF file #
############################################################

# Read in the annotation information
miRNA.info <- read.table(file = "../miRNA_Btaurus.txt", header = TRUE,
                        sep = "\t", quote = "")

# Determine which miRNA have identical mature sequence
miRNA.duplicate <- aggregate(gene_id ~ sequence, FUN = "as.vector",
                             data = miRNA.info, na.action = "as.vector")

# Combine new info with miRNA annotation
miRNA.info <- merge(x = miRNA.info, y = miRNA.duplicate, by.x = "sequence",
                   by.y = "sequence", all = TRUE)
miRNA.info <- cbind(miRNA.info[, 2:7], miRNA.info[, 1], 
                    miRNA.info[, 8:ncol(miRNA.info)])
colnames(miRNA.info)[c(1,7,12)] <- c("gene_id", "sequence",
                                     "identical_sequence")
head(miRNA.info)
dim(miRNA.info)

# Merge the annotation information with the count table
Count <- merge(x = miRNA.info, y = Count$counts, 
                         by.x = "gene_id", by.y = "row.names", all = TRUE)
head(Count)

# Check the merged matrix in terms of size and content
dim(Count)
table(duplicated(Count[, "gene_id"]))
length(unique(Count[, "gene_id"], incomparables = FALSE))

# Ouptut samples data
write.matrix(x = Count, file = "Novo-feature_counts.txt",  sep = "\t")

###############################
# Create groups and a DGElist #
###############################

# Create a target matrix and a experimental group vector and animal block
target <- colnames(Count[, (ncol(Count)-69):ncol(Count)])
target <- cbind(target, matrix(data = unlist(strsplit(x = target, split = "_",
                                                      fixed = TRUE)),
                               nrow = 70, ncol = 2, byrow = TRUE))
colnames(target) = c("sample", "animal", "time_point")
group <- factor(x = target[, "time_point"])
animal <- factor(x = target[, "animal"])

# Create a DGElist containing the group information
dgelist <- DGEList(counts = Count[, (ncol(Count)-69):ncol(Count)], 
                         lib.size = NULL, norm.factors = NULL, group = group, 
                         genes = Count[, 1:(ncol(Count)-70)], 
                         remove.zeros = FALSE)
rownames(dgelist$counts) <- rownames(
  dgelist$genes) <- dgelist$genes[, "gene_id"]
dgelist$genes["gene_id"] <- NULL
names(dgelist)
head(dgelist$samples)
head(dgelist$counts)
head(dgelist$genes)

# Ouptut samples data
write.table(x = dgelist$samples, file = "Novo-feature_samples.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

###########################################################
# Quality check of libraries by plotting density of count #
###########################################################

# Log2 transform the count data for better visualization
count.log2 <- log(x = (dgelist$counts[,]+1), base = 2)

# Plot density of count for all libraries
png(filename = "Novo-feature_Dens_all_miRNA.png", width = 1500, height = 1300,
    units = "px")
plot(x = density(count.log2[, 1]), main = "Density plot of count per gene", 
     lty =  1, xlab = "Log2 of count per gene", ylab = "Density", col = 1, 
     xlim = c(-0.5, 6), ylim = c(0.0, 0.8), lwd = 2 , cex.axis = 1.5, 
     cex.lab = 1.5, cex.main = 2)
for (i in 2:ncol(count.log2)) {
  lines(x = density(count.log2[, i]), lty = 1, col = i)
}
dev.off()

#####################################
# Filtering of lowly expressed tags #
#####################################

# Identify genes with zero counts across all samples
dim(dgelist[rowSums(dgelist$counts) == 0, ])
head(dgelist[rowSums(dgelist$counts) == 0, ])

# Filter lowly expressed tags, retaining only tags with at least 50 counts per 
# million in 10 or more libraries (10 libraries correspond to one time point)
dgelist.filt <- dgelist[rowSums(cpm(dgelist$counts) > 50) >= 10, ]
dim(dgelist.filt$counts)

# Compute the library size
dgelist.filt$samples$lib.size <- colSums(dgelist.filt$counts)
head(dgelist.filt$samples)
head(dgelist$samples)

###########################################################################
# Quality check of libraries by plotting density of count after filtering #
###########################################################################

# Log2 transform the count data for better visualization
count.log2 <- log(x = (dgelist.filt$counts[,]+1), base = 2)

# Plot density of count for all libraries
png(filename = "Novo-feature_Dens_filt_miRNA.png", width = 1366, height = 768,
    units = "px")
plot(x = density(count.log2[, 1]), main = "Density plot of count per gene", 
     lty =  1, xlab = "Log2 of count per gene", ylab = "Density", col = 1, 
     xlim = c(-0.5, 15), ylim = c(0.0, 0.2), lwd = 2 , cex.axis = 1.5, 
     cex.lab = 1.5, cex.main = 2)
for (i in 2:ncol(count.log2)) {
  lines(x = density(count.log2[, i]), lty = 1, col = i)
}
dev.off()

############################################
# Histogram of CPM per expression category #
############################################

# Prepare the histogram picture output
png(filename = "Novo-feature_Hist_miRNA.png", width = 1500, height = 1500,
    units = "px")
par(mfrow = c(2,2), cex = 3, cex.axis = 0.7)

# Plot histogram of number of gene per CPM category for all microRNA
hist_plot <- hist(x = rowMeans(x = cpm(dgelist$counts), na.rm = FALSE), 
                  breaks = c(
                    0, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000), 
                  plot = FALSE)
plot(x = hist_plot$counts, log = "x", type = "h", axes = FALSE, 
     xlab = "Average CPM", ylab = "Number of gene", lwd = 10, 
     col = "lightblue", main = "All miRNA")
axis(side = 1, at = c(1:10), labels = c(
  "0", "0.01", "0.1", "1", "10", "100", "1000", "10000", "100000", "1000000"))
axis(side = 2)

# Plot histogram of number of gene per CPM category for microRNA with 
# at least one count
hist_plot <- hist(x = rowMeans(x = cpm(dgelist$counts[rowSums(
  dgelist$counts) > 0, ]), na.rm = FALSE), breaks = c(
    0, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000), plot = FALSE)
plot(hist_plot$counts, log = "x", type = "h", axes = FALSE, 
     xlab = "Average CPM", ylab = "Number of gene", lwd = 10, 
     col = "lightgreen", main = "Min one count miRNA")
axis(side = 1, at = c(1:10), labels = c(
  "0", "0.01", "0.1", "1", "10", "100", "1000", "10000", "100000", "1000000"))
axis(side = 2)

# Plot histogram of number of gene per CPM category for microRNA 
# filtered for low expression
hist_plot <- hist(x = rowMeans(x = cpm(dgelist.filt$counts), na.rm = FALSE), 
                  breaks = c(
                    0, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000), 
                  plot = FALSE)
plot(hist_plot$counts, log = "x", type = "h", axes = FALSE, 
     xlab = "Average CPM", ylab = "Number of gene", lwd = 10, 
     col = "lightcoral", main = "Low expression filtered miRNA")
axis(side = 1, at = c(1:10), labels = c(
  "0", "0.01", "0.1", "1", "10", "100", "1000", "10000", "100000", "1000000"))
axis(side = 2)

# Close graphic device
dev.off()

#################################################
# Gene expression correlation between libraries #
#################################################

# Note that CPM values are not required with the use of Spearman rank test

# Perform gene expression correlation between libraries at pre 2 weeks to 
# identify outlier library
png(filename = "Novo-feature_Pre2_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_pre2", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Pre 2 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at pre 1 weeks to 
# identify outlier library
png(filename = "Novo-feature_Pre1_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_pre1", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Pre 1 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 1 weeks to 
# identify outlier library
png(filename = "Novo-feature_Post1_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_1$", x = colnames(dgelist.filt$counts), perl = TRUE)] + 1),
  base = 2), histogram = TRUE, method = "spearman", 
  main = "Post 1 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 2 weeks to 
# identify outlier library
png(filename = "Novo-feature_Post2_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_2", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Post 2 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 6 weeks to 
# identify outlier library
png(filename = "Novo-feature_Post6_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_6", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Post 6 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 10 weeks to 
# identify outlier library
png(filename = "Novo-feature_Post10_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_10", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Post 10 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 12 weeks to 
# identify outlier library
png(filename = "Novo-feature_Post12_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_12", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Post 12 weeks time point expression correlation")
dev.off()

########################################################
# Normalization of data using trimmed mean of M-values #
########################################################

# Calculate normalisation factor for our DGElist, note that with edgeR 
# the counts are not transformed in any way after normalization, instead 
# normalization will modify library size
dgelist.norm <- calcNormFactors(dgelist.filt)
dgelist.norm$samples

################################################
# Multidimensional scaling plot of all samples #
################################################

# Output value for MDS plot (dimension 1 and 2 in this case)
MDS <- plotMDS(x = dgelist.norm, top = 1000000, gene.selection = "pairwise", 
               xlab = "Dimension 1", ylab = "Dimension 2", dim.plot = c(1, 2), 
               cex = 1)

# MDS values are then plotted with ggplot2
MDS.ggplot <- data.frame(target, MDS$x, MDS$y)
tiff(filename = "Novo-feature_MDS_miRNA.tif", width = (
  (range(MDS.ggplot$MDS.x)[2]-range(MDS.ggplot$MDS.x)[1])*11.5), height = (
    (range(MDS.ggplot$MDS.y)[2]-range(MDS.ggplot$MDS.y)[1])*10), units = "in",
     res = 600, compression = "lzw")
ggplot(data = MDS.ggplot, aes(x = MDS.ggplot$MDS.x, y = MDS.ggplot$MDS.y, 
                              colour = MDS.ggplot$time_point)) + 
  geom_point(size = 4.5) + theme(panel.background = element_rect(
    fill = 'wheat'), legend.title = element_text(size = 12, face = "bold"), 
    legend.text = element_text(size = 10, face = "bold"), 
    axis.title = element_text(face = "bold", size = 12), 
    axis.text = element_text(face = "bold", size = 10), 
    plot.title = element_text(face = "bold", size = 15)) + 
  ggtitle("MDS plot") + xlab("Dimension 1") + ylab("Dimension 2") + 
  scale_colour_discrete(name = "Time point\npost-infection", breaks = c(
    "pre2", "pre1", "1", "2", "6", "10", "12"), 
    labels = c("pre 2w", "pre 1w", "1w", "2w", "6w", "10w", "12w"))
dev.off()

# Write into a table the coordinates of each library for the MDS plot
write.table(x = MDS.ggplot, file = "Novo-feature_MDS_xy.txt", sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = TRUE)

#########################################################################
# Determine if there has been inversion between 6522_pre2 and 6526_pre2 #
#########################################################################

# Output value for MDS plot (dimension 1 and 2) for animals 6522 and 6526
MDS <- plotMDS(x = dgelist.norm[, grep(pattern = "652(2|6)",
                                       x = colnames(dgelist.norm),
                                       perl = TRUE)], top = 1000000,
               gene.selection = "pairwise", xlab = "Dimension 1",
               ylab = "Dimension 2", dim.plot = c(1, 2), cex = 1)

# MDS values are then plotted with ggplot2
MDS.ggplot <- data.frame(target[grep(pattern = "652(2|6)",
                                     x = target[, "sample"], perl = TRUE),],
                         MDS$x, MDS$y)
ggplot(data = MDS.ggplot, aes(x = MDS.ggplot$MDS.x, y = MDS.ggplot$MDS.y, 
                              colour = MDS.ggplot$time_point,
                              shape = MDS.ggplot$animal)) + 
  geom_point(size = 8) + theme(panel.background = element_rect(
    fill = 'wheat'), legend.title = element_text(size = 20, face = "bold"), 
    legend.text = element_text(size = 15, face = "bold"), 
    axis.title = element_text(face = "bold", size = 30), 
    axis.text = element_text(face = "bold", size = 20), 
    plot.title = element_text(face = "bold", size = 40)) + 
  ggtitle("MDS plot") + xlab("Dimension 1") + ylab("Dimension 2")

##############################################
# Create a design matrix for paired analysis #
##############################################

# Create a design matrix
design <- model.matrix(~ animal + group)
rownames(design) <- rownames(dgelist.norm$samples)
colnames(design) <- gsub(pattern = "(animal)|(group)", replacement = "",
                         x = colnames(design), perl = TRUE)
design

########################################################################
# Estimate the dispersion parameter for each tag using Cox-Reid method #
########################################################################

# Calculate the dispersion (common, trended and tagwise)
dgelist.disp <- estimateGLMCommonDisp(y = dgelist.norm, design = design,
                                      verbose = TRUE)
dgelist.disp <- estimateGLMTrendedDisp(y = dgelist.disp, design = design)
dgelist.disp <- estimateGLMTagwiseDisp(y = dgelist.disp, design = design)
names(dgelist.disp)

# Plot the dispersion
png(filename = "Novo-feature_BCV.png", width = 1366, height = 768,
    units = "px")
plotBCV(dgelist.disp)
dev.off()

# Show the calculated dispersion and the coefficient of biological variation
dgelist.disp$common.dispersion
sqrt(dgelist.disp$common.dispersion)

##################################################################
# Determine differential expression using negative binomial GLMs #
##################################################################

# Fit a negative binomial generalized linear model for each tag using
# the design matrix and calculated dispersion
dgeglm.fit <- glmFit(y = dgelist.disp, design = design)
names(dgeglm.fit)

################################
# Differential expression call #
################################

# Test for differential expression for -1 week versus -2 week
miR.DE("pre1", 1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "Novofeature.dgelrt.pre1w", dedata = "Novofeature.de.pre1w",
       smearfile = "Novo-feature_smear_pre1w")

# Test for differential expression for 1 week versus -2 week and -1 week
miR.DE(1, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "Novofeature.dgelrt.1w", dedata = "Novofeature.de.1w",
       smearfile = "Novo-feature_smear_1w")

# Test for differential expression for 2 week versus -2 week and -1 week
miR.DE(2, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "Novofeature.dgelrt.2w", dedata = "Novofeature.de.2w",
       smearfile = "Novo-feature_smear_2w")

# Test for differential expression for 6 week versus -2 week and -1 week
miR.DE(6, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "Novofeature.dgelrt.6w", dedata = "Novofeature.de.6w",
       smearfile = "Novo-feature_smear_6w")

# Test for differential expression for 10 week versus -2 week and -1 week
miR.DE(10, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "Novofeature.dgelrt.10w", dedata = "Novofeature.de.10w",
       smearfile = "Novo-feature_smear_10w")

# Test for differential expression for 12 week versus -2 week and -1 week
miR.DE(12, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "Novofeature.dgelrt.12w", dedata = "Novofeature.de.12w",
       smearfile = "Novo-feature_smear_12w")

#########################################################
# Merge all DE call data from the different time points #
#########################################################

# Merge all DE table for the different time points into a single dataframe
DEtable.merge("Novofeature.de.1w", "Novofeature.de.2w", "Novofeature.de.6w",
              "Novofeature.de.10w", "Novofeature.de.12w",
              output = "Novofeature.DE", pattern = "^Novofeature.de.")
head(Novofeature.DE)

# Write into a table the full DE call data
write.matrix(x = Novofeature.DE, file = "Novo-feature_full_DE.txt", sep = "\t")

##############################################
# Comparison of DE genes between time points #
##############################################

# Identify as a vector list the significant DE genes per time point
sig.1w <- as.character(Novofeature.DE[!is.na(Novofeature.DE$FDR_1w) & (
  Novofeature.DE$FDR_1w < 0.05), "gene_id"])
sig.2w <- as.character(Novofeature.DE[!is.na(Novofeature.DE$FDR_2w) & (
  Novofeature.DE$FDR_2w < 0.05), "gene_id"])
sig.6w <- as.character(Novofeature.DE[!is.na(Novofeature.DE$FDR_6w) & (
  Novofeature.DE$FDR_6w < 0.05), "gene_id"])
sig.10w <- as.character(Novofeature.DE[!is.na(Novofeature.DE$FDR_10w) & (
  Novofeature.DE$FDR_10w < 0.05), "gene_id"])
sig.12w <- as.character(Novofeature.DE[!is.na(Novofeature.DE$FDR_12w) & (
  Novofeature.DE$FDR_12w < 0.05), "gene_id"])

# Create the Venn diagram of overlapping DE genes between time points
venn.diagram(x = list("1w" = sig.1w, "2w" = sig.2w, "6w" = sig.6w,
                      "10w" = sig.10w, "12w" = sig.12w),
             filename = "Novo-feature_Venn_All.tiff", na="remove")

# Identify the DE miRNA common to all post-infection time points
Novofeature.overlap <- Reduce(intersect, list(sig.1w, sig.2w, sig.6w, sig.10w, sig.12w))
Novofeature.overlap

############################
# Clean up the R workspace #
############################

# Rename variables to be kept according to the dataset used for generation
Novofeature.dgelist.disp <- dgelist.disp

# Remove temporary and unrequired variables
rm(miRNA.duplicate, target, animal, group, Count, i, count.log2, hist_plot,
   MDS, MDS.ggplot, dgelist, dgelist.filt, design, dgelist.norm, dgeglm.fit, 
   Novofeature.dgelrt.pre1w, Novofeature.de.pre1w, sig.1w, sig.2w, sig.6w,
   sig.10w, sig.12w, dgelist.disp)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
##########################################################
# Analysis of counts data obtained via miRdeep2 pipeline #
##########################################################

################################################
# Read in and concatenate input files within R #
################################################

# Read in the count files via the read.table function (cannot use readDGE)
files <- list.files(path = paste("C:/Users/nnalpas/Dropbox/Home_work_sync/",
                                 "Work/TIDA/miRNA-seq/Results/Counts/mirdeep2",
                                 sep = ""), pattern = "6")
for (i in 1:length(files)) {
  Dat <- read.table(file = paste("C:/Users/nnalpas/Dropbox/Home_work_sync/",
                                  "Work/TIDA/miRNA-seq/Results/Counts/",
                                  "mirdeep2", files[i], sep = "/"), quote = "")
  Dat <- Dat[, c(1:3)]
  sample <- gsub(pattern = "_expressed.csv", replacement = "", x = files[i])
  colnames(Dat) <- c("gene_name", sample, "precursor_name")
  Dat <- merge(x = Dat, y = miRNA.info, by = c("gene_name", "precursor_name"))
  Dat <- Dat[, c("gene_id", sample)]
  if (i == 1) {
    Count <- Dat
  }
  else {
    Count <- merge(x = Count, y = Dat, by = "gene_id")
  }
}
dim(Count)
head(Count)
tail(Count)

######################################################
# Merge the annotation information with counts table #
######################################################

# Merge previously imported annotation information with the count table
Count <- merge(x = miRNA.info, y = Count, by = "gene_id", all = TRUE)
head(Count)

# Check the merged matrix in terms of size and content
dim(Count)
table(duplicated(Count[, "gene_id"]))
length(unique(Count[, "gene_id"], incomparables = FALSE))

# Ouptut samples data
write.matrix(x = Count, file = "miRdeep2_counts.txt",  sep = "\t")

###############################
# Create groups and a DGElist #
###############################

# Create a target matrix and a experimental group vector and animal block
target <- colnames(Count[, (ncol(Count)-69):ncol(Count)])
target <- cbind(target, matrix(data = unlist(strsplit(x = target, split = "_",
                                                      fixed = TRUE)),
                               nrow = 70, ncol = 2, byrow = TRUE))
colnames(target) = c("sample", "animal", "time_point")
group <- factor(x = target[, "time_point"])
animal <- factor(x = target[, "animal"])

# Create a DGElist containing the group information
dgelist <- DGEList(counts = Count[, (ncol(Count)-69):ncol(Count)], 
                   lib.size = NULL, norm.factors = NULL, group = group, 
                   genes = Count[, 1:(ncol(Count)-70)], 
                   remove.zeros = FALSE)
rownames(dgelist$counts) <- rownames(
  dgelist$genes) <- dgelist$genes[, "gene_id"]
dgelist$genes["gene_id"] <- NULL
names(dgelist)
head(dgelist$samples)
head(dgelist$counts)
head(dgelist$genes)

# Ouptut samples data
write.table(x = dgelist$samples, file = "miRdeep2_samples.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

###########################################################
# Quality check of libraries by plotting density of count #
###########################################################

# Log2 transform the count data for better visualization
count.log2 <- log(x = (dgelist$counts[,]+1), base = 2)

# Plot density of count for all libraries
png(filename = "miRdeep2_Dens_all_miRNA.png", width = 1500, height = 1300,
    units = "px")
plot(x = density(count.log2[, 1]), main = "Density plot of count per gene", 
     lty =  1, xlab = "Log2 of count per gene", ylab = "Density", col = 1, 
     xlim = c(-0.5, 5), ylim = c(0.0, 1), lwd = 2 , cex.axis = 1.5, 
     cex.lab = 1.5, cex.main = 2)
for (i in 2:ncol(count.log2)) {
  lines(x = density(count.log2[, i]), lty = 1, col = i)
}
dev.off()

#####################################
# Filtering of lowly expressed tags #
#####################################

# Identify genes with zero counts across all samples
dim(dgelist[rowSums(dgelist$counts) == 0, ])
head(dgelist[rowSums(dgelist$counts) == 0, ])

# Filter lowly expressed tags, retaining only tags with at least 50 counts per 
# million in 10 or more libraries (10 libraries correspond to one time point)
dgelist.filt <- dgelist[rowSums(cpm(dgelist$counts) > 50) >= 10, ]
dim(dgelist.filt$counts)

# Compute the library size
dgelist.filt$samples$lib.size <- colSums(dgelist.filt$counts)
head(dgelist.filt$samples)
head(dgelist$samples)

###########################################################################
# Quality check of libraries by plotting density of count after filtering #
###########################################################################

# Log2 transform the count data for better visualization
count.log2 <- log(x = (dgelist.filt$counts[,]+1), base = 2)

# Plot density of count for all libraries
png(filename = "miRdeep2_Dens_filt_miRNA.png", width = 1366, height = 768,
    units = "px")
plot(x = density(count.log2[, 1]), main = "Density plot of count per gene", 
     lty =  1, xlab = "Log2 of count per gene", ylab = "Density", col = 1, 
     xlim = c(-0.5, 15), ylim = c(0.0, 0.2), lwd = 2 , cex.axis = 1.5, 
     cex.lab = 1.5, cex.main = 2)
for (i in 2:ncol(count.log2)) {
  lines(x = density(count.log2[, i]), lty = 1, col = i)
}
dev.off()

############################################
# Histogram of CPM per expression category #
############################################

# Prepare the histogram picture output
png(filename = "miRdeep2_Hist_miRNA.png", width = 1500, height = 1500,
    units = "px")
par(mfrow = c(2,2), cex = 3, cex.axis = 0.7)

# Plot histogram of number of gene per CPM category for all microRNA
hist_plot <- hist(x = rowMeans(x = cpm(dgelist$counts), na.rm = FALSE), 
                  breaks = c(
                    0, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000), 
                  plot = FALSE)
plot(x = hist_plot$counts, log = "x", type = "h", axes = FALSE, 
     xlab = "Average CPM", ylab = "Number of gene", lwd = 10, 
     col = "lightblue", main = "All miRNA")
axis(side = 1, at = c(1:10), labels = c(
  "0", "0.01", "0.1", "1", "10", "100", "1000", "10000", "100000", "1000000"))
axis(side = 2)

# Plot histogram of number of gene per CPM category for microRNA with 
# at least one count
hist_plot <- hist(x = rowMeans(x = cpm(dgelist$counts[rowSums(
  dgelist$counts) > 0, ]), na.rm = FALSE), breaks = c(
    0, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000), plot = FALSE)
plot(hist_plot$counts, log = "x", type = "h", axes = FALSE, 
     xlab = "Average CPM", ylab = "Number of gene", lwd = 10, 
     col = "lightgreen", main = "Min one count miRNA")
axis(side = 1, at = c(1:10), labels = c(
  "0", "0.01", "0.1", "1", "10", "100", "1000", "10000", "100000", "1000000"))
axis(side = 2)

# Plot histogram of number of gene per CPM category for microRNA 
# filtered for low expression
hist_plot <- hist(x = rowMeans(x = cpm(dgelist.filt$counts), na.rm = FALSE), 
                  breaks = c(
                    0, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000), 
                  plot = FALSE)
plot(hist_plot$counts, log = "x", type = "h", axes = FALSE, 
     xlab = "Average CPM", ylab = "Number of gene", lwd = 10, 
     col = "lightcoral", main = "Low expression filtered miRNA")
axis(side = 1, at = c(1:10), labels = c(
  "0", "0.01", "0.1", "1", "10", "100", "1000", "10000", "100000", "1000000"))
axis(side = 2)

# Close graphic device
dev.off()

#################################################
# Gene expression correlation between libraries #
#################################################

# Note that CPM values are not required with the use of Spearman rank test

# Perform gene expression correlation between libraries at pre 2 weeks to 
# identify outlier library
png(filename = "miRdeep2_Pre2_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_pre2", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Pre 2 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at pre 1 weeks to 
# identify outlier library
png(filename = "miRdeep2_Pre1_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_pre1", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Pre 1 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 1 weeks to 
# identify outlier library
png(filename = "miRdeep2_Post1_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_1$", x = colnames(dgelist.filt$counts), perl = TRUE)] + 1),
  base = 2), histogram = TRUE, method = "spearman", 
  main = "Post 1 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 2 weeks to 
# identify outlier library
png(filename = "miRdeep2_Post2_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_2", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Post 2 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 6 weeks to 
# identify outlier library
png(filename = "miRdeep2_Post6_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_6", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Post 6 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 10 weeks to 
# identify outlier library
png(filename = "miRdeep2_Post10_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_10", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Post 10 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 12 weeks to 
# identify outlier library
png(filename = "miRdeep2_Post12_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_12", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Post 12 weeks time point expression correlation")
dev.off()

########################################################
# Normalization of data using trimmed mean of M-values #
########################################################

# Calculate normalisation factor for our DGElist, note that with edgeR 
# the counts are not transformed in any way after normalization, instead 
# normalization will modify library size
dgelist.norm <- calcNormFactors(dgelist.filt)
dgelist.norm$samples

################################################
# Multidimensional scaling plot of all samples #
################################################

# Output value for MDS plot (dimension 1 and 2 in this case)
MDS <- plotMDS(x = dgelist.norm, top = 1000000, gene.selection = "pairwise", 
               xlab = "Dimension 1", ylab = "Dimension 2", dim.plot = c(1, 2), 
               cex = 1)

# MDS values are then plotted with ggplot2
MDS.ggplot <- data.frame(target, MDS$x, MDS$y)
tiff(filename = "miRdeep2_MDS_miRNA.tif", width = (
  (range(MDS.ggplot$MDS.x)[2]-range(MDS.ggplot$MDS.x)[1])*11.5), height = (
    (range(MDS.ggplot$MDS.y)[2]-range(MDS.ggplot$MDS.y)[1])*10), units = "in",
  res = 600, compression = "lzw")
ggplot(data = MDS.ggplot, aes(x = MDS.ggplot$MDS.x, y = MDS.ggplot$MDS.y, 
                              colour = MDS.ggplot$time_point)) + 
  geom_point(size = 4.5) + theme(panel.background = element_rect(
    fill = 'wheat'), legend.title = element_text(size = 12, face = "bold"), 
    legend.text = element_text(size = 10, face = "bold"), 
    axis.title = element_text(face = "bold", size = 12), 
    axis.text = element_text(face = "bold", size = 10), 
    plot.title = element_text(face = "bold", size = 15)) + 
  ggtitle("MDS plot") + xlab("Dimension 1") + ylab("Dimension 2") + 
  scale_colour_discrete(name = "Time point\npost-infection", breaks = c(
    "pre2", "pre1", "1", "2", "6", "10", "12"), 
    labels = c("pre 2w", "pre 1w", "1w", "2w", "6w", "10w", "12w"))
dev.off()

# Write into a table the coordinates of each library for the MDS plot
write.table(x = MDS.ggplot, file = "miRdeep2_MDS_xy.txt", sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = TRUE)

#########################################################################
# Determine if there has been inversion between 6522_pre2 and 6526_pre2 #
#########################################################################

# Output value for MDS plot (dimension 1 and 2) for animals 6522 and 6526
MDS <- plotMDS(x = dgelist.norm[, grep(pattern = "652(2|6)",
                                       x = colnames(dgelist.norm),
                                       perl = TRUE)], top = 1000000,
               gene.selection = "pairwise", xlab = "Dimension 1",
               ylab = "Dimension 2", dim.plot = c(1, 2), cex = 1)

# MDS values are then plotted with ggplot2
MDS.ggplot <- data.frame(target[grep(pattern = "652(2|6)",
                                     x = target[, "sample"], perl = TRUE),],
                         MDS$x, MDS$y)
ggplot(data = MDS.ggplot, aes(x = MDS.ggplot$MDS.x, y = MDS.ggplot$MDS.y, 
                              colour = MDS.ggplot$time_point,
                              shape = MDS.ggplot$animal)) + 
  geom_point(size = 8) + theme(panel.background = element_rect(
    fill = 'wheat'), legend.title = element_text(size = 20, face = "bold"), 
    legend.text = element_text(size = 15, face = "bold"), 
    axis.title = element_text(face = "bold", size = 30), 
    axis.text = element_text(face = "bold", size = 20), 
    plot.title = element_text(face = "bold", size = 40)) + 
  ggtitle("MDS plot") + xlab("Dimension 1") + ylab("Dimension 2")

##############################################
# Create a design matrix for paired analysis #
##############################################

# Create a design matrix
design <- model.matrix(~ animal + group)
rownames(design) <- rownames(dgelist.norm$samples)
colnames(design) <- gsub(pattern = "(animal)|(group)", replacement = "",
                         x = colnames(design), perl = TRUE)
design

########################################################################
# Estimate the dispersion parameter for each tag using Cox-Reid method #
########################################################################

# Calculate the dispersion (common, trended and tagwise)
dgelist.disp <- estimateGLMCommonDisp(y = dgelist.norm, design = design,
                                      verbose = TRUE)
dgelist.disp <- estimateGLMTrendedDisp(y = dgelist.disp, design = design)
dgelist.disp <- estimateGLMTagwiseDisp(y = dgelist.disp, design = design)
names(dgelist.disp)

# Plot the dispersion
png(filename = "miRdeep2_BCV.png", width = 1366, height = 768,
    units = "px")
plotBCV(dgelist.disp)
dev.off()

# Show the calculated dispersion and the coefficient of biological variation
dgelist.disp$common.dispersion
sqrt(dgelist.disp$common.dispersion)

##################################################################
# Determine differential expression using negative binomial GLMs #
##################################################################

# Fit a negative binomial generalized linear model for each tag using
# the design matrix and calculated dispersion
dgeglm.fit <- glmFit(y = dgelist.disp, design = design)
names(dgeglm.fit)

################################
# Differential expression call #
################################

# Test for differential expression for -1 week versus -2 week
miR.DE("pre1", 1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "miRdeep2.dgelrt.pre1w", dedata = "miRdeep2.de.pre1w",
       smearfile = "miRdeep2_smear_pre1w")

# Test for differential expression for 1 week versus -2 week and -1 week
miR.DE(1, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "miRdeep2.dgelrt.1w", dedata = "miRdeep2.de.1w",
       smearfile = "miRdeep2_smear_1w")

# Test for differential expression for 2 week versus -2 week and -1 week
miR.DE(2, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "miRdeep2.dgelrt.2w", dedata = "miRdeep2.de.2w",
       smearfile = "miRdeep2_smear_2w")

# Test for differential expression for 6 week versus -2 week and -1 week
miR.DE(6, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "miRdeep2.dgelrt.6w", dedata = "miRdeep2.de.6w",
       smearfile = "miRdeep2_smear_6w")

# Test for differential expression for 10 week versus -2 week and -1 week
miR.DE(10, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "miRdeep2.dgelrt.10w", dedata = "miRdeep2.de.10w",
       smearfile = "miRdeep2_smear_10w")

# Test for differential expression for 12 week versus -2 week and -1 week
miR.DE(12, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "miRdeep2.dgelrt.12w", dedata = "miRdeep2.de.12w",
       smearfile = "miRdeep2_smear_12w")

#########################################################
# Merge all DE call data from the different time points #
#########################################################

# Merge all DE table for the different time points into a single dataframe
DEtable.merge("miRdeep2.de.1w", "miRdeep2.de.2w", "miRdeep2.de.6w",
              "miRdeep2.de.10w", "miRdeep2.de.12w",
              output = "miRdeep2.DE", pattern = "^miRdeep2.de.")
head(miRdeep2.DE)

# Write into a table the full DE call data
write.matrix(x = miRdeep2.DE, file = "miRdeep2_full_DE.txt", sep = "\t")

##############################################
# Comparison of DE genes between time points #
##############################################

# Identify as a vector list the significant DE genes per time point
sig.1w <- as.character(miRdeep2.DE[!is.na(miRdeep2.DE$FDR_1w) & (
  miRdeep2.DE$FDR_1w < 0.05), "gene_id"])
sig.2w <- as.character(miRdeep2.DE[!is.na(miRdeep2.DE$FDR_2w) & (
  miRdeep2.DE$FDR_2w < 0.05), "gene_id"])
sig.6w <- as.character(miRdeep2.DE[!is.na(miRdeep2.DE$FDR_6w) & (
  miRdeep2.DE$FDR_6w < 0.05), "gene_id"])
sig.10w <- as.character(miRdeep2.DE[!is.na(miRdeep2.DE$FDR_10w) & (
  miRdeep2.DE$FDR_10w < 0.05), "gene_id"])
sig.12w <- as.character(miRdeep2.DE[!is.na(miRdeep2.DE$FDR_12w) & (
  miRdeep2.DE$FDR_12w < 0.05), "gene_id"])

# Create the Venn diagram of overlapping DE genes between time points
venn.diagram(x = list("1w" = sig.1w, "2w" = sig.2w, "6w" = sig.6w,
                      "10w" = sig.10w, "12w" = sig.12w),
             filename = "miRdeep2_Venn_All.tiff", na="remove")

# Identify the DE miRNA common to all post-infection time points
miRdeep2.overlap <- Reduce(intersect, list(sig.1w, sig.2w, sig.6w, sig.10w, sig.12w))
miRdeep2.overlap

############################
# Clean up the R workspace #
############################

# Rename variables to be kept according to the dataset used for generation
miRdeep2.dgelist.disp <- dgelist.disp

# Remove temporary and unrequired variables
rm(target, animal, group, Count, i, count.log2, hist_plot, MDS, MDS.ggplot,
   dgelist, dgelist.filt, design, dgelist.norm, dgeglm.fit,
   miRdeep2.dgelrt.pre1w, miRdeep2.de.pre1w, sig.1w, sig.2w, sig.6w,
   sig.10w, sig.12w, dgelist.disp, files, Dat, sample)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
##########################################################
# Analysis of counts data obtained via miRdeep* pipeline #
##########################################################

################################################
# Read in and concatenate input files within R #
################################################

# Read in the count files via the read.table function (cannot use readDGE)
files <- list.files(path = paste("C:/Users/nnalpas/Dropbox/Home_work_sync/",
                                 "Work/TIDA/miRNA-seq/Results/Counts/",
                                 "miRdeep-star", sep = ""), pattern = "6")
for (i in 1:length(files)) {
  Dat <- read.table(file = paste("C:/Users/nnalpas/Dropbox/Home_work_sync/",
                                 "Work/TIDA/miRNA-seq/Results/Counts/",
                                 "miRdeep-star", files[i], sep = "/"),
                    quote = "", sep = "\t", header = TRUE)
  Dat <- Dat[, c(1,3,4,6,7)]
  sample <- gsub(pattern = "_result.txt", replacement = "", x = files[i])
  colnames(Dat)[1:4] <- c("precursor_name", "chromosome", "strand", sample)
  Dat <- cbind(Dat[,1:4], colsplit(string = as.character(
    Dat[,5]), pattern = "-", names = c("start_position", "end_position")))
  Dat <- merge(x = Dat, y = miRNA.info, by = c(
    "precursor_name", "chromosome", "strand", "start_position"), all = TRUE)
  Dat <- Dat[, c("gene_id", sample)]
  Dat[is.na(x = Dat)] <- 0
  if (i == 1) {
    Count <- Dat
  }
  else {
    Count <- merge(x = Count, y = Dat, by = "gene_id", all = TRUE)
  }
}
dim(Count)
head(Count)
tail(Count)

######################################################
# Merge the annotation information with counts table #
######################################################

# Merge previously imported annotation information with the count table
Count <- merge(x = miRNA.info, y = Count, by = "gene_id", all = TRUE)
head(Count)

# Check the merged matrix in terms of size and content
dim(Count)
table(duplicated(Count[, "gene_id"]))
length(unique(Count[, "gene_id"], incomparables = FALSE))

# Ouptut samples data
write.matrix(x = Count, file = "miRdeepstar_counts.txt",  sep = "\t")

###############################
# Create groups and a DGElist #
###############################

# Create a target matrix and a experimental group vector and animal block
target <- colnames(Count[, (ncol(Count)-69):ncol(Count)])
target <- cbind(target, matrix(data = unlist(strsplit(x = target, split = "_",
                                                      fixed = TRUE)),
                               nrow = 70, ncol = 2, byrow = TRUE))
colnames(target) = c("sample", "animal", "time_point")
group <- factor(x = target[, "time_point"])
animal <- factor(x = target[, "animal"])

# Create a DGElist containing the group information
dgelist <- DGEList(counts = Count[, (ncol(Count)-69):ncol(Count)], 
                   lib.size = NULL, norm.factors = NULL, group = group, 
                   genes = Count[, 1:(ncol(Count)-70)], 
                   remove.zeros = FALSE)
rownames(dgelist$counts) <- rownames(
  dgelist$genes) <- dgelist$genes[, "gene_id"]
dgelist$genes["gene_id"] <- NULL
names(dgelist)
head(dgelist$samples)
head(dgelist$counts)
head(dgelist$genes)

# Ouptut samples data
write.table(x = dgelist$samples, file = "miRdeepstar_samples.txt", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

###########################################################
# Quality check of libraries by plotting density of count #
###########################################################

# Log2 transform the count data for better visualization
count.log2 <- log(x = (dgelist$counts[,]+1), base = 2)

# Plot density of count for all libraries
png(filename = "miRdeepstar_Dens_all_miRNA.png", width = 1500, height = 1300,
    units = "px")
plot(x = density(count.log2[, 1]), main = "Density plot of count per gene", 
     lty =  1, xlab = "Log2 of count per gene", ylab = "Density", col = 1, 
     xlim = c(-0.5, 5), ylim = c(0.0, 0.9), lwd = 2 , cex.axis = 1.5, 
     cex.lab = 1.5, cex.main = 2)
for (i in 2:ncol(count.log2)) {
  lines(x = density(count.log2[, i]), lty = 1, col = i)
}
dev.off()

#####################################
# Filtering of lowly expressed tags #
#####################################

# Identify genes with zero counts across all samples
dim(dgelist[rowSums(dgelist$counts) == 0, ])
head(dgelist[rowSums(dgelist$counts) == 0, ])

# Filter lowly expressed tags, retaining only tags with at least 50 counts per 
# million in 10 or more libraries (10 libraries correspond to one time point)
dgelist.filt <- dgelist[rowSums(cpm(dgelist$counts) > 50) >= 10, ]
dim(dgelist.filt$counts)

# Compute the library size
dgelist.filt$samples$lib.size <- colSums(dgelist.filt$counts)
head(dgelist.filt$samples)
head(dgelist$samples)

###########################################################################
# Quality check of libraries by plotting density of count after filtering #
###########################################################################

# Log2 transform the count data for better visualization
count.log2 <- log(x = (dgelist.filt$counts[,]+1), base = 2)

# Plot density of count for all libraries
png(filename = "miRdeepstar_Dens_filt_miRNA.png", width = 1366, height = 768,
    units = "px")
plot(x = density(count.log2[, 1]), main = "Density plot of count per gene", 
     lty =  1, xlab = "Log2 of count per gene", ylab = "Density", col = 1, 
     xlim = c(-0.5, 15), ylim = c(0.0, 0.2), lwd = 2 , cex.axis = 1.5, 
     cex.lab = 1.5, cex.main = 2)
for (i in 2:ncol(count.log2)) {
  lines(x = density(count.log2[, i]), lty = 1, col = i)
}
dev.off()

############################################
# Histogram of CPM per expression category #
############################################

# Prepare the histogram picture output
png(filename = "miRdeepstar_Hist_miRNA.png", width = 1500, height = 1500,
    units = "px")
par(mfrow = c(2,2), cex = 3, cex.axis = 0.7)

# Plot histogram of number of gene per CPM category for all microRNA
hist_plot <- hist(x = rowMeans(x = cpm(dgelist$counts), na.rm = FALSE), 
                  breaks = c(
                    0, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000), 
                  plot = FALSE)
plot(x = hist_plot$counts, log = "x", type = "h", axes = FALSE, 
     xlab = "Average CPM", ylab = "Number of gene", lwd = 10, 
     col = "lightblue", main = "All miRNA")
axis(side = 1, at = c(1:10), labels = c(
  "0", "0.01", "0.1", "1", "10", "100", "1000", "10000", "100000", "1000000"))
axis(side = 2)

# Plot histogram of number of gene per CPM category for microRNA with 
# at least one count
hist_plot <- hist(x = rowMeans(x = cpm(dgelist$counts[rowSums(
  dgelist$counts) > 0, ]), na.rm = FALSE), breaks = c(
    0, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000), plot = FALSE)
plot(hist_plot$counts, log = "x", type = "h", axes = FALSE, 
     xlab = "Average CPM", ylab = "Number of gene", lwd = 10, 
     col = "lightgreen", main = "Min one count miRNA")
axis(side = 1, at = c(1:10), labels = c(
  "0", "0.01", "0.1", "1", "10", "100", "1000", "10000", "100000", "1000000"))
axis(side = 2)

# Plot histogram of number of gene per CPM category for microRNA 
# filtered for low expression
hist_plot <- hist(x = rowMeans(x = cpm(dgelist.filt$counts), na.rm = FALSE), 
                  breaks = c(
                    0, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000), 
                  plot = FALSE)
plot(hist_plot$counts, log = "x", type = "h", axes = FALSE, 
     xlab = "Average CPM", ylab = "Number of gene", lwd = 10, 
     col = "lightcoral", main = "Low expression filtered miRNA")
axis(side = 1, at = c(1:10), labels = c(
  "0", "0.01", "0.1", "1", "10", "100", "1000", "10000", "100000", "1000000"))
axis(side = 2)

# Close graphic device
dev.off()

#################################################
# Gene expression correlation between libraries #
#################################################

# Note that CPM values are not required with the use of Spearman rank test

# Perform gene expression correlation between libraries at pre 2 weeks to 
# identify outlier library
png(filename = "miRdeepstar_Pre2_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_pre2", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Pre 2 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at pre 1 weeks to 
# identify outlier library
png(filename = "miRdeepstar_Pre1_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_pre1", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Pre 1 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 1 weeks to 
# identify outlier library
png(filename = "miRdeepstar_Post1_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_1$", x = colnames(dgelist.filt$counts), perl = TRUE)] + 1),
  base = 2), histogram = TRUE, method = "spearman", 
  main = "Post 1 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 2 weeks to 
# identify outlier library
png(filename = "miRdeepstar_Post2_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_2", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Post 2 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 6 weeks to 
# identify outlier library
png(filename = "miRdeepstar_Post6_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_6", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Post 6 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 10 weeks to 
# identify outlier library
png(filename = "miRdeepstar_Post10_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_10", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Post 10 weeks time point expression correlation")
dev.off()

# Perform gene expression correlation between libraries at post 12 weeks to 
# identify outlier library
png(filename = "miRdeepstar_Post12_cor.png", width = 1500, height = 1500,
    units = "px")
chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
  pattern = "_12", x = colnames(dgelist.filt$counts))] + 1), base = 2), 
  histogram = TRUE, method = "spearman", 
  main = "Post 12 weeks time point expression correlation")
dev.off()

########################################################
# Normalization of data using trimmed mean of M-values #
########################################################

# Calculate normalisation factor for our DGElist, note that with edgeR 
# the counts are not transformed in any way after normalization, instead 
# normalization will modify library size
dgelist.norm <- calcNormFactors(dgelist.filt)
dgelist.norm$samples

################################################
# Multidimensional scaling plot of all samples #
################################################

# Output value for MDS plot (dimension 1 and 2 in this case)
MDS <- plotMDS(x = dgelist.norm, top = 1000000, gene.selection = "pairwise", 
               xlab = "Dimension 1", ylab = "Dimension 2", dim.plot = c(1, 2), 
               cex = 1)

# MDS values are then plotted with ggplot2
MDS.ggplot <- data.frame(target, MDS$x, MDS$y)
tiff(filename = "miRdeepstar_MDS_miRNA.tif", width = (
  (range(MDS.ggplot$MDS.x)[2]-range(MDS.ggplot$MDS.x)[1])*11.5), height = (
    (range(MDS.ggplot$MDS.y)[2]-range(MDS.ggplot$MDS.y)[1])*10), units = "in",
  res = 600, compression = "lzw")
ggplot(data = MDS.ggplot, aes(x = MDS.ggplot$MDS.x, y = MDS.ggplot$MDS.y, 
                              colour = MDS.ggplot$time_point)) + 
  geom_point(size = 4.5) + theme(panel.background = element_rect(
    fill = 'wheat'), legend.title = element_text(size = 12, face = "bold"), 
    legend.text = element_text(size = 10, face = "bold"), 
    axis.title = element_text(face = "bold", size = 12), 
    axis.text = element_text(face = "bold", size = 10), 
    plot.title = element_text(face = "bold", size = 15)) + 
  ggtitle("MDS plot") + xlab("Dimension 1") + ylab("Dimension 2") + 
  scale_colour_discrete(name = "Time point\npost-infection", breaks = c(
    "pre2", "pre1", "1", "2", "6", "10", "12"), 
    labels = c("pre 2w", "pre 1w", "1w", "2w", "6w", "10w", "12w"))
dev.off()

# Write into a table the coordinates of each library for the MDS plot
write.table(x = MDS.ggplot, file = "miRdeepstar_MDS_xy.txt", sep = "\t",
            quote = FALSE, row.names = TRUE, col.names = TRUE)

#########################################################################
# Determine if there has been inversion between 6522_pre2 and 6526_pre2 #
#########################################################################

# Output value for MDS plot (dimension 1 and 2) for animals 6522 and 6526
MDS <- plotMDS(x = dgelist.norm[, grep(pattern = "652(2|6)",
                                       x = colnames(dgelist.norm),
                                       perl = TRUE)], top = 1000000,
               gene.selection = "pairwise", xlab = "Dimension 1",
               ylab = "Dimension 2", dim.plot = c(1, 2), cex = 1)

# MDS values are then plotted with ggplot2
MDS.ggplot <- data.frame(target[grep(pattern = "652(2|6)",
                                     x = target[, "sample"], perl = TRUE),],
                         MDS$x, MDS$y)
ggplot(data = MDS.ggplot, aes(x = MDS.ggplot$MDS.x, y = MDS.ggplot$MDS.y, 
                              colour = MDS.ggplot$time_point,
                              shape = MDS.ggplot$animal)) + 
  geom_point(size = 8) + theme(panel.background = element_rect(
    fill = 'wheat'), legend.title = element_text(size = 20, face = "bold"), 
    legend.text = element_text(size = 15, face = "bold"), 
    axis.title = element_text(face = "bold", size = 30), 
    axis.text = element_text(face = "bold", size = 20), 
    plot.title = element_text(face = "bold", size = 40)) + 
  ggtitle("MDS plot") + xlab("Dimension 1") + ylab("Dimension 2")

##############################################
# Create a design matrix for paired analysis #
##############################################

# Create a design matrix
design <- model.matrix(~ animal + group)
rownames(design) <- rownames(dgelist.norm$samples)
colnames(design) <- gsub(pattern = "(animal)|(group)", replacement = "",
                         x = colnames(design), perl = TRUE)
design

########################################################################
# Estimate the dispersion parameter for each tag using Cox-Reid method #
########################################################################

# Calculate the dispersion (common, trended and tagwise)
dgelist.disp <- estimateGLMCommonDisp(y = dgelist.norm, design = design,
                                      verbose = TRUE)
dgelist.disp <- estimateGLMTrendedDisp(y = dgelist.disp, design = design)
dgelist.disp <- estimateGLMTagwiseDisp(y = dgelist.disp, design = design)
names(dgelist.disp)

# Plot the dispersion
png(filename = "miRdeepstar_BCV.png", width = 1366, height = 768,
    units = "px")
plotBCV(dgelist.disp)
dev.off()

# Show the calculated dispersion and the coefficient of biological variation
dgelist.disp$common.dispersion
sqrt(dgelist.disp$common.dispersion)

##################################################################
# Determine differential expression using negative binomial GLMs #
##################################################################

# Fit a negative binomial generalized linear model for each tag using
# the design matrix and calculated dispersion
dgeglm.fit <- glmFit(y = dgelist.disp, design = design)
names(dgeglm.fit)

################################
# Differential expression call #
################################

# Test for differential expression for -1 week versus -2 week
miR.DE("pre1", 1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "miRdeepstar.dgelrt.pre1w", dedata = "miRdeepstar.de.pre1w",
       smearfile = "miRdeepstar_smear_pre1w")

# Test for differential expression for 1 week versus -2 week and -1 week
miR.DE(1, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "miRdeepstar.dgelrt.1w", dedata = "miRdeepstar.de.1w",
       smearfile = "miRdeepstar_smear_1w")

# Test for differential expression for 2 week versus -2 week and -1 week
miR.DE(2, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "miRdeepstar.dgelrt.2w", dedata = "miRdeepstar.de.2w",
       smearfile = "miRdeepstar_smear_2w")

# Test for differential expression for 6 week versus -2 week and -1 week
miR.DE(6, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "miRdeepstar.dgelrt.6w", dedata = "miRdeepstar.de.6w",
       smearfile = "miRdeepstar_smear_6w")

# Test for differential expression for 10 week versus -2 week and -1 week
miR.DE(10, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "miRdeepstar.dgelrt.10w", dedata = "miRdeepstar.de.10w",
       smearfile = "miRdeepstar_smear_10w")

# Test for differential expression for 12 week versus -2 week and -1 week
miR.DE(12, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
       group = group, adjpvalue = 0.05, method = "BH",
       lrtdata = "miRdeepstar.dgelrt.12w", dedata = "miRdeepstar.de.12w",
       smearfile = "miRdeepstar_smear_12w")

#########################################################
# Merge all DE call data from the different time points #
#########################################################

# Merge all DE table for the different time points into a single dataframe
DEtable.merge("miRdeepstar.de.1w", "miRdeepstar.de.2w", "miRdeepstar.de.6w",
              "miRdeepstar.de.10w", "miRdeepstar.de.12w",
              output = "miRdeepstar.DE", pattern = "^miRdeepstar.de.")
head(miRdeepstar.DE)

# Write into a table the full DE call data
write.matrix(x = miRdeepstar.DE, file = "miRdeepstar_full_DE.txt", sep = "\t")

##############################################
# Comparison of DE genes between time points #
##############################################

# Identify as a vector list the significant DE genes per time point
sig.1w <- as.character(miRdeepstar.DE[!is.na(miRdeepstar.DE$FDR_1w) & (
  miRdeepstar.DE$FDR_1w < 0.05), "gene_id"])
sig.2w <- as.character(miRdeepstar.DE[!is.na(miRdeepstar.DE$FDR_2w) & (
  miRdeepstar.DE$FDR_2w < 0.05), "gene_id"])
sig.6w <- as.character(miRdeepstar.DE[!is.na(miRdeepstar.DE$FDR_6w) & (
  miRdeepstar.DE$FDR_6w < 0.05), "gene_id"])
sig.10w <- as.character(miRdeepstar.DE[!is.na(miRdeepstar.DE$FDR_10w) & (
  miRdeepstar.DE$FDR_10w < 0.05), "gene_id"])
sig.12w <- as.character(miRdeepstar.DE[!is.na(miRdeepstar.DE$FDR_12w) & (
  miRdeepstar.DE$FDR_12w < 0.05), "gene_id"])

# Create the Venn diagram of overlapping DE genes between time points
venn.diagram(x = list("1w" = sig.1w, "2w" = sig.2w, "6w" = sig.6w,
                      "10w" = sig.10w, "12w" = sig.12w),
             filename = "miRdeepstar_Venn_All.tiff", na="remove")

# Identify the DE miRNA common to all post-infection time points
miRdeepstar.overlap <- Reduce(intersect, list(sig.1w, sig.2w, sig.6w, sig.10w, sig.12w))
miRdeepstar.overlap

############################
# Clean up the R workspace #
############################

# Rename variables to be kept according to the dataset used for generation
miRdeepstar.dgelist.disp <- dgelist.disp

# Remove temporary and unrequired variables
rm(target, animal, group, Count, i, count.log2, hist_plot, MDS, MDS.ggplot,
   dgelist, dgelist.filt, design, dgelist.norm, dgeglm.fit,
   miRdeepstar.dgelrt.pre1w, miRdeepstar.de.pre1w, sig.1w, sig.2w, sig.6w,
   sig.10w, sig.12w, miR.DE, DEtable.merge, DESeq.merge, diff_expr_edgeR,
   loadpackage, sig_label, dgelist.disp, files, Dat, sample)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#############################################################################
# Comparison of time overlapping DE genes between various analysis pipeline #
#############################################################################

# Create the Venn diagram of overlapping DE genes between time points
venn.diagram(x = list("Novo-feature" = Novofeature.overlap,
                      "miRdeep2" = miRdeep2.overlap,
                      "miRdeep-star" = miRdeepstar.overlap),
             filename = "Venn_3pipelines.tiff", na="remove")

# Get all the gene expression information for the candidate biomarkers miRNA
candidates <- merge(x = miRNA.info[(miRNA.info$gene_id %in% unique(c(
  Novofeature.overlap, miRdeep2.overlap, miRdeepstar.overlap))),],
  y = Novofeature.DE[, (ncol(Novofeature.DE)-24):ncol(Novofeature.DE)],
  by.x = "gene_id", by.y = "row.names", all.x = TRUE)
candidates <- merge(x = candidates, y = miRdeep2.DE[, (ncol(miRdeep2.DE)-24):ncol(
  miRdeep2.DE)], by.x = "gene_id", by.y = "row.names", all.x = TRUE)
candidates <- merge(x = candidates, y = miRdeepstar.DE[, (ncol(
  miRdeepstar.DE)-24):ncol(miRdeepstar.DE)], by.x = "gene_id",
  by.y = "row.names", all.x = TRUE)

# Write into a table all the time overlapping DE genes
write.matrix(x = candidates, file = "Candidate_biomarkers.txt", sep = "\t")

# Remove temporary and unrequired variables
rm(Novofeature.overlap, miRdeep2.overlap, miRdeepstar.overlap)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##########################################################
# Identification of reference genes for RT-qPCR analysis #
##########################################################

# Calculate average CPM and Standard deviation for each gene
mean <- as.matrix(apply(X = cpm(x = Novofeature.dgelist.disp), MARGIN = 1,
                        FUN = function(x) mean(x)))
std.dev <- as.matrix(apply(X = cpm(x = Novofeature.dgelist.disp), MARGIN = 1,
                           FUN = function(x) sd(x)))

# Merge results together
gene.stable <- merge(x = mean, y = std.dev, by = "row.names", all = TRUE)
gene.stable <- merge(x = miRNA.info[, 1:2], y = gene.stable, by.x = "gene_id",
                     by.y = "Row.names")
rownames(gene.stable) <- gene.stable[, "gene_id"]
gene.stable <- gene.stable[, -1]
colnames(gene.stable)[2:3] <- c("meanCPM", "SdevCPM")
head(gene.stable)

# Calculate percentage of Standard deviation in order to compare all genes
gene.stable <- cbind(gene.stable, (
  gene.stable[, "SdevCPM"]*100/gene.stable[, "meanCPM"]))
colnames(gene.stable)[4] <- "Perc_Sdev"
gene.stable <- gene.stable[order(gene.stable$Perc_Sdev),]

# Plot the mean CPM and Standard deviation CPM of all genes
ggplot(data = gene.stable, aes(x = log(x = gene.stable$meanCPM, base = 2),
                               y = log(x = gene.stable$SdevCPM,
                                       base = 2))) + geom_point(size = 4.5)

# Write into a table the analysis for reference gene discovery
write.matrix(x = gene.stable, file = "Potential_reference_gene.txt", sep = "\t")

# Remove temporary and unrequired variables
rm(mean, std.dev, gene.stable)

#######
# END #
#######