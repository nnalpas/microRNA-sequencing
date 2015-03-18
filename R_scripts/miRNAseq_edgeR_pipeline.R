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

# Extract User path to Dropbox
user <- gsub(pattern = "(^.*/Dropbox).*$", replacement = "\\1", x = getwd(),
             perl = TRUE)

# Source the common functions used across this script
source(file = paste(user, "/Home_work_sync/Work/Bioinformatics/R",
                    "/General_function.R", sep = ""))

# Load the required packages
loadpackage(package = edgeR)
loadpackage(package = MASS)
loadpackage(package = ggplot2)
loadpackage(package = PerformanceAnalytics)
loadpackage(package = VennDiagram)
loadpackage(package = reshape2)
loadpackage(package = randomForest)
loadpackage(package = gplots)
loadpackage(package = RColorBrewer)
loadpackage(package = magrittr)
loadpackage(package = gridExtra)

#######################################################################
# Analysis of counts data obtained via different pipelines (pick one) #
#######################################################################

# Define the method used for count generation ("Novoalign", "miRdeep2"
# or "miRdeepstar")
#method <- "Novoalign"
#method <- "miRdeep2"
method <- "miRdeepstar"

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
miRNA.info <- miRNA.info[, c(2:7, 1, 8:ncol(miRNA.info))]
colnames(miRNA.info)[c(1,7,12)] <- c("gene_id", "sequence",
                                     "identical_sequence")
head(miRNA.info)
dim(miRNA.info)

#####################################################################
# Read in and concatenate input files within R for method Novoalign #
#####################################################################

# Reads and merges a set of files containing counts according to the method
{
  if (method == "Novoalign") {
    Count <- readDGE(files = list.files(
      path = paste(user, "/Home_work_sync/Work/TIDA/miRNA-seq",
                   "/Results/Counts/Novo-feature/mature_miRNA",
                   sep = ""), pattern = "6"),
      path = paste(user, "/Home_work_sync/Work/TIDA/miRNA-seq",
                   "/Results/Counts/Novo-feature/mature_miRNA",
                   sep = ""), columns = c(1,7), 
      skip = 1, header = TRUE)
    Count <- data.frame(gene_id = rownames(Count), Count$counts)
    colnames(Count) <- gsub(pattern = "^X", replacement = "", x = colnames(Count),
                            perl = TRUE)
  }
  else if (method == "miRdeep2") {
    files <- list.files(path = paste(user, "/Home_work_sync/Work/TIDA/miRNA-seq",
                                     "/Results/Counts/mirdeep2", sep = ""),
                        pattern = "6")
    for (i in 1:length(files)) {
      Dat <- read.table(file = paste(user, "Home_work_sync/Work/TIDA/miRNA-seq",
                                     "/Results/Counts/mirdeep2", files[i],
                                     sep = "/"), quote = "")
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
  }
  else if (method == "miRdeepstar") {
    files <- list.files(path = paste(user, "/Home_work_sync/Work/TIDA/miRNA-seq/",
                                     "Results/Counts/miRdeep-star", sep = ""),
                        pattern = "6")
    for (i in 1:length(files)) {
      Dat <- read.table(file = paste(user, "Home_work_sync/Work/TIDA/miRNA-seq/",
                                     "Results/Counts/miRdeep-star", files[i],
                                     sep = "/"), quote = "", sep = "\t",
                        header = TRUE)
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
  }
  else {
    print("Error: The provided method for count generation is unknown!")
  }
}

# Merge the annotation information with the count table
Count <- merge(x = miRNA.info, y = Count, by = "gene_id", all = TRUE)
head(Count)

# Check the merged matrix in terms of size and content
table(duplicated(Count[, "gene_id"]))
length(unique(Count[, "gene_id"], incomparables = FALSE))

# Check the raw data
head(Count)
tail(Count)
dim(Count)

# Ouptut samples data
write.matrix(x = Count, file = paste(method, "_counts.txt", sep = ""),
             sep = "\t")

###############################
# Create groups and a DGElist #
###############################

# Create a target matrix and a experimental group vector and animal block
target <- colnames(Count[, (ncol(Count)-69):ncol(Count)])
target <- cbind(target, matrix(data = unlist(strsplit(
  x = target, split = "_", fixed = TRUE)), nrow = 70, ncol = 2, byrow = TRUE))
colnames(target) = c("sample", "animal", "time_point")
group <- factor(x = target[, "time_point"])
animal <- factor(x = target[, "animal"])

# Create a DGElist containing the group information
dgelist <- DGEList(counts = Count[, (ncol(Count)-69):ncol(Count)],
                   lib.size = NULL, norm.factors = NULL, group = group,
                   genes = Count[, 1:(ncol(Count)-70)], remove.zeros = FALSE)
rownames(dgelist$counts) <- rownames(
  dgelist$genes) <- dgelist$genes[, "gene_id"]
dgelist$genes["gene_id"] <- NULL
names(dgelist)
head(dgelist$samples)
head(dgelist$counts)
head(dgelist$genes)

# Create a variable containing the current DGElist
assign(x = paste(method, ".dgelist", sep = ""), value = dgelist)

# Ouptut samples data
write.table(x = dgelist$samples, file = paste(method, "_sample.txt", sep = ""),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

###########################################################
# Quality check of libraries by plotting density of count #
###########################################################

# Log2 transform the count data for better visualization
count.log2 <- log(x = (dgelist$counts[,]+1), base = 2)

# Plot density of count for all libraries
png(filename = paste(method, "_Dens_all_miRNA.png", sep = ""), width = 1500,
    height = 1300, units = "px")
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
dgelist.filt <- dgelist[rowSums(
  cpm(dgelist$counts) > 50) >= median(summary(group)), ]
dim(dgelist.filt$counts)

# Compute the library size
dgelist.filt$samples$lib.size <- colSums(dgelist.filt$counts)
head(dgelist.filt$samples)
head(dgelist$samples)

###########################################################################
# Quality check of libraries by plotting density of count after filtering #
###########################################################################

# Log2 transform the count data for better visualization
count.log2 <- log(x = (dgelist.filt$counts[,] + 1), base = 2)

# Plot density of count for all libraries
png(filename = paste(method, "_Dens_filt_miRNA.png", sep = ""),
    width = 1366, height = 768, units = "px")
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
png(filename = paste(method, "_Hist_miRNA.png", sep = ""),
    width = 1500, height = 1500, units = "px")
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

# Perform gene expression correlation between libraries at each time points to 
# identify potential outlier library (CPM values not required with Spearman)
for (val in levels(group)) {
  png(filename = paste(method, val, "cor.png", sep = "_"),
      width = 1500, height = 1500, units = "px")
  chart.Correlation(R = log(x = (dgelist.filt$counts[, grep(
    pattern = paste("_", val, "$", sep = ""), x = colnames(dgelist.filt$counts),
    perl = TRUE)] + 1), base = 2), histogram = TRUE, method = "spearman",
    main = paste(val, " time point expression correlation", sep = ""))
  dev.off()
}

########################################################
# Normalization of data using trimmed mean of M-values #
########################################################

# Calculate normalisation factor for our DGElist, note that with edgeR 
# the counts are not transformed in any way after normalization
dgelist.norm <- calcNormFactors(dgelist.filt)
dgelist.norm$samples

################################################
# Multidimensional scaling plot of all samples #
################################################

# Define a set of colors
colours <- brewer.pal(n = length(levels(factor(
  target[, "time_point"]))), name = "Set1")

# Use custom function to output a MDS plot of all samples
multi.MDS(data = dgelist.norm, target = target, prefix = method,
          suffix = "all", plotmds = list(
            top = 1000000, gene.selection = "pairwise", dim.plot = c(1, 2)),
          aes.colour = "time_point", fill = "lightgrey", size = 3,
          legend.pos = "right", manual.colour = colours,
          breaks = c("pre2", "pre1", "1", "2", "6", "10", "12"))

#################################################
# Multidimensional scaling plot per time points #
#################################################

# Use custom function to output a MDS plot per post-infection time points
multi.MDS(pattern = c("_pre2|_pre1|_1$", "_pre2|_pre1|_2$", "_pre2|_pre1|_6$", "_pre2|_pre1|_10$",
                      "_pre2|_pre1|_12$"),
          data = dgelist.norm, target = target, prefix = method,
          suffix = c("W1", "W2", "W6", "W10", "W12"),
          plotmds = list(top = 1000000, gene.selection = "pairwise",
                         dim.plot = c(1, 2)),
          aes.colour = "time_point", fill = "wheat", size = 3,
          legend.pos = "right", combine = TRUE,
          manual.colour = colours[c(2, 3, rep(x = 1, times = 5))],
          breaks = c("pre2", "pre1", "1", "2", "6", "10", "12"))

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
png(filename = paste(method, "_BCV.png", sep = ""),
    width = 1366, height = 768, units = "px")
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
multi.DE("pre1", 1, "pre2", -1, data = dgeglm.fit, design = design,
         group = group, adjpvalue = 0.05, method = "BH",
         lrtdata = paste(method, ".dgelrt.pre1w", sep = ""),
         dedata = paste(method, ".de.pre1w", sep = ""),
         smearfile = paste(method, "_smear_pre1w", sep = ""))

# Test for differential expression for 1 week versus -2 week and -1 week
multi.DE(1, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
         group = group, adjpvalue = 0.05, method = "BH",
         lrtdata = paste(method, ".dgelrt.1w", sep = ""),
         dedata = paste(method, ".de.1w", sep = ""),
         smearfile = paste(method, "_smear_1w", sep = ""))

# Test for differential expression for 2 week versus -2 week and -1 week
multi.DE(2, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
         group = group, adjpvalue = 0.05, method = "BH",
         lrtdata = paste(method, ".dgelrt.2w", sep = ""),
         dedata = paste(method, ".de.2w", sep = ""),
         smearfile = paste(method, "_smear_2w", sep = ""))

# Test for differential expression for 6 week versus -2 week and -1 week
multi.DE(6, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
         group = group, adjpvalue = 0.05, method = "BH",
         lrtdata = paste(method, ".dgelrt.6w", sep = ""),
         dedata = paste(method, ".de.6w", sep = ""),
         smearfile = paste(method, "_smear_6w", sep = ""))

# Test for differential expression for 10 week versus -2 week and -1 week
multi.DE(10, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
         group = group, adjpvalue = 0.05, method = "BH",
         lrtdata = paste(method, ".dgelrt.10w", sep = ""),
         dedata = paste(method, ".de.10w", sep = ""),
         smearfile = paste(method, "_smear_10w", sep = ""))

# Test for differential expression for 12 week versus -2 week and -1 week
multi.DE(12, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
         group = group, adjpvalue = 0.05, method = "BH",
         lrtdata = paste(method, ".dgelrt.12w", sep = ""),
         dedata = paste(method, ".de.12w", sep = ""),
         smearfile = paste(method, "_smear_12w", sep = ""))

#########################################################
# Merge all DE call data from the different time points #
#########################################################

# Create a vector of DE table to merge
DEtable.name <- c("1w", "2w", "6w", "10w", "12w") %>% lapply(
  X = ., FUN = function(x) paste(method, "de", x, sep = '.')) %>% unlist(.)

# Merge all DE table for the different time points into a single dataframe
DEtable.merge(DEtable.name, output = paste(method, ".DE", sep = ""),
              pattern = paste("^", method, '.de.', sep = ""))

# Create a variable containing the current DEtable
DEtable <- eval(parse(text = paste(method, '.DE', sep = "")))
head(DEtable)

# Write into a table the full DE call data
write.matrix(x = DEtable, file = paste(method, "_full_DE.txt", sep = ""),
             sep = "\t")

##############################################
# Comparison of DE genes between time points #
##############################################

# Identify as a vector list the significant DE genes per time point
sig.1w <- as.character(DEtable[!is.na(DEtable$FDR_1w) & (
  DEtable$FDR_1w < 0.05), "gene_id"])
sig.2w <- as.character(DEtable[!is.na(DEtable$FDR_2w) & (
  DEtable$FDR_2w < 0.05), "gene_id"])
sig.6w <- as.character(DEtable[!is.na(DEtable$FDR_6w) & (
  DEtable$FDR_6w < 0.05), "gene_id"])
sig.10w <- as.character(DEtable[!is.na(DEtable$FDR_10w) & (
  DEtable$FDR_10w < 0.05), "gene_id"])
sig.12w <- as.character(DEtable[!is.na(DEtable$FDR_12w) & (
  DEtable$FDR_12w < 0.05), "gene_id"])

# Define a set of colors
colours <- brewer.pal(n = 5, name = "Set1")

# Create the Venn diagram of overlapping DE genes between time points
venn.diagram(x = list("1w" = sig.1w, "2w" = sig.2w, "6w" = sig.6w,
                      "10w" = sig.10w, "12w" = sig.12w),
             filename = paste(method, "_Venn_All.tiff", sep = ""),
             na = "remove", res = 600, fill = colours, lwd = 0.7,
             cex = 1, cat.cex = 1, cat.col = colours)

# Identify the DE miRNA common to all post-infection time points
overlap <- assign(x = paste(method, '.overlap', sep = ""),
                  value = Reduce(intersect, list(
                    sig.1w, sig.2w, sig.6w, sig.10w, sig.12w)))
overlap

##################################################
# Clustering of samples based on common DE genes #
##################################################

# Get the CPM values for each filtered miRNA
CPM <- cpm(x = dgelist.disp, normalized.lib.sizes = TRUE)
dim(CPM)
head(CPM)

# Create a phenotypic data table with color coding and infection status
pheno <- colnames(CPM)
pheno <- cbind(sample = pheno,
               status = unlist(as.character(lapply(
                 X = pheno, FUN = function(x) {if(
                   length(grep(pattern = "_pre(1|2)$", x = x,
                               perl = TRUE)) == 1) {"healthy"} else {
                                 "infected"}}))),
               col.side = unlist(as.character(lapply(
                 X = pheno, FUN = function(x) {if(
                   length(grep(pattern = "_pre(1|2)$", x = x,
                               perl = TRUE)) == 1) {colours[2]} else {
                                 colours[1]}}))))
rownames(pheno) <- pheno[, "sample"]
col.side <- pheno[, "col.side"]
pheno <- factor(x = pheno[, "status"])
pheno

# Define breaks and colours for the heatmap
colors <- c(seq(from = 0, to = 1, length = 10), seq(
  from = 1, to = 10, length = 10), seq(from = 10, to = 100, length = 10),
  seq(from = 100, to = 1000, length = 10), seq(
    from = 1000, to = 10000, length = 10), seq(
      from = 10000, to = 100000, length = 10), seq(
        from = 100000, to = 1000000, length = 10))
my.palette <- colorRampPalette(c("blue", "red"))(n = 69)

# Represent the data in a heatmap based on the overlapping miRNA
tiff(filename = paste(method, "_overlap_Heatmap.tif", sep = ""),
     width = 15, height = 15, units = "in", res = 600, compression = "lzw")
heatmap.2(x = CPM[overlap,], labRow = rownames(CPM), labCol = colnames(CPM),
          cexRow = .5, cexCol = .5, breaks = colors, col = my.palette,
          symm = F, symkey = F, symbreaks = T, scale = "none",
          ColSideColors = col.side)
dev.off()

###################################################
# Biomarkers discovery using Random Forest method #
###################################################

# Use randomForest on the CPM values per miRNA
rf <- randomForest(x = t(CPM), y = pheno, importance = TRUE,
                   do.trace = 1, ntree = 5000, mtry = sqrt(nrow(
                     CPM)), keep.forest = TRUE, proximity = TRUE)

# Order the miRNA by their importance value
rf.ImportanceOrdered <- importance(rf)[order(importance(
  rf)[, "MeanDecreaseGini"], decreasing = TRUE),]
head(rf.ImportanceOrdered, n = 20)

# Create a MDS plot
MDSplot(rf = rf, fac = pheno, palette = c(1,3))

# Represent the data in a heatmap based on the top miRNA (by importance)
for (n in 2:15) {
  tiff(filename = paste(method, "_Heatmap_n", n, ".tif", sep = ""),
       width = 15, height = 15, units = "in", res = 600, compression = "lzw")
  pat <- grep(pattern = "_", x = rownames(rf.ImportanceOrdered), value = TRUE,
       invert = TRUE)
  dat <- CPM[head(pat, n = n),]
  heatmap.2(x = dat, labRow = rownames(dat), labCol = colnames(dat),
            cexRow = .5, cexCol = .5, breaks = colors, col = my.palette,
            symm = F, symkey = F, symbreaks = T, scale = "none",
            ColSideColors = col.side)
  dev.off()
}

##########################################################
# Identification of reference genes for RT-qPCR analysis #
##########################################################

# Merge verage CPM and Standard deviation for each gene together
gene.stable <- merge(
  x = as.matrix(apply(X = CPM, MARGIN = 1, FUN = function(x) mean(x))),
  y = as.matrix(apply(X = CPM, MARGIN = 1, FUN = function(x) sd(x))),
  by = "row.names", all = TRUE)
gene.stable <- merge(x = miRNA.info[, 1:2], y = gene.stable, by.x = "gene_id",
                     by.y = "Row.names")
rownames(gene.stable) <- gene.stable[, "gene_id"]
gene.stable <- gene.stable[, -1]
colnames(gene.stable)[2:3] <- c("meanCPM", "SdevCPM")

# Calculate percentage of Standard deviation in order to compare all genes
gene.stable <- cbind(gene.stable, Perc_Sdev = (
  gene.stable[, "SdevCPM"]*100/gene.stable[, "meanCPM"]))
gene.stable <- gene.stable[order(gene.stable$Perc_Sdev),]
head(gene.stable)

# Plot the mean CPM and Standard deviation CPM of all genes
ggplot(data = gene.stable, aes(x = log(x = gene.stable$meanCPM, base = 2),
                               y = log(x = gene.stable$SdevCPM,
                                       base = 2))) + geom_point(size = 4.5)

# Write into a table the analysis for reference gene discovery
write.matrix(x = gene.stable, file = paste(
  method, "_potential_refgene.txt", sep = ""), sep = "\t")

###########################################################
# Comparison of Exiqon PCR results with miRNA-seq results #
###########################################################

# Read in the csv file containing PCR results from Kirsten
PCR.Kirsten <- read.table(file = paste(
  user, "/Home_work_sync/Work/Colleagues shared work/Kirsten McLoughlin/",
  "miRNA_work/PCR_data/PCR_to_compare.txt", sep = ""), header = TRUE,
  sep = "\t", quote = "", na.strings = '\'N/A\'')

# Reformat Kirsten's PCR data
PCR.Kirsten <- PCR.Kirsten[, c(1,2,4,7,12)]
PCR.Kirsten <- cbind(PCR.Kirsten, as.matrix(as.numeric(lapply(X = PCR.Kirsten[
  , 'Geometric.mean.fold.change..relative.to..time.point.1.'],
  FUN = function(i) {if (i < 0){-log(x = abs(i), base = 2)} else {log(x = abs(
    i), base = 2)}}))))
colnames(PCR.Kirsten) <- c("hsapiens_gene_name", "sequence", "time_point",
                   "fold_change", "FDR", "logFC")
PCR.Kirsten$time_point %<>% gsub(pattern = "\\+(\\d*) wk", replacement = "\\1",
                                 x = .) %>% gsub(pattern = "(-\\d*) wk",
                                                 replacement = "\\1", x = .)
PCR.Kirsten <- cbind(PCR.Kirsten, Methods = rep(
  x = "Kirsten_PCR", times = length(rownames(PCR.Kirsten))))
head(PCR.Kirsten)

# Read the file containing biomarkers RT-qPCR results
RT.qPCR <- read.table(file = paste(user, "/Home_work_sync/Work/TIDA/miRNA_PCR",
                                   "/Diff_expr/RT-qPCR_DE.txt", sep = ""),
                      sep = "\t", header = TRUE)
rownames(RT.qPCR) <- RT.qPCR$gene_name %>% paste(
  ., RT.qPCR$time_point, sep = "_") %>% gsub(
    pattern = "-(\\d*)$", replacement = "pre\\1", x = .)
RT.qPCR <- RT.qPCR[, -c(1)]
head(RT.qPCR)

# Merge the miRNA-seq, RT-qPCR and Kirsten PCR data
Comp <- NULL
for (time in c(1, 2, 6, 10, 12)) {
  data1 <- RT.qPCR[RT.qPCR$time_point == time, c(
    "sequence", "meanlogFC", "se", "final.Pvalue", "Methods")]
  colnames(data1) <- c("sequence", "logFC_N.C.PCR", "se_N.C.PCR",
                       "Pvalue_N.C.PCR", "Methods_N.C.PCR")
  data2 <- PCR.Kirsten[PCR.Kirsten$time_point == time & !is.na(
    PCR.Kirsten$sequence), c("sequence", "hsapiens_gene_name", "logFC", "FDR",
                             "Methods")]
  colnames(data2) <- c("sequence", "hsapiens_gene_name", "logFC_K.PCR",
                       "FDR_K.PCR", "Methods_K.PCR")
  method.DE <- eval(parse(text = paste(method, ".DE", sep = "")))
  pattern.time <- colnames(method.DE) %>% grep(pattern = paste(
    time, "w", sep = ""), x = ., value = TRUE) %>% grep(
      pattern = "logFC|FDR", x = ., value = TRUE)
  pattern.id <- grep(pattern = "_", x = rownames(method.DE),
                     invert = TRUE)
  data3 <- method.DE[pattern.id, c("gene_id", "gene_name", "sequence",
                                        pattern.time[1], pattern.time[2])]
  colnames(data3) <- c("gene_id", "gene_name", "sequence", "logFC_RNAseq",
                       "FDR_RNAseq")
  if (is.null(Comp)) {
    Comp <- merge(x = data1, y = data2, by = "sequence", all = TRUE) %>% merge(
      x = data3, y = ., by = "sequence", all.y = TRUE)
    Comp <- cbind(Comp, time_point = rep(x = time, times = length(
      Comp$sequence)))
  }
  else {
    tmp <- merge(x = data1, y = data2, by = "sequence", all = TRUE) %>% merge(
      x = data3, y = ., by = "sequence", all.y = TRUE)
    tmp <- cbind(tmp, time_point = rep(x = time, times = length(tmp$sequence)))
    Comp <- rbind(Comp, tmp)
  }
}
Comp <- Comp[!is.na(Comp$sequence) & !is.na(Comp$gene_id),]
Comp <- Comp[!(is.na(Comp$logFC_N.C.PCR & is.na(Comp$logFC_RNAseq))), ]
Comp

# Perform overall correlation of log fold-change between PCR and RNA-seq data
cor.gene <- list()
cor.gene["logFC_RNAseq_vs_logFC_N.C.PCR"] <- list(cor.test(
  x = Comp$logFC_RNAseq, y = Comp$logFC_N.C.PCR, alternative = "two.sided",
  method = "spearman", na.action = na.omit()))
cor.gene["logFC_RNAseq_vs_logFC_K.PCR"] <- list(cor.test(
  x = Comp$logFC_RNAseq, y = Comp$logFC_K.PCR, alternative = "two.sided",
  method = "spearman", na.action = na.omit()))
cor.gene["logFC_N.C.PCR_vs_logFC_K.PCR"] <- list(cor.test(
  x = Comp$logFC_N.C.PCR, y = Comp$logFC_K.PCR, alternative = "two.sided",
  method = "spearman", na.action = na.omit()))

# Perform correlation of log fold-change between PCR and RNA-seq data per gene
for (g in unique(Comp[!is.na(Comp$logFC_N.C.PCR), "gene_id"])) {
  dat1 <- Comp[Comp$gene_id == g, "logFC_N.C.PCR"]
  dat2 <- Comp[Comp$gene_id == g, "logFC_RNAseq"]
  res <- cor.test(x = dat1, y = dat2, alternative = "two.sided",
                  method = "spearman", na.action = na.omit())
  cor.gene[g] <- list(res)
}

# Create a variable containing the current correlation results
assign(paste(method, ".cor.gene", sep = ""), cor.gene)

# Determine concordance of direction of expression between techniques
# and concordance of P-value significance between techniques
concordance <- list()
concor.fc <- 0
concor.Pval <- 0
counter <- 0
for (x in 1:nrow(Comp)) {
  if (!is.na(Comp$logFC_N.C.PCR[x])) {
    counter <- counter+1
    if (Comp[x, "logFC_N.C.PCR"] > 0 & Comp[x, "logFC_RNAseq"] > 0) {
      concor.fc <- concor.fc+1
    }
    else if (Comp[x, "logFC_N.C.PCR"] < 0 & Comp[x, "logFC_RNAseq"] < 0) {
      concor.fc <- concor.fc+1
    }
    if (Comp[x, "Pvalue_N.C.PCR"] > 0.05 & Comp[x, "FDR_RNAseq"] > 0.05) {
      concor.Pval <- concor.Pval+1
    }
    else if (Comp[x, "Pvalue_N.C.PCR"] < 0 & Comp[x, "FDR_RNAseq"] < 0) {
      concor.Pval <- concor.Pval+1
    }
  }
  if (x == nrow(Comp)) {
    concor.fc <- (concor.fc*100/counter)
    concor.Pval <- (concor.Pval*100/length(Comp$gene_id))
  }
}
concordance["N.C.PCR-RNAseq"] <- list(list("logFC_concordance" = concor.fc, "Pvalue_concordance" = concor.Pval))

concor.fc <- 0
concor.Pval <- 0
counter <- 0
for (x in 1:nrow(Comp)) {
  if (!is.na(Comp$logFC_K.PCR[x])) {
    counter <- counter+1
    if (Comp[x, "logFC_K.PCR"] > 0 & Comp[x, "logFC_RNAseq"] > 0) {
      concor.fc <- concor.fc+1
    }
    else if (Comp[x, "logFC_K.PCR"] < 0 & Comp[x, "logFC_RNAseq"] < 0) {
      concor.fc <- concor.fc+1
    }
    if (Comp[x, "FDR_K.PCR"] > 0.05 & Comp[x, "FDR_RNAseq"] > 0.05) {
      concor.Pval <- concor.Pval+1
    }
    else if (Comp[x, "FDR_K.PCR"] < 0 & Comp[x, "FDR_RNAseq"] < 0) {
      concor.Pval <- concor.Pval+1
    }
  }
  if (x == nrow(Comp)) {
    concor.fc <- (concor.fc*100/counter)
    concor.Pval <- (concor.Pval*100/length(Comp$gene_id))
  }
}
concordance["K.PCR-RNAseq"] <- list(list("logFC_concordance" = concor.fc, "Pvalue_concordance" = concor.Pval))

# Create a variable containing the current concordance results
assign(paste(method, ".concordance", sep = ""), concordance)

############################
# Clean up the R workspace #
############################

# Remove temporary and unrequired variables
rm(CPM, Comp, Count, DEtable, MDS.ggplot, PCR.Kirsten, RT.qPCR, count.log2,
   Dat, data1, data2, data3, design, gene.stable, method.DE, miRNA.duplicate,
   rf.ImportanceOrdered, target, tmp, DEtable.name, MDS, animal, col.side, dat,
   colors, colours, concor.Pval, concor.fc, concordance, cor.gene, counter,
   dat1, dat2, dgeglm.fit, dgelist, dgelist.disp, dgelist.filt, dgelist.norm,
   g, group, hist_plot, i, method, my.palette, n, overlap, pat, pattern.id,
   pattern.time, pheno, res, rf, sig.10w, sig.12w, sig.1w, sig.2w, sig.6w,
   time, user, val, x, DESeq.merge, DEtable.merge, diff_expr_edgeR, g.legend,
   loadpackage, miR.DE, multi.DE, multi.MDS, sig_label, sample, files)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
##############################################################################
# Comparison of time overlapping DE genes between various analysis pipelines #
##############################################################################

# Create the Venn diagram of overlapping DE genes between time points
venn.diagram(x = list("Novoalign" = Novoalign.overlap,
                      "miRdeep2" = miRdeep2.overlap,
                      "miRdeepstar" = miRdeepstar.overlap),
             filename = "Venn_3pipelines.tiff", na="remove")

# Get all the gene expression information for the candidate biomarkers miRNA
biomarker <- merge(x = miRNA.info[(miRNA.info$gene_id %in% unique(c(
  Novoalign.overlap, miRdeep2.overlap, miRdeepstar.overlap))),],
  y = Novoalign.DE[, (ncol(Novoalign.DE)-24):ncol(Novoalign.DE)],
  by.x = "gene_id", by.y = "row.names", all.x = TRUE)
biomarker <- merge(x = biomarker,
                   y = miRdeep2.DE[, (ncol(miRdeep2.DE)-24):ncol(miRdeep2.DE)],
                   by.x = "gene_id", by.y = "row.names", all.x = TRUE)
biomarker <- merge(x = biomarker, y = miRdeepstar.DE[, (ncol(
  miRdeepstar.DE)-24):ncol(miRdeepstar.DE)], by.x = "gene_id",
  by.y = "row.names", all.x = TRUE)

# Format the variable columns
colnames(biomarker) %<>% gsub(pattern = "w$", replacement = "_miRdeepstar",
                              x = ., perl = TRUE) %>% gsub(
                                pattern = "\\.y$", replacement = "_miRdeep2",
                                x = ., perl = TRUE) %>% gsub(
                                  pattern = "\\.x$", replacement = "_novoalign",
                                  x = ., perl = TRUE)

# Write into a table all the time overlapping DE genes
write.matrix(x = biomarker, file = "Candidate_biomarkers.txt", sep = "\t")

# Remove the temporary variables
rm(miRNA.info, biomarker)

#######
# END #
#######