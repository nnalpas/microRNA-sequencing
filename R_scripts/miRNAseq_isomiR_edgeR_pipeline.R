

### miRNA-seq data analysis of isomiR counts (paired data) ------------------


# Analysis of 10 animals infected over a 15 weeks time course for which 
# miRNA-seq libraries where prepared from serum samples, and this analysis
# focuses on isomiRs expression
# Information regarding samples preparation can be found in LabBook number 68


### List of required packages -----------------------------------------------

# Extract the User path to Dropbox
user <- gsub(
    pattern = "(^.*/Dropbox)/.*", replacement = "\\1",
    x = getwd(), perl = TRUE)

# Source the user's general functions used across this script
source(
    file = paste(
        user, "/Home_work_sync/Work/Bioinformatics/R/General_function.R",
        sep = ""))

# Load the required packages (or install if not already in library)
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
loadpackage(package = plyr)


### Gene annotation using information obtained from GTF file ----------------

# Read in the annotation information
miRNA.info <- read.table(
    file = "../miRNA_Btaurus.txt", header = TRUE, sep = "\t", quote = "")

# Determine which miRNA have identical mature sequence
miRNA.duplicate <- aggregate(
    gene_id ~ sequence, FUN = "as.vector",
    data = miRNA.info, na.action = "as.vector")

# Combine new info with miRNA annotation
miRNA.info <- merge(
    x = miRNA.info, y = miRNA.duplicate, by = "sequence", all = TRUE)

# Reorganise and rename columns for miRNA annotation
miRNA.info <- miRNA.info[, c(2:7, 1, 8:ncol(miRNA.info))]
colnames(miRNA.info)[c(1,7,12)] <- c(
    "gene_id", "sequence", "identical_sequence")

# Check the miRNA annotation variable
head(miRNA.info)
tail(miRNA.info)
dim(miRNA.info)

# Save the data to keep in R objects
save(miRNA.info, file = "miRNA_info.RData", compress = TRUE)


### Read in and concatenate count input files within R -------------------------

# Define a set of files containing isomiR counts to process
files <- list.files(
    path = paste(
        user, "/Home_work_sync/Work/TIDA/miRNA-seq/Results",
        "/Sequencing_depth/isomiR_counts",
        sep = ""),
    pattern = "6",
    full.names = TRUE)

# Use a loop to read all the isomiR count file
for (i in 1:length(files)) {
    
    # Read in the current file and keep only the sequence and count columns
    Dat <- read.table(
        file = files[i], quote = "\"", sep = "\t",
        header = TRUE, stringsAsFactors = FALSE)
    
    # Rename the matching pri-miRNA information column
    colnames(Dat)[3] <- c("primiRNA")
    
    # Collect the library name from the file name
    lib <- gsub(
        pattern = "^.*\\/(6\\d+.*)_isomiR.txt$", replacement = "\\1_",
        x = files[i], perl = TRUE)
    
    # Rename the pri-miRNA information and count columns with library name
    colnames(Dat)[3:4] <- gsub(
        pattern = "^", replacement = lib, x = colnames(Dat)[3:4], perl = TRUE)
    
    # Transform all sequence to uppercase (for merging compatibility)
    Dat$Sequence <- toupper(x = Dat$Sequence)
    
    # For the first ocurence of i define the count variable that will
    # store all files isomiR counts
    if (i == 1) {
        Count <- Dat[, -1]
    }
    # Otherwise merge the count variable with the isomiR count from the
    # current file being processed
    else {
        Count <- merge(x = Count, y = Dat[, -1], by = "Sequence", all = TRUE)
    }
    
}

# Add a rownames containing unique ID to the count variable
rownames(Count) <- paste("iso", rownames(Count), sep = "_")

# Extract the annotation information columns
primirna <- Count[, c(
    "Sequence", grep(pattern = "_primiRNA", x = colnames(Count), value = TRUE))
    ]

# Check the raw pri-miRNA information data
dim(primirna)
head(primirna, n = 3)
tail(primirna, n = 3)

# Extract the sequence and count columns
Count <- Count[, c(
    "Sequence", grep(pattern = "_Count", x = colnames(Count), value = TRUE))
    ]

# Create an artificial 0 count for the sequence not found in a library
Count[is.na(x = Count)] <- 0

# Include the rownames of Count dataframe within the dataframe
Count <- cbind(rownames(Count), Count)

# Rename the column with artificial isomiR ID
colnames(Count)[1] <- c("isomiR_ID")
colnames(Count) <- gsub(
    pattern = "_Count", replacement = "", x = colnames(Count))

# Check the raw count data
dim(Count)
head(Count, n = 3)
tail(Count, n = 3)

# Save the data to keep in R objects
save(Count, file = "All_count.RData", compress = TRUE)
save(primirna, file = "All_primiRNA.RData", compress = TRUE)


### Prepare miRNA annotation information ---------------------------------

# Reformat pri-miRNA info into a single column, remove all row with NA value
# and remove duplicated sequence values
annot <- melt(data = primirna, id = "Sequence") %>%
    .[!is.na(.$value), -2] %>%
    .[!duplicated(.$Sequence),]

# Rename the pri-miRNA info dataframe column name
colnames(annot)[2] <- c("primiRNA_info")

# Extract the pri-miRNA info that are consistent across libraries
annot$primiRNA_info <- gsub(
        pattern = "(MI\\d+ \\(.+?),.+?(, mm.+?\\[.+?), #.+?(\\]\\))",
        replacement = "\\1\\2\\3",
        x = annot$primiRNA_info,
        perl = TRUE)

# Check the pri-miRNA info data
dim(annot)
head(annot)
tail(annot)

# Merge the pri-miRNA info data with the count data
Count.annot <- merge(x = annot, y = Count, by = "Sequence", all = TRUE)

# Reformat the annotated count dataframe and include rownames
Count.annot <- Count.annot[, c(3, 1:2, 4:ncol(Count.annot))]
rownames(Count.annot) <- Count.annot$isomiR_ID

# Check the annotated count data
dim(Count.annot)
head(Count.annot, n = 3)
tail(Count.annot, n = 3)

# Save the data to keep in R objects
save(Count.annot, file = "All_count_annotated.RData", compress = TRUE)





# To be checked at later stage


# Merge the annotation information with the count table
Count <- merge(x = miRNA.info, y = Count, by = "gene_id", all = TRUE)
head(Count)

# Check the merged matrix in terms of size and content
table(duplicated(Count[, "gene_id"]))
length(unique(Count[, "gene_id"], incomparables = FALSE))


### Create groups and a DGElist ------------------------------------------

# Create a target matrix and a experimental group vector and animal block
target <- colnames(Count.annot[, 4:ncol(Count.annot)])
target %<>%
    strsplit2(x = ., split = "_", fixed = TRUE) %>%
    cbind(target, .)

# Add column name to the target matrix
colnames(target) <- c("sample", "animal", "time_point", "depth")

# Define the different factors fo time point, animal ID and sequencing depth
time.point <- factor(
    x = target[, "time_point"], levels = c("pre2", "pre1", 1, 2, 6, 10, 12),
    ordered = TRUE)
animal <- factor(x = target[, "animal"])
depth <- factor(
    x = target[, "depth"],
    levels = c(
        "500000", "1000000", "1500000", "2000000", "2500000", "3000000",
        "3500000", "4000000", "4500000", "5000000", "5500000", "6000000",
        "full"),
    ordered = TRUE)

# Create a DGElist containing the group information
dgelist.all <- DGEList(
    counts = Count.annot[, 4:ncol(Count.annot)],
    lib.size = NULL,
    norm.factors = NULL,
    group = time.point,
    genes = Count.annot[, 1:3],
    remove.zeros = FALSE)

# Check the DGElist data
names(dgelist.all)
head(dgelist.all$samples)
head(dgelist.all$counts)
head(dgelist.all$genes)

# Ouptut samples data
write.table(
    x = dgelist.all$samples, file = "All_sample.txt", sep = "\t",
    quote = FALSE, row.names = TRUE, col.names = TRUE)

# Save the data to keep in R objects
save(dgelist.all, file = "All_dgelist.RData", compress = TRUE)
save(
    target, time.point, animal, depth,
    file = "All_phenotype.RData", compress = TRUE)


### Data subsetting ------------------------------------------------------

# Create a variable to hold the method value for subsetting
method <- levels(depth)[13]

# Create a variable containing the subsetted DGElist
dgelist <- dgelist.all[, grep(
    pattern = paste("^.*_", method, "$", sep = ""),
    x = colnames(dgelist.all),
    perl = TRUE)
    ]


### Quality check of libraries by plotting density of count --------------

# Log2 transform the count data for better visualization
count.log2 <- log(x = (dgelist$counts[,]+1), base = 2)

# Plot density of count for all libraries
png(filename = paste(method, "_Dens_all_miRNA.png", sep = ""),
    width = 1500, height = 1300, units = "px")
plot(
    x = density(count.log2[, 1]), main = "Density plot of count per isomiR", 
    lty =  1, xlab = "Log2 of count per isomiR", ylab = "Density",
    col = 1, xlim = c(-0.5, 6), ylim = c(0.0, 1.0), lwd = 2 ,
    cex.axis = 1.5, cex.lab = 1.5, cex.main = 2
    )
for (i in 2:ncol(count.log2)) {
    lines(x = density(count.log2[, i]), lty = 1, col = i)
}
dev.off()


### Filtering of lowly expressed tags ------------------------------------

# Identify genes with zero counts across all samples
dim(dgelist[rowSums(dgelist$counts) == 0, ])
head(dgelist[rowSums(dgelist$counts) == 0, ])

# Filter lowly expressed tags, retaining only tags with at least 50 counts per 
# million in 10 or more libraries (10 libraries correspond to one time point)
dgelist.filt <- dgelist[rowSums(
    cpm(dgelist$counts) > 50) >= median(summary(dgelist$samples$group)), ]
dim(dgelist.filt$counts)

# Compute the library size
dgelist.filt$samples$lib.size <- colSums(dgelist.filt$counts)
head(dgelist.filt$samples)
head(dgelist$samples)


### Quality check of libraries by plotting density of count after --------

# Log2 transform the count data for better visualization
count.log2 <- log(x = (dgelist.filt$counts[,] + 1), base = 2)

# Plot density of count for all libraries
png(filename = paste(method, "_Dens_filt_miRNA.png", sep = ""),
    width = 1366, height = 768, units = "px")
plot(
    x = density(count.log2[, 1]), main = "Density plot of count per isomiR", 
    lty =  1, xlab = "Log2 of count per isomiR", ylab = "Density", col = 1, 
    xlim = c(-0.5, 15), ylim = c(0.0, 1.0), lwd = 2 , cex.axis = 1.5, 
    cex.lab = 1.5, cex.main = 2
    )
for (i in 2:ncol(count.log2)) {
    lines(x = density(count.log2[, i]), lty = 1, col = i)
}
dev.off()


### Histogram of CPM per expression category -----------------------------

# Prepare the histogram picture output
png(filename = paste(method, "_Hist_miRNA.png", sep = ""),
    width = 1500, height = 1500, units = "px")
par(mfrow = c(2,2), cex = 3, cex.axis = 0.7)

# Plot histogram of number of gene per CPM category for all microRNA
hist_plot <- hist(
    x = rowMeans(x = cpm(dgelist$counts), na.rm = FALSE),
    breaks = c(0, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000),
    plot = FALSE)
plot(
    x = hist_plot$counts, log = "x", type = "h", axes = FALSE, 
    xlab = "Average CPM", ylab = "Number of gene", lwd = 10, 
    col = "lightblue", main = "All miRNA")
axis(
    side = 1, at = c(1:10),
    labels = c(
        "0", "0.01", "0.1", "1", "10", "100", "1000",
        "10000", "100000", "1000000"))
axis(side = 2)

# Plot histogram of number of gene per CPM category for microRNA with 
# at least one count
hist_plot <- hist(
    x = rowMeans(
        x = cpm(dgelist$counts[rowSums(dgelist$counts) > 0, ]), na.rm = FALSE),
    breaks = c(0, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000),
    plot = FALSE)
plot(
    x = hist_plot$counts, log = "x", type = "h", axes = FALSE, 
    xlab = "Average CPM", ylab = "Number of gene", lwd = 10, 
    col = "lightgreen", main = "Min one count miRNA")
axis(
    side = 1, at = c(1:10),
    labels = c(
        "0", "0.01", "0.1", "1", "10", "100", "1000",
        "10000", "100000", "1000000"))
axis(side = 2)

# Plot histogram of number of gene per CPM category for microRNA 
# filtered for low expression
hist_plot <- hist(
    x = rowMeans(x = cpm(dgelist.filt$counts), na.rm = FALSE),
    breaks = c(0, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000), 
    plot = FALSE)
plot(
    x = hist_plot$counts, log = "x", type = "h", axes = FALSE, 
    xlab = "Average CPM", ylab = "Number of gene", lwd = 10, 
    col = "lightcoral", main = "Low expression filtered miRNA")
axis(
    side = 1, at = c(1:10), 
    labels = c(
        "0", "0.01", "0.1", "1", "10", "100", "1000",
        "10000", "100000", "1000000"))
axis(side = 2)

# Close graphic device
dev.off()


### Gene expression correlation between libraries ------------------------

# Perform gene expression correlation between libraries at each time points to 
# identify potential outlier library (CPM values not required with Spearman)
for (val in levels(time.point)) {
    png(filename = paste(method, val, "cor.png", sep = "_"),
        width = 1500, height = 1500, units = "px")
    chart.Correlation(
        R = log(x = (dgelist.filt$counts[, grep(
            pattern = paste("_", val, "_", sep = ""),
            x = colnames(dgelist.filt$counts),
            perl = TRUE)] + 1), base = 2),
        histogram = TRUE, method = "spearman",
        main = paste(val, " time point expression correlation", sep = ""))
    dev.off()
}


### Normalization of data using trimmed mean of M-values -----------------

# Calculate normalisation factor for our DGElist, note that with edgeR 
# the counts are not transformed in any way after normalization
dgelist.norm <- calcNormFactors(dgelist.filt)
dgelist.norm$samples


### Multidimensional scaling plot of all samples -------------------------

# Define a set of colors
colours <- brewer.pal(
    n = length(levels(time.point)),
    name = "Set1")

# Use custom function to output a MDS plot of all samples
multi.MDS(
    data = dgelist.norm,
    target = target[grep(
        pattern = paste("^.*_", method, "$", sep = ""),
        x = target[, c("sample")], perl = TRUE),],
    prefix = method, suffix = "all",
    plotmds = list(
        top = 1000000, gene.selection = "pairwise",
        dim.plot = c(1, 2)),
    aes.colour = "time_point", fill = "lightgrey", size = 3,
    legend.pos = "right", manual.colour = colours, breaks = levels(time.point))


### Multidimensional scaling plot per time points ------------------------

# Use custom function to output a MDS plot per post-infection time points
multi.MDS(
    pattern = c(
        "_pre2_|_pre1_|_1_", "_pre2_|_pre1_|_2_", "_pre2_|_pre1_|_6_",
        "_pre2_|_pre1_|_10_", "_pre2_|_pre1_|_12_"),
    data = dgelist.norm,
    target = target[grep(
        pattern = paste("^.*_", method, "$", sep = ""),
        x = target[, c("sample")], perl = TRUE),],
    prefix = method,
    suffix = c("W1", "W2", "W6", "W10", "W12"), 
    plotmds = list(
        top = 1000000, gene.selection = "pairwise", dim.plot = c(1, 2)),
    aes.colour = "time_point", fill = "wheat", size = 3, legend.pos = "right",
    combine = TRUE, manual.colour = colours[c(2, 3, rep(x = 1, times = 5))], 
    breaks = levels(time.point))


### Create a design matrix for paired analysis ---------------------------

# Define the different factors for time point and animal ID after subsetting
time.point.sub <- factor(x = gsub(
    pattern = "^.*?_(.*?)_.*$", replacement = "\\1",
    x = rownames(dgelist.norm$samples), perl = TRUE))
animal.sub <- factor(x = gsub(
    pattern = "(6.*?)_.*$", replacement = "\\1",
    x = rownames(dgelist.norm$samples), perl = TRUE))

# Create a design matrix
design <- model.matrix(~ animal.sub + time.point.sub)
rownames(design) <- rownames(dgelist.norm$samples)
colnames(design) <- gsub(
    pattern = "(animal.sub)|(time.point.sub)", replacement = "",
    x = colnames(design), perl = TRUE)
design


### Estimate the dispersion parameter for each tag using Cox-Reid --------

# Calculate the dispersion (common, trended and tagwise)
dgelist.disp <- estimateGLMCommonDisp(
    y = dgelist.norm, design = design, verbose = TRUE)
dgelist.disp <- estimateGLMTrendedDisp(y = dgelist.disp, design = design)
dgelist.disp <- estimateGLMTagwiseDisp(y = dgelist.disp, design = design)
names(dgelist.disp)

# Plot the dispersion
png(
    filename = paste(method, "_BCV.png", sep = ""),
    width = 1366, height = 768, units = "px")
plotBCV(dgelist.disp)
dev.off()

# Show the calculated dispersion and the coefficient of biological variation
dgelist.disp$common.dispersion
sqrt(dgelist.disp$common.dispersion)


### Determine differential expression using negative binomial GLMs -------

# Fit a negative binomial generalized linear model for each tag using
# the design matrix and calculated dispersion
dgeglm.fit <- glmFit(y = dgelist.disp, design = design)
names(dgeglm.fit)


### Main differential expression call ------------------------------------

# Test for differential expression for -1 week versus -2 week
multi.DE(
    "pre1", 1, "pre2", -1, data = dgeglm.fit, design = design,
    group = time.point.sub, adjpvalue = 0.05, method = "BH",
    dedata = paste("de", method, "pre1w", sep = "."),
    smearfile = paste(method, "_smear_pre1w", sep = ""))

# Test for differential expression for 1 week versus -2 week and -1 week
multi.DE(
    1, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
    group = time.point.sub, adjpvalue = 0.05, method = "BH",
    dedata = paste("de", method, "1w", sep = "."),
    smearfile = paste(method, "_smear_1w", sep = ""))

# Test for differential expression for 2 week versus -2 week and -1 week
multi.DE(
    2, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
    group = time.point.sub, adjpvalue = 0.05, method = "BH",
    dedata = paste("de", method, "2w", sep = "."),
    smearfile = paste(method, "_smear_2w", sep = ""))

# Test for differential expression for 6 week versus -2 week and -1 week
multi.DE(
    6, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
    group = time.point.sub, adjpvalue = 0.05, method = "BH",
    dedata = paste("de", method, "6w", sep = "."),
    smearfile = paste(method, "_smear_6w", sep = ""))

# Test for differential expression for 10 week versus -2 week and -1 week
multi.DE(
    10, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
    group = time.point.sub, adjpvalue = 0.05, method = "BH",
    dedata = paste("de", method, "10w", sep = "."),
    smearfile = paste(method, "_smear_10w", sep = ""))

# Test for differential expression for 12 week versus -2 week and -1 week
multi.DE(
    12, 1, "pre1", -1, "pre2", -1, data = dgeglm.fit, design = design,
    group = time.point.sub, adjpvalue = 0.05, method = "BH",
    dedata = paste("de", method, "12w", sep = "."),
    smearfile = paste(method, "_smear_12w", sep = ""))


### Merge all DE call data from the different time points ----------------

# Create a vector of DE table to merge
DEtable.name <- c("1w", "2w", "6w", "10w", "12w") %>%
    lapply(X = ., FUN = function(x) paste("de", method, x, sep = ".")) %>%
    unlist(.)

# Merge all DE table for the different time points into a single dataframe
DEtable.merge(
    DEtable.name, output = paste("Full.de", method, sep = "."),
    pattern = paste("^de.", method, '.', sep = ""))

# Create a variable containing the current DEtable
DEtable <- eval(parse(text = paste("Full.de", method, sep = ".")))
head(DEtable)

# Write into a table the full DE call data
write.matrix(
    x = DEtable, file = paste(method, "_full_DE.txt", sep = ""), sep = "\t")


### Plot the number of differentially expressed genes --------------------

# Plot the number of differentially expressed genes per comparisons
plot.numb.DE(
    data = DEtable, comparison = c("1w", "2w", "6w", "10w", "12w"),
    filename = paste(method, "_DE.tif", sep = ""))


### Comparison of DE genes between time points ---------------------------

# Identify as a vector list the significant DE genes per time point
venn.de(
    data = DEtable,
    comparison = c("1w", "2w", "6w", "10w", "12w"),
    picname = paste(method, "_Venn_All", sep = ""),
    overlapname = paste("overlap", method, sep = "."),
    lwd = 0.7, cex = 1, cat.cex = 1)

# Create a variable containing the current overlap genes
overlap <- eval(parse(text = paste("overlap", method, sep = ".")))


### Clustering of samples based on common DE genes -----------------------

# Get the CPM values for each filtered miRNA
CPM <- cpm(x = dgelist.disp, normalized.lib.sizes = TRUE)
dim(CPM)
head(CPM)

# Create a phenotypic data table with color coding and infection status
pheno <- colnames(CPM)
pheno <- cbind(
    sample = pheno,
    status = unlist(as.character(lapply(
        X = pheno,
        FUN = function(x) {
            if(length(grep(
                pattern = "_pre(1|2)_", x = x, perl = TRUE)
                ) == 1) { "healthy" }
            else { "infected" }
            }))),
    col.side = unlist(as.character(lapply(
        X = pheno, FUN = function(x) {
            if(length(grep(
                pattern = "_pre(1|2)_", x = x, perl = TRUE)
                ) == 1) { colours[2] }
            else { colours[1] }
            }))))
rownames(pheno) <- pheno[, "sample"]
col.side <- pheno[, "col.side"]
pheno <- factor(x = pheno[, "status"])
pheno

# Define breaks and colours for the heatmap
colors <- c(
    seq(from = 0, to = 1, length = 10),
    seq(from = 1, to = 10, length = 10),
    seq(from = 10, to = 100, length = 10),
    seq(from = 100, to = 1000, length = 10),
    seq(from = 1000, to = 10000, length = 10),
    seq(from = 10000, to = 100000, length = 10),
    seq(from = 100000, to = 1000000, length = 10))
my.palette <- colorRampPalette(c("blue", "red"))(n = 69)

# Represent the data in a heatmap based on the overlapping miRNA
tiff(
    filename = paste(method, "_overlap_Heatmap.tif", sep = ""),
    width = 15, height = 15, units = "in", res = 600, compression = "lzw")
heatmap.2(
    x = CPM[overlap,], labRow = rownames(CPM), labCol = colnames(CPM),
    cexRow = .5, cexCol = .5, breaks = colors, col = my.palette,
    symm = F, symkey = F, symbreaks = T, scale = "none",
    ColSideColors = col.side)
dev.off()


### Biomarkers discovery using Random Forest method ----------------------

# Use randomForest on the CPM values per miRNA
rf <- randomForest(
    x = t(CPM), y = pheno, importance = TRUE, replace = TRUE, do.trace = FALSE,
    ntree = 5000, mtry = sqrt(nrow(CPM)), keep.forest = TRUE, proximity = TRUE)

# Order the miRNA by their importance value
rf.ImportanceOrdered <- importance(rf)[order(
    importance(rf)[, "MeanDecreaseGini"], decreasing = TRUE),]
head(rf.ImportanceOrdered, n = 20)

# Create a MDS plot
MDSplot(rf = rf, fac = pheno, palette = c(1,3))

# Represent the data in a heatmap based on the top miRNA (by importance)
for (n in 2:15) {
    tiff(
        filename = paste(method, "_Heatmap_n", n, ".tif", sep = ""),
        width = 15, height = 15, units = "in", res = 600, compression = "lzw")
    dat <- CPM[head(rownames(rf.ImportanceOrdered), n = n),]
    heatmap.2(
        x = dat, labRow = rownames(dat), labCol = colnames(dat),
        cexRow = .5, cexCol = .5, breaks = colors, col = my.palette,
        symm = F, symkey = F, symbreaks = T, scale = "none",
        ColSideColors = col.side)
    dev.off()
}


### Additional differential expression calls -----------------------------

# Perform all differential expression analysis within a for loop
time <- levels(time.point)
for (i in time) {
    print("Performing comparison versus:")
    print(i)
    DEtable.name <- c()
    for (j in 1:length(time)) {
        if (time[j] != i) {
            # Test for differential expression for comparison
            name <- paste(time[j], "wvs", i, "w", sep = "")
            DEtable.name <- c(DEtable.name, name)
            multi.DE(
                time[j], 1, i, -1, data = dgeglm.fit, design = design,
                group = time.point.sub, adjpvalue = 0.05, method = "BH",
                dedata = paste("de", method, name, sep = '.'))
        }
        else {
            next
        }
    }
    # Create a vector of DE table to merge
    DEtable.method <- DEtable.name %>%
        lapply(X = ., FUN = function(x) paste("de", method, x, sep = '.')) %>%
        unlist(.)
    # Merge all DE table for the different time points into a single dataframe
    DEtable.merge(
        DEtable.method, output = paste("DE.vs.", i, "w.", method, sep = ""),
        pattern = paste("^", "de.", method, '.', sep = ""))
    # Create a variable containing the current DEtable
    DEtable <- eval(
        expr = parse(text = paste("DE.vs.", i, "w.", method, sep = "")),
        envir = .GlobalEnv)
    colnames(DEtable)
    # Write into a table the full DE call data
    write.matrix(
        x = DEtable, file = paste(method, "_DE_vs", i, "w.txt", sep = ""),
        sep = "\t")
    # Plot the number of differentially expressed genes per comparisons
    plot.numb.DE(
        data = DEtable, comparison = DEtable.name,
        pattern = "vs", replace = " vs ",
        filename = paste(method, "_DE_vs", i, "w.tif", sep = ""))
    # Comparison of DE genes between time points
    # Identify as a vector list the significant DE genes per time point
    venn.de(
        data = DEtable, comparison = DEtable.name[-1],
        pattern = "vs", replace = " vs ",
        picname = paste(method, "_Venn_vs", i, "w", sep = ""),
        overlapname = paste("overlap.vs", i, "w.", method, sep = ""),
        lwd = 0.7, cex = 1, cat.cex = 1)
}

# Rename the files containing overlapping list of isomiRs between time points
list.files(pattern = "^overlap.*.txt$") %>%
    lapply(
        X = .,
        FUN = function(y) {
            file.rename(
                from = y,
                to = gsub(
                    pattern = "\\.(vs.+?)", replacement = "_\\1", x = y) %>%
                    gsub(
                        pattern = "^(.+)\\.(.+?)(\\.txt)$",
                        replacement = "\\2_\\1\\3", x = .))
            }
        )


### Comparison of Exiqon PCR results with miRNA-seq results -----------------

# Read in the file containing PCR results from Kirsten
PCR.Kirsten <- read.table(file = paste(
    user, "/Home_work_sync/Work/Colleagues shared work/Kirsten McLoughlin/",
    "miRNA_work/PCR_data/PCR_to_compare.txt", sep = ""),
    header = TRUE, sep = "\t", quote = "", na.strings = '\'N/A\'')

# Reformat Kirsten's PCR data
PCR.Kirsten <- PCR.Kirsten[, c(1,2,4,7,12)]
colnames(PCR.Kirsten) <- c(
    "hsapiens_gene_name", "sequence", "time_point",
    "fold_change_Kirsten", "FDR_Kirsten")
PCR.Kirsten <- cbind(
    PCR.Kirsten, 
    logFC_Kirsten = as.matrix(as.numeric(lapply(
        X = PCR.Kirsten[, "fold_change_Kirsten"],
        FUN = function(i) {
            if (i < 0) { -log(x = abs(i), base = 2) }
            else { log(x = abs(i), base = 2) } }))))
PCR.Kirsten$time_point %<>%
    gsub(pattern = "\\+(\\d*) wk", replacement = "\\1", x = .)
PCR.Kirsten <- PCR.Kirsten[, -4]
head(PCR.Kirsten)

# Read in the file containing biomarkers RT-qPCR results from Carol
PCR.Carol <- read.table(
    file = paste(
        user, "/Home_work_sync/Work/TIDA/miRNA_PCR/Diff_expr/RT-qPCR_DE.txt",
        sep = ""),
    sep = "\t", header = TRUE)
PCR.Carol <- PCR.Carol[, c(3:5, 8, 16)]
colnames(PCR.Carol) <- c(
    "gene_name", "sequence", "time_point", "logFC_Carol", "FDR_Carol")
head(PCR.Carol)

# Merge the miRNA-seq, Carol PCR and Kirsten PCR data
Comp <- NULL
for (time in c(1, 2, 6, 10, 12)) {
    data1 <- PCR.Carol[
        PCR.Carol$time_point == time & !is.na(PCR.Carol$sequence),
        c("sequence", "gene_name", "time_point", "logFC_Carol", "FDR_Carol")]
    data2 <- PCR.Kirsten[
        PCR.Kirsten$time_point == time & !is.na(PCR.Kirsten$sequence),
        c("sequence", "hsapiens_gene_name", "time_point",
          "logFC_Kirsten", "FDR_Kirsten")]
    pattern.time <- colnames(DEtable) %>%
        grep(
            pattern = paste("(logFC|FDR)_", time, "w", sep = ""),
            x = ., value = TRUE)
    data3 <- DEtable[, c("feature_id", "Sequence", pattern.time)]
    colnames(data3) <- c(
        "feature_id", "sequence", "logFC_RNAseq", "FDR_RNAseq")
    if (is.null(Comp)) {
        Comp <- merge(x = data1, y = data2, by = "sequence", all = TRUE) %>%
            merge(x = data3, y = ., by = "sequence", all.y = TRUE)
        Comp <- cbind(
            Comp,
            time_point = rep(x = time, times = length(Comp$sequence)))
    }
    else {
        tmp <- merge(x = data1, y = data2, by = "sequence", all = TRUE) %>%
            merge(x = data3, y = ., by = "sequence", all.y = TRUE)
        tmp <- cbind(
            tmp,
            time_point = rep(x = time, times = length(tmp$sequence)))
        Comp <- rbind(Comp, tmp)
    }
}
Comp <- Comp[
    (!is.na(Comp$logFC_Carol) || !is.na(Comp$logFC_Kirsten)) &
        !is.na(Comp$logFC_RNAseq), -c(6, 10)]
head(Comp)

# Create a loop for comparison versus Carol PCR or Kirsten PCR
cor.gene <- list()
for (val in c("Carol", "Kirsten")) {
    title <- paste("vs_PCR_", val, sep = "")
    logfc <- paste("logFC_", val, sep = "")
    fdr <- paste("FDR_", val, sep = "")
    # Perform overall correlation of log fold-change between PCR
    # and RNA-seq data
    cor.gene[paste("Correlation", title, "overall", sep = "_")] <- list(
        cor.test(
            x = Comp$logFC_RNAseq, y = Comp[, logfc],
            alternative = "two.sided", method = "spearman",
            na.action = na.omit()))
    # Perform correlation of log fold-change between PCR and RNA-seq
    # data per gene
    for (g in unique(Comp[!is.na(Comp[, logfc]), "sequence"])) {
        dat1 <- Comp[Comp$sequence == g, logfc]
        dat2 <- Comp[Comp$sequence == g, "logFC_RNAseq"]
        res <- cor.test(
            x = dat1, y = dat2, alternative = "two.sided",
            method = "spearman", na.action = na.omit())
        cor.gene[paste("Correlation", title, g, sep = "_")] <- list(res)
    }
    # Determine concordance of fold-change and P-value between techniques
    concor.val <- 0
    for (x in 1:nrow(Comp)) {
        if (!is.na(Comp[x, logfc])) {
            if ((Comp[x, "FDR_RNAseq"] < 0.05) && (Comp[x, fdr] < 0.05)) {
                if ((
                    Comp[x, logfc] > 0 && Comp[x, "logFC_RNAseq"] > 0) || (
                        Comp[x, logfc] < 0 && Comp[x, "logFC_RNAseq"] < 0)) {
                    concor.val <- concor.val + 1
                }
            }
        }
    }
    cor.gene[paste("Concordance", title, sep = "_")] <- (
        concor.val*100/length(Comp[Comp$FDR_RNAseq < 0.05, "sequence"]))
}

# Create a variable containing the current concordance results
assign(x = paste("cor.", method, sep = ""), value = cor.gene)

# Save the data to keep in R objects
save(
    list = ls(pattern = paste(".*", method, "$", sep = "")),
    file = paste(method, "_data.RData", sep = ""), compress = TRUE)


### Clean up the R workspace ---------------------------------------------

# Remove temporary and unrequired variables
rm(list = setdiff(
    x = ls(),
    y = c("target", "animal", "depth", "dgelist.all", "method", "time.point")))


### END ------------------------------------------------------------------

