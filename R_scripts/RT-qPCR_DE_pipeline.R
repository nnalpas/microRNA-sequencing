##########################################################
# Differential expression analysis of miRNA RT-qPCR data #
##########################################################

# Analysis of 10 animals infected over a 15 weeks time course for which 
# miRNA-seq libraries where prepared from PBMC samples and following which 
# RT-qPCR technical validation was performed

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
loadpackage(package = ggplot2)
loadpackage(package = grid)
loadpackage(package = VennDiagram)
loadpackage(package = RColorBrewer)
loadpackage(package = magrittr)
loadpackage(package = gdata)
loadpackage(package = psych)

################################
# Read in input files within R #
################################

# Read in the excel RT-qPCR expression file
PCR <- read.xls(xls = "cnrq_mean.xlsx", sheet = 1, row.names = 1,
                header = TRUE, na.strings = "NaN")

# Clean up the column name in the qPCR dataset
colnames(PCR) %<>% gsub(
  pattern = "(^X|^miR\\.|^miR)(\\d)", replacement = "miR\\2", x = .,
  perl = TRUE) %>% gsub(pattern = "\\.SE.*", replacement = "_SE", x = ., 
                        perl = TRUE) %>% gsub(
                          pattern = "\\.CNRQ.*", replacement = "_CNRQ", x = ., 
                          perl = TRUE) %>% gsub(
                            pattern = "^(hs|bt)a\\.", replacement = "", x = .,
                            perl = TRUE) %>% gsub(
                              pattern = "\\.", replacement = "_", x = .,
                              perl = TRUE)
head(PCR)

# Read in and format the excel sample code name file
sample <- read.xls(xls = "Sample_code.xlsx", sheet = 1, header = TRUE)
colnames(sample) <- c("animal", "time_point", "sample_id")
sample$time_point <- gsub(pattern = "-", replacement = "pre",
                          x = sample$time_point)

# Merge the sample information with the PCR data
PCR.data <- as.data.frame(merge(x = sample, y = PCR, by.x = "sample_id",
                                by.y = "row.names", all.x = TRUE))
rownames(PCR.data) <- paste(PCR.data$animal, PCR.data$time_point, sep = "_")

###############################################
# Calculate the log fold-change in expression #
###############################################

# Log 2 transform the CNRQ value
PCR.data <- data.frame(PCR.data[, 1:(ncol(sample))], log(
  x = PCR.data[, grep(pattern = "_CNRQ", x = colnames(PCR.data))],
  base = 2), row.names = row.names(PCR.data))
colnames(PCR.data) <- gsub(pattern = "_(CNRQ)", replacement = "_log\\1",
                           x = colnames(PCR.data))
head(PCR.data)

# Define the variables required to compute fold-change in expression
targets.animal <- unique(PCR.data$animal)
targets.time <- unique(PCR.data$time_point)
animal.rownames <- lapply(X = unique(PCR.data$animal), FUN = function(x) rep(
  x = x, length(unique(PCR.data$time_point)))) %>% unlist()
time.rownames <- rep(targets.time, length(targets.animal))
full.rownames <- paste(animal.rownames, time.rownames, sep = "_")
logFC <- data.frame(row.names = full.rownames)
logFC <- cbind(logFC, matrix(data = unlist(x = strsplit(
  x = row.names(logFC), split = "_", fixed = TRUE)), nrow = length(rownames(
    logFC)), ncol = 2, byrow = TRUE))
colnames(logFC) <- c("animal", "time_point")

# Compute the log fold-change in expression versus pre1 week
for (gene in colnames(PCR.data)[(ncol(sample)+1):ncol(PCR.data)]){
  col.data <- c()
  for (animal in targets.animal){
    for (t.target in targets.time){
      data.ref <- PCR.data[PCR.data$time_point == 'pre1' &
                             PCR.data$animal == animal, gene]
      data.target <- PCR.data[PCR.data$time_point == t.target &
                                PCR.data$animal == animal, gene]
      col.data <- c(col.data, data.target-data.ref)
    }
  }
  logFC[, gsub(pattern = "CNRQ", replacement = "FC", x = gene)] <- col.data
}
head(logFC)

########################################################
# Assess normal distribution of PCR data for each gene #
########################################################

# Use the Shapiro-Wilk test on overall data
shapiro <- apply(X = logFC[, ncol(sample):ncol(logFC)],
                 MARGIN = 2, FUN = function(x) shapiro.test(x = x))

# Plot the Q-Q plots for the overall data
qqplot <- apply(X = logFC[, ncol(sample):ncol(logFC)], MARGIN = 2,
      FUN = function(x) qqnorm(y = x, main = colnames(x)))

############################################################
# Compute significance values of fold-change in expression #
############################################################

# Prepare dataframe to include all fold-change and significance evaluation
gene.list <- colnames(logFC)[ncol(sample):ncol(logFC)] %>% 
lapply(X = ., FUN = function(x) rep(
  x = x, length(targets.time))) %>% unlist()
full.rownames <- paste(gene.list, time.rownames, sep = ".") %>% gsub(
  pattern = "_logFC.", replacement = ".", x = .)
full.colnames <- c("gene", "time_point", "shapiro", "n", "meanlogFC", "median",
                   "sd", "se", "min", "max", "t.Pvalue", "w.Pvalue",
                   "final.Pvalue")
sig <- matrix(nrow = length(full.rownames), ncol = length(full.colnames))
sig <- data.frame(x = sig, row.names = full.rownames)
colnames(sig) <- full.colnames
sig[, c("gene", "time_point")] <- matrix(data = unlist(x = strsplit(
  x = row.names(sig), split = ".", fixed = TRUE)), nrow = length(rownames(
    sig)), ncol = 2, byrow = TRUE)

# Compute the descriptive statistics and Pvalues of the log Fold-changes
shapiro.time <- list()
t.value <- list()
w.value <- list()
for (gene in gene.list){
  for (t.target in targets.time){
    x.target <- logFC[logFC$time_point == t.target  & 
                        logFC$animal %in% targets.animal, gene]
    ref <- logFC[logFC$time_point == "pre1"  & 
                   logFC$animal %in% targets.animal, gene]
    stat.value <- describe(x.target)
    if (t.target == "pre1") {
      shap.eval$p.value <- list("NaN")
    }
    else {
      shap.eval <- shapiro.test(x = x.target)
    }
    shapiro.time[t.target] <- list(shap.eval)
    t.eval <- t.test(x = x.target, y = ref, alternative = "two.sided",
                     paired = TRUE, conf.level = 0.95, na.action = omit())
    t.value[t.target] <- list(t.eval)
    w.eval <- wilcox.test(x = x.target, y = ref,alternative = "two.sided",
                          paired=TRUE, na.action = omit())
    w.value[t.target] <- list(w.eval)
    if (shap.eval == "NaN") {
      final.pvalue <- "NaN"      
    }
    else if (shap.eval < 0.1) {
      final.pvalue <- w.eval$p.value
    }
    else {
      final.pvalue <- t.eval$p.value
    }
    gene.val <- paste(gene, t.target, sep = ".") %>% gsub(
      pattern = "_logFC.", replacement = ".", x = .)
    sig[gene.val, c("shapiro", "n", "meanlogFC", "median", "sd", "se", "min",
                    "max", "t.Pvalue", "w.Pvalue", "final.Pvalue")] <- c(
                      shap.eval$p.value, stat.value$n, stat.value$mean,
                      stat.value$median, stat.value$sd, stat.value$se,
                      stat.value$min, stat.value$max, t.eval$p.value,
                      w.eval$p.value, final.pvalue)
  }
}
sig

##################################
# Plot fold-change in expression #
##################################

# Add significance label
sig <- sig_label(arg1 = sig, arg2 = "final.Pvalue")
sig <- cbind(sig, Methods = rep(x = "RT-qPCR", times = length(rownames(sig))))
sig$time_point <- as.numeric(gsub(pattern = "pre", replacement = "-",
                                  x = sig$time_point))
sig$meanlogFC <- as.numeric(sig$meanlogFC)
head(sig)

# Plot the expression data
for (i in unique(sig$gene)) {
  file <- paste(i, "png", sep=".")
  dat <- sig[sig$gene == i, c("time_point", "meanlogFC", "se",
                              "Significance_label", "Methods")]
  plot1 <- ggplot(data = dat, aes(x = dat$time_point, y = dat$meanlogFC,
                                  colour = dat$Methods)) + 
    geom_point(size = 15) + geom_line(size = 5) + geom_text(
      aes(x = dat$time_point, y = (dat$meanlogFC + dat$se),
          label = dat$Significance_label), size = 22) + geom_errorbar(aes(
        x = dat$time_point, ymin = (dat$meanlogFC - dat$se),
        ymax = (dat$meanlogFC + dat$se)),
        width = 1, size = 2.5) + theme(
          panel.background = element_rect(fill = 'wheat'), 
          title = element_text(face = "bold", size = 45),
          text = element_text(face = "bold", size = 35),
          plot.title = element_text(face = "bold", size = 55),
          legend.position = "right") + ggtitle(label = i) +
    xlab("Time point (weeks)") + ylab("log2 fold-change") +
    scale_colour_discrete(name = "Methods")
  png(filename = file, width = 1366, height = 1366, units = "px")
  print(plot1)
  dev.off()
}

# Read in the excel file containing miRNA gene information
gene.info <- read.xls(xls = "miRNA_RTqPCR_gene.xlsx", sheet = 1,
                header = TRUE, na.strings = "NaN")

# Include the gene information with the differential expression results
sig <- merge(x = gene.info, y = sig, by.x = "Gene_performed", by.y = "gene",
             all.y = TRUE)
sig

# Output the RT-qPCR differential expression results
write.table(x = sig, file = "RT-qPCR_DE.txt", quote = FALSE, sep = "\t",
            na = "NaN", row.names = FALSE, col.names = TRUE)

############################
# Clean up the R workspace #
############################

# Remove temporary and unrequired variables
rm(animal, col.data, dat, stat.value, animal.rownames, data.ref, data.target,
   file, final.pvalue, full.colnames, full.rownames, gene, gene.list, gene.val,
   i, plot1, ref, shap.eval, t.eval, t.target, targets.animal, targets.time,
   time.rownames, user, w.eval, x.target, DESeq.merge, DEtable.merge,
   diff_expr_edgeR, g.legend, loadpackage, miR.DE, multi.DE, multi.MDS,
   sig_label, gene.info)

#######
# END #
#######