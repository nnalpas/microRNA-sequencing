#################################################
# Define functions used accross several scripts #
#################################################

#############################
# List of required packages #
#############################

# Create a function to load or install (then load) the required packages
loadpackage <- function(package) {
  if (suppressWarnings(require(package = deparse(substitute(package)),
                               character.only = TRUE, quietly = TRUE))) {
    print(paste(deparse(substitute(package)), " is loaded correctly!",
                sep = ""))
  }
  else {
    print(paste("Trying to install ", deparse(substitute(package)), sep=""))
    source(file = "http://bioconductor.org/biocLite.R", verbose = FALSE)
    biocLite(pkgs = deparse(substitute(package)), suppressUpdates = TRUE)
    if(require(package = deparse(substitute(package)), character.only = TRUE,
               quietly = TRUE)) {
      print(paste(deparse(substitute(package)),
                  " is correctly installed and loaded!", sep = ""))
    }
    else {
      stop(paste('"', "Could not install ", deparse(substitute(package)), '"',
                 sep = ""))
    }
  }
  print(paste(deparse(substitute(package)), " version: ",
              packageVersion(pkg = deparse(substitute(package))), sep = ""))
}

################################
# Differential expression call #
################################

# Create a function to perform the differential expression within edgeR according to provided parameters
diff_expr_edgeR <- function(treat1, treat2, data, design, group, adjpvalue, method, LRTdata, DEdata, DEfile, Smearfile) {
  print("This function is being replaced! Please use the function 'multi.DE()' for your edgeR analysis!")
#  if(length(grep(pattern=treat1, x=colnames(design)))==1 && length(grep(pattern=treat2, x=colnames(design)))==1) {
#    contr <- rep(x=0, times=length(colnames(design)))
#    contr[c(grep(pattern=treat1, x=colnames(design)), grep(pattern=treat2, x=colnames(design)))] <- c(1, -1)
#    lrt <- glmLRT(glmfit=data, contrast=contr)
#  }
#  else if(length(grep(pattern=treat1, x=colnames(design)))==0 && length(grep(pattern=treat1, x=levels(group)))==1) {
#    contr <- grep(pattern=treat2, x=colnames(design))
#    lrt <- glmLRT(glmfit=data, coef=contr)
#  }
#  else if(length(grep(pattern=treat2, x=colnames(design)))==0 && length(grep(pattern=treat2, x=levels(group)))==1) {
#    contr <- grep(pattern=treat1, x=colnames(design))
#    lrt <- glmLRT(glmfit=data, coef=contr)
#  }
#  else {
#    stop("Error: Check that the treatments provided are in group table!")
#  }
#  de <- topTags(object=lrt, n="inf", adjust.method=method)
#  print("Names of the edgeR likeli-hood ratio test dataframe:")
#  print(names(lrt))
#  print("Comparison perfomed in the edgeR likeli-hood ratio test:")
#  print(lrt$comparison)
#  print("Heading of the edgeR likeli-hood ratio test dataframe:")
#  print(head(lrt$table))
#  print("Summary of the number of edgeR DEG:")
#  print(summary(decideTestsDGE(lrt, p.value=adjpvalue)))
#  print("Names of the edgeR multiple correction test dataframe:")
#  print(names(de))
#  print("Heading of the edgeR multiple correction test dataframe:")
#  print(head(de$table))
#  if(length(grep(pattern="antisense", x=DEfile))==1) {
#    write.table(x=de$table[,c("sense_ensembl_gene_id", "external_gene_id","description","logFC", "logCPM", "LR", "PValue", "FDR")], file=paste(DEfile, "txt", sep="."), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
#  }
#  else if(length(grep(pattern="novel", x=DEfile))==1) {
#    write.table(x=de$table[,c("Hsapiens_ensembl_gene_id", "external_gene_id","description","logFC", "logCPM", "LR", "PValue", "FDR")], file=paste(DEfile, "txt", sep="."), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
#  }
#  else {
#    write.table(x=de$table[,c("external_gene_id","description","logFC", "logCPM", "LR", "PValue", "FDR")], file=paste(DEfile, "txt", sep="."), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
#  }
#  png(filename=paste(Smearfile, "png", sep="."), width=1366, height=768, units="px")
#  plotSmear(object=lrt, de.tags=(rownames(lrt$table)[as.logical(decideTestsDGE(lrt, p.value = 0.05))]))
#  abline(h=c(-1, 1), col="blue")
#  dev.off()
#  assign(x=LRTdata, value=lrt, envir=.GlobalEnv)
#  assign(x=DEdata, value=de, envir=.GlobalEnv)
}

################################
# Significance label for graph #
################################

# Create a function to add significance label
sig_label <- function(arg1, arg2) {
  Significance_label <- vector()
  for (j in 1:length(arg1[,arg2])) {
    if (is.na(arg1[j,arg2])) {
      Significance_label <- c(Significance_label, "")
    }
    else if (arg1[j,arg2] < 0.001) {
      Significance_label <- c(Significance_label, "***")
    }
    else if (arg1[j,arg2] < 0.01) {
      Significance_label <- c(Significance_label, "**")
    }
    else if (arg1[j,arg2] < 0.05) {
      Significance_label <- c(Significance_label, "*")
    }
    else {
      Significance_label <- c(Significance_label, "")
    }
  }
  arg1 <- cbind(arg1, Significance_label)
  return(arg1)
}

##########################################
# Differential expression call for miRNA #
##########################################

# Create a function to perform differential expression of miRNA genes
miR.DE <- function(..., data, design, group, adjpvalue, method, lrtdata, dedata,
                   smearfile) {
  print("This function is being replaced! Please use the function 'multi.DE()' for your edgeR analysis!")
#  comparison <- matrix(data = unlist(list(...)), ncol = 2, byrow = TRUE)
#  contr <- rep(x = 0, times = length(colnames(design)))
#  for (i in 1:length(comparison[, 1])) {
#    if (length(grep(pattern = paste('^', comparison[i, 1], '$', sep = ""),
#                    x = colnames(design), perl = TRUE)) == 1) {
#      contr[grep(pattern = paste('^', comparison[i, 1], '$', sep = ""),
#                 x = colnames(design), perl = TRUE)] <-
#        as.numeric(comparison[i, 2])
#    }
#    else if (length(grep(pattern = paste('^', comparison[i, 1], '$', sep = ""),
#                         x = colnames(design), perl = TRUE)) == 0 && length(
#                           grep(pattern = paste('^', comparison[i, 1], '$',
#                                                sep = ""), x = levels(group),
#                                perl = TRUE)) == 1) {
#      print(x = "Intersect value used!")
#    }
#    else {
#      stop("Error: Check that the comparison values are in group table!")
#    }
#  }
#  lrt <- glmLRT(glmfit = data, contrast = contr)
#  print("Names in the edgeR likeli-hood ratio test dataframe:")
#  print(names(lrt))
#  print("Comparison perfomed in the edgeR likeli-hood ratio test:")
#  print(lrt$comparison)
#  print(x = contr)
#  print("Heading of the edgeR likeli-hood ratio test table:")
#  print(head(lrt$table))
#  print("Summary of the number of DE gene:")
#  print(summary(decideTestsDGE(object = lrt, p.value = adjpvalue,
#                               adjust.method = "BH")))
#  de <- topTags(object = lrt, n = "inf", adjust.method = method)
#  print("Names in the edgeR multiple correction test dataframe:")
#  print(names(de))
#  print("Heading of the edgeR multiple correction test table:")
#  print(head(de$table))
#  png(filename = paste(smearfile, "png", sep="."), width = 1366, height = 768,
#      units = "px")
#  plotSmear(object = lrt, de.tags = (rownames(lrt$table)[as.logical(
#    decideTestsDGE(object = lrt, p.value = 0.05, adjust.method = "BH"))]))
#  abline(h = c(-0.5, 0.5), col = "blue")
#  dev.off()
#  assign(x = lrtdata, value = lrt, envir = .GlobalEnv)
#  assign(x = dedata, value = de, envir = .GlobalEnv)
}

############################################################
# Merge differential expression tables obtained from edgeR #
############################################################

# Create a function to merge the differential expression tables
# obtained in edgeR
DEtable.merge <- function(..., output, pattern) {
  var <- c(...)
  suffix <- gsub(pattern = pattern, replacement = '_', x = var, perl = TRUE)
  for (i in 1:length(var)) {
    tomerge <- eval(parse(text = var[i]))$table
    if (i == 1) {
      data <- tomerge[,]
      colnames(data)[(ncol(data) - 4):ncol(data)] <- c(gsub(
        pattern = "$", replacement = suffix[i], x = colnames(data)[(
          ncol(data) - 4):ncol(data)], perl = TRUE))
    }
    else {
      data <- merge(x = data, y = tomerge[, (ncol(tomerge) - 4):ncol(tomerge)],
                    by = "row.names")
      rownames(data) <- data[, "Row.names"]
      data <- data[, -1]
      colnames(data)[(ncol(data) - 4):ncol(data)] <- c(gsub(
        pattern = "$", replacement = suffix[i], x = colnames(data)[(
          ncol(data) - 4):ncol(data)], perl = TRUE))
    }
  }
  data <- cbind(rownames(data), data)
  colnames(data)[1] <- "gene_id"
  assign(x = output, value = data, envir = .GlobalEnv)
}

############################################################
# Merge differential expression tables obtained from DESeq #
############################################################

# Create a function to merge the differential expression tables
# obtained in DESeq
DESeq.merge <- function(..., output, pattern) {
  var <- c(...)
  suffix <- gsub(pattern = pattern, replacement = '_', x = var, perl = TRUE)
  for (i in 1:length(var)) {
    tomerge <- as.data.frame(eval(parse(text = var[i])))
    if (i == 1) {
      data <- tomerge[,]
      colnames(data)<- c(gsub(pattern = "$", replacement = suffix[i],
                              x = colnames(data), perl = TRUE))
    }
    else {
      data <- merge(x = data, y = tomerge, by = "row.names")
      rownames(data) <- data[, "Row.names"]
      data <- data[, -1]
      colnames(data)[(ncol(data) - 5):ncol(data)] <- c(gsub(
        pattern = "$", replacement = suffix[i], x = colnames(data)[(
          ncol(data) - 5):ncol(data)], perl = TRUE))
    }
  }
  data <- cbind(rownames(data), data)
  colnames(data)[1] <- "gene_id"
  assign(x = output, value = data, envir = .GlobalEnv)
}

####################################################
# Multiple Differential expression call with edgeR #
####################################################

# Create a function to perform differential expression of genes using edgeR
multi.DE <- function(..., data, design, group, adjpvalue, method, lrtdata, dedata,
                   smearfile) {
  comparison <- matrix(data = unlist(list(...)), ncol = 2, byrow = TRUE)
  contr <- rep(x = 0, times = length(colnames(design)))
  for (i in 1:length(comparison[, 1])) {
    eval <- grep(pattern = paste('^', comparison[i, 1], '$', sep = ""),
                x = colnames(design), perl = TRUE)
    if (length(eval) == 1) {
      contr[eval] <- as.numeric(comparison[i, 2])
    }
    else if (length(eval) == 0 && length(grep(pattern = paste('^', comparison[
      i, 1], '$', sep = ""), x = levels(group), perl = TRUE)) == 1) {
      print(x = "Intersect value used!")
    }
    else {
      stop("Error: Check that the comparison values are in group table!")
    }
  }
  lrt <- glmLRT(glmfit = data, contrast = contr)
  print("Names in the edgeR likeli-hood ratio test dataframe:")
  print(names(lrt))
  print("Comparison perfomed in the edgeR likeli-hood ratio test:")
  print(lrt$comparison)
  print(x = contr)
  print("Heading of the edgeR likeli-hood ratio test table:")
  print(head(lrt$table))
  print("Summary of the number of DE gene:")
  print(summary(decideTestsDGE(object = lrt, p.value = adjpvalue,
                               adjust.method = "BH")))
  de <- topTags(object = lrt, n = "inf", adjust.method = method)
  print("Names in the edgeR multiple correction test dataframe:")
  print(names(de))
  print("Heading of the edgeR multiple correction test table:")
  print(head(de$table))
  png(filename = paste(smearfile, "png", sep="."), width = 1366, height = 768,
      units = "px")
  plotSmear(object = lrt, de.tags = (rownames(lrt$table)[as.logical(
    decideTestsDGE(object = lrt, p.value = 0.05, adjust.method = "BH"))]))
  abline(h = c(-0.5, 0.5), col = "blue")
  dev.off()
  assign(x = lrtdata, value = lrt, envir = .GlobalEnv)
  assign(x = dedata, value = de, envir = .GlobalEnv)
}

################################################
# Function to collect the legend from a ggplot #
################################################

# Function to collect ggplot legend obtained from website: 
# 'http://stackoverflow.com/questions/11883844/'
# 'inserting-a-table-under-the-legend-in-a-ggplot2-histogram'
g.legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]] 
  return(legend)
}

################################################
# Function to create MDS plot from data subset #
################################################

# Version in development
multi.MDS <- function(pattern = NULL, data, target, prefix = NULL,
                      suffix = NULL, plotmds = list(), aes.colour, size = 3,
                      manual.colour = NULL, fill = "grey", legend = list(),
                      combine = FALSE, legend.pos = "right", breaks = NULL,
                      tiff.picture = list()) {
  .e <- environment()
  legend.name <- gsub(pattern = "_", replacement = " ", x = aes.colour,
                      perl = TRUE) %>% gsub(
                        pattern = "^(.)", replacement = "\\U\\1", x = .,
                        perl = TRUE)
  myplot <- list(NULL)
  if (is.null(pattern)) {
    pattern <- ".*"
  }
  if (is.null(prefix)) {
    prefix = ""
  }
  else {
    prefix <- gsub(pattern = "$", replacement = "_", x = prefix, perl = TRUE)
  }
  if (is.null(suffix)) {
    if (pattern == ".*") {
      suffix <- ""
    }
    else {
      stop("Error: A suffix for name output is required when using a pattern!")
    }
  }
  else {
    suffix.origin <- suffix
    for (j in 1:length(suffix)) {
      suffix[j] <- gsub(pattern = "^", replacement = "_", x = suffix[j],
                        perl = TRUE)
    }
  }
  if (length(suffix) != length(pattern)) {
    stop(paste("Error: Number of values in suffix does not match number of",
               " values in pattern!", sep = ""))
  }
  if (is.null(row.names(target))) {
    rownames(target) <- target[, 1]
  }
  if (is.null(manual.colour)) {
    manual.colour <- brewer.pal(n = length(levels(factor(
      target[, aes.colour]))), name = "Set1")
  }
  if (is.null(breaks)) {
    breaks <- unique(target[, aes.colour])
  }
  else if (length(breaks) != length(manual.colour)) {
    stop("Length of 'breaks' needs to match length of 'manual.colour'!")
  }
  if (length(tiff.picture) == 0) {
    if (combine == FALSE) {
      tiff.picture <- list(units = "in", res = 600, compression = "lzw")
    }
    else {
      number.combine <- (length(pattern) + 1)
      col <- 1
      for (row in 1:number.combine) {
        if ((row*col) >= number.combine) {
          break
        }
        else {
          col <- row
          if ((row*col) >= number.combine) {
            break
          }
        }
      }
      tiff.picture <- list(width = 5*col, height = 5*row, units = "in",
                           res = 600, compression = "lzw")
    }
  }
  if (combine == TRUE) {
    if (is.null(pattern)) {
      print("Note: No plots combining is possible without pattern specified!")
      combine <- FALSE
    }
    if (!is.null(suffix)) {
      print("Note: No plots combining is possible without pattern specified!")
      if (length(suffix) > 1) {
        suffix.combine <- paste(suffix.origin, collapse = "-") %>% gsub(
          pattern = "^", replacement = "_", x = ., perl = TRUE)
      }
      else {
        suffix.combine <- suffix.origin %>% gsub(
          pattern = "^", replacement = "_", x = ., perl = TRUE)
      }
    }
    else {
      suffix.combine <- ""
    }
    if (!is.null(legend.pos)) {
      print(paste("Note: the argument 'legend.pos' will be overwritten",
                  " because plot combining is TRUE!", sep = ""))
      legend.pos <- "none"
    }
    full.pat <- pattern %>% paste(., collapse = "|")
    leg.dat <- target[grep(pattern = full.pat, x = rownames(target)),] %>% 
      cbind(., as.matrix(runif(length(rownames(.)))), as.matrix(
        runif(length(rownames(.))))) %>% data.frame(.)
    colnames(leg.dat)[(ncol(leg.dat)-1):ncol(leg.dat)] <- c("x", "y")
    plot.leg <- ggplot(data = leg.dat, aes(x = x, y = y), environment = .e) +
      geom_point(size = size, aes(colour = leg.dat[, aes.colour])) + theme(
        legend.position = "right",
        legend.title = element_text(size = size*4,  face = "bold"),
        legend.text = element_text(size = size*3.5, face = "bold")) +
      scale_colour_manual(name = legend.name, breaks = breaks,
                          limits  = breaks, values = manual.colour)
    plot.leg <- g.legend(plot.leg)
    x.pos <- ((row*col)-(length(pattern)))
    plot.leg$vp$x <- unit(x = (1-(1/((col*2)/x.pos))), units = "npc")
    plot.leg$vp$y <- unit(x = (1/(row*2)), units = "npc")
  }
  else {
    suffix.origin <- ""
  }
  for (i in 1:length(pattern)) {
    MDS <- do.call(what = "plotMDS", args = c(list(data[, grep(
      pattern = pattern[i], x = colnames(data), perl = TRUE)]), plotmds)) %>%
      .[c("x", "y")] %>% data.frame(target[grep(
        pattern = pattern[i], x = rownames(target), perl = TRUE),], .) %T>%
      write.table(file = paste(prefix, "MDS_xy", suffix[i], ".txt", sep = ""),
                  sep = "\t", quote = FALSE, row.names = TRUE,
                  col.names = TRUE)
    if (is.null(MDS)) {
      stop("Error: The MDS was not generated by function plotMDS!")
    }
    MDS.ggplot <- ggplot(data = MDS, aes(x = x, y = y), environment = .e) + 
      geom_point(size = size, aes(colour = MDS[, aes.colour])) + theme(
          legend.position = legend.pos,
          panel.background = element_rect(fill = fill),
          title = element_text(face = "bold", size = (size*4)),
          text = element_text(face = "bold", size = (size*3.5)),
          plot.title = element_text(face = "bold", size = (size*5))) + ggtitle(
              suffix.origin[i]) + xlab(
                paste("Dimension ", plotmds$dim.plot[1], sep = "")) + ylab(
                  paste("Dimension ", plotmds$dim.plot[2], sep = "")) + 
      scale_colour_manual(name = legend.name, breaks = breaks, limits = breaks,
                          values = manual.colour)
    myplot[i] <- list(ggplotGrob(MDS.ggplot))
    if (combine == FALSE) {
      do.call(what = "tiff", args = c(list(filename = paste(
        prefix, "MDS", suffix[i], ".tif", sep = ""),
        width = ((range(MDS$x)[2]-range(MDS$x)[1])*11),
        height = ((range(MDS$y)[2]-range(MDS$y)[1])*10)),
        tiff.picture))
      print(MDS.ggplot)
      dev.off()
    }
  }
  if (length(myplot) != length(pattern)) {
    stop("Error: The number of plots does not match the number of pattern!")
  }
  if (combine == TRUE) {
    do.call(what = "tiff", args = c(list(filename = paste(
      prefix, "MDS", suffix.combine, ".tif", sep = "")), tiff.picture))
    grid.newpage()
    do.call(what = "grid.arrange", args = c(myplot, list(nrow = row,
                                                         ncol = col)))
    grid.draw(plot.leg)
    dev.off()
  }
print("It looks like a successfull run!")
}

#######
# END #
#######