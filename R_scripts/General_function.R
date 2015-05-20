

### Define functions used accross several scripts ------------------------

# Several functions are used across several pipelines and are therefore
# sourced into the various pippelines



### List of required packages --------------------------------------------

# Create a function to load or install (then load) the required packages
loadpackage <- function(package) {
    if (
        suppressWarnings(require(
            package = deparse(substitute(package)),
            character.only = TRUE, quietly = TRUE))) {
        print(paste(
            deparse(substitute(package)), " is loaded correctly!", sep = ""))
    }
    else {
        print(paste(
            "Trying to install ", deparse(substitute(package)), sep=""))
        source(file = "http://bioconductor.org/biocLite.R", verbose = FALSE)
        biocLite(pkgs = deparse(substitute(package)), suppressUpdates = TRUE)
        if(
            require(
                package = deparse(substitute(package)),
                character.only = TRUE, quietly = TRUE)) {
            print(paste(
                deparse(substitute(package)),
                " is correctly installed and loaded!", sep = ""))
        }
        else {
            stop(paste(
                '"', "Could not install ",
                deparse(substitute(package)), '"', sep = ""))
        }
    }
    print(paste(
        deparse(substitute(package)), " version: ",
        packageVersion(pkg = deparse(substitute(package))), sep = ""))
}



### Old differential expression call -----------------------------------------

# An old version of a function to perform the differential expression
# within edgeR according to provided parameters
diff_expr_edgeR <- function(
    treat1,
    treat2,
    data,
    design,
    group,
    adjpvalue,
    method,
    LRTdata,
    DEdata,
    DEfile,
    Smearfile) {
    print(paste(
        "This function is being replaced! Please use the function",
        " 'multi.DE()' for your edgeR analysis!", sep = ""))
}



### Significance label for graph -----------------------------------------

# Create a function to add significance label
sig_label <- function(
    arg1,
    arg2) {
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



### Old differential expression call for miRNA -------------------------------

# An old version of a function to perform differential expression of
# miRNA genes within edgeR package
miR.DE <- function(
    ...,
    data,
    design,
    group,
    adjpvalue,
    method,
    lrtdata,
    dedata,
    smearfile) {
    print(paste(
        "This function is being replaced! Please use the function",
        " 'multi.DE()' for your edgeR analysis!", sep = ""))
}



### Merge differential expression tables obtained from edgeR -------------

# Create a function to merge the differential expression tables
# obtained in edgeR
DEtable.merge <- function(
    ...,
    output,
    pattern) {
    var <- c(...)
    suffix <- gsub(pattern = pattern, replacement = '_', x = var, perl = TRUE)
    for (i in 1:length(var)) {
        tomerge <- eval(parse(text = var[i]))$table
        if (i == 1) {
            data <- tomerge[,]
        }
        else {
            data <- merge(
                x = data,
                y = tomerge[, (ncol(tomerge) - 4):ncol(tomerge)],
                by = "row.names")
            rownames(data) <- data[, "Row.names"]
            data <- data[, -1]
        }
        colnames(data)[(ncol(data) - 4):ncol(data)] <- c(
            gsub(
                pattern = "$", replacement = suffix[i],
                x = colnames(data)[(ncol(data) - 4):ncol(data)], perl = TRUE))
    }
    data <- cbind(rownames(data), data)
    colnames(data)[1] <- "feature_id"
    assign(x = output, value = data, envir = .GlobalEnv)
}



### Merge differential expression tables obtained from DESeq -------------

# Create a function to merge the differential expression tables
# obtained in DESeq
DESeq.merge <- function(
    ...,
    output,
    pattern) {
    var <- c(...)
    suffix <- gsub(pattern = pattern, replacement = '_', x = var, perl = TRUE)
    for (i in 1:length(var)) {
        tomerge <- as.data.frame(eval(parse(text = var[i])))
        if (i == 1) {
            data <- tomerge[,]
            colnames(data)<- c(gsub(
                pattern = "$", replacement = suffix[i],
                x = colnames(data), perl = TRUE))
        }
        else {
            data <- merge(x = data, y = tomerge, by = "row.names")
            rownames(data) <- data[, "Row.names"]
            data <- data[, -1]
            colnames(data)[(ncol(data) - 5):ncol(data)] <- c(gsub(
                pattern = "$", replacement = suffix[i],
                x = colnames(data)[(ncol(data) - 5):ncol(data)], perl = TRUE))
        }
    }
    data <- cbind(rownames(data), data)
    colnames(data)[1] <- "feature_id"
    assign(x = output, value = data, envir = .GlobalEnv)
}



### Multiple Differential expression call with edgeR ---------------------

# Create a function to perform differential expression of genes using edgeR
multi.DE <- function(
    ...,
    data,
    design,
    group,
    adjpvalue,
    method,
    lrtdata = NULL,
    dedata = NULL,
    smearfile = NULL) {
    comparison <- matrix(data = unlist(list(...)), ncol = 2, byrow = TRUE)
    contr <- rep(x = 0, times = length(colnames(design)))
    for (i in 1:length(comparison[, 1])) {
        eval <- grep(
            pattern = paste('^', comparison[i, 1], '$', sep = ""),
            x = colnames(design), perl = TRUE)
        if (length(eval) == 1) {
            contr[eval] <- as.numeric(comparison[i, 2])
        }
        else if (length(eval) == 0 && length(grep(
            pattern = paste('^', comparison[i, 1], '$', sep = ""),
            x = levels(group), perl = TRUE)) == 1) {
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
    print(summary(decideTestsDGE(
        object = lrt, p.value = adjpvalue, adjust.method = "BH")))
    de <- topTags(object = lrt, n = "inf", adjust.method = method)
    print("Names in the edgeR multiple correction test dataframe:")
    print(names(de))
    print("Heading of the edgeR multiple correction test table:")
    print(head(de$table))
    if (!is.null(smearfile)) {
        png(
            filename = paste(smearfile, "png", sep="."),
            width = 1366, height = 768, units = "px")
        plotSmear(
            object = lrt,
            de.tags = (rownames(lrt$table)[as.logical(decideTestsDGE(
                object = lrt, p.value = 0.05, adjust.method = "BH"))]))
        abline(h = c(-0.5, 0.5), col = "blue")
        dev.off()
    }
    if (!is.null(lrtdata)) {
        assign(x = lrtdata, value = lrt, envir = .GlobalEnv)
    }
    if (!is.null(dedata)) {
        assign(x = dedata, value = de, envir = .GlobalEnv)
    }
}



### Function to collect the legend from a ggplot -------------------------

# Function to collect ggplot legend obtained from website: 
# 'http://stackoverflow.com/questions/11883844/'
# 'inserting-a-table-under-the-legend-in-a-ggplot2-histogram'
g.legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}



### Function to create MDS plot from data subset -------------------------

# Define the function to plot MDS
multi.MDS <- function(
    pattern = NULL,
    data,
    target,
    prefix = NULL,
    suffix = NULL,
    plotmds = list(),
    aes.colour,
    size = 3,
    manual.colour = NULL,
    fill = "grey",
    legend = list(),
    combine = FALSE,
    legend.pos = "right",
    breaks = NULL,
    tiff.picture = list()
    ) {
    .e <- environment()
    legend.name <- gsub(
        pattern = "_", replacement = " ", x = aes.colour, perl = TRUE) %>%
        gsub(pattern = "^(.)", replacement = "\\U\\1", x = ., perl = TRUE)
    myplot <- list(NULL)
    if (is.null(pattern)) {
        pattern <- ".*"
    }
    if (is.null(prefix)) {
        prefix = ""
    }
    else {
        prefix <- gsub(
            pattern = "$", replacement = "_", x = prefix, perl = TRUE)
    }
    if (is.null(suffix)) {
        if (pattern == ".*") {
            suffix <- ""
        }
        else {
            stop(paste(
                "Error: A suffix for name output is required",
                "when using a pattern!",
                sep = " "))
        }
    }
    else {
        suffix.origin <- suffix
        for (j in 1:length(suffix)) {
            suffix[j] <- gsub(
                pattern = "^", replacement = "_", x = suffix[j], perl = TRUE)
        }
    }
    if (length(suffix) != length(pattern)) {
        stop(paste(
            "Error: Number of values in suffix does not match number",
            "of values in pattern!",
            sep = " "))
    }
    if (is.null(row.names(target))) {
        rownames(target) <- target[, 1]
    }
    if (is.null(manual.colour)) {
        manual.colour <- brewer.pal(
            n = length(levels(factor(target[, aes.colour]))),
            name = "Set1")
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
            tiff.picture <- list(
                width = 5*col, height = 5*row, units = "in",
                res = 600, compression = "lzw")
        }
    }
    if (combine == TRUE) {
        if (is.null(pattern)) {
            print(paste(
                "Note: No plots combining is possible without",
                "pattern specified!",
                sep = " "))
            combine <- FALSE
        }
        if (!is.null(suffix)) {
            print(paste(
                "Note: No plots combining is possible without",
                "pattern specified!",
                sep = " "))
            if (length(suffix) > 1) {
                suffix.combine <- paste(suffix.origin, collapse = "-") %>%
                    gsub(pattern = "^", replacement = "_", x = ., perl = TRUE)
            }
            else {
                suffix.combine <- suffix.origin %>%
                    gsub(pattern = "^", replacement = "_", x = ., perl = TRUE)
            }
        }
        else {
            suffix.combine <- ""
        }
        if (!is.null(legend.pos)) {
            print(paste(
                "Note: the argument 'legend.pos' will be overwritten",
                "because plot combining is TRUE!",
                sep = " "))
            legend.pos <- "none"
        }
        full.pat <- pattern %>% paste(., collapse = "|")
        leg.dat <- target[grep(pattern = full.pat, x = rownames(target)),] %>%
            cbind(
                .,
                as.matrix(runif(length(rownames(.)))),
                as.matrix(runif(length(rownames(.))))) %>%
            data.frame(.)
        colnames(leg.dat)[(ncol(leg.dat)-1):ncol(leg.dat)] <- c("x", "y")
        plot.leg <- ggplot(
            data = leg.dat, aes(x = x, y = y), environment = .e) +
            geom_point(size = size, aes(colour = leg.dat[, aes.colour])) +
            theme(
                legend.position = "right",
                legend.title = element_text(size = size*4,  face = "bold"),
                legend.text = element_text(size = size*3.5, face = "bold")) +
            scale_colour_manual(
                name = legend.name, breaks = breaks, limits  = breaks,
                values = manual.colour)
        plot.leg <- g.legend(plot.leg)
        x.pos <- ((row*col)-(length(pattern)))
        plot.leg$vp$x <- unit(x = (1-(1/((col*2)/x.pos))), units = "npc")
        plot.leg$vp$y <- unit(x = (1/(row*2)), units = "npc")
    }
    else {
        suffix.origin <- ""
    }
    for (i in 1:length(pattern)) {
        MDS <- do.call(
            what = "plotMDS",
            args = c(list(data[, grep(
                pattern = pattern[i], x = colnames(data), perl = TRUE)]),
                plotmds)) %>%
            .[c("x", "y")] %>%
            data.frame(target[grep(
                pattern = pattern[i], x = rownames(target), perl = TRUE),
                ], .) %T>%
            write.table(
                file = paste(prefix, "MDS_xy", suffix[i], ".txt", sep = ""),
                sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
        if (is.null(MDS)) {
            stop("Error: The MDS was not generated by function plotMDS!")
        }
        MDS.ggplot <- ggplot(data = MDS, aes(x = x, y = y), environment = .e) +
            geom_point(size = size, aes(colour = MDS[, aes.colour])) +
            theme(
                legend.position = legend.pos,
                panel.background = element_rect(fill = fill),
                title = element_text(face = "bold", size = (size*4)),
                text = element_text(face = "bold", size = (size*3.5)),
                plot.title = element_text(face = "bold", size = (size*5))) +
            ggtitle(suffix.origin[i]) +
            xlab(paste("Dimension ", plotmds$dim.plot[1], sep = "")) +
            ylab(paste("Dimension ", plotmds$dim.plot[2], sep = "")) + 
            scale_colour_manual(
                name = legend.name, breaks = breaks,
                limits = breaks, values = manual.colour)
        myplot[i] <- list(ggplotGrob(MDS.ggplot))
        if (combine == FALSE) {
            do.call(
                what = "tiff",
                args = c(
                    list(
                        filename = paste(
                            prefix, "MDS", suffix[i], ".tif", sep = ""),
                        width = ((range(MDS$x)[2]-range(MDS$x)[1])*11),
                        height = ((range(MDS$y)[2]-range(MDS$y)[1])*10)),
                    tiff.picture))
            print(MDS.ggplot)
            dev.off()
        }
    }
    if (length(myplot) != length(pattern)) {
        stop("Error: Number of plots does not match the number of pattern!")
    }
    if (combine == TRUE) {
        do.call(
            what = "tiff",
            args = c(
                list(
                    filename = paste(
                        prefix, "MDS", suffix.combine, ".tif", sep = "")),
                tiff.picture))
        grid.newpage()
        do.call(
            what = "grid.arrange",
            args = c(myplot, list(nrow = row, ncol = col)))
        grid.draw(plot.leg)
        dev.off()
    }
    print("It looks like a successfull run!")
}



### Function to plot the number of DE ------------------------------------

# Define the dataframe that will contain the number of genes info
plot.numb.DE <- function(
    data,
    comparison = NULL,
    pattern = NULL,
    replace = NULL,
    filename) {
    if (is.null(comparison)) {
        Stop("Error: Comparison values are required for data selection!")
    }
    if (length(grep(pattern = "feature_id", x = colnames(data))) != 1) {
        stop("Error: The ID column in data must be named 'feature_id'!")
    }
    number.DE <- matrix(nrow = (length(comparison) * 8), ncol = 6)
    number.DE <- data.frame(number.DE)
    colnames(number.DE) <- c(
        "time", "expression", "FC_range", "FC_number",
        "FDR_range", "FDR_number")
    if (!is.null(pattern) && !is.null(replace) &&
            length(pattern == length(replace))) {
        comp.val <- gsub(
            pattern = pattern, replacement = replace,
            x = comparison, perl = TRUE)
    }
    else if (is.null(pattern) && is.null(replace)) {
        comp.val <- comparison
    }
    else {
        stop("Error: A single pattern and replace values are to be provided!")
    }
    number.DE[, "time"] <- rep(x = comp.val, each = 8)
    number.DE$time <- factor(number.DE$time, levels = comp.val)
    number.DE[, "expression"] <- rep(
        x = c("Upregulated", "Downregulated"),
        each = 4, times = length(comp.val))
    # Define the matrix containing the range of logFC and Pvalue
    range.val <- matrix(
        data = c(seq(
            from = 0, to = 6, by = 2), 2, 4, 6, Inf,
            0.05, 0.01, 0.001, 0.0001, 0.01, 0.001, 0.0001, 0),
        nrow = 4, ncol = 4, byrow = FALSE)
    colnames(range.val) <- c("FCmin", "FCmax", "FDRmin", "FDRmax")
    rownames(range.val) <- c("low", "mid_low", "mid_high", "high")
    # Calculate the number of DE gene in each logFC and Pvalue range
    counter <- 1
    for (time.val in comparison) {
        fc.var <- paste("logFC_", time.val, sep = "")
        fdr.var <- paste("FDR_", time.val, sep = "")
        for (x in rownames(range.val)) {
            FCmin <- range.val[x, "FCmin"]
            FCmax <- range.val[x, "FCmax"]
            FCmin.invert <- (-1 * FCmin)
            FCmax.invert <- (-1 * FCmax)
            FDRmin <- range.val[x, "FDRmin"]
            FDRmax <- range.val[x, "FDRmax"]
            number.DE[counter, "FC_range"] <- paste(
                FCmin, "to", FCmax, "- Upregulated", sep = " ")
            number.DE[(counter + 4), "FC_range"] <- paste(
                FCmin, "to", FCmax, "- Downregulated", sep = " ")
            number.DE[counter, "FC_number"] <- length(
                data[(data[, fc.var] >= FCmin & data[
                    , fc.var] < FCmax & data[, fdr.var] <= 0.05), "feature_id"])
            number.DE[(counter + 4), "FC_number"] <- length(
                data[(data[, fc.var] <= FCmin.invert & data[
                    , fc.var] > FCmax.invert & data[
                        , fdr.var] <= 0.05), "feature_id"])
            number.DE[counter, "FDR_range"] <- paste(
                FDRmin, "to", FDRmax, "- Upregulated", sep = " ")
            number.DE[(counter + 4), "FDR_range"] <- paste(
                FDRmin, "to", FDRmax, "- Downregulated", sep = " ")
            number.DE[counter, "FDR_number"] <- length(
                data[(data[, fdr.var] <= FDRmin & data[
                    , fdr.var] > FDRmax & data[, fc.var] > 0), "feature_id"])
            number.DE[(counter + 4), "FDR_number"] <- length(
                data[(data[, fdr.var] <= FDRmin & data[
                    , fdr.var] > FDRmax & data[, fc.var] < 0), "feature_id"])
            counter <- counter + 1
        }
        counter <- counter + 4
    }
    # Define the ranges as factor, define a set of colours and the max y limit
    FC.breaks <- number.DE$FC_range <- factor(
        x = number.DE$FC_range, levels = unique(number.DE$FC_range))
    FDR.breaks <- number.DE$FDR_range <- factor(
        x = number.DE$FDR_range, levels = unique(number.DE$FDR_range))
    colours <- c(
        "indianred1", "indianred2", "indianred3", "indianred4",
        "skyblue1", "skyblue2", "skyblue3", "skyblue4")
    max <- (((max(aggregate(
        FC_number ~ time + expression,
        data = number.DE,
        FUN = sum)$FC_number) %/% 10) + 1) * 10)
    # Plot the number of DE gene per range category and time points
    FC.plot <- ggplot(number.DE, aes(x = time)) +
        geom_bar(
            stat = "identity", subset = .(expression == "Upregulated"),
            aes(y = FC_number, fill = FC_range)) +
        geom_bar(
            stat = "identity", subset = .(expression == "Downregulated"),
            aes(y = (FC_number*(-1)), fill = FC_range)) +
        xlab(label = "Comparisons") +
        ylab(label = "Number of DE genes") +
        coord_cartesian(ylim = c(-max, max)) +
        geom_abline(intercept = 0, slope = 0) +
        scale_fill_manual(
            name = "Range of logFC", limits = levels(FC.breaks),
            values = colours) +
        theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
    FDR.plot <- ggplot(number.DE, aes(x = time)) +
        geom_bar(stat = "identity", subset = .(expression == "Upregulated"),
                 aes(y = FDR_number, fill = FDR_range)) +
        geom_bar(stat = "identity", subset = .(expression == "Downregulated"),
                 aes(y = (FDR_number*(-1)), fill = FDR_range)) +
        xlab(label = "Comparisons") +
        ylab(label = "Number of DE genes") +
        coord_cartesian(ylim = c(-max, max)) +
        geom_abline(intercept = 0, slope = 0) +
        scale_fill_manual(
            name = "Range of FDR", limits = levels(FDR.breaks),
            values = colours) +
        theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0))
    # Output a figure of the number of DE per time point
    tiff(
        filename = filename, width = 7, height = 7,
        units = "in", res = 600, compression = "lzw")
    grid.newpage()
    grid.arrange(FC.plot, FDR.plot, nrow = 2, ncol = 1)
    dev.off()
}



### Function to plot Venn diagram of differentially expressed gen --------

# Define the function to plot Venn diagram
venn.de <- function(
    data,
    comparison = NULL,
    pattern = NULL,
    replace = NULL,
    picname = NULL,
    overlapname = NULL,
    ...) {
    if (is.null(comparison)) {
        Stop("Error: Comparison values are required for data selection!")
    }
    else if (length(comparison) > 5) {
        Stop("Error: This function cannot handle more than 5 sets comparison!")
    }
    if (is.null(picname)) {
        Stop("Error: Picname value is required to output the venn plot!")
    }
    if (!is.null(pattern) && !is.null(replace) &&
            length(pattern == length(replace))) {
        comp.val <- gsub(
            pattern = pattern, replacement = replace,
            x = comparison, perl = TRUE)
    }
    else if (is.null(pattern) && is.null(replace)) {
        comp.val <- comparison
    }
    else {
        stop("Error: A single pattern and replace values are to be provided!")
    }
    if (length(grep(pattern = "feature_id", x = colnames(data))) != 1) {
        stop("Error: The ID column in data must be named 'feature_id'!")
    }
    # Identify as a vector list the significant DE genes per time point
    venn.list <- list()
    for (i in 1:length(comparison)) {
        fdr.var <- paste("FDR_", comparison[i], sep = "")
        venn.list[comp.val[i]] <- list(
            as.character(data[!is.na(data[, fdr.var]) & (
                data[, fdr.var] < 0.05), "feature_id"]))
    }
    # Define a set of colors
    colours <- brewer.pal(n = length(comparison), name = "Set1")
    # Create the Venn diagram of overlapping DE genes between time points
    venn.diagram(
        x = venn.list, filename = paste(picname, ".tiff", sep = ""),
        na = "remove", res = 600, fill = colours, cat.col = colours, ...)
    if (!is.null(overlapname)) {
        # Identify and output the DE common to all comparisons
        overlap <- assign(
            x = overlapname, value = Reduce(intersect, venn.list),
            envir=.GlobalEnv)
        write.matrix(
            x = data[overlap,], file = paste(overlapname, ".txt", sep = ""),
            sep = "\t")
    }
}



### END ------------------------------------------------------------------

