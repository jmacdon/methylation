##' Create a plot for a region showing differential methylation.
##'
##' This function is used to plot a genomic region where there is differential methylation
##' that includes the methylation data, a smoothed line, and any genes in the region. This
##' function is in general not called by the end user. Instead use \code{methByRegion}.
##' @title Create methylation region plot
##' @param bumpsObj The output from bumphunter
##' @param eset Usually a \code{GenomicRatioSet} created by the minfi package
##' @param row Which row of the bumpsObj are we plotting?
##' @param txdb A transcript database object (e.g., TxDb.Hsapiens.UCSC.hg19.knownGene)
##' @param orgpkg A species-level annotation package (e.g., Homo.sapiens)
##' @param samps A data.frame containing sample information. This data.frame must contain two columns; one called Gender, listing the gender of
##' each sex, and one called Category, listing the category for each sample
##' @param dontuse A character vector of Category levels that won't be plotted. If using all levels, this argument is "".
##' @param sexfirst Which sex should be plotted first? In general this is the sex for which the bump being shown was significant.
##' @param use.symbols Should we convert transcript IDs to gene symbols?
##' @param fitobj An MArrayLM object, after fitting the same model as bumphunter. This is the source of the methylation data in the plot
##' @return This function doesn't return anything. It is only called for the side effect of creating a plot.
##' @author James W. MacDonald (\email{jmacdon@@u.washington.edu})
##' @export
makeMethPlot <- function(bumpsObj, eset, row, txdb, orgpkg, samps, dontuse = "", sexfirst = "Male", use.symbols = TRUE){
    if(!"Gender" %in% colnames(samps))
        stop("There must be a 'Gender' column in the samps data.frame!\n", call. = FALSE)
    if(!all(c("Male","Female") %in% samps$Gender))
        stop("There must be both Male and Female in the Gender column of the samps data.frame!\n", call. = FALSE)
    if(is.character(txdb)) txdb <- get(txdb)
    if(is.character(orgpkg)) orgpkg <- get(orgpkg)
    genome <- genome(eset)[1]
    chr <- bumpsObj$table[row,1]
    start <- bumpsObj$table[row,2]
    end <- bumpsObj$table[row,3]
    reg <- GRanges(chr, IRanges(start, end))
    reg <- resize(reg, width(reg) * 10, "center")
    iTrack <- IdeogramTrack(genome = genome, chromosome = chr)
    gTrack <- GenomeAxisTrack()
    grTrack <- GeneRegionTrack(txdb, genome = genome, chromosome = chr, start = start(reg), end = end(reg),
                               showId = TRUE, name = "Transcripts")
    tmp <- grTrack@range
    if(use.symbols && length(tmp) > 0){
        mapper <- select(orgpkg, elementMetadata(tmp)$symbol, "SYMBOL","TXNAME")
        mp <- mapper[,2]
        names(mp) <- mapper[,1]
        elementMetadata(tmp)$symbol <- mp[elementMetadata(tmp)$symbol]
        grTrack@range <- tmp
    }
    dtm <- dtf <- rowData(eset)[rowData(eset) %over% reg,]
    elementMetadata(dtm) <- getM(eset)[rowData(eset) %over% reg, samps$Gender == "Male" & !samps$Category %in% dontuse, drop = FALSE]
    elementMetadata(dtf) <- getM(eset)[rowData(eset) %over% reg, samps$Gender == "Female" & !samps$Category %in% dontuse, drop = FALSE]
    if(sexfirst == "Male"){
        dTrackM <- DataTrack(dtm, name = "Male methylation", 
                             groups = as.character(samps$Category[samps$Gender == "Male" & !samps$Category %in% dontuse]),
                             type = c("a","p"))
        dTrackF <- DataTrack(dtf, name = "Female methylation", 
                             groups = as.character(samps$Category[samps$Gender == "Female" & !samps$Category %in% dontuse]),
                             type = c("a","p"), legend = TRUE)
        plotTracks(list(iTrack, gTrack, grTrack, dTrackM, dTrackF), from = start(reg), to = end(reg), background.title = "darkblue")
    } else {
         dTrackM <- DataTrack(dtm, name = "Male methylation", 
                             groups = as.character(samps$Category[samps$Gender == "Male" & !samps$Category %in% dontuse]),
                             type = c("a","p"), legend = TRUE)
        dTrackF <- DataTrack(dtf, name = "Female methylation", 
                             groups = as.character(samps$Category[samps$Gender == "Female" & !samps$Category %in% dontuse]),
                             type = c("a","p"))
        plotTracks(list(iTrack, gTrack, grTrack, dTrackF, dTrackM), from = start(reg), to = end(reg), background.title = "darkblue")
     }
}

##' Generate sex-stratified dotplot of methylation for a given region.
##'
##' This is an internal function and not intended for direct use. This is intended to create a dotplot
##' stratified by sex, to show the mean methylation for a given region of the genome, presumably because
##' there is differential methylation for at least on sex.
##' @title Sex stratified mean methylation
##' @param samps A data.frame that maps samples to phenotype. There must be both a Category and Gender column.
##' @param bumpavg Mean methylation data, usually from a call to \code{getMeans}
##' @param dontuse Category levels that are not to be used in the dotplot. If there are only two levels, use "".
##' @return Nothing is returned. Only called for the side effect of creating a dotplot.
##' @author James W. MacDonald (\email{jmacdon@@u.washington.edu})
bwplotfun <- function(samps, bumpavg, dontuse = "AGA"){
    for(i in seq_len(ncol(bumpavg))){
        tmp <- data.frame(meth = bumpavg[!samps$Category %in% dontuse,i],
                          cat = factor(samps$Category[!samps$Category %in% dontuse]),
                          gend = factor(samps$Gender[!samps$Category %in% dontuse], labels = c("Female","Male")))
        print(dotplot(meth~cat|gend, tmp, ylab = paste0("Methylation of genomic region ", colnames(bumpavg)[i])))
    }
}



##' A function to create an HTML page with links to plots showing methylation status over each differentially methylated region,
##' a dotplot showing sex-stratified mean methylation, and an HTML table showing correlation between methylation status and gene
##' expression for all genes within 1 Mb of the methylation site. This is the main function for this package, and is likely the only
##' one that an end user should use.
##'
##' This function is intended to create plots corresponding to a region of the genome that is considered
##' to be differentially methylated, usually after running \code{bumphunter} to detect differentially methylated
##' regions. The methylation region plot will consist of the probe-wise methylation, with a smoothed line to indicate the
##' portion that is differentially methylated. The dotplot will show sex-stratified differential methylation (as it is usually
##' safer to do differential methylation for each sex separately). If there are expression data available, there will be links to
##' tables showing the correlation between methylation and gene expression.
##' @title Create plots showing differential methylation and correlation to expression data.
##' @param bmpsObj The output from \code{bumphunter}
##' @param eset Usually a GenomicRatioSet, created by a call to \code{preprocessQuantile} from the minfi package
##' @param samps A data.frame that maps samples to phenotype. This data.frame MUST contain columns named Category and Gender!
##' @param contname A contrast name, used to name the directory where these data will be stored. Usually of the form 'this_vs_that'.
##' @param longname A long form of the contrast name, usually of the form 'This versus that'
##' @param txdb A transcript database (e.g., TxDb.Hsapiens.UCSC.hg19.knownGene)
##' @param gene.data Default is NULL. Otherwise, a data.frame or matrix containing gene expression data. The columns of this
##' data.frame must correspond exactly to columns of the methylation data. Row names must be array IDs or Entrez Gene IDs.
##' @param chip.db The chip-level array data package corresponding to the gene expression data. If NULL, the assumption will
##' be made that the row.names are ENTREZ GENE IDs.
##' @param fitcol Which column of the MArrayLM object corresponds to the coefficient tested by \code{bumphunter}?
##' @param cut The p-value cutoff used to select significant 'bumps'.
##' @param cutcol Which column of the bumpsObj table item to use for defining the p-value cutoff?
##' @param dontuse Which Categories from the samps data.frame should we NOT use? If only two Category levels, use "".
##' @param orgpkg The organism-level annotation package (e.g., Homo.sapiens)
##' @param sexfirst Which sex should be plotted first? Defaults to Male.
##' @param use.symbols Should transcript IDs be converted to gene symbols when plotting gene regions?
##' @param fit An MArrayLM object, created by fitting the same model as used by \code{bumphunter}, but probe-wise using the limma package
##' @return This returns an HTMLReportRef that can be used to create an index.html page.
##' @export An organism level annotation package (e.g., org.Hs.eg.db)
##' @author James W. MacDonald (\email{jmacdon@@u.washington.edu})
methByRegion <- function(bmpsObj, eset, samps, contname, longname, txdb, gene.data = NULL, chip.db = NULL, fitcol, cut = 0.001,
                         cutcol = c("p.value", "fwer","p.valueArea", "fwerArea"), dontuse = "", orgpkg, sexfirst = "Male", use.symbols = TRUE){
    cutcol <- match.arg(cutcol,  c("p.value", "fwer","p.valueArea", "fwerArea"))
    cutcol <- bmpsObj$table[, cutcol]
    tab <- bmpsObj$table[cutcol <= cut,]
    bmpavg <- getMeans(tab[,1:3], eset)
    if(!file.exists("reports")) dir.create("reports")
    if(!file.exists(paste0("reports/", contname))) dir.create(paste0("reports/", contname))
    for(i in seq_len(ncol(bmpavg))){
        png(paste0("reports/", contname, "/", colnames(bmpavg)[i], "methplot.png"))
        makeMethPlot(bumpsObj = bmpsObj, eset = eset, row = i, txdb = txdb, samps = samps, dontuse = dontuse, sexfirst = sexfirst, orgpkg = orgpkg,
                     use.symbols = use.symbols)
        dev.off()
        png(paste0("reports/", contname, "/", colnames(bmpavg)[i], "bwplot.png"))
        bwplotfun(samps, bmpavg[,i,drop = FALSE], dontuse)
        dev.off()
    }
    uri.meth <- sapply(colnames(bmpavg), function(x) paste0(contname, "/", x, "methplot.png"))
    uri.bw <- sapply(colnames(bmpavg), function(x) paste0(contname, "/", x, "bwplot.png"))
    uri.meth <- paste0("<a href=\"", uri.meth, "\">Methylation region plot</a>")
    uri.bw <- paste0("<a href=\"", uri.bw, "\">Methylation dotplot</a>")
    if(!is.null(gene.data)){
        if(is.character(txdb)) txdb <- get(txdb)
        if(is.character(orgpkg)) orgpkg <- get(orgpkg)
        genes <- genes(txdb)
        gbm <- geneByMeth(tab = tab[,1:3], genes = genes, eset = eset, samps = samps,
                          gene.data = gene.data, chip.db = chip.db, contname = contname, dontuse = dontuse, orgpkg = orgpkg)
        uris2 <- gsub("reports//", "", sapply(gbm, function(x) path(x[[2]])))
        uris2 <- paste0("<a href=\"", uris2,"\">",gsub("\\.html","", uris2), "</a>")
        out <- data.frame(Regions = uris2, p.value = tab$p.value,  Gene.regions = uri.meth, Dotplots = uri.bw)
    }else{
        out <- data.frame(Regions = apply(tab, 1, paste, collapse = "_"),  Gene.regions = uri.meth, Dotplots = uri.bw)
    }
    htmlFile <- HTMLReport(contname, longname, "reports/")
    publish(out, htmlFile)
    finish(htmlFile)
    return(htmlFile)
}

##' A function to create tables showing correlation between methylation status and gene expression.
##'
##' 
##' This function is not intended to be called by the end user. Instead it is intended to create the HTML table showing correlation
##' between gene gene expression and methylation status, as well as xyplots that show the correlation graphically.
##' @title Create tables and plots showing correlation between methylation and gene expression data.
##' @param lstitm A data.frame containing gene expression data for genes that are within the 1 Mb region centered on the methylation 'bump'.
##' @param eset Usually a GenomicRatioSet, created by a call to \code{preprocessQuantile} from the mimfi package
##' @param prb A vector of probe IDs that correspond to the differently methylated (Illumina 450K) probe IDs.
##' @param samps A data.frame that maps samples to phenotype. This data.frame MUST contain columns named Category and Gender!
##' @param file A filename. In general this is the genomic region, separated by underscores (e.g., chr1_12345_23456)
##' @param contname A contrast name. Usually this is lowercase and separated by underscores (e.g., this_vs_that)
##' @param orgpkg An organism level annotation package (e.g., Homo.sapiens)
##' @return This returns an HTMLReportRef.
##' @author James W. MacDonald (\email{jmacdon@@u.washington.edu})
plotAndOut <- function(lstitm, eset, prb, samps, file, contname, orgpkg){
    if(nrow(lstitm) == 0) return(NULL)
    tmp <- as.vector(colMeans(getM(eset)[prb,,drop = FALSE]))
    samps.x <- samps
    lstitm <- as.matrix(t(lstitm))
    if(is.character(orgpkg)) orgpkg <- get(orgpkg)
    cn <- select(orgpkg, colnames(lstitm), "SYMBOL", "ENTREZID")
    colnames(lstitm) <- gsub("-", "_", cn[,2])
    naind <- apply(lstitm, 1, is.na)
    if(is.vector(naind)) dim(naind) <- c(1, length(naind))
    samps.x <- cbind(samps.x, tmp, lstitm)
    out <- lapply(seq_len(ncol(lstitm)), function(x) {
        dep <- colnames(samps.x)[x + ncol(samps) + 1]
        mod <- lm(as.formula(paste(dep, "~tmp")), samps.x, subset = Gender == "Male")
        mal <- summary(mod)$coefficients
        fem <- summary(update(mod, subset = Gender == "Female"))$coefficients 
        return(c(fem[2,c(1,4)], mal[2,c(1,4)]))
    })
    out <- do.call("rbind", out)
    colnames(out) <- paste(rep(c("Female","Male"), each = 2), rep(c("methylation effect", "p-value"), 2))
    for(i in colnames(lstitm)){
        png(paste0("reports/",contname, "/", file, "_", i, ".png"))
        print(xyplot(as.formula(paste0(i, "~tmp|Gender")), samps.x, 
                     panel = function(x, y, ...){panel.xyplot(x, y, ...); 
                                                 panel.lmline(x, y, ...)},
                     xlab = paste("Methylation status at", file),
                     ylab = paste("Expression data for gene", i)))
        dev.off()
    }
    uris <- sapply(colnames(lstitm), function(x) paste0(contname, "/",  file, "_", x,".png"))
    uris <- paste0("<a href=\"", uris, "\">", colnames(lstitm), "</a>")
    row.names(out) <-  colnames(lstitm)
    out2 <-data.frame(GeneSymbol = uris, out)
    htmlFile <- HTMLReport(file, paste("Gene expression as a function of methylation at", file), "reports/")
    publish(out2, htmlFile)
    finish(htmlFile)
    return(list(df = out, thelink = htmlFile))
}


##' Get genes within a 1 Mb region centered on a methylation region and fit linear models with methylation status as independent variable
##' and gene expression as dependent variable
##'
##' This function takes a table of differentially methylated regions, selects the genes that are within 1 Mb of a given region, and then fits models
##' that test for a relationship between methylation status and gene expression. In order for this to work correctly the gene expression and methylation data
##' must be in corresponding columns (e.g., sample 1 must be in the first column of the methylation data as well as the expression data). In addition, the row.names of the
##' gene expression data have to be the manufacturer IDs, so we can figure out what genes are being interrogated. Otherwise, the row.names can be ENTREZ GENE IDs.
##' @title Correlate methylation and gene expression data.
##' @param tab A data.frame containing the chromosome, start, stop for a methylation region. Usually extracted from the table list item of a 'bumps' object.
##' @param genes A GRanges object that lists known genes. Usually generated by e.g. genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
##' @param eset Usually a GenomicRatioSet, created by a call to \code{preprocessQuantile} from the mimfi package
##' @param samps A data.frame that maps samples to phenotype. This data.frame MUST contain columns named Category and Gender!
##' @param gene.data Default is NULL. Otherwise, a data.frame or matrix containing gene expression data. The columns of this data.frame must correspond
##' exactly to columns of the methylation data. Row names must be array IDs or Entrez gene IDs.
##' @param chip.db The chip-level array data package corresponding to the gene expression data. If NULL, the assumption will be made that the row.names
##' are ENTREZ GENE IDs.
##' @param contname A contrast name, used to name the directory where these data will be stored. Usually of the form 'this_vs_that'.
##' @param dontuse Which Categories from the samps data.frame should we NOT use? If only two Category levels, use "".
##' @param orgpkg An organism level annotation package (e.g., Homo.sapiens)
##' @return This function returns a list of HTMLReportRef items that can be used to create links.
##' @author James W. MacDonald (\email{jmacdon@@u.washington.edu})
geneByMeth <- function(tab,  genes, eset, samps, gene.data, chip.db, contname, dontuse = "", orgpkg){
    methranges <- GRanges(tab[,1], IRanges(tab[,2], tab[,3]))
    probes <- apply(tab, 1, paste, collapse = "_")
    probes <- gsub("\\s+", "", probes, perl = TRUE)
    probelst <- lapply(1:length(methranges), function(x) names(rowData(eset))[rowData(eset) %over% methranges[x,]])
    methranges <- resize(methranges, 1e6, "center")
    if(!is.null(chip.db)){
        if(is.character(chip.db)) chip.db <- get(chip.db)
        annot <- select(chip.db, row.names(gene.data), "ENTREZID")
    }else{
        annot <- data.frame(ENTREZID = row.names(gene.data))
    }
    genlst <- lapply(seq_len(length(methranges)), function(x) names(genes[subjectHits(findOverlaps(methranges[x,], genes)),]))
    gendatlst <- lapply(genlst, function(x) { gd <- gene.data[annot$ENTREZID %in% x, ]
                                              gd <- gd[,!samps$Category %in% dontuse]
                                              nam <- annot$ENTREZID[annot$ENTREZID %in% x]
                                              gd <- gd[!duplicated(nam),]
                                              nam <- nam[!duplicated(nam)]
                                              row.names(gd) <- nam
                                              return(gd)})
    methdat.out <- lapply(seq_len(length(gendatlst)), function(x) plotAndOut(gendatlst[[x]], eset[,!samps$Category %in% dontuse], probelst[[x]], 
                                                                             samps[!samps$Category %in% dontuse,], 
                                                                             probes[x], contname, orgpkg))
    methdat.out
}

##' Compute means from methylation data
##'
##' This function is an internal function and not intended for use by end users. The purpose is to
##' compute mean methylation data that will be used for creating dotplots
##' @title Get mean methylation data
##' @param bmptab The table item from a bumps object created from running \code{bumphunter} on methylation data
##' @param eset A GenomicRatioSet or MethylSet.
##' @return A data.frame containing the mean expression for all probes within all differentially methylated regions.
##' @author James W. MacDonald (\email{jmacdon@@u.washington.edu})
getMeans <- function(bmptab, eset){
    gr <- GRanges(bmptab$chr, IRanges(start = bmptab$start, end = bmptab$end))
    dat <- sapply(1:nrow(bmptab), function(x) colMeans(getM(eset[rowData(eset) %over% gr[x,],])))
    colnames(dat) <- apply(bmptab, 1, function(x) paste(gsub("\\s+", "", x, perl=TRUE), collapse = "_"))
    dat
}
