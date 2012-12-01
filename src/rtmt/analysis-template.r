#
# Analysis code template for TMT clustering experiments
#
# Copyright (c) 2011 Marc Kirchner 
#

# get libs
#source ("rtmt.r")
require(compositions)
require(mclust)

tmt.analysis <- function(prefix, filename, tmtcols, 
  fileIdentifiers, nClusters=2:9)
{
    ################################################################
    # user parameters, libs, etc
    ################################################################
    #filename <- "BB/H export filename"
    #prefix <- "project-analysis-prefix"
    #tmtcols <- 3:8 # TMT sixplex
    #nClusters <- 2:9
    
    ################################################################
    # load and process data
    ################################################################
    # load data
    x <- tmt.loadBigBangExport(filename)
    # construct the file IDs
    # yields "X1"  "1X1" "X3"  "1X3" "X5"  "1X5"
    #fileIdentifiers <- apply(expand.grid(list(c("_X", "_1X"), c("1_", "3_",
    #  "5_"))), 1, paste, collapse="")
    # loading normalization
    x <- tmt.filePatternWiseLoadingNormalization(x, tmtcols, fileIdentifiers, 
      prefix=prefix, doPlot=TRUE)
    # get rid of all-zero traces
    x <- tmt.filterZeroTraces(x)
    # mark shared peptides (if available)
    x <- tmt.markSharedPeptides(x)
    
    ################################################################
    # go to the protein level
    ################################################################
    # throw away all shared peptides
    x <- x[x$shared == 0,]
    # save the total reporter ion intensities
    a <- apply(x[,tmtcols], 1, sum)
    # enable projection approaches and go compsitional
    x[,tmtcols][x[,tmtcols] == 0] <- NaN;
    xc <- acomp(x[,tmtcols])
    # protein level inference
    f <- x$AccNumber
    p <- tmt.peptide2protein(xc, a, f, k=5)
    # number of peptide per protein
    nPeptides <- tapply(matrix(1, nrow(x), 1), f, sum)
    ## don't plot
    png(file=paste(prefix, "-nPeptides_ecdf.png", sep=""))
    plot(ecdf(nPeptides), lty=1, pch=20, bty="l", col=4, xlab="number of peptides")
    abline(v=3, col=4)
    dev.off()
    # throw out all proteins with less than three peptides
    good <- rownames(nPeptides)[nPeptides > 2]
    p <- p[rownames(p) %in% good,]
    #p <- data.frame(EnsemblID=rownames(p), p, row.names=NULL)
    
    ################################################################
    # clustering
    ################################################################
    mc <- Mclust(idt(acomp(p))[,], G=nClusters)
    png(file=paste(prefix, "-clustering-Aitchison.png", sep=""))
    plot(acomp(p), col=mc$class)
    dev.off()
    ################################################################
    # generate excel output
    ################################################################
    descs <- as.character(x[match(rownames(p), x$AccNumber),"Description"])
    xls <- data.frame(as.character(rownames(p)), as.character(descs),
      nPeptides[rownames(p)], mc$class, as.matrix(p[,]))
    colnames(xls) <- c("AccNum", "Description", "nPeptides", "Cluster",
      paste("TMT", 1:6))
    write.table(file=paste(prefix, "-clustering.xls", sep=""), row.names=FALSE,
      xls)
}
