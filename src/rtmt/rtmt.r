#
# TMT analysis helper functions.
# 
# Copyright (C) 2011 Marc Kirchner 
#

#
# Load data from a tab-delimited Steen Lab BB/H export
# @param[in] filename The filename of the BB/H export file to load.
#
tmt.loadBigBangExport <- function(filename)
{
    # load data
    x <- read.delim(file=filename, header=T, na.strings="Null")
    # adjust the peptide sequences
    x[,"Sequence"] <- gsub(".[A-Z-]$", "", gsub("^[^.].", "", x[,"Sequence"]))
    x
}

# 
# Loads the combined set output from the Steen lab BB/H pipeline.
# Tailored to TMT 6-plex loading. Can be extended to work with iTRAQ etc.
#
# @param[in] filename The filename of the combined set file to load.
# @return A data frame holding most of the export data (leaving out some of 
#         the BB/H internal variables and cross-ref IDs). Notably, this
#         contains a column of NAs, termed "id"; depending on the analysis
#         the user should replace the column using
#               dat$id <- dat$accession
#         for using e.g. IPI numbers as ID; or
#               dat$id <- dat$geneName
#         in order to use gene names as IDs.
# 
tmt.loadCombinedSet <- function(filename)
{
	dat <- read.table(file=filename, colClasses="character", header=T)
	dat <- data.frame(
			id=rep(NA, nrow(dat)),
			accession=gsub(".[0-9]+$", "", gsub("^IPI:", "", dat[["AccNumber"]])),
			geneName=gsub("Tax_Id=9606 Gene_Symbol=([A-Z0-9]+).*", "\\1", dat[["Description"]]),
			description=gsub("Tax_Id=9606 Gene_Symbol=[^ ]+ ", "", dat[["Description"]]),
			sequence=dat[["Seq"]], 
			modifications=dat[["Var_Mod"]],
			TMT1=as.numeric(dat[["Ch1"]]),
			TMT2=as.numeric(dat[["Ch2"]]),
			TMT3=as.numeric(dat[["Ch3"]]),
			TMT4=as.numeric(dat[["Ch4"]]),
			TMT5=as.numeric(dat[["Ch5"]]),
			TMT6=as.numeric(dat[["Ch6"]]),
			score=as.numeric(dat[["Score"]]),
			rt=as.numeric(dat[["Elution"]]),
			exp_mz=as.numeric(dat[["Exp_MZ"]]),
			delta=as.numeric(dat[["Delta"]]),
			miss=as.numeric(dat[["Miss"]]),
			charge=as.numeric(dat[["Charge"]]),
			sample=dat[["Sample"]],
			filename=dat[["FileName"]]
	)
}

#
# make use of the Filename info to get a run-wise TMT loading correction.
# 'fileIdentifiers' should contain a list of substrings that are used to 
# find groups of filenames in x$FileNames.
#
tmt.filePatternWiseLoadingNormalization <- function(x, tmtcols, fileIdentifiers, prefix="", doPlot=TRUE)
{
    if (doPlot) {
        png(file=paste(prefix, "-loadings.png", sep=""), h=480/2, w=480/2*length(fileIdentifiers))
        par(mfrow=c(1,length(fileIdentifiers)))
    }
    for (i in 1:length(fileIdentifiers)) {
        print(fileIdentifiers[[i]])
        # figure out which filenames are in the experiment
        patternIdxs <- unique(as.numeric(unlist(
          sapply(fileIdentifiers[[i]], grep, as.character(x$FileName)))))
        experimentFiles <- unique(as.character(x$FileName)[patternIdxs])
        # find the peptides in question in the complete dataset
        idx <- as.character(x$FileName) %in% experimentFiles
        # channel loading correction (median)
        med <- apply(x[idx, tmtcols], 2, median)
        if (doPlot) {
            plot(med, type="l", bty="l", main=paste(fileIdentifiers[[i]],
            sep="/"), col=4)
        }
        x[idx, tmtcols] <- x[idx, tmtcols] / (matrix(1,nrow(x),1) %*% matrix(med, 1, length(med)))
    }
    if (doPlot) {
        dev.off()
    }
    x
}

#
# filter out all-zero traces
#
tmt.filterZeroTraces <- function(x) {
    # eliminate traces with all zeros
    x[is.na(x)] <- 0
    z <- apply(x != 0, 1, sum)
    x <- x[!(z==0),]
    x
}

# find all shared peptides
tmt.markSharedPeptides <- function(x) {
    ## count how often a PeptideHitID occurs
    nOccurences <- tapply(matrix(1, nrow(x), 1), x$PeptideHitID, sum) 
    ## determine indexes of peptides whose PeptideHitID occurs more than once
    sharedIdxs <- which(x$PeptideHitID %in% as.numeric(rownames(nOccurences[nOccurences > 1])))
    ## create indicator variable column for excel
    x$shared <- matrix(0, nrow(x), 1)
    x$shared[sharedIdxs] <- 1
    x
}

#
# protein level inference
#

# (a) protein similarity screening code
rpssGeomeanTopKAbundance <- function(x, a, f, k=5)
{
    stopifnot(class(x) == "acomp")
    sp <- split(cbind(x,a), f)
    topKMean <- function(x, nc, k=5) {
        # get back the matrix
        x <- unlist(x)
        nr <- length(x) / nc;
        x <- matrix(x, nr, nc);
        l <- nrow(x)
        if (l >= k) {
            idx <- sort(x[,ncol(x)], decreasing=T, index.return=T)$ix[1:k]
        } else {
            idx <- 1:nrow(x)
        }
        res <- (mean(acomp(x[,1:(ncol(x)-1),drop=F])))
        res <- matrix(unclass(res), 1, ncol(x)-1)
    }
    gm <- t(sapply(sp, topKMean, nc=ncol(x)+1))
    rownames(gm) <- names(sp)
    gm
}

# (b) wrap it into something more usable
tmt.peptide2protein <- function(xc, a, f=x$AccNumber, k=5)
{
    p <- rpssGeomeanTopKAbundance(xc, a, f, k=k)
    incomplete <- is.nan(apply(p, 1, sum))
    p <- acomp(p[!incomplete,])
    p
}

#
# plots a trace plot
#
traceplot <- function(trace, id, desc, class, scheme=rainbow, legend.location="topright", ...)
{
    f <- as.factor(class)
    if (nlevels(f) > 1) {
        colors <- scheme(nlevels(f))[f]
    } else {
        colors <- scheme(length(f))
    }
    matplot(t(trace), type="l", col=colors, bty="l", lty=1, ...)
    #legend(legend.location, bty="n", lty=1, col=colors, legend=id)
    # return the color vector
    invisible(colors)
}

# hugely convenient wrapper function for p-values
pvalue.or.na <- function(x, test=t.test) {
    t <- try(test(na.omit(x)), TRUE)
    if (is.null(t) || inherits(t, "try-error")) {
        p <- NA
    } else {
        p <- t$p.value
    }
    p
}

#
# timepointRatioTest
#
# The function takes a set of peptide-level TMT/iTRAQ measurements and 
# calculates the up-/down-regulation q-values for two specified time points.
#
# x: the peptide measurement, one per row
# tmtcols: the TMT measurement column numbers in x
# timepoints: a pair of TMT measurement indexes (1..6), specifying which TMT
#             labels should be compared
#
tmt.timepointRatioTest <- function(x, tmtcols, timepoints, mediancorrection=FALSE, test=wilcox.test)
{
    # calculate the logratios
    peptide.logratios <- log(x[,tmtcols[timepoints[1]]] / x[,tmtcols[timepoints[2]]])
    if (mediancorrection) {
        peptide.logratios <- peptide.logratios - median(peptide.logratios, na.rm=TRUE)
        peptide.logratios[abs(peptide.logratios) == Inf] <- NA
        peptide.logratios[is.nan(peptide.logratios)] <- NA
    }
    # have a look at how the peptide-level logratios are distributed; no subset
    # of peptides will be closer to a normal than this
    par(mfrow=c(1,2))
    qqnorm(na.omit(peptide.logratios))
    qqline(na.omit(peptide.logratios))
    plot(density(na.omit(peptide.logratios)), bty="l", col=4)
    abline(v=0, lty=2)
    # We have pretty big deviances; using a t-test is not a good idea.
    # Consequently, we make use of a non-parametric wilcoxon test.
    # The test is wrapped into a function that allows it to fail if the number
    # of available samples (i.e. number of peptides) for a protein becomes
    # too small.
    # run the test for each protein
    pvals <- tapply(peptide.logratios, x["AccNumber"], pvalue.or.na, test=test)
    # Also calculate the median peptide logratio for each protein.
    # Without the variance estimate, this is a pretty meaningless number, but it
    # allows the biologist to quickly see how much of a change we are talking about
    # for the proteins that are called significant.
    protein.logratios <- tapply(peptide.logratios, x["AccNumber"], median, na.rm=TRUE)
    # count how many peptides we have for each protein
    protein.npeptides <- tapply(!is.na(peptide.logratios), x["AccNumber"], sum)
    # Wrap everything into a data frame; this includes a multiple testing correction
    # to control the FDR, according to Benjamini/Yekutieli.
    proteins <- data.frame(
        AccNumber = names(protein.logratios),
        nPeptides = protein.npeptides,
        MedianLogratio = protein.logratios,
        pvalue = pvals,
        qvalue = p.adjust(pvals, method="fdr"))
    # Filter out all proteins for which we were unable to determine a p-value (i.e.
    # where the Wilcoxon test could not be carried out).
    bad <- is.na(pvals)
    proteins <- proteins[!bad,]
    # get the protein description
    m <- match(as.character(proteins[,"AccNumber"]), x[,"AccNumber"])
    proteins <- cbind(proteins, x[m, "Description"])
    names(proteins)[6] <- "Description"
    # sort by p-value (which implicitly sorts by q-value)
    proteins <- proteins[sort(proteins$pvalue, index.return=T)$ix,]
}


