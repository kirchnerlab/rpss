#
# rpssGeomeanTopKAbundance.R
#
# Copyright (C) 2009-2012 Marc Kirchner
#
# This function calculates the geometric mean of the top k abundance
# peptides in a table, following the idea that a projection into the
# simplicial domain will effectively decorrelate the sum-normalized
# TMT measurements.
#
# Parameters:
#   x   -   An m x n table of m peptide measurements (i.e. TMT intensities)
#           across a set of n conditions.
#   a   -   An m x 1 abundance vector (this often corresponds to the row
#           sums of x, but sometimes x is already row-normalized; in these
#           cases, the vector a is necessary to select the k top-abundance
#           peptides for a factor f (see below).
#   f   -   A m x 1 factor vector, in most cases this is a vector that
#           holds the protein names/numbers/IDs for every column in the
#           matrix x. The "top k" approach is carried out for each
#           distinct value in this vector.
#   k   -   The maximum (!) number of peptides that contribute to the
#           estimated protein trace.
#
# Return value:
#   A simplicial matrix gm with protein-level relative abundance traces.
#   Use idt(gm) to regain the protein traces in the raw data domain.
#
# Example:
#   library(compositions)
#   x <- ... read some abundance data from disk ...
#   f <- ... extract a vector with corresponding protein names ...
#   a <- apply(x, 1, sum)
#   xc <- acomp(x)
#   p <- rpssGeomeanTopKAbundance(x, a, f, 10)
#
rpssGeomeanTopKAbundance <- function(x, a, f, k=5)
{
    stopifnot(class(x) == "acomp")
    sp <- split(cbind(x,a), f)
    topKMean <- function(x, nc, k=5) {
        #browser()
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
