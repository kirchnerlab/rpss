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
