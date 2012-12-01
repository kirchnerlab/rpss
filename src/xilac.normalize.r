#
# median-normalize the columns of a matrix
#   x - a matrix
#   reference - if not NULL, the function will extract the
#               rows whose rownames are in 'reference' and
#               determine the median on this subset.
#
xilac.normalize <- function(x, reference=NULL)
{
    ret <- list();
    if (is.null(reference)) {
        ret$median <- apply(x, 2, median)
    } else {
        idx <- rownames(x) %in% reference
        ret$median <- apply(x[idx], 2, median)
    }
    ret$x <- x / (matrix(1,nrow(x),1) %*% matrix(ret$median,
      1, length(ret$median)))
    ret
}

