# load an orbiter results.txt file
#   filename - filename as a string
#   removeZeros - toggle removal of rows with a zero entry
xilac.load <- function(filename, removeZeros=FALSE) 
{
    x <- read.delim(filename, header=F)
    # no zeros
    if (removeZeros) {
        z <- apply(x[,8:12]==0, 1, sum)
        x <- x[z==0,]
    }
    # clean up IPI representation
    x[,1] <- gsub(".[0-9]+$", "", gsub("^IPI:", "", x[,1]))
    x
}

