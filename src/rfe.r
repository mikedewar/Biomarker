library(varSelRF)
library(Biobase)
library(log4r)

logger <- create.logger()
logfile(logger) <- file.path('../logs/rfe.log')
level(logger) <- log4r:::DEBUG
debug(logger,"--rfe log begins--")


rfe <- function(exprset){
    # apply the Random Forest variable selection measure
    data = t(exprs(exprset))
    ptype = factor(pData(exprset)$phenotype)
    #ptype = pData(exprset)$day
    debug(logger,"running varSelRF...")
    probenames <- varSelRF(data, ptype, ntreeIterat=5000)$selected.vars
    debug(logger,"chosen the following probenames")
    debug(logger,probenames)
    # extract the indices of each probe
    indices <- which(rownames(fData(exprset))%in%probenames)
    debug(logger,"with the following indices")
    debug(logger,indices)
    debug(logger,"which have the following symbols:")
    debug(logger,unlist(fData(exprset[indices,])$symbol))
    # form a signature set
    exprset[indices,]
}