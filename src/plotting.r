library(ggplot2)
source('/Users/mike/Dropbox/Projects/BioCUtils/exprset_utils.r')


plot_bar_expression <- function(exprset,
    raise=T, filename=NULL ,max_per_page=16, ...){
    # exponentiate the data for better plotting
    if (raise){
        df <- as.data.frame(2**exprs(exprset))
    } else {
        df <- as.data.frame(exprs(exprset))
    }
    colnames(df) <- make.unique(as.character(pData(exprset)$phenotype))
    df$ID <- as.character(fData(exprset[rownames(df),])$GeneSymbol)
    df.m <- melt(df, id = "ID")
    df.m$variable <- sapply(
        as.character(df.m$variable),
        function(name){
            strsplit(name,'\\.')[[1]][1]
        }
    )
    m <- melt(cast(df.m, ID ~ variable, fun = mean))
    colnames(m) <- c('ID','mean','phenotype')
    s <- melt(cast(df.m, ID ~ variable, fun = sd))
    colnames(s) <- c('ID','sd','phenotype')
    df.merged <- merge(m,s)
    if (exists("filename")){
        num_phenotypes = length(unique(m$phenotype))
        num_graphs = nrow(df.merged)/num_phenotypes
        num_pages = ceiling(num_graphs/max_per_page)
        num_rows = max_per_page*num_phenotypes
        n_facets = ceiling(sqrt(max_per_page))
        pdf(filename)
        for (i in seq(num_pages-1)){
            slice = seq(((i-1)*num_rows)+1,(i*num_rows))
            slice = slice[!(slice > nrow(df.merged))]
            df.segment=df.merged[slice,]
            p <- ggplot(df.segment, aes(y=mean,x=phenotype))
            p <- p + geom_bar(stat="identity",fill="white",colour="black") 
            print(p + facet_wrap("ID",ncol=n_facets,nrow=n_facets,...))
        }
        dev.off()
    } else {
        limits <- aes(
            ymax = df.merged$mean+df.merged$sd, 
            ymin=df.merged$mean-df.merged$sd
        )
        p <- ggplot(df.merged, aes(y=mean,x=phenotype))
        p <- p + geom_bar(stat="identity",fill="white",colour="black")
        p <- p + geom_errorbar(limits,width=0.25)
        p + facet_wrap("ID",scales="free_y",...)
    }
}

plot_pheno_dists <- function(exprset, sorted_indices=NULL, n=NULL, ...){
    #  this plots smoothed distributions of the data, coloured by phenotype
    if (is.null(sorted_indices)) sorted_indices = seq(nrow(exprset))
    if (is.null(n)) n=length(sorted_indices)
    # extract the relevant values
    relevant <- exprs(exprset)[sorted_indices,][1:n,]
    # get the probes and their corresponding gene symbols
    probes <- rownames(relevant)
    genenames = as.character(fData(exprset[probes,])$Gene.Symbol)
    if (length(genenames)==0){
        genenames = as.character(fData(exprset[probes,])$GeneSymbol)
    }
    if (length(genenames)==0){
        stop()
    }
    # uniqueify the genenames
    rownames(relevant) = make.unique(genenames)
    # build a data frame of the relevant genes and their phenotype
    top_n <- data.frame(
        ptype = factor(pData(exprset)$phenotype), 
        t(relevant)
    )
    # melt top_tfs for happy plotting
    top_n <- melt(top_n, id="ptype")
    colnames(top_n)[2] <- "probe"
    colnames(top_n)[3] <- "expression"
    # plot the top tfs
    p <- ggplot(top_n,aes(x=expression,fill=ptype))
    p <- p + geom_histogram(binwidth=0.1,alpha=0.3)
    p <- p + geom_density(alpha = 0.2)
    p + facet_wrap(~probe,...)
}

plot_gene_ts <- function(exprset,
    sorted_indices=NULL, n=NULL, zero_basal=T, 
    title="", ylabel="fold change", xlabel="time (days)", grid=F, ...
){
    if (is.null(sorted_indices)) sorted_indices = seq(nrow(exprset))
    if (is.null(n)) n=length(sorted_indices)
    relevant <- exprs(exprset[sorted_indices[1:n],])
    if (zero_basal){
        basal <- get_basal_expression(exprset[sorted_indices[1:n],])
        relevant <- relevant - basal
    }
    # get the probes and their corresponding gene symbols
    probes <- rownames(relevant)
    genenames = as.character(fData(exprset)[probes,1])
    # uniqueify the genenames
    rownames(relevant) = make.unique(genenames)
    titles <- colnames(relevant)
    # we need to pull out the time for each array
    time <- sub(".*\\.D([0-9]+).*","\\1",titles)
    # call the memory cells that aren't marked as d=106
    time[grep("MEM",time)] = 106
    # call the naive cells as d = 0
    time[grep("NVE",time)] = 0
    time <- as.numeric(time)
    relevant = as.data.frame(cbind(time=time,t(relevant)))
    data <- melt(relevant,id.vars="time")
    p <- ggplot(data, aes(x=time, y=value))
    p <- p + geom_point() + geom_smooth()
    p <- p + opts(title = title)
    if (grid){
        p <- p + facet_grid(~variable,...)
    } else{
        p <- p + facet_wrap(~variable,...)
    }
    
    p <- p + ylab(ylabel) 
    p <- p + xlab(xlabel)
}
