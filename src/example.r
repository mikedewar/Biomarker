source('../conventions/io.r')
source('filter.r')
source('plotting.r')
source('rfe.r')
library(Biobase)

#load('../cache/tcell.data') #LIS

#LIS.filtered = filter_probesets(LIS,by=has_human_homologue)
#LIS.filtered = filter_probesets(LIS.filtered,by=is_differentially_expressed)

#save(LIS.filtered,file='../cache/LIS.filtered.out')

load('../cache/LIS.filtered.out') #LIS.filtered

names = c()

for (i in seq(1000)){
    signature = rfe(LIS.filtered)
    names = c(names,fData(signature)$symbol)
}
    

#cache.save(signature,'signature.out')

#p <- plot_gene_ts(signature,scales="free_y")
#p