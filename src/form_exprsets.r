source('../conventions/io.r')
library(Biobase)

load('../data/immgen.data')
# remove anything from not-spleen
immgen <- immgen[,grepl("SP",colnames(exprs(immgen)))]
# only use Listeria plus all naive CD8 Tcells
i = grepl("LIS",colnames(exprs(immgen))) | grepl("NVE",colnames(exprs(immgen)))
LIS <- immgen[,i]
pData(LIS)$inclusion <- pData(LIS)$phenotype
# form the phenotype label
phenotype = rep("",ncol(LIS))
phenotype[grepl("Nve",pData(LIS)$title)]="Naive"
phenotype[grepl("Eff",pData(LIS)$title)]="Effector"
phenotype[grepl("Mem",pData(LIS)$title)]="Memory"
pData(LIS)$phenotype = factor(phenotype)
# save
outfile="../cache/tcell.data"
save(LIS, file=outfile)