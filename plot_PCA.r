pdf("PCA_plots.pdf")

# DATASET A

datasetA="A_datasetPCA_filt.eigs"
resA = read.table(datasetA, col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Pop"))

popGroups=read.table("indiv_pop.tab", col.names=c("Sample", "Population","Color","Point"))
mergedA = merge(resA, popGroups, by="Sample") # notice that it will only merge individuals that are present in both files, so you can actually exclude some indiv using this.
plot(mergedA$PC1, mergedA$PC2, col=mergedA$Color, pch=mergedA$Point)
legend("topleft", legend=c("YRI","CEU","CHB","PER","SH", "NC","NC_lc"), col=c("yellow","darkorchid1","green","black","red", "blue","dodgerblue4"),pch=19)

# DATASET B

datasetB="B_datasetPCA_filt.eigs"
resB = read.table(datasetB, col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Pop"))

popGroups=read.table("indiv_pop.tab", col.names=c("Sample", "Population","Color","Point"))
mergedB = merge(resB, popGroups, by="Sample") # notice that it will only merge individuals that are present in both files, so you can actually exclude some indiv using this.
plot(mergedB$PC1, mergedB$PC2, col=mergedB$Color, pch=mergedB$Point)
legend("topleft", legend=c("YRI","CEU","CHB","PER","SH", "NC","NC_lc"), col=c("yellow","darkorchid1","green","black","red", "blue","dodgerblue4"),pch=19)

# DATASET C

datasetC="C_datasetPCA_all.eigs"
resC = read.table(datasetC, col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Pop"))

popGroups=read.table("indiv_pop.tab", col.names=c("Sample", "Population","Color","Point"))
mergedC = merge(resC, popGroups, by="Sample") # notice that it will only merge individuals that are present in both files, so you can actually exclude some indiv using this.
plot(mergedC$PC1, mergedC$PC2, col=mergedC$Color, pch=mergedC$Point)
legend("topleft", legend=c("YRI","CEU","CHB","PER","SH", "NC","NC_lc"), col=c("yellow","darkorchid1","green","black","red", "blue","dodgerblue4"),pch=19)

#DATASET D

datasetD="D_datasetPCA_all.eigs"
resD = read.table(datasetD, col.names=c("Sample", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Pop"))

popGroups=read.table("indiv_pop.tab", col.names=c("Sample", "Population","Color","Point"))
mergedD = merge(resD, popGroups, by="Sample") # notice that it will only merge individuals that are present in both files, so you can actually exclude some indiv using this.
plot(mergedD$PC1, mergedD$PC2, col=mergedD$Color, pch=mergedD$Point)
legend("topleft", legend=c("YRI","CEU","CHB","PER","SH", "NC","NC_lc"), col=c("yellow","darkorchid1","green","black","red", "blue","dodgerblue4"),pch=19)

dev.off()