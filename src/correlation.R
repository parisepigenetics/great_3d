geneCorr <- read.table("../data/for_corr.csv")
geneCorr <- cor(t(geneCorr), method="spearman")
write.table(geneCorr, "../data/correlation.csv", sep="\t")