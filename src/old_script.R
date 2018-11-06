#' Compute the 3-dimensional trancriptome map.
#'
#' Calulate the 3d transcription correlation score. Gets two matrices, one for gene expression and one for gene position.
#'
#' @param genePos a data frame with X, Y and Z column names which specify the 3D coordinates of the genes. The row names much be common gene nakes.
#' @param geneExpr a data frame with the gene expression values in different conditions (columns) of genes (row names).
transcriptomeMap3D <- function(genePos, geneExpr, n=10, ...) {
  # Initialise the two report vectors.
  geneNames <- vector();
  gene3Dcors <- vector();
  # Get data that is overlaping.
  overlapGenes <- intersect(rownames(geneExpr), rownames(genePos));
  # Calulate the distance and correlation matrices.
  geneDist <- as.matrix(dist(genePos[overlapGenes, c("X", "Y", "Z")], diag=TRUE, upper=TRUE));
  geneCorr <- cor(t(geneExpr[overlapGenes, ]), method="spearman");
  # Calculate teh 3D gene expresion score.
  for (gene in rownames(geneDist)) {
    closeGenesNames <- names(sort(geneDist[gene,])[2:(n+1)]);
    gene3Dcorr <- sum(abs(geneCorr[gene,][closeGenesNames]));
    geneNames <- append(geneNames, gene);
    gene3Dcors <- append(gene3Dcors, gene3Dcorr;)
  }
  transMap3D <- data.frame(corr3D=gene3Dcors);
  rownames(transMap3D) <- geneNames;
  return(transMap3D)
}

# Computes an one dimensional transcriptome map.
#'
#' Calulate the 1d transcription correlation score. Gets two matrices, one for gene expression and one for gene position.
#'
#' @param geneStart a data frame with a column which specifies the genome start site of the genes and one the chromosome number. The row names much be common gene names.
#' @param geneExpr a data frame with the gene expression values in different conditions (columns) of genes (row names).
transcriptomeMap1D <- function(geneStart, geneExpr, n=10, ...) {
  # Initilise report daata frame.
  transMap1D <- data.frame(corr1D = NULL, chrom = NULL);
  # Get data that is overlaping.
  overlapGenes <- intersect(rownames(geneExpr), rownames(geneStart));
  geneStart <- geneStart[overlapGenes,, drop=FALSE];
  geneExpr <- geneExpr[overlapGenes, ];
  # Do the 1D analysis per chomosome.
  for (chr in unique(geneStart$chrom)) {
    # Initialise the report vector
    gene1Dcors <- vector();
    gsChr <- subset(geneStart, chrom == chr);
    genesChr <- rownames(gsChr);
    # Calculate the distance and correlation matrices.
    geneDist <- as.matrix(dist(gsChr[,"geneStart", drop=FALSE], diag=TRUE, upper=TRUE));
    geneCorr <- cor(t(geneExpr[genesChr, ]), method="spearman");
    # Calculate teh 3D gene expresion score.
    for (gene in genesChr) {
      closeGenesNames <- names(sort(geneDist[gene,])[2:(n+1)]);
      gene1Dcorr <- sum(abs(geneCorr[gene,][closeGenesNames]));
      gene1Dcors <- append(gene1Dcors, gene1Dcorr);
    }
    transMapChrom <- data.frame(corr1D = gene1Dcors, chrom = chr);
    rownames(transMapChrom) <- genesChr;
    transMap1D <- rbind(transMap1D, transMapChrom);
  }
  return(transMap1D)
}
