"""Module to deal with the translation of genomic coordinates.

Originating from a standard annotation file (i.e. GTF) into Hi-C fragments/bins cordinates for
integration to the GREAT framework analysis.
"""

import pandas as pd


def compute_gene_coordinates(pos_file_Genes, pos_file_Genome3D, proximity=5000):
    """Get two tab files, one with the genenames andtheir respective start sites and one with the 3D genome coordinates (as midpoints of fragments of fixed length i.e. resoluton).

    Return a dataFrame with the 3D coordinates of genes.
    """
    # TODO in the first run the closest known midpoint will be considered, however n later versions amore proper geometric clculation of the 3D position of a gene must be considered.
    dfGenome = pd.read_table(pos_file_Genome3D)
    dfGenes = pd.read_table(pos_file_Genes)
    # First simple approach, the gene is assigned to its closest bin median point.
    # Check chromosome consistancy
    chromsGenome = list(set(dfGenome.loc[:, "chr"].tolist()))
    chromsGenes = list(set(dfGenes.loc[:, "chr"].tolist()))
    if chromsGenes != chromsGenome:
        raise Exception("Chromosome consistency between genes and the genome is not met. Check the input files.")
    # Prepare the data frame lists
    names = []
    chrs = []
    Xs = []
    Ys = []
    Zs = []
    # Loop through all chromosomes in the files
    for ch in chromsGenes:
        dfChromGene = dfGenes[dfGenes["chr"] == ch]
        dfChromGenome = dfGenome[dfGenome["chr"] == ch]
        for row1 in dfChromGene.itertuples():
            for row2 in dfChromGenome.itertuples():
                if abs(row1.start - row2.midpoint) < proximity:
                    break
            names.append(row1.Index)
            chrs.append(ch)
            Xs.append(row2.X)
            Ys.append(row2.Y)
            Zs.append(row2.Z)
    # Generate the new data frame and reurn it
    df = pd.DataFrame({"chr": chr, "X": Xs, "Y": Ys, "Z": Zs})
    df.index = names
    return df
