#!/usr/bin/env python3

import pandas
import argparse
import FileManager
import numpy

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="3D transcription map, main file")
    parser.add_argument('gene3Dpos', type=str,
                        help='path to the file contains 3D positions of genes')
    parser.add_argument('geneTranscriptionTable', type=str,
                        help='path to the file contains expression of genes')

    args = parser.parse_args()


    GENE_3D_POS = FileManager.checkFile(args.gene3Dpos)
    GENE_EXPRESSION = FileManager.checkFile(args.geneTranscriptionTable)
    DF_3DPOS = pandas.read_table(GENE_3D_POS, header=0)
    DF_GENE_EXPR = pandas.read_table(GENE_EXPRESSION, header=0)
    gene_in_3DPOS = list(DF_3DPOS.index)
    gene_in_EXPR = list(DF_GENE_EXPR.index)
    overlap_gene = set(gene_in_EXPR).intersection(gene_in_3DPOS)
    #print(list(overlap_gene))
    DF_3DPOS = DF_3DPOS.loc[list(overlap_gene)]
    DF_GENE_EXPR = DF_GENE_EXPR.loc[list(overlap_gene)]
    # DF_GENE_EXPR = DF_GENE_EXPR.transpose()
    test = numpy.corrcoef(DF_GENE_EXPR.values, rowvar=False)
    print(DF_GENE_EXPR.values)
    # DF_GENE_EXPR = DF_GENE_EXPR.transpose().corr(method="spearman")
    # print(DF_GENE_EXPR)
    # c=expr.loc[list(inter)[:10]].transpose().corr(method="spearman")
    # pandas.DataFrame.as_matrix??
    # numpy.corrcoeff???
