 #!/usr/bin/env python3

import pandas
import argparse
from pathlib import Path

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="3D transcription map, main file")
    parser.add_argument('gene3Dpos', type=str,
                        help='path to the file contains 3D positions of genes')
    parser.add_argument('geneTranscriptionTable', type=str,
                        help='path to the file contains expression of genes')

    args = parser.parse_args()

    GENE_3D_POS = Path(args.gene3Dpos)
    GENE_EXPRESSION = Path(args.geneTrasncriptionTable)
    if GENE_3D_POS.is_file():
        data_frame_3d_pos = pandas.read_table(GENE_3D_POS, header=0)
    if GENE_EXPRESSION.is_file():
        data_frame_gene_expression = pandas.read_table(GENE_EXPRESSION, header=0)
        print(data_frame_gene_expression.head())
