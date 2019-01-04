#! /usr/bin/python3
"""
Genes position file Parsing
"""
import argparse

def gene_position_parser (position_file):
    """
    function that proceeds the parsing returns a dictionary of
    which the keys are the gene names and the values are the 
    x y z position coordinates
    """
    POS_DIC = {}
    with open(position_file, "r") as pf:
        for line in pf :
            line=line.split()
            if len(line) > 4 :
                POS_DIC[line[0]]=line[1:]
    return POS_DIC

if __name__ == "__main__":

    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("pos_file", help="the file containing the x y z coordinates of the genes", type=str)
    ARGS = PARSER.parse_args()
    POS_FILE = ARGS.pos_file
    
    GENE_POS = gene_position_parser(POS_FILE)
    for gene in GENE_POS:
        print("{} : {}".format(gene,GENE_POS[gene]))
