#!/usr/bin/python3
"""selecting a file to give the distance Matrix and
the dictionnary of genes Distance as an output """

import argparse
import multiprocessing
import multiepigenomics_3d as me3d

parser = argparse.ArgumentParser()
parser.add_argument("infile" , metavar = "input_file" , type = str , help = 'Path of the input file of the gene ID file ')
parser.add_argument("matrice" , nargs='?' ,  default='-')


args = parser.parse_args()
matrice_dist = me3d.distance_matrice(args.infile)

from multiprocessing import Pool
with Pool(5) as p :
    x= p.map(me3d.dico_matrice,[matrice_dist])

if args.matrice:
    print(matrice_dist)


print('le nombre de processeur est de : ', multiprocessing.cpu_count())
