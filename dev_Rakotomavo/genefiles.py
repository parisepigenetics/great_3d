#!/usr/bin/python3
"""
gene expression and position files management programm to create
genes Dctionnary file of the N closest Genes
"""

'''
usage: genefiles.py [-h] [-e EXP_FILE] [-n NB_GENES] [output_file]

3D transcription map, main file

positional arguments:
  output_file           Path to the outputfile of genes Dictionnary . (or
                        STDOUT)

optional arguments:
  -h, --help            show this help message and exit
  -e EXP_FILE, --exp_file EXP_FILE
                        path to the file containing the gene expression
                        --DEFAULT :test_GE_Dist_Miara.tab
  -n NB_GENES, --nb_genes NB_GENES
                        Select the number of the closest genes-- DEFAULT = 10


'''
import argparse
import multiprocessing
import multiepigenomics_3d as me3d

parser = argparse.ArgumentParser(description="3D transcription map, main file")

#Positionnal argument containing the 3D coordinate of genes
'''parser.add_argument("pos_file" , type = str ,
help = 'Path to the file containing the x, y , z coordinates of genes ')'''

parser.add_argument("outfile", nargs='?', default="-" , type = argparse.FileType('w') ,
metavar='output_file', help="Path to the outputfile of genes Dictionnary . (or STDOUT)")


parser.add_argument("-e", "--exp_file" , help = 'path to the file containing the gene expression --DEFAULT :test_GE_Dist_Miara.tab ',
type = str, default = "test_GE_Dist_Miara.tab")
parser.add_argument("-n" , "--nb_genes" , help = "Select the number of the closest genes-- DEFAULT = 10" ,
type = int , default = 10)

args = parser.parse_args()

#CALCUL OF THE DISTANCE MATRIX  :
exp_file = args.exp_file
nb_genes = args.nb_genes
matrice_dist = me3d.distance_matrice(exp_file)
N_closest_genes = me3d.dico_matrice(matrice_dist , nb_genes)
#WRITTING and FILLIN IN THE OUTPUT FILE :
outfile = args.outfile
for i in N_closest_genes :
    outfile.write("ID__{} : \n{} ".format(i, N_closest_genes[i]))


print('\nle nombre de processeur utilisé est de : ', multiprocessing.cpu_count())
