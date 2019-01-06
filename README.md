# 3DtranscriptionMap

## Author

Adam Bellaïche

adam.bellaiche@gmail.com

https://github.com/toontun/3DTranscriptionMap

University of Paris 7
M2BI
2018-2019

# Short Description

This is a program made for a project at the university of Paris 7 in master 2 of bio-informatics.
This program permits to display a 3D transcription map of a genome.

# Needed Files 

* file_manager.py
* map3d.py
* main.py
* correlation.R

# How to use it

* make sure you have python3
* make sure you are in the repository "src" with:
    * file_manager.py
    * map3d.py
    * main.py
    * correlation.R

* commands to run main.py:
    * python3 main.py -h
    * python3 main.py path_to_3Dpos_file path_to_gene_expression_profile_file N_closer_genes_to_take colormap maxDist

* command to have help on class file_manager.py or map3d.py:
    * python3 map3d.py
    * python3 file_manager.py

# Some examples

* python3 main.py ../data/SCHIZONTS.genes_pos.txt ../data/profiles_Otto2010.min 7 magma 0.5
    * path_to_3Dpos_file = ../data/SCHIZONTS.genes_pos.txt
    * path_to_gene_expression_profile_file = ../data/profiles_Otto2010.min
    * N_closer_genes_to_take = 7
    * colormap = magma
    * maxDist = 0.5

* you will have a scatter plot of matplotlib, there are some key event, mouse event:
    
    * "clt+alt+s" = save your current image
    * "mouse click" on a point = display gene name (thanks to Hocine Meraouna who shows me https://stackoverflow.com/questions/7908636/possible-to-make-labels-appear-when-hovering-over-a-point-in-matplotlib?fbclid=IwAR0XLYQ8h4wJaSBPzA5OGMmJLLcbZjdtQrC8KtVUX6FSTHsIxVW-HQF2zbA )
    * "mouse click left" = rotate the figure
    * "mouse click right" = zoom in terms of mouse mouvement

# Documentation

* you can read doc with command line: see in "How to use it" with help function. 
* go to 3DtranscriptionMap/doc/pythonDoc/Doxygen_doc/index.html
* you have the rules for the project in 3DtranscriptionMap/doc/projectInfo