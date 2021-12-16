# GREAT 3d

## The Genome Regulatory Epigenetics Tools 3D.
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/great_3d/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/great_3d/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/great_3D/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/great_3D/branch/master)

A software suite for multi-3D (epi)genomics

Suite of tools to analyse and visualise 3D genome architecture, transcription and epigenetics.

## Installation

- [ ] Add a proper setup file.
- [ ] Generate a proper python package structure.

## Usage

To test the software with 100 plasmodium genes run:
```
./trans3Dmap -p test/genePositions100.tab -e test/geneExpression100.tab <outfile.html>
```
Then an html file with the plotly graphic will be generated in the current folder.



### Copyright

Copyright (c) 2015-2021, Costas Bouyioukos


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.


Several students have contributed code in the initial stages of the project.

- First of all Miara Rakotomavo's M1 internship. Code into dev_Miara directory on dev_students branch.

- We have cloned-merged three repositories from 3 projet_longues student coding projects into respective directories on the dev_students branch.

1. Helene Kabbech [3D_transcription_map](https://github.com/kabhel/3D_transcription_map)
2. Hocine Meraouna [3D_transmap_Meraouna](https://github.com/hocinebib/3D_transmap_Meraouna)
3. Adam Belaiche [3DtranscriptionMap](https://github.com/toontun/3DTranscriptionMap)
