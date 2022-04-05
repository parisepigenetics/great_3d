# GREAT 3D
[![GitHub Actions Build Status](https://github.com/cbouyio/parisepigenetics/workflows/CI/badge.svg)](https://github.com/parisepigenetics/great_3d/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/great_3D/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/great_3D/branch/master)

## The Genome Regulatory Epigenetics Tools on 3D

A software suite for multi-3D (epi)genomics

Suite of tools to analyse and visualise 3D genome architecture, transcription and epigenetics.

## Installation

```python setup.py install --user```

## Usage
To run the transcriptome 3D map you must have the 2 3D coordinates files ready (for the genome and for the genes) and then run:

```great3D -g <genome_coordinates_table> -s <genes_coordinates_table> -e <gene_expression_table>```


Then a Dash application will launch on your local machine a you can access it by pointing your browser at http://127.0.0.1:8050/


### Copyright

Copyright (c) 2016-2022, Costas Bouyioukos, Universite Paris Cite et UMR7216


#### Acknowledgements

Project based on the
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.

Several students have contributed code in the initial stages of the project.

- First of all Miara Rakotomavo's M1 internship. Code into dev_Miara directory on dev_students branch.

- We have cloned-merged three repositories from 3 projet_longues student coding projects into respective directories on the dev_students branch.

1. Helene Kabbech [3D_transcription_map](https://github.com/kabhel/3D_transcription_map)
2. Hocine Meraouna [3D_transmap_Meraouna](https://github.com/hocinebib/3D_transmap_Meraouna)
3. Adam Belaiche [3DtranscriptionMap](https://github.com/toontun/3DTranscriptionMap)
