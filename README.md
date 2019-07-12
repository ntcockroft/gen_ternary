# Generating PROTAC Mediated Ternary Complex Structure Predictions

This is a collection of python scripts that can be used to generate protac
mediated ternary structure predictions. The bash scripts that can be used to
connect the inputs and ouputs of the individual python scripts are also
provided. See the commented section at the beginning of the "run.sh" script
for additional information.


### Prerequisites

ZDOCK is used to perform the protein-protein docking and ZRANK (optional) is
used for scoring of the complexes. Both software packages can be obtained from:

http://zdock.umassmed.edu/software/


A list of the python libraries used is contained in the "environment.yml" file.
The python3 version of miniconda can be used to automatically set up an
environment with this libraries if desired. Miniconda can be obtained from:

https://docs.conda.io/en/latest/miniconda.html

After downloading and navigating to this repository the enivronment can be
built by using the command line or the now installed Anaconda Prompt. 

In the command line or prompt type:
```
conda env create -f environment.yml
```
The creation of the conda environment may take several minutes.


