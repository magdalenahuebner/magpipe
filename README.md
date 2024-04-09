# MagPipe

MagPipe is a comprehensive analysis pipeline developed specifically for performing differential phosphosite occupancy analysis on phosphoproteomics data. It has been tailored to process the dataset published in the following paper, which explores phosphorylation patterns across various cell lines and kinase perturbations:

*Hijazi M, Smith R, Rajeeve V, Bessant C, Cutillas PR. Reconstructing kinase network topologies from phosphoproteomics data reveals cancer-associated rewiring. [Nat Biotechnol. 2020 Apr](https://www.nature.com/articles/s41587-019-0391-9)*

However, MagPipe can also be used to analyse other phosphoproteomics datasets, as it includes functionalities for harmonising protein IDs, quality control, quantile normalisation, batch correction, and differential phosphosite occupancy analysis. Unique features of the MagPipe package are its method for estimating missing fold changes and the calculation of signal-intensity dependent confidence scores.

MagPipe was created with [cookiecutter-bioinformatics-project](https://github.com/maxplanck-ie/cookiecutter-bioinformatics-project), a bioinformatics pipeline template provided by the Max Planck Institute of Immunobiology and Epigenetics, which builds on [cookiecutter-fair-data-science](https://github.com/FAIR4HEP/cookiecutter4fair) and [Snakemake](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#distribution-and-reproducibility). This ensures a streamlined, reproducible workflow for phosphoproteomics analysis.

## Getting Started

### Dependencies
- Python 3.8+
- R 4.0+
- Required Python packages: snakemake, uniprot-id-mapper (see `environment.yaml` for a complete list)
- Required R packages: limma, reshape2, plyr

### Cloning and Running the Snakemake Project

1. **Clone the MagPipe repository:**
```bash
git clone https://github.com/magdalenahuebner/magpipe.git
cd magpipe
```

2. **Activate the Conda environment:**
```bash
conda env create -f environment.yaml
conda activate magpipe
```

3. **Customise:**
Add you data to folder. Make sure formatted correctly (link to next section). And that paths in the snakemake file are correct and name of replcate name of your data. You can also run multiple samples at the same time (link)

4. **Execute the pipeline:**
```bash
snakemake --cores all
```

### Alternative: Using MagPipe as a Module
MagPipe can also be integrated into your Python projects as a module:

1. **Install MagPipe:**
```bash
pip install git+https://github.com/magdalenahuebner/magpipe.git
```

## Input File Formatting

Under construction...

## Workflow Overview

The MagPipe workflow encompasses several critical steps: 

![Workflow Diagram](img/dag.png "Workflow Overview")

Under construction...
