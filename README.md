# CF-FBA, a pipeline for creating cell-free metabolic models based on experimental data
CF-FBA generates cell-free metabolic models from the appropriate organism's GEM.
These models are reduced versions of the GEM based on experimental data.
To reduce the GEM, CF-FBA uses a VAE to perform dimensionality reduction on a set of FBA models.
CF-FBA uses a particular type of VAE called a Corr-VAE which uses a loss function that includes a correlation term with the experimental data.
Though CF-FBA has only been tested on E. coli models, it was written with generalizability in mind and should work with any COBRA compatible GEM.

## Project structure
- bio_models: where the FBA models live. Files will be written to here and the organism GEMs should be kept here.
- data: contains information about the cell-free systems. This includes experimental data as well as the experimental setups.
- genes: sequences of genes used in the cell-free systems.
- models: where computational models (primarily VAEs) are stored.
- notebooks: scratch work used in the development of this thesis. Contains the beginning of the RL system for model reduction.
- paper: contains the Master's thesis writeup of this work.
- results: contains results used to make comparison table in the thesis.
- scripts: scripts that are used to run individual parts of the pipeline. Scripts should be run from within the src/ directory. run_pipeline.sh provides an example of a full CF-FBA pipeline.
- src: the most important directory. Contains all of the code used to run the CF-FBA pipeline. Most scripts can be run as standalone python scripts using python scripts or imported as a module to use the functions.

## Current work
Changing from a threshold system of reduction to using reinforcement learning.
Standardization of command-line argument usage.

## Code location and Author
This code was written by Nicholas Larus-Stone and lives at https://github.com/nlarusstone/cf_fba
