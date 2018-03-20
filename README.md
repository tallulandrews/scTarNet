# scTarNet

This package provides the network inference method used in:
Bergiers I, Andrews T, et. al. (2018) Single-cell transcriptomics reveals a new dynamical function of transcription factors during embryonic hematopoiesis. eLife. doi: 10.7554/eLife.29312

It builds networks based on an input set of transcription factors and tests for indirect relationships and TF-TF interactions using distance correlations from the "energy" package.

The current version is highly parallelized but not yet optimized for large datasets. We are currently working on optimizing all parts of the code, and regular updates will be made as this progresses.

Installation :

>require("devtools")
>install_github('tallulandrews/scTarNet')

