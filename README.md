# B_ALL_differentiation
Original code developed for [Geßner, Jolly et al](https://doi.org/10.1101/2025.09.10.674954)

- [ABCLeukemogenesisModelselection.ipynb](ABCLeukemogenesisModelselection.ipynb) is a jupyter notebook containing the whole dynamical modeling procedure (ODE model, clones*parameters matrix construction,  Approximate Bayesian Computation parameter inference with pyABC).

  
- [BarcodeExtraction.r](BarcodeExtraction.r): pipeline for the extraction and filtering of lineage barcode sequences generated as described in the manuscript, and then add the barcode information to a Bioconductor SingleCellExperiment object. The Barcodes are extracted from a BAM file containing exclusively lineage barcode sequences assigned to single-cells. 


- [XGboostTrainingandModelTesting.r](XGboostTrainingandModelTesting.r) : pipeline for learning the clone classification model with XGBoost as described in the manuscript and then testing the statistical significance of the model predictions by label permutations.
