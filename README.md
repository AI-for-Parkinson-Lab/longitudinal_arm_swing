# Longitudinal arm swing
This repository contains code related to the scientific publication: <b>'Longitudinal progression of digital arm swing measures during free-living gait in early Parkinson disease'</b>.

## Overview
The repository contains code to transform weekly aggregates of arm swing into results presented in the scientific publication. The notebooks are enumerated to ensure proper serial execution
and improve readability:
* `0.data_preparation.ipynb`: Applying inclusion/exclusion criteria, preparing clinical and measurement files;
* `1.descriptives.ipynb`: Group characteristics, between-group differences in characteristics and gait quantity, and longitudinal changes in gait quantity;
* `2.cross-sectional.ipynb`: Cross-sectional between-group differences in digital measures, correlations with clinical ratings, and reliability assessments;
* `3.l1tf_results.ipynb`: Evaluation of l1 trend filter applied to biweekly measurements, applying interpolation to extracted trends, and evaluating trends;
* `4.survival_ipynb`: Application and evaluation of survival analysis to account for possible censoring bias;
* `5.srm.ipynb`: Longitudinal evaluation of extracted digital measures;
* `6.regression.ipynb`: Assessment of the effects of medication (and its adverse effects such as dyskinesia) on (un)observed changes in digital measures;

Note that, between notebooks `2.cross-sectional.ipynb` and `3.l1tf_results.ipynb`, we computed the l1 trend filter in Matlab using code available in subdirectory `matlab`.
The l1 trend filter is run by executing `l1tf.m`.

## Running the notebooks
Biweekly digital measures required to run the notebooks will be made available to researchers applying for data access after publication of the paper.

This repository uses Poetry for dependency management, which can be installed using pip: 

```
pip install poetry 
```

After installation, optionally make a new virtual environment, and install the dependencies using `poetry install`. 
Also update the path to the data directory in `.env_sample.txt` and rename the file to `.env` (this is parsed in `src/constants.py`, where further pathing can be updated).

## Contact
For questions, please reach out to erik.post@radboudumc.nl.
