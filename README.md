# Understanding relationships between chlamydia infection, symptoms and testing behavior: an analysis of data from Natsal-3

_Joanna Lewis and Peter J. White_

_NIHR Health Protection Research Unit in Modelling Methodology, Department of Infectious Disease Epidemiology, School of Public Health, Imperial College London_

This repository contains code for analysing data on chlamydia infection, testing and diagnosis from the Natsa-l3 study: http://www.natsal.ac.uk/natsal-3.aspx The data is available from the UK Data Archive: http://doi.org/10.5255/UKDA-SN-8178-1

The code for the descriptive analysis is organised in two files:

* `reason_setting.R` is used to examine reason for test, test setting, number of new partners in the last year, and correlations between them.
* `positivity_risk_newpart.R` compares positivity of tests carried out for different reasons, and by men and women reporting different numbers of partners.

The model-based analysis was carried out using the Stan software. There are four Stan model files: using priors for natural history parameters for men and women, and with and without partner notification (PN):
* `3-comp-model-stratified-men.stan` (no PN; men),
* `3-comp-model-stratified-women.stan` (no PN; women),
* `3-comp-model-stratified-PN-men.stan` (PN; men),
* `3-comp-model-stratified-PN-women.stan` (PN; women)

Each model is run using an `R` file named `run-<name of model file>`:
* `run-3-comp-model-stratified-men.R` (no PN; men),
* `run-3-comp-model-stratified-women.R` (no PN; women),
* `run-3-comp-model-stratified-PN-men.R` (PN; men),
* `run-3-comp-model-stratified-PN-women.R` (PN; women)

For the models without partner notification, choose to run the models with or without stratification by adjusting the element `str` in the Stan data input list `dt0`: either set `str` to 1 for all respondents, or use it to indicate the number of new partners in the last year reported by each person. See the `R` files for further details.