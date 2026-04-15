# Microplastics-in-South-African-Wildlife

## Overview

This repository contains all data processing, statistical analyses, and visualization scripts used in the study investigating microplastic (MP) concentration, composition, and size across wildlife samples in South Africa.

## Data

Sample-level metadata are derived from a relational SQL database associated with the Aleph∞One platform (project code: A1A1).

Due to data sensitivity, raw data are not publicly distributed, while processed datasets are generated within the workflow. Data may be available upon reasonable request.

## Methods

All analyses were conducted in R (v4.5.1) using the RStudio Integrated Development Environment.

Statistical modelling was performed using the glmmTMB package, allowing flexible specification of non-Gaussian distributions:

- Gamma models (log link) for strictly positive continuous responses  
- Tweedie models (log link) for zero-inflated continuous responses  

The log link function was used throughout to allow multiplicative interpretation of model coefficients.

## Repository structure

data/ → raw and processed datasets  
scripts/ → fully reproducible analysis scripts  
docs/ → supplementary material (optional)

## Analysis pipeline

The scripts should be executed in the following order:

01_data_cleaning.R  
02_exploratory_analysis.R  
03_gamma_model_concentration.R  
04_tweedie_model_concentration.R  
05_gamma_model_concentration_mpstype.R  
06_gamma_model_feret.R  

## Reproducibility

To reproduce the analyses:

Install required packages:
install.packages(c("data.table", "dplyr", "ggplot2", "glmmTMB", "broom", "broom.mixed", "leaflet", "sf", "ggtext"))

Run scripts in order:

01_data_cleaning.R  
02_exploratory_analysis.R  
03_gamma_model_concentration.R  
04_tweedie_model_concentration.R  
05_gamma_model_concentration_mpstype.R  
06_gamma_model_feret.R  

## Key outputs

Model summaries are generated within each script.

Figures and visual outputs are produced during the exploratory analysis workflow.

## Citation

If you use this repository, please cite the associated publication (to be added upon acceptance).

## Contact

For questions, collaborations, or data access requests, please contact: [your email]
