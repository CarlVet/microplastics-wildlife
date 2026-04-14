# Microplastics-in-South-African-Wildlife
## Overview

This repository contains all data processing, statistical analyses, and visualization scripts used in the study investigating microplastic (MP) concentration, composition, and size across wildlife samples in South Africa.

## Data

Sample-level metadata are derived from a relational SQL database associated with the Aleph∞One platform (project code: A1A1).

Due to data sensitivity, raw data are either provided in anonymized form in /data/processed or available upon reasonable request

## Methods

All analyses were conducted in R (v4.5.1).

Statistical models were implemented using the glmmTMB package, allowing flexible specification of non-Gaussian distributions:

- Gamma models (log link) for strictly positive continuous data

- Tweedie models (log link) for zero-inflated continuous data

## Repository structure

data/ → raw and processed datasets

scripts/ → fully reproducible analysis scripts

outputs/ → tables and figures used in the manuscript

## Reproducibility

To reproduce the analyses:

  Install required packages: install.packages(c("data.table", "glmmTMB", "broom", "dplyr"))

  Run scripts in order:
   01_data_cleaning.R
   02_exploratory_analysis.R
   03_gamma_model_concentration.R
   04_tweedie_model.R

## Outputs

Model summaries are saved in /outputs/tables

Figures are saved in /outputs/figures
