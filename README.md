# Hybrid Deep Metabolic Network

**Embedding Genome-Scale Model within Neural Networks as a Layer**

This repository contains the code and data for demostration of the method proposed in the paper:
> **Hybrid Deep Metabolic Network: Embedding Genome-Scale Model within Neural Networks as a Layer**
> *João R. C. Ramos, José Pinto, Ludovic Peeters, Patrick Dumas, Rui Oliveira*


## Overview
This toolbox implements a deep hybrid modeling approach for bioreactor cell culture data. It combines neural networks (FFNN or LSTM) with a genome-scale metabolic model (GEM) and first-principles equations.

## Getting Started

### Prerequisites
- MATLAB (The code files are `.m` files)

### Main Execution
To run the hybrid model training or simulate using previously trained models, run the main script in MATLAB:
`hybnet_train_main.m`

When you run the script, a dialog box will appear asking you to choose between **Simulate** and **New training**.

#### 1. Simulation Mode
Select **"Simulate"** to load the pre-trained models (`hybrid_FFNN_1.mat` and `hybrid_LSTM_1.mat`) and generate results.
This mode will:
- Simulate the saved LSTM and FFNN models.
- Generate comparison plots for concentrations (e.g., Br1 and Br6).
- Plot normalized predictions vs. experimental concentrations for all experiments (Br1–Br9).
- Plot predicted rates.

#### 2. Training Mode
Select **"New training"** to train the hybrid models from scratch.
This mode will:
- Define the hybrid model architectures (FFNN and LSTM) with the metabolic layer.
- Use the parameters defined in the script (e.g., `niter=150000`, `nruns=1`).
- Train the models using the specified training (Br1-Br4, Br9), cross-validation (Br7), and test (Br5, Br6, Br8) datasets.
- Save the trained models to `.mat` files.
- Export a summary of statistical results (RMSE, MSE, AIC, etc.) to `structures_fit_results.xlsx`.

## Repository Structure

### Root Directory
- `hybnet_train_main.m`: Main entry point for training and simulation.
- `hybrid_FFNN_1.mat` / `hybrid_LSTM_1.mat`: Pre-trained model files.
- `runs.mat`: Saved run data.
- `plot_*.m`: Various plotting scripts for paper figures.

### /data
Contains the dataset and processing scripts.
- `data.xlsx`: Feed DoE details and concentration simulations (synthetic dataset based on Robitaille et al., 2015).
- `data.mat`: Processed data structure containing:
    - `data(i).accum`: Total amount of metabolite in bioreactor (accumulated).
    - `data(i).m_r`: Reacted amounts over time.
    - `data(i).val`: Estimated rates at different growth phases.
- `main_data_processing.m`: Script used to generate `data.mat`.
- `RatesEstimation.m`: Script for estimating rates (`data(i).val`).
- `read_me_data_processing.doc`: Supplementary material describing data processing.

### /GEM
Contains the Genome-Scale Metabolic Model.
- `model.mat`: The GEM file.
- `main_least_square_regression.m`: Performs regression of estimated rates vs. model predictions to validate compatibility.

### /hybnet
Contains the core functions for the hybrid network toolbox (initialization, training, layers, etc.).

## Experimental Design (DoE)

The dataset covers experiments Br1 - Br9 (see `data.xlsx`).

```text
           Br5      
            |
            |
  Br1 ----- Br2     
   |        |
   |        |
  Br7 ----| Br9 |---- Br8
   |        |
   |        |
  Br3 ----- Br4     
            |       
            |       
           Br6      
```

## Contact
**Rui Oliveira**
LAQV-REQUIMTE, Department of Chemistry, NOVA School of Science and Technology, NOVA University Lisbon
Email: rmo@fct.unl.pt
