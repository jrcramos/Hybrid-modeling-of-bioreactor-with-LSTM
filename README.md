# Hybrid Modeling of Bioreactor with LSTM

A MATLAB implementation of deep hybrid modeling for HEK293 cell culture processes, combining Long Short-Term Memory (LSTM) networks with first-principles equations for bioprocess prediction and optimization.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Data Description](#data-description)
- [Model Architecture](#model-architecture)
- [Usage Examples](#usage-examples)
- [Output Files](#output-files)
- [Customization](#customization)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

## Overview

This repository contains the implementation of hybrid neural network models that combine:
- **Deep learning**: LSTM and Feed-Forward Neural Networks (FFNN)
- **First-principles**: Mass balance equations and biokinetic relationships
- **Bioprocess data**: Synthetic HEK293 cell culture datasets

The hybrid approach leverages the strengths of both mechanistic knowledge and data-driven learning to predict bioprocess dynamics with improved accuracy and interpretability.

## Features

- 🧬 **Bioprocess-specific**: Optimized for mammalian cell culture (HEK293)
- 🤖 **Hybrid AI**: Combines LSTM/FFNN with mechanistic equations
- 📊 **Pre-trained models**: Ready-to-use LSTM and FFNN models included
- 📈 **Comprehensive analysis**: Automated plotting and statistical evaluation
- 🔄 **Flexible training**: Easy model retraining with custom parameters
- 📁 **Structured data**: Well-organized synthetic datasets with DoE design

## System Requirements

### Software Requirements
- **MATLAB** R2019a or later
- **Required MATLAB Toolboxes:**
  - Deep Learning Toolbox
  - Statistics and Machine Learning Toolbox
  - Optimization Toolbox
  - Curve Fitting Toolbox (recommended)

### Hardware Requirements
- **RAM**: Minimum 8GB (16GB+ recommended for large-scale training)
- **Storage**: At least 500MB free space
- **CPU**: Multi-core processor recommended for faster training

## Installation

1. **Clone the repository:**
   ```bash
   git clone https://github.com/jrcramos/Hybrid-modeling-of-bioreactor-with-LSTM.git
   cd Hybrid-modeling-of-bioreactor-with-LSTM
   ```

2. **Open MATLAB** and navigate to the repository directory:
   ```matlab
   cd('/path/to/Hybrid-modeling-of-bioreactor-with-LSTM')
   ```

3. **Add all subdirectories to MATLAB path:**
   ```matlab
   addpath(genpath(pwd))
   ```

## Quick Start

### Option 1: Use Pre-trained Models (Recommended for first-time users)
```matlab
% Run the main script
hybnet_train_main

% When prompted, select "Simulate" to use pre-trained models
% This will generate plots and analysis using existing LSTM and FFNN models
```

### Option 2: Train New Models
```matlab
% Run the main script
hybnet_train_main

% When prompted, select "New training" to train models from scratch
% Note: Training may take 10-30 minutes depending on your system
```

## Data Description

### Dataset Overview
- **Source**: Synthetic HEK293 cell culture data
- **Base model**: Robitaille et al. (2015) metabolic model
- **Design**: 3³ factorial Design of Experiments (DoE)
- **Experiments**: 9 bioreactor runs (Br1-Br9)
- **Duration**: 240 hours per run

### DoE Design Structure
```
         Br5
          |
          |
    Br1 --+-- Br2
     |         |
     |         |
Br7 --+   Br9  +-- Br8
     |         |
     |         |
    Br3 --+-- Br4
          |
          |
         Br6
```

### Key Data Files
- **`data/data.xlsx`**: Raw experimental data with feed compositions
- **`data/data.mat`**: Processed data in MATLAB format
- **Pre-trained models**: `hybrid_LSTM_1.mat`, `hybrid_FFNN_1.mat`

### Data Structure
Each experiment contains:
- **`data(i).time`**: Time points (hours)
- **`data(i).conc`**: Metabolite concentrations (mM)
- **`data(i).accum`**: Accumulated masses in bioreactor
- **`data(i).m_r`**: Reacted amounts (metabolic rates)
- **`data(i).vol`**: Reactor volume (L)

## Model Architecture

### LSTM Hybrid Model
- **Input layer**: 27 features (concentrations + process variables)
- **Hidden layers**: 
  - Dense layer (27 → 10 neurons, ReLU activation)
  - Dense layer (10 → 10 neurons, ReLU activation)  
  - LSTM layer (10 → 4 neurons)
- **Output**: 4 principal components of reaction rates

### FFNN Hybrid Model
- **Input layer**: 27 features
- **Hidden layers**:
  - Dense layer (27 → 10 neurons, ReLU activation)
  - Dense layer (10 → 10 neurons, ReLU activation)
  - Dense layer (10 → 4 neurons, ReLU activation)
- **Output**: 4 principal components of reaction rates

### Training Parameters
- **Algorithm**: ADAM optimizer
- **Learning rate**: Adaptive (starts at 1e-3, decays to 1e-6)
- **Iterations**: 5000 (default, configurable)
- **Training runs**: 5 repetitions with different initializations
- **Validation**: Cross-validation with held-out experiments

## Usage Examples

### Basic Simulation
```matlab
% Load and simulate pre-trained LSTM model
hybnet_train_fun('hybrid_LSTM_1');

% Load and simulate pre-trained FFNN model  
hybnet_train_fun('hybrid_FFNN_1');
```

### Custom Training
```matlab
% Define custom model architecture
layers = {
    hnetfflayer(27, 20, 'relu'),    % Input layer
    hnetfflayer(20, 15, 'relu'),    % Hidden layer  
    hnetLSTMlayer(15, 4)            % LSTM output layer
};

% Training parameters
niter = 3000;           % Number of iterations
nruns = 10;             % Number of training repetitions
npcs = 4;               % Principal components

% Define training/validation splits
Indtr = [1 2 3 4];      % Training experiments
Indcr = [7 9];          % Validation experiments  

% Run training
hybnet_train_fun('my_custom_model', layers, Indtr, Indcr, {niter, nruns, npcs});
```

### Data Processing
```matlab
% Regenerate processed data from raw Excel files
cd data
main_data_processing
cd ..
```

## Output Files

After running the models, you'll get:

### Plots
- **Concentration profiles**: Time-series of all metabolites
- **Parity plots**: Predicted vs. experimental values
- **Training curves**: Loss evolution during optimization
- **Statistical analysis**: R², RMSE, AIC metrics

### Data Files
- **`structures_fit_results.xlsx`**: Comprehensive model performance metrics
- **Model files**: Saved trained models (`.mat` format)
- **Figure files**: All generated plots

### Performance Metrics
- **RMSE**: Root Mean Square Error
- **R²**: Coefficient of determination
- **AIC**: Akaike Information Criterion
- **Cross-validation scores**: Validation set performance

## Customization

### Adding New Experiments
1. Add data to `data/data.xlsx` following the existing format
2. Run `data/main_data_processing.m` to generate `data.mat`
3. Update experiment indices in `hybnet_train_main.m`

### Modifying Model Architecture
```matlab
% Example: Deeper LSTM network
layers = {
    hnetfflayer(27, 30, 'relu'),
    hnetfflayer(30, 20, 'relu'),
    hnetLSTMlayer(20, 10),
    hnetfflayer(10, 4, 'linear')
};
```

### Adjusting Training Parameters
```matlab
% In hybnet_train_main.m, modify:
niter = 10000;          % More iterations for complex models
nruns = 20;             % More repetitions for robustness
npcs = 6;               % More principal components
```

## Troubleshooting

### Common Issues

**Problem**: "Undefined function or variable" errors
**Solution**: Ensure all subdirectories are added to MATLAB path:
```matlab
addpath(genpath(pwd))
```

**Problem**: Out of memory during training
**Solution**: 
- Reduce batch size in model parameters
- Use fewer training repetitions (`nruns`)
- Close other MATLAB applications

**Problem**: Poor model convergence
**Solution**:
- Increase number of iterations (`niter`)
- Adjust learning rate in `hnet` parameters
- Try different random initializations

**Problem**: Missing toolbox functions
**Solution**: Install required MATLAB toolboxes or contact your MATLAB administrator

### Getting Help
1. Check the troubleshooting section above
2. Review the original paper for methodological details
3. Examine the `data/read_me_data_processing.doc` file
4. Contact the authors (see [Contact](#contact) section)

## Citation

**Please cite this work as:**

```bibtex
@article{ramos2024deep,
  title={Deep hybrid modeling of a HEK293 process: Combining long short-term memory networks with first principles equations},
  author={Ramos, João RC and Pinto, José and Poiares-Oliveira, Gil and Peeters, Ludovic and Dumas, Patrick and Oliveira, Rui},
  journal={Biotechnology and Bioengineering},
  pages={1--15},
  year={2024},
  publisher={Wiley},
  doi={10.1002/bit.28668}
}
```

**Paper DOI**: https://doi.org/10.1002/bit.28668

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Contact

**Corresponding Author:**  
Rui Oliveira  
LAQV-REQUIMTE, Department of Chemistry  
NOVA School of Science and Technology  
NOVA University Lisbon, Portugal  
📧 Email: rmo@fct.unl.pt

**Authors:**
- João R. C. Ramos¹ 
- José Pinto¹  
- Gil Poiares-Oliveira¹
- Ludovic Peeters²
- Patrick Dumas²
- Rui Oliveira¹

¹ LAQV-REQUIMTE, Department of Chemistry, NOVA School of Science and Technology, NOVA University Lisbon, 2829-516 Caparica, Portugal  
² GSK, 89 rue de l'Institut, 1330 Rixensart, Belgium

---
*This README was designed to help scientists easily understand and reuse this hybrid modeling framework for bioprocess applications.*
