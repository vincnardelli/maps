# Spatial Econometrics Analysis

> MAPS Summer School – Spatial Analysis | Instructor: Vincenzo Nardelli

Spatial econometric analyses with R and Python: LISA, Spatial Error Models (SEM), and Spatial Lag Models (SLM).

## Quick Start

### 🐍 Python (Google Colab - Recommended)

**Click to open - Zero setup required!** Data loads automatically from GitHub.

1. **LISA Analysis** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/vincnardelli/maps/blob/main/Python/notebooks/1_lisa_colab.ipynb)
2. **Spatial Error Model** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/vincnardelli/maps/blob/main/Python/notebooks/2_spatial_error_model_colab.ipynb)
3. **Spatial Lag Model** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/vincnardelli/maps/blob/main/Python/notebooks/3_spatial_lag_model_colab.ipynb)

### 📊 R

```r
# Load required packages
library(sf)
library(spdep)
library(spatialreg)

# Run scripts
source("R/1_lisa.R")
source("R/2_spatial_error_model.R")
source("R/3_spatial_lag_model.R")
source("R/4_h3_lisa.R")
```

## Data

- **Shapefile**: `data/sdr_subnational_data_migration/shps/` - DHS subnational data (Tanzania, 2022)
- **Raster**: `data/TZ2022DHS_EDLITRWLIT_MS_v01/` - Literacy rate raster for H3 analysis

## Code Structure

```
.
├── Python/
│   └── notebooks/          # Jupyter notebooks for Google Colab
│       ├── 1_lisa_colab.ipynb
│       ├── 2_spatial_error_model_colab.ipynb
│       └── 3_spatial_lag_model_colab.ipynb
├── R/
│   ├── 1_lisa.R
│   ├── 2_spatial_error_model.R
│   ├── 3_spatial_lag_model.R
│   └── 4_h3_lisa.R
└── data/
    ├── sdr_subnational_data_migration/shps/
    └── TZ2022DHS_EDLITRWLIT_MS_v01/
