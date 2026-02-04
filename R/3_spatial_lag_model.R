# -------------------------------------------------------------
#  MAPS Summer School – Spatial Analysis
#  Module: Spatial Analysis – Teacher: Vincenzo Nardelli
# -------------------------------------------------------------

# 1) Load packages
# -------------------------------------------------------------
library(sf)          # vector geometry
library(spdep)       # spatial dependence objects & tests
library(spatialreg)  # SAR, SEM, SLM, etc.
library(dplyr)       # data wrangling
library(ggplot2)     # (optional) visualisation

set.seed(123)        # reproducibility
# -------------------------------------------------------------
# 2) Read the shapefile
# -------------------------------------------------------------
shp_path <- "../data/sdr_subnational_data_migration/shps/sdr_subnational_data_dhs_2022_lvl_2.shp"
tz <- st_read(shp_path, quiet = TRUE)   # keep geometry for later plotting

# -------------------------------------------------------------
# 3) Variable names & quick summaries
# -------------------------------------------------------------
vars_y <- "AHMIGRWEMP"                     # migration reason: employment
vars_x <- c("EDMDIAWN3M", "WSSRCEHIMP",
            "AHBTHPWIMG", "EDEDUCWSEH")    # explanatory variables

# Short summary statistics
summary(tz$AHMIGRWEMP)
summary(tz$EDMDIAWN3M)
summary(tz$WSSRCEHIMP)
summary(tz$AHBTHPWIMG)
summary(tz$EDEDUCWSEH)

# -------------------------------------------------------------
# 4) Drop geometry for modelling
# -------------------------------------------------------------
tz_df <- tz %>% st_drop_geometry()

# -------------------------------------------------------------
# 5) Baseline OLS regression
#    y ~ x1 + x2 + x3 + x4
# -------------------------------------------------------------
ols <- lm(AHMIGRWEMP ~ EDMDIAWN3M + WSSRCEHIMP + AHBTHPWIMG + EDEDUCWSEH,
          data = tz_df)
print(summary(ols))

# -------------------------------------------------------------
# 6) Spatial weights
#    Queen contiguity (share an edge or a vertex)
# -------------------------------------------------------------
nb_q  <- poly2nb(tz, queen = TRUE)
lw_q  <- nb2listw(nb_q, style = "W", zero.policy = TRUE)

# -------------------------------------------------------------
# 7) Moran's I on OLS residuals
# -------------------------------------------------------------
ols_resid_moran <- moran.test(residuals(ols), lw_q, zero.policy = TRUE)
print(ols_resid_moran)

# -------------------------------------------------------------
# 8) LM tests to decide between spatial lag and error
# -------------------------------------------------------------
lm_tests <- lm.RStests(ols, lw_q,
                       test = c("LMlag", "LMerr", "RLMlag", "RLMerr"),
                       zero.policy = TRUE)
print(lm_tests)
# Interpretation: compare the p‑values for LMlag, LMerr, etc.

# -------------------------------------------------------------
# 9) Fit a Spatial Lag Model (SAR / SLM)
# -------------------------------------------------------------
sar_lag <- lagsarlm(AHMIGRWEMP ~ EDMDIAWN3M + WSSRCEHIMP + AHBTHPWIMG + EDEDUCWSEH,
                    data = tz,
                    listw = lw_q,
                    zero.policy = TRUE,
                    method = "eigen")   # fast eigenvalue algorithm
print(summary(sar_lag, Nagelkerke = TRUE))

# -------------------------------------------------------------
# 10) Residual Moran's I after SAR lag
# -------------------------------------------------------------
cat("\nResidual Moran's I (SAR lag):\n")
print(moran.test(residuals(sar_lag), lw_q, zero.policy = TRUE))

# -------------------------------------------------------------
# 11) (Optional) Compute direct, indirect and total impacts
# -------------------------------------------------------------
imp <- impacts(sar_lag, listw = lw_q, R = 1000, zero.policy = TRUE)
print(imp)
