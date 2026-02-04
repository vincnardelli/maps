# -------------------------------------------------------------
#  MAPS Summer School – Spatial Analysis
#  Module: Spatial Analysis – Teacher: Vincenzo Nardelli
# -------------------------------------------------------------

# --------------------------------------------------------------------
# 0) Load required packages
# --------------------------------------------------------------------
library(sf)          # vector data handling
library(spdep)       # spatial dependence tools
library(spatialreg)  # spatial regression models
library(dplyr)       # data manipulation
library(ggplot2)     # plotting
library(patchwork)   # combine ggplots
set.seed(123)        # reproducibility

# --------------------------------------------------------------------
# 1) Read polygons and quick data checks
# --------------------------------------------------------------------
tz <- st_read("../data/sdr_subnational_data_migration/shps/sdr_subnational_data_dhs_2022_lvl_2.shp",
              quiet = TRUE)

# Summaries of key variables
summary(tz$AHMIGRWEMP)   # migration rate
summary(tz$EDEDUCWSEH)   # education level

# --------------------------------------------------------------------
# 2) Define colour scheme for choropleth maps
# --------------------------------------------------------------------
fertility_fill_low  <- "#fff7bc"
fertility_fill_high <- "#d7301f"

# --------------------------------------------------------------------
# 3) Exploratory plots
# --------------------------------------------------------------------
## 3a) Scatter of education vs migration
scatter <- ggplot(data = tz) +
  geom_point(aes(x = EDEDUCWSEH, y = AHMIGRWEMP)) +
  labs(x = "Education", y = "Migration") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))

## 3b) Choropleth of education
map_education <- ggplot(tz) +
  geom_sf(aes(fill = EDEDUCWSEH), color = "grey35", linewidth = 0.1) +
  scale_fill_gradient(name = "Education",
                      low = fertility_fill_low, high = fertility_fill_high) +
  labs(title = "Education") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

## 3c) Choropleth of migration
map_migration <- ggplot(tz) +
  geom_sf(aes(fill = AHMIGRWEMP), color = "grey35", linewidth = 0.1) +
  scale_fill_gradient(name = "Migration",
                      low = fertility_fill_low, high = fertility_fill_high) +
  labs(title = "Migration") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

## 3d) Combine scatter and maps into one figure
model_error_map <- scatter | map_migration | map_education
ggsave("model_error_map.pdf", plot = model_error_map,
       path = "plots", width = 6, height = 3)

# Convert sf object to plain dataframe for modelling
tz_df <- tz %>% st_drop_geometry()

# --------------------------------------------------------------------
# 4) Baseline OLS regression & residual spatial autocorrelation
# --------------------------------------------------------------------
## 4a) Fit simple OLS
ols <- lm(AHMIGRWEMP ~ EDEDUCWSEH, data = tz_df)
print(summary(ols))

## 4b) Plot OLS line on scatter
linear_regression <- ggplot(data = tz) +
  geom_point(aes(x = EDEDUCWSEH, y = AHMIGRWEMP)) +
  geom_abline(intercept = 11.95, slope = 0.001) +  # fitted line
  labs(x = "Education", y = "Migration") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))
linear_regression

## 4c) Build queen contiguity weights
nb_q <- poly2nb(tz, queen = TRUE)
lw_q <- nb2listw(nb_q, style = "W", zero.policy = TRUE)

## 4d) Moran's I on OLS residuals
ols_resid_moran <- moran.test(residuals(ols), lw_q, zero.policy = TRUE)
print(ols_resid_moran)

## 4e) Store residuals in sf object
tz$residual_ols <- residuals(ols)

## 4f) Plot residuals on scatter with colour scale
linear_regression <- ggplot(data = tz) +
  geom_point(aes(x = EDEDUCWSEH, y = AHMIGRWEMP, fill = residual_ols),
             shape = 21, color = "black", stroke = 0.5, size = 3) +
  geom_abline(intercept = 11.95, slope = 0.001) +
  labs(x = "Education", y = "Migration") +
  scale_fill_gradient2(name = "Residual OLS",
                       low = "#2b83ba", mid = "white", high = "#d7191c",
                       midpoint = 0) +
  theme_minimal() +
  labs(title = "Linear Regression") +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5))
linear_regression

## 4g) Residual map
residual_map <- ggplot(tz) +
  geom_sf(aes(fill = residual_ols), color = "grey35", linewidth = 0.1) +
  scale_fill_gradient2(name = "Residual Map",
                       low = "#2b83ba", mid = "white", high = "#d7191c",
                       midpoint = 0) +
  labs(title = "Residual OLS") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))
residual_map

## 4h) Combine regression and residual map
model_error_ols <- linear_regression | residual_map
ggsave("model_error_ols.pdf", plot = model_error_ols,
       path = "plots", width = 6, height = 3)

# --------------------------------------------------------------------
# 5) Spatial error model (SEM)
# --------------------------------------------------------------------
# y = Xβ + u ; u = λWu + ε
sar_err <- errorsarlm(AHMIGRWEMP ~ EDEDUCWSEH,
                      data = tz,
                      listw = lw_q,
                      zero.policy = TRUE,
                      method = "eigen")
print(summary(sar_err, Nagelkerke = TRUE))