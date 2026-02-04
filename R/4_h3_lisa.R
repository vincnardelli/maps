# -------------------------------------------------------------
#  MAPS Summer School – Spatial Analysis
#  Module: H3 Grid + LISA Analysis from TIF data
#  Teacher: Vincenzo Nardelli
# -------------------------------------------------------------

# -------------------------------------------------------------
#  0) Load required packages
# -------------------------------------------------------------
library(terra)       # raster data handling
library(sf)          # spatial vector handling
library(h3jsr)       # H3 hexagonal grid system
library(spdep)       # spatial dependence functions
library(dplyr)       # data manipulation
library(ggplot2)     # plotting
library(patchwork)   # combine plots

set.seed(123)

# -------------------------------------------------------------
#  1) Load TIF raster and create H3 grid
# -------------------------------------------------------------
raster_data <- rast("../data/TZ2022DHS_EDLITRWLIT_MS_v01/TZ2022DHS_EDLITRWLIT_MS_MEAN_v01.tif")

# Create H3 grid (resolution 5 = ~252 km² per hexagon)
bbox_wgs84 <- st_bbox(st_transform(st_as_sfc(st_bbox(raster_data)), 4326))
h3_indices <- polygon_to_cells(st_as_sfc(bbox_wgs84), res = 5, simple = TRUE)
h3_sf <- cell_to_polygon(h3_indices, simple = FALSE)

# Extract raster values to H3 hexagons
h3_sf_proj <- st_transform(h3_sf, st_crs(raster_data))
h3_sf$value <- extract(raster_data, vect(h3_sf_proj), fun = mean, na.rm = TRUE)[, 2]
h3 <- h3_sf %>% filter(!is.na(value) & is.finite(value))

# -------------------------------------------------------------
#  2) Exploratory visualisations (map + histogram)
# -------------------------------------------------------------
p_map <- ggplot(h3) +
  geom_sf(aes(fill = value), color = "grey35", linewidth = 0.1) +
  scale_fill_gradient(name = "Literacy rate",
                      low = "#fff7bc", high = "#d7301f") +
  labs(title = "Map") +
  theme_void() + theme(plot.title = element_text(hjust = 0.5))

p_dist <- ggplot(st_drop_geometry(h3), aes(x = value)) +
  geom_histogram(aes(fill = after_stat(x)),
                 bins = 30, alpha = 0.9, na.rm = TRUE) +
  scale_fill_gradient(low = "#fff7bc", high = "#d7301f", guide = "none") +
  labs(title = "Distribution",
       x = "Literacy rate", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank())

fig_dist_map <- p_dist + p_map + plot_layout(widths = c(0.7, 1.3))
ggsave("h3_map.pdf", plot = fig_dist_map, path = "plots", width = 6, height = 3)

# -------------------------------------------------------------
#  2b) Comparison: Original raster vs H3 aggregation
# -------------------------------------------------------------
# Convert raster to data frame for plotting
raster_df <- as.data.frame(raster_data, xy = TRUE)
colnames(raster_df) <- c("x", "y", "value")

# Original raster map
p_raster <- ggplot(raster_df) +
  geom_raster(aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(name = "Literacy rate",
                      low = "#fff7bc", high = "#d7301f",
                      na.value = "grey90") +
  labs(title = "Original raster") +
  coord_sf() +
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

# H3 grid map
p_h3 <- ggplot(h3) +
  geom_sf(aes(fill = value), color = "grey35", linewidth = 0.1) +
  scale_fill_gradient(name = "Literacy rate",
                      low = "#fff7bc", high = "#d7301f") +
  labs(title = "H3 aggregation") +
  theme_void() + 
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

fig_comparison <- p_raster + p_h3 + plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
ggsave("h3_comparison.pdf", plot = fig_comparison, path = "plots", width = 8, height = 4)

# -------------------------------------------------------------
#  3) Construct spatial weights – Queen contiguity
# -------------------------------------------------------------
nb_q  <- poly2nb(h3, queen = TRUE)
lw_q  <- nb2listw(nb_q, style = "W", zero.policy = TRUE)

# Filter to largest connected component
comp <- n.comp.nb(nb_q)
if(comp$nc > 1) {
  largest_comp <- as.integer(names(which.max(table(comp$comp.id))))
  h3 <- h3[comp$comp.id == largest_comp, ]
  nb_q <- poly2nb(h3, queen = TRUE)
  lw_q <- nb2listw(nb_q, style = "W", zero.policy = TRUE)
}

# -------------------------------------------------------------
#  4) Visualise the spatial network
# -------------------------------------------------------------
coords <- st_coordinates(st_centroid(st_geometry(h3)))

edges <- do.call(rbind, lapply(seq_along(nb_q), function(i){
  js <- nb_q[[i]]; js <- js[js > i]
  if(length(js) == 0) return(NULL)
  data.frame(from = i, to = js)
}))

if(!is.null(edges) && nrow(edges) > 0){
  edges_geom <- st_sfc(lapply(seq_len(nrow(edges)), function(k){
    i <- edges$from[k]; j <- edges$to[k]
    st_linestring(rbind(coords[i,], coords[j,]))
  }), crs = st_crs(h3))
  edges_sf <- st_sf(edges, geometry = edges_geom)
} else {
  edges_sf <- st_sf(from = integer(0), to = integer(0),
                    geometry = st_sfc(crs = st_crs(h3)))
}

nodes_sf <- st_as_sf(data.frame(id = seq_len(nrow(coords)),
                                x = coords[,1], y = coords[,2]),
                     coords = c("x","y"), crs = st_crs(h3))

p_net <- ggplot() +
  geom_sf(data = h3, fill = "grey95", color = "grey75", linewidth = 0.15) +
  geom_sf(data = edges_sf, color = "grey25", linewidth = 0.45, alpha = 0.85) +
  geom_sf(data = nodes_sf, color = "#d7301f", size = 1.7, alpha = 0.95) +
  labs(title = paste0("H3 spatial network (", nrow(h3), " hexagons)")) +
  theme_void() + theme(plot.title = element_text(hjust = 0.5))

ggsave("h3_network.pdf", plot = p_net, path = "plots", width = 6, height = 3)

# -------------------------------------------------------------
#  5) Global Moran's I
# -------------------------------------------------------------
global_moran_mc <- moran.mc(h3$value, lw_q, nsim = 999, zero.policy = TRUE)
print(global_moran_mc)

# -------------------------------------------------------------
#  6) Local Moran's I (LISA) and cluster classification
# -------------------------------------------------------------
lisa <- localmoran(h3$value, lw_q, zero.policy = TRUE)

h3$Ii   <- lisa[,"Ii"]
h3$Zi   <- lisa[,"Z.Ii"]
h3$Pr_z <- lisa[,"Pr(z != E(Ii))"]

# Standardise variables for LISA classification
z_y     <- as.numeric(scale(h3$value))
lag_z_y <- lag.listw(lw_q, z_y, zero.policy = TRUE)

# Define significant clusters (α = 0.05)
alpha <- 0.05
sig   <- is.finite(lisa[,"Pr(z != E(Ii))"]) & lisa[,"Pr(z != E(Ii))"] <= alpha

h3$lisa_cluster <- "Not significant"
h3$lisa_cluster[sig & z_y >= 0 & lag_z_y >= 0] <- "High-High"
h3$lisa_cluster[sig & z_y <= 0 & lag_z_y <= 0] <- "Low-Low"
h3$lisa_cluster[sig & z_y >= 0 & lag_z_y <= 0] <- "High-Low"
h3$lisa_cluster[sig & z_y <= 0 & lag_z_y >= 0] <- "Low-High"

# Plot LISA clusters
p2 <- ggplot(h3) +
  geom_sf(aes(fill = lisa_cluster), color = "grey35", linewidth = 0.2) +
  scale_fill_manual(values = c("High-High" = "#B2182B",
                               "Low-Low"   = "#2166AC",
                               "High-Low"  = "#EF8A62",
                               "Low-High"  = "#67A9CF",
                               "Not significant" = "grey85")) +
  labs(title = "Local Moran", fill = "LISA cluster") +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# -------------------------------------------------------------
#  7) Moran scatterplot
# -------------------------------------------------------------
moran_i <- global_moran_mc$statistic

df_scatter <- data.frame(
  z_y      = z_y,
  lag_z_y  = lag_z_y,
  cluster  = h3$lisa_cluster
)

p_scatter <- ggplot(df_scatter, aes(x = z_y, y = lag_z_y)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_point(aes(color = cluster), size = 2.4, alpha = 0.9) +
  geom_abline(intercept = 0, slope = moran_i, color = "grey25", linewidth = 0.7) +
  scale_color_manual(values = c("High-High" = "#B2182B",
                                "Low-Low"   = "#2166AC",
                                "High-Low"  = "#EF8A62",
                                "Low-High"  = "#67A9CF",
                                "Not significant" = "grey85")) +
  labs(title = "Moran scatterplot",
       x = "z(Literacy rate)",
       y = "W z(Literacy rate)",
       color = "LISA cluster") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# -------------------------------------------------------------
#  8) Combined visualisations (scatter + LISA map)
# -------------------------------------------------------------
fig_scatter_lisa <- p_scatter + p2 + plot_layout(widths = c(0.8, 1.2))
ggsave("h3_lisa.pdf", plot = fig_scatter_lisa, path = "plots", width = 6, height = 3)

# Print plots
print(fig_dist_map)
print(fig_comparison)
print(fig_scatter_lisa)
print(p_net)
