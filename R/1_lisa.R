# -------------------------------------------------------------
#  MAPS Summer School – Spatial Analysis
#  Module: Spatial Analysis – Teacher: Vincenzo Nardelli
# -------------------------------------------------------------

# -------------------------------------------------------------
#  0) Load required packages
# -------------------------------------------------------------
library(sf)          # spatial vector handling
library(spdep)       # spatial dependence functions
library(dplyr)       # data manipulation
library(ggplot2)     # plotting
library(patchwork)   # combine plots

set.seed(123)  # reproducibility

# -------------------------------------------------------------
#  1) Read polygons and quick data check
# -------------------------------------------------------------
tz <- st_read("../data/sdr_subnational_data_migration/shps/sdr_subnational_data_dhs_2022_lvl_2.shp",
              quiet = TRUE)

# -------------------------------------------------------------
#  2) Exploratory visualisations (map + histogram)
# -------------------------------------------------------------
# 2a) Choropleth of fertility rate
p_map <- ggplot(tz) +
  geom_sf(aes(fill = FEFRTRWTFR), color = "grey35", linewidth = 0.1) +
  scale_fill_gradient(name = "Fertility rate",
                      low = "#fff7bc", high = "#d7301f") +
  labs(title = "Map") +
  theme_void() + theme(plot.title = element_text(hjust = 0.5))

# 2b) Histogram of fertility rate
p_dist <- ggplot(st_drop_geometry(tz), aes(x = FEFRTRWTFR)) +
  geom_histogram(aes(fill = after_stat(x)),
                 bins = 5, alpha = 0.9, na.rm = TRUE) +
  scale_fill_gradient(low = "#fff7bc", high = "#d7301f", guide = "none") +
  labs(title = "Distribution",
       x = "Fertility rate", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank())

# Combine histogram and map side‑by‑side
fig_dist_map <- p_dist + p_map + plot_layout(widths = c(0.7, 1.3))
ggsave("map.pdf", plot = fig_dist_map, path = "plots", width = 6, height = 3)

# -------------------------------------------------------------
#  3) Construct spatial weights – Queen contiguity
# -------------------------------------------------------------
# Queen adjacency (share edge or vertex)
nb_q  <- poly2nb(tz, queen = TRUE)
lw_q  <- nb2listw(nb_q, style = "W", zero.policy = TRUE)

# -------------------------------------------------------------
#  3b) Visualise the spatial network (nodes + edges)
# -------------------------------------------------------------
coords <- st_coordinates(st_centroid(st_geometry(tz)))

# Build edge list from neighbour object
edges <- do.call(rbind, lapply(seq_along(nb_q), function(i){
  js <- nb_q[[i]]; js <- js[js > i]   # avoid duplicate edges
  if(length(js) == 0) return(NULL)
  data.frame(from = i, to = js)
}))

# Create LINESTRING geometries for edges
if(!is.null(edges) && nrow(edges) > 0){
  edges_geom <- st_sfc(lapply(seq_len(nrow(edges)), function(k){
    i <- edges$from[k]; j <- edges$to[k]
    st_linestring(rbind(coords[i,], coords[j,]))
  }), crs = st_crs(tz))
  edges_sf <- st_sf(edges, geometry = edges_geom)
} else {
  edges_sf <- st_sf(from = integer(0), to = integer(0),
                    geometry = st_sfc(crs = st_crs(tz)))
}

# Node geometries
nodes_sf <- st_as_sf(data.frame(id = seq_len(nrow(coords)),
                                x = coords[,1], y = coords[,2]),
                     coords = c("x","y"), crs = st_crs(tz))

# Number of connected components
comp <- n.comp.nb(nb_q)

# Plot network overlaying the polygons
p_net <- ggplot() +
  geom_sf(data = tz, fill = "grey95", color = "grey75", linewidth = 0.15) +
  geom_sf(data = edges_sf, color = "grey25", linewidth = 0.45, alpha = 0.85) +
  geom_sf(data = nodes_sf, color = "#d7301f", size = 1.7, alpha = 0.95) +
  labs(title = paste0("Spatial weights network (queen contiguity) — components: ",
                      comp$nc)) +
  theme_void() + theme(plot.title = element_text(hjust = 0.5))

ggsave("network.pdf", plot = p_net, path = "plots", width = 6, height = 3)

# -------------------------------------------------------------
#  4) Global Moran's I (standard & Monte‑Carlo)
# -------------------------------------------------------------
y   <- tz$FEFRTRWTFR
idx <- is.finite(y)           # indices of non‑missing values
y_use <- y[idx]
lw_use <- if(all(idx)) lw_q else subset.listw(lw_q, idx, zero.policy = TRUE)

# Analytical test
global_moran <- moran.test(y_use, lw_use, zero.policy = TRUE,
                           alternative = "two.sided")
print(global_moran)

# Monte‑Carlo approximation
global_moran_mc <- moran.mc(y_use, lw_use, nsim = 999, zero.policy = TRUE)
print(global_moran_mc)

# -------------------------------------------------------------
#  5) Local Moran's I (LISA) and cluster classification
# -------------------------------------------------------------
lisa <- localmoran(y_use, lw_use, zero.policy = TRUE)

# Attach LISA statistics to the sf object
tz$Ii   <- NA; tz$Zi   <- NA; tz$Pr_z <- NA
tz$Ii[idx]   <- lisa[,"Ii"]
tz$Zi[idx]   <- lisa[,"Z.Ii"]
tz$Pr_z[idx] <- lisa[,"Pr(z != E(Ii))"]

# Standardise variables for LISA classification
z_y       <- as.numeric(scale(y_use))
lag_z_y   <- lag.listw(lw_use, z_y, zero.policy = TRUE)

# Define significant clusters (α = 0.05)
alpha <- 0.05
tz$lisa_cluster <- "Missing"
tz$lisa_cluster[idx] <- "Not significant"

sig       <- is.finite(lisa[,"Pr(z != E(Ii))"]) & lisa[,"Pr(z != E(Ii))"] <= alpha
cluster_use <- rep("Not significant", length(z_y))
cluster_use[sig & z_y >= 0 & lag_z_y >= 0] <- "High-High"
cluster_use[sig & z_y <= 0 & lag_z_y <= 0] <- "Low-Low"
cluster_use[sig & z_y >= 0 & lag_z_y <= 0] <- "High-Low"
cluster_use[sig & z_y <= 0 & lag_z_y >= 0] <- "Low-High"
tz$lisa_cluster[idx] <- cluster_use

# Plot LISA clusters
p2 <- ggplot(tz) +
  geom_sf(aes(fill = lisa_cluster), color = "grey35", linewidth = 0.2) +
  scale_fill_manual(values = c("High-High" = "#B2182B",
                               "Low-Low"   = "#2166AC",
                               "High-Low"  = "#EF8A62",
                               "Low-High"  = "#67A9CF",
                               "Not significant" = "grey85",
                               "Missing" = "grey92")) +
  labs(title = "Local Moran", fill = "LISA cluster") +
  theme_void() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# -------------------------------------------------------------
#  6) Moran scatterplot (showing significant points)
# -------------------------------------------------------------
moran_i <- unname(global_moran$estimate[["Moran I statistic"]])

# Data for scatterplot (standardised)
df_scatter <- data.frame(
  y        = z_y,
  lag_z_y  = lag_z_y,
  p_local  = lisa[,"Pr(z != E(Ii))"],
  sig      = sig,
  cluster  = cluster_use,
  stringsAsFactors = FALSE
)

# Raw‑scale scatterplot for reference
df_scatter_first <- data.frame(
  y        = y_use,
  lag_z_y  = lag.listw(lw_use, y_use, zero.policy = TRUE)
)

p_scatter_first <- ggplot(df_scatter_first, aes(x = y, y = lag_z_y)) +
  geom_point(size = 2.4, alpha = 0.9) +
  labs(title = "Variable vs Lag Variable",
       x = "Fertility rate",
       y = "W Fertility rate") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# Map coloured by fertility rate (for layout purposes)
p_net_color <- ggplot() +
  geom_sf(data = tz, aes(fill = FEFRTRWTFR),
          color = "grey35", linewidth = 0.1) +
  scale_fill_gradient(name = "Fertility rate",
                      low = "#fff7bc", high = "#d7301f") +
  geom_sf(data = edges_sf, color = "grey25", linewidth = 0.45, alpha = 0.85) +
  geom_sf(data = nodes_sf, color = "grey25", size = 1.7, alpha = 0.95) +
  labs(title = "Map") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5))

# Combine network map and raw scatterplot
fig_lag <- p_net_color + p_scatter_first + plot_layout(widths = c(1,1))
ggsave("lag.pdf", plot = fig_lag, path = "plots", width = 6, height = 3)

# Standardised Moran scatterplot with LISA colours
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
       x = "z(Fertility rate)",
       y = "W z(Fertility rate)",
       color = "LISA cluster") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# -------------------------------------------------------------
#  7) Combined visualisations (scatter + LISA map)
# -------------------------------------------------------------
fig_scatter_lisa <- p_scatter + p2 + plot_layout(widths = c(0.8, 1.2))
ggsave("lisa.pdf", plot = fig_scatter_lisa, path = "plots", width = 6, height = 3)

# Print plots to the interactive device
print(fig_dist_map)
print(fig_scatter_lisa)
print(p_net)
