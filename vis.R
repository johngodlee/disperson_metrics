# Visualise and compare various metrics of dispersion
# John Godlee (johngodlee@gmail.com)
# 2022-03-17

# Packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scico)
library(sf)
library(GGally)
library(deldir)
library(ggnewscale)

source("functions.R")

# Import data
even_csr_metrics <- readRDS("./dat/even_csr_metrics.rds")
csr_clust_metrics <- readRDS("./dat/csr_clust_metrics.rds")
bicuar_metrics <- readRDS("./dat/bicuar_metrics.rds")
even_csr_polys <- readRDS("dat/even_csr_polys.rds")
even_csr_points <- readRDS("dat/even_csr_points.rds")
csr_clust_polys <- readRDS("dat/csr_clust_polys.rds")
csr_clust_points <- readRDS("dat/csr_clust_points.rds")

# Compare metrics for evenness to csr
even_csr_gather <- even_csr_metrics %>% 
  gather(key, value, -repl, -adj)

## Line plots of metrics vs. adj
even_csr_comp <- ggplot() + 
  geom_path(data = even_csr_gather, 
    aes(x = adj, y = value, group = repl)) + 
  facet_wrap(~key, ncol = 1, scales = "free_y") +
  mytheme()
ggsave(even_csr_comp, width = 5, height = 15, file = "img/even_csr_comp.png")

## Pairwise plots
even_csr_pairs <- ggpairs(even_csr_metrics,
  columns = c("cell_area_dev", "cell_area_cv", "pp_cv", "elong_cv", 
    "wi_mean", "wi_cv", "neighb_cv", "centre_gravity_cv"))
ggsave(even_csr_pairs, width = 10, height = 10, file = "img/even_csr_pairs.png")

# Compare metrics for csr to clustered
csr_clust_gather <- csr_clust_metrics %>% 
  gather(key, value, -repl, -adj)

## Line plots of metrics vs. adj
csr_clust_comp <- ggplot() + 
  geom_path(data = csr_clust_gather, 
    aes(x = adj, y = value, group = repl)) + 
  facet_wrap(~key, ncol = 1, scales = "free_y") + 
  mytheme()
ggsave(csr_clust_comp, width = 5, height = 15, file = "img/csr_clust_comp.png")

## Pairwise plots
csr_clust_pairs <- ggpairs(csr_clust_metrics,
  columns = c("cell_area_dev", "cell_area_cv", "pp_cv", "elong_cv", 
    "wi_mean", "wi_cv", "neighb_cv", "centre_gravity_cv"))
ggsave(csr_clust_pairs, width = 10, height = 10, file = "img/csr_clust_pairs.png")

# Make pretty plots showing variation in spatial distribution of stems
adj_seq <- seq(1,101,25)
polys_fil <- even_csr_polys[[1]][adj_seq]
points_fil <- even_csr_points[[1]][adj_seq]

polys_fil_sf <- do.call(rbind, lapply(seq_along(polys_fil), function(x) {
  polys_sf <- st_sf(st_sfc(lapply(polys_fil[[x]], function(y) { 
    st_polygon(list(rbind(y, y[1,])))
  })))
  polys_sf$cell_area <- st_area(polys_sf)
  polys_sf$adj <- adj_seq[x] - 1
  return(polys_sf)
}))

points_fil_comb <- do.call(rbind, lapply(seq_along(points_fil), function(x) {
  points_fil[[x]]$adj <- adj_seq[x] - 1
  return(points_fil[[x]])
    }))

voronoi_maps <- ggplot() + 
  geom_sf(data = polys_fil_sf, aes(fill = cell_area), 
    colour = "black", size = 0.2) + 
  scale_fill_scico(palette = "bamako", name = "Cell area") + 
  geom_point(data = points_fil_comb, aes(x = x, y = y),
    fill = "darkgrey", shape = 21, stroke = 0.2, size = 0.5) + 
  facet_wrap(~adj, nrow = 1) + 
  theme_bw() + 
  mytheme() + 
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")

cell_area_plot <- ggplot() + 
  geom_line(data = even_csr_metrics, 
    aes(x = adj, y = cell_area_cv, group = repl)) + 
  geom_vline(xintercept = adj_seq - 1, colour = "red", linetype = 2) +
  theme_bw() +
  labs(x = "N substitutions", y = expression(CV~sqrt("A")~(m^2)))

cell_vor_plot <- cell_area_plot + voronoi_maps + 
  plot_layout(ncol = 1, heights = c(2,1))
ggsave(cell_vor_plot, width = 8, height = 6, file = "img/cell_vor_plot.png")

# Investigate important transitions
point_distrib_norm_anom_repl <- even_csr_metrics %>% 
  filter(point_distrib_norm > 60) %>% 
  pull(repl) %>%
  unique()

inf_point <- even_csr_metrics %>%
  filter(repl == point_distrib_norm_anom_repl) %>% 
  dplyr::select(adj, point_distrib_norm) %>% 
  mutate(
    val_lag = lag(point_distrib_norm),
    lag_sum = point_distrib_norm - val_lag) %>% 
  arrange(desc(abs(lag_sum))) %>% 
  pull(adj) %>% 
  head(1)
inf_vec <- c(inf_point, inf_point+1)

point_distrib_norm_anom_polys <- even_csr_polys[[point_distrib_norm_anom_repl]][inf_vec]
point_distrib_norm_anom_points <- lapply(
  even_csr_points[[point_distrib_norm_anom_repl]][inf_vec], function(x) { 
    x$pt <- seq_len(nrow(x))
    return(x)
  })

point_distrib_norm_anom_polys_sf <- lapply(
  seq_along(point_distrib_norm_anom_polys), function(x) { 
    out <- st_sf(st_sfc(lapply(point_distrib_norm_anom_polys[[x]], function(y) {
      st_polygon(list(rbind(y, y[1,])))
    })))
    out$cell_area <- st_area(out)
    out$adj <- inf_vec[x]
    out$poly <- seq_len(nrow(out))
    return(out)
  })

pt_move <- which(!unlist(lapply(seq_len(nrow(point_distrib_norm_anom_points[[1]])), function(x) {
  all(point_distrib_norm_anom_points[[1]][x,] == point_distrib_norm_anom_points[[2]][x,])
  })))

point_distrib_norm_anom_points_move <- lapply(point_distrib_norm_anom_points, function(x) { 
  x$move <- ifelse(x$pt == pt_move, TRUE, FALSE)
  return(x)
  })

point_move_plot <- ggplot() + 
  geom_sf(data = point_distrib_norm_anom_polys_sf[[1]], 
    colour = "black", size = 2, fill = NA) + 
  geom_sf(data = point_distrib_norm_anom_polys_sf[[2]], 
    colour = "orange", fill = NA) + 
  geom_point(data = point_distrib_norm_anom_points_move[[2]], 
    aes(x = x, y = y, shape = move, fill = move), size = 4) +
  geom_point(data = point_distrib_norm_anom_points_move[[1]], 
    aes(x = x, y = y, shape = move, fill = move), size = 2) + 
  scale_shape_manual(values = c(21,22)) + 
  scale_fill_manual(values = c("lightgrey", "cyan3")) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  labs(x = "", y = "")
ggsave(point_move_plot, width = 6, height = 6, file = "./img/point_move_plot.png")

# Compare metrics in Bicuar 
bicuar_pairs <- ggpairs(bicuar_metrics,
  columns = c("cell_area_dev", "cell_area_cv", "pp_cv", "elong_cv", "wi_mean", "wi_cv", "neighb_cv", "centre_gravity_cv"),
  lower = list(continuous = "smooth"))
ggsave(bicuar_pairs, width = 10, height = 10, file = "img/bicuar_pairs.png")

# Three panel image of Evenness, CSR and clustering
even_ex <- even_csr_points[[1]][[1]]
even_ex$group <- "1: Even"
csr_ex <- even_csr_points[[1]][[101]]
csr_ex$group <- "2: Random"
clst_ex <- csr_clust_points[[1]][[50]]
clst_ex$group <- "3: Clustered"

pan_ex <- rbind(even_ex, csr_ex, clst_ex)

csr_diag <- ggplot() + 
  geom_point(data = pan_ex, aes(x = x, y = y )) + 
  facet_wrap(~group) + 
  coord_equal() + 
  mytheme() + 
  theme(axis.title = element_blank())
ggsave(csr_diag, width = 8, height = 3.5, file = "img/csr_diag.png")

# Changing centre of gravity of Voronoi cells with deviation from evenness
g1 <- data.frame(
  x = rep(c(10,30,50,70,90), times = 5),
  y = rep(c(10,30,50,70,90), each = 5))

g2 <- g1 %>% 
  mutate(
    x = case_when(
      x == 30 ~ 40,
      x == 70 ~ 60,
      TRUE ~ x),
    y = case_when(
      y == 30 ~ 40,
      y == 70 ~ 60,
      TRUE ~ y))

g1_vor <- tile.list(deldir(g1, rw = c(0,100,0,100)))
g2_vor <- tile.list(deldir(g2, rw = c(0,100,0,100)))

g1_polys <- lapply(g1_vor, function(z) { cbind(z$x, z$y) })
g2_polys <- lapply(g2_vor, function(z) { cbind(z$x, z$y) })

g1_polys_sf <- st_sf(st_sfc(lapply(g1_polys, function(x) {
      st_polygon(list(rbind(x, x[1,])))
})))
g1_centroid <- st_centroid(g1_polys_sf)

g2_polys_sf <- st_sf(st_sfc(lapply(g2_polys, function(x) {
      st_polygon(list(rbind(x, x[1,])))
})))
g2_centroid <- st_centroid(g2_polys_sf)

g1_plot <- ggplot() + 
  geom_sf(data = g1_polys_sf) + 
  geom_sf(data = g1_centroid, colour = "red", size = 4) + 
  geom_point(data = g1, aes(x = x, y = y)) + 
  coord_equal() + 
  mytheme() + 
  theme(axis.title = element_blank())
g2_plot <- ggplot() + 
  geom_sf(data = g2_polys_sf) + 
  geom_sf(data = g2_centroid, colour = "red", size = 4) + 
  geom_point(data = g2, aes(x = x, y = y)) + 
  mytheme() + 
  coord_equal() + 
  theme(axis.title = element_blank())

g_plot <- g1_plot + g2_plot + 
  plot_layout(ncol = 2)
ggsave(g_plot, width = 7.5, height = 3.5, file = "img/g_plot.png")

# Minimum bounding box 
polyg_ex <- even_csr_polys[[1]][[101]][[92]]
min_box <- minBox(polyg_ex)
min_box_sf <- st_sf(st_sfc(st_polygon(list(rbind(min_box$pts, min_box$pts[1,])))))
polyg_ex_sf <- st_sf(st_sfc(st_polygon(list(rbind(polyg_ex, polyg_ex[1,])))))

min_box_plot <- ggplot() + 
  geom_sf(data = polyg_ex_sf, fill = NA) + 
  geom_sf(data = min_box_sf, fill = NA, colour = "red", linetype = 2) + 
  coord_sf() + 
  mytheme() + 
  theme(axis.title = element_blank())
ggsave(min_box_plot, width = 5, height = 5, file = "img/min_box_plot.png")
