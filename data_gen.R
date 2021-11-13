# Create data for Voronoi tesselation test
# John Godlee (johngodlee@gmail.com)
# 2021-08-30

# Packages
library(parallel)
library(sf)
library(nngeo)
library(dplyr)
library(ggplot2)
library(scico)
library(patchwork)
library(maptools)
library(spatstat)

source("functions.R")

# Create a grid of points
# 2.5 m buffer inside plot
# 5 m between points
# 100x100 m grid
xy_vec <- seq(2.5, 100, 5)
dat <- expand.grid(xy_vec, xy_vec)
names(dat) <- c("x", "y")
dat_sf <- st_as_sf(dat, coords = c("x", "y"))
p <- st_polygon(list(matrix(c(0,0,100,100,0, 0,100,100,0,0), ncol = 2)))

saveRDS(p, "./dat/plot_poly.rds")

# Increase from evenness towards complete spatial randomness

# Number of replicates
reps <- 100

# Replicate grid of points
dat_list <- replicate(reps, dat_sf, simplify = FALSE)

# Possible coordinates for movement
c_repls <- seq(0,100,0.1)

# # For each plot, sequentially randomly move a point
# dat_even_csr_points <- mclapply(seq_along(dat_list), function(x) {
#   dat <- dat_list[[x]]
#   dat$adj <- 0
#   dat <- list(dat)
#   ord <- sample(seq_len(nrow(dat[[1]])), nrow(dat[[1]]))
#   for (i in seq_len(nrow(dat[[1]]))) {
#     message(x, ":", i)
#     dat[[i + 1]] <- dat[[i]]
#     st_geometry(dat[[i + 1]][ord[i],]) <- st_sfc(st_point(sample(c_repls, 2)))
#     dat[[i + 1]]$adj <- i
#   }
#   return(dat)
# }, mc.cores = 4)
# 
# # Write data
# saveRDS(dat_even_csr_points, "dat/dat_even_csr_points.rds")
dat_even_csr_points <- readRDS("dat/dat_even_csr_points.rds")
    
# # Voronoi tessellation of points
# dat_even_csr_polys <- mclapply(seq_along(dat_even_csr_points), function(x) {
#   lapply(seq_along(dat_even_csr_points[[x]]), function(y) {
#     message(x, ":", y)
#     cds <- st_coordinates(dat_even_csr_points[[x]][[y]])
#     cds_ppp <- ppp(
#       x = cds[,1],
#       y = cds[,2],
#       window = owin(c(0,100), c(0,100)))
#     vor <- dirichlet(cds_ppp) # Dirichlet tesselation
#     vor_sf <- st_as_sf(as(vor, "SpatialPolygons"))
#     vor_sf$poly_id <- seq_len(nrow(vor_sf))
#     return(vor_sf)
#   })
# }, mc.cores = 4)
# 
# # Write data
# saveRDS(dat_even_csr_polys, "dat/dat_even_csr_polys.rds")
dat_even_csr_polys <- readRDS("dat/dat_even_csr_polys.rds")

# Test plot
polys_plot <- dat_even_csr_polys[[1]][[400]] %>%
  mutate(poly_area = st_area(.))
points_plot <- dat_even_csr_points[[1]][[400]]
pdf(file = "img/dat_even_plot.pdf", width = 8, height = 8)
ggplot() + 
  geom_sf(data = polys_plot, 
    aes(fill = poly_area), colour = "black") +
  scale_fill_scico(name = "Area", palette = "bamako") + 
  geom_sf(data = points_plot, 
    shape = 21, fill = "darkgrey") + 
  mytheme()
dev.off()

# Increase from complete spatial randomness to clustered
# Move random individual closer to nearest cluster point by half its distance
xy_vec <- seq(0, 100, 1)
dat <- expand.grid(xy_vec, xy_vec)
names(dat) <- c("x", "y")
clust_points <- st_as_sf(dat, coords = c("x", "y"))

# Number of replicates
reps <- 100

# Replicate plots
csr_dat_list <- replicate(reps, {
  st_as_sf(dat[sample(seq_len(nrow(dat)), 100),], coords = c("x", "y"))
    }, simplify = FALSE)

# Number of points to cluster around
n_clust_points <- 5

# # For each plot, sequentially move points closer to nearest cluster points
# dat_csr_clust_points <- mclapply(seq_along(csr_dat_list), function(x) {
#   dat <- csr_dat_list[[x]]
#   dat$adj <- 0
#   dat <- list(dat)
#   cluster_points <- clust_points[sample(seq_len(nrow(clust_points)), n_clust_points),]
#   for (i in seq_len(500)) {
#     message(x, ":", i)
#     dat[[i + 1]] <- dat[[i]]
#     chosen_row <- sample(seq_len(nrow(dat[[i]])), 1)
#     chosen <- dat[[i + 1]][chosen_row,]
#     nn <- cluster_points[st_nn(chosen, cluster_points, k = 1)[[1]],]
#     chosen_nn_line <- st_cast(st_combine(c(st_geometry(chosen), st_geometry(nn))),
#       "LINESTRING")
#     midpoint <- st_cast(st_line_sample(chosen_nn_line, 1), "POINT")
#     st_geometry(dat[[i + 1]][chosen_row,]) <- st_geometry(midpoint)
#     dat[[i + 1]]$adj <- i
#   }
#   return(dat)
# }, mc.cores = 4)
# 
# # Write data
# saveRDS(dat_csr_clust_points, "dat/dat_csr_clust_points.rds")
dat_csr_clust_points <- readRDS("dat/dat_csr_clust_points.rds")
    
# # Voronoi tessellation
# dat_csr_clust_polys <- mclapply(seq_along(dat_csr_clust_points), function(x) {
#   lapply(seq_along(dat_csr_clust_points[[x]]), function(y) {
#     message(x, ":", y)
#     cds <- st_coordinates(dat_csr_clust_points[[x]][[y]])
#     cds_ppp <- ppp(
#       x = cds[,1],
#       y = cds[,2],
#       window = owin(c(0,100), c(0,100)))
#     vor <- dirichlet(cds_ppp) # Dirichlet tesselation
#     vor_sf <- st_as_sf(as(vor, "SpatialPolygons"))
#     vor_sf$poly_id <- seq_len(nrow(vor_sf))
#     return(vor_sf)
#   })
# }, mc.cores = 4)
# 
# # Write data
# saveRDS(dat_csr_clust_polys, "dat/dat_csr_clust_polys.rds")
dat_csr_clust_polys <- readRDS("dat/dat_csr_clust_polys.rds")

# Test plot
pdf(file = "img/dat_clust_plot.pdf", width = 10, height = 8)
wrap_plots(lapply(seq(1,length(dat_csr_clust_points[[1]]), 100), function(x) {
  polys_plot <- dat_csr_clust_polys[[1]][[x]] %>%
    mutate(poly_area = st_area(.))
  points_plot <- dat_csr_clust_points[[1]][[x]]
  ggplot() + 
    geom_sf(data = polys_plot, 
      aes(fill = poly_area), colour = "black") +
    scale_fill_scico(name = "Area", palette = "bamako", limits = c(0,1200)) + 
    geom_sf(data = points_plot, 
      shape = 21, fill = "darkgrey") + 
    theme_bw()
})) +
  plot_layout(guides = "collect")
dev.off()

