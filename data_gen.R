# Create data for Voronoi tesselation test
# John Godlee (johngodlee@gmail.com)
# 2021-08-30

# Packages
library(parallel)
library(sf)
library(dplyr)
library(ggplot2)
library(maptools)
library(spatstat)

# Create a grid of points
# 2.5 m buffer inside plot
# 5 m between points
# 100x100 m grid
xy_vec <- seq(2.5, 100, 5)
dat <- expand.grid(xy_vec, xy_vec)
names(dat) <- c("x", "y")
dat_sf <- st_as_sf(dat, coords = c("x", "y"))
p <- st_polygon(list(matrix(c(0,0,100,100,0, 0,100,100,0,0), ncol = 2)))

# Voronoi tessellation test
cds <- st_coordinates(dat_sf)
cds_ppp <- ppp(
  x = cds[,1],
  y = cds[,2],
  window = owin(c(0,100), c(0,100)))
vor <- dirichlet(cds_ppp) # Dirichlet tesselation
vor_sf <- st_as_sf(as(vor, "SpatialPolygons"))
vor_sf$poly_id <- seq_len(nrow(vor_sf))

# Test plot 
ggplot() + 
  geom_sf(data = p) + 
  geom_sf(data = vor_sf, aes(fill = poly_id), colour = "black") + 
  geom_sf(data = dat_sf) +
  theme_bw()

# Increase from evenness towards complete spatial randomness
reps <- 100
dat_list <- replicate(reps, dat_sf, simplify = FALSE)

c_repls <- seq(0,100,0.1)

dat_even_csr_points <- mclapply(seq_along(dat_list), function(x) {
  dat <- dat_list[[x]]
  dat$adj <- 0
  dat <- list(dat)
  ord <- sample(seq_len(nrow(dat[[1]])), nrow(dat[[1]]))
  for (i in seq_len(nrow(dat[[1]]))) {
    message(x, ":", i)
    dat[[i + 1]] <- dat[[i]]
    st_geometry(dat[[i + 1]][ord[i],]) <- st_sfc(st_point(sample(c_repls, 2)))
    dat[[i + 1]]$adj <- i
  }
  return(dat)
}, mc.cores = 4)

saveRDS(dat_even_csr_points, "dat/dat_even_csr_points.rds")
dat_even_csr_points <- readRDS("dat/dat_even_csr_points.rds")
    
dat_even_csr_polys <- mclapply(seq_along(dat_even_csr_points), function(x) {
  lapply(seq_along(dat_even_csr_points[[x]]), function(y) {
    message(x, ":", y)
    cds <- st_coordinates(dat_even_csr_points[[x]][[y]])
    cds_ppp <- ppp(
      x = cds[,1],
      y = cds[,2],
      window = owin(c(0,100), c(0,100)))
    vor <- dirichlet(cds_ppp) # Dirichlet tesselation
    vor_sf <- st_as_sf(as(vor, "SpatialPolygons"))
    vor_sf$poly_id <- seq_len(nrow(vor_sf))
    return(vor_sf)
  })
}, mc.cores = 4)

saveRDS(dat_even_csr_polys, "dat/dat_even_csr_polys.rds")
dat_even_csr_polys <- readRDS("dat/dat_even_csr_polys.rds")

# Test plot
ggplot() + 
  geom_sf(data = dat_even_csr_polys[[1]][[400]], 
    aes(fill = poly_id), colour = "black") +
  geom_sf(data = dat_even_csr_points[[1]][[400]]) 

# Increase from complete spatial randomness to clustered

