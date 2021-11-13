# Generate frames for animated random sampling of even grid
# John Godlee (johngodlee@gmail.com)
# 2021-10-31

# Packages
library(ggplot2)
library(sf)
library(scico)
library(mclapply)

# Even-ness to CSR 
## Import data
dat_even_csr_points_list <- readRDS("dat/dat_even_csr_points.rds")
dat_even_csr_polys_list <- readRDS("dat/dat_even_csr_polys.rds")

## Isolate one example plot
dat_even_csr_polys <- dat_even_csr_polys_list[[1]]
dat_even_csr_points <- dat_even_csr_points_list[[1]]

## Calculate areas
dat_even_csr_polys <- lapply(dat_even_csr_polys, function(x) {
  x$area <- st_area(x)
  return(x) 
})

col_range <- range(unlist(lapply(dat_even_csr_polys, "[[", "area")))

## Create all plots
mclapply(seq_along(dat_even_csr_polys), function(x) {
  png(file = file.path("frames", "even_csr", 
      paste0(sprintf("%03d", x), "_even_csr", ".png")), 
    width = 500, height = 500)
    print(ggplot() + 
      geom_sf(data = dat_even_csr_polys[[x]], aes(fill = area), 
        colour = "black") + 
      scale_fill_scico(palette = "bamako", limits = col_range) + 
      geom_sf(data = dat_even_csr_points[[x]], colour = "black") +
      theme_void() + 
      theme(legend.position = "none"))
  dev.off()
}, mc.cores = 4)

# CSR to clustered
## Import data
dat_csr_clust_points_list <- readRDS("dat/dat_csr_clust_points.rds")
dat_csr_clust_polys_list <- readRDS("dat/dat_csr_clust_polys.rds")

## Isolate one example plot
dat_csr_clust_polys <- dat_csr_clust_polys_list[[1]]
dat_csr_clust_points <- dat_csr_clust_points_list[[1]]

## Calculate areas
dat_csr_clust_polys <- lapply(dat_csr_clust_polys, function(x) {
  x$area <- st_area(x)
  return(x) 
})

col_range <- range(unlist(lapply(dat_csr_clust_polys, "[[", "area")))

## Create all plots
mclapply(seq_along(dat_csr_clust_polys), function(x) {
  png(file = file.path("frames", "csr_clust", 
      paste0(sprintf("%03d", x), "_csr_clust", ".png")), 
    width = 500, height = 500)
    print(ggplot() + 
      geom_sf(data = dat_csr_clust_polys[[x]], aes(fill = area), 
        colour = "black") + 
      scale_fill_scico(palette = "bamako", limits = col_range) + 
      geom_sf(data = dat_csr_clust_points[[x]], colour = "black") +
      theme_void() + 
      theme(legend.position = "none"))
  dev.off()
}, mc.cores = 4)


