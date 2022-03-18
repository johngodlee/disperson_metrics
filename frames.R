# Generate frames for animated random sampling of even grid
# John Godlee (johngodlee@gmail.com)
# 2021-10-31

# Packages
library(ggplot2)
library(sf)
library(scico)
library(parallel)

# Even-ness to CSR 
## Import data
even_csr_points <- readRDS("dat/even_csr_points.rds")
even_csr_polys <- readRDS("dat/even_csr_polys.rds")

## Isolate one example plot
even_csr_polys_fil <- even_csr_polys[[1]]
even_csr_points_fil <- even_csr_points[[1]]

## Calculate areas
even_csr_polys_fil_sf <- lapply(seq_along(even_csr_polys_fil), function(x) {
  polys_sf <- st_sf(st_sfc(lapply(even_csr_polys_fil[[x]], function(y) { 
    st_polygon(list(rbind(y, y[1,])))
  })))
  polys_sf$cell_area <- st_area(polys_sf)
  return(polys_sf)
})

col_range <- range(unlist(lapply(even_csr_polys_fil_sf, "[[", "cell_area")))

## Create all plots
lapply(seq_along(even_csr_polys_fil_sf), function(x) {
    out <- ggplot() + 
      geom_sf(data = even_csr_polys_fil_sf[[x]], aes(fill = cell_area), 
        colour = "black") + 
      scale_fill_scico(palette = "bamako", limits = col_range) + 
      geom_point(data = even_csr_points_fil[[x]], aes(x = x, y = y), 
        colour = "black") +
      theme_void() + 
      theme(legend.position = "none")
    ggsave(out, width = 5, height = 5, 
      file = file.path("frames", "even_csr", 
        paste0(sprintf("%03d", x), "_even_csr", ".png")))
})

# CSR to clustered
## Import data
csr_clust_points <- readRDS("dat/csr_clust_points.rds")
csr_clust_polys <- readRDS("dat/csr_clust_polys.rds")

## Isolate one example plot
csr_clust_polys_fil <- csr_clust_polys[[1]]
csr_clust_points_fil <- csr_clust_points[[1]]

## Calculate areas
csr_clust_polys_fil_sf <- lapply(seq_along(csr_clust_polys_fil), function(x) {
  polys_sf <- st_sf(st_sfc(lapply(csr_clust_polys_fil[[x]], function(y) { 
    st_polygon(list(rbind(y, y[1,])))
  })))
  polys_sf$cell_area <- st_area(polys_sf)
  return(polys_sf)
})

col_range <- range(unlist(lapply(csr_clust_polys_fil_sf, "[[", "cell_area")))

## Create all plots
lapply(seq_along(csr_clust_polys_fil_sf), function(x) {
    out <- ggplot() + 
      geom_sf(data = csr_clust_polys_fil_sf[[x]], aes(fill = cell_area), 
        colour = "black") + 
      scale_fill_scico(palette = "bamako", limits = col_range) + 
      geom_point(data = csr_clust_points_fil[[x]], aes(x = x, y = y), 
        colour = "black") +
      theme_void() + 
      theme(legend.position = "none")
    ggsave(out, width = 5, height = 5, 
      file = file.path("frames", "csr_clust", 
        paste0(sprintf("%03d", x), "_csr_clust", ".png")))
})
