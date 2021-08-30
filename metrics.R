# Calculate various metrics on voronoi cells
# John Godlee (johngodlee@gmail.com)
# 2021-08-23

# Packages
library(sf)
library(nngeo)
library(dplyr)

# Import data
dat_even_csr_polys <- readRDS("dat/dat_even_csr_polys.rds")
dat_even_csr_points <- readRDS("dat/dat_even_csr_points.rds")

mclapply(seq_along(dat_even_csr_polys), function(x) {
  lapply(seq_along(dat_even_csr_polys[[x]]), function(y) {
    apply(seq_len(nrow(dat_even_csr_polys[[x]][[y]])), function(z) {
      polyg <- dat_even_csr_polys[[x]][[y]][z,]
      vertex <- st_cast(polyg, "POINT")
      point <- dat_even_csr_points[[x]][[y]][unlist(st_intersects(polyg, dat_even_csr_points[[x]][[y]])),]

      # Distance of point to furthest vertex
      vertex_dists <- st_distance(point, vertex)

      # Distance between point and centre of gravity of cell
      centre_gravity_dist <- st_distance(st_centroid(polyg), point)

      # Cell area
      cell_area <- st_area(polyg)

      # Cell elongation (Polsby-Popper)
      pp <- (4*pi*cell_area)/(st_length(st_cast(polyg, "LINESTRING"))^2)
      
      # Cell elongation (minimum bounding rectangle ratio)

    })

    # Nearest neighbour distances
    neighb_dist <- unlist(lapply(st_nn(dat_even_csr_points[[x]][[y]], dat_even_csr_points[[x]][[y]], 
          k = 2, returnDist = TRUE)[[2]], 
        "[[", 2))
  })
}, mc.cores = 4)
