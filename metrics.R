# Calculate various metrics on voronoi cells
# John Godlee (johngodlee@gmail.com)
# 2021-08-23

# Packages
library(sf)
library(nngeo)
library(parallel)

source("functions.R")

# Import data
dat_even_csr_polys <- readRDS("dat/dat_even_csr_polys.rds")
dat_even_csr_points <- readRDS("dat/dat_even_csr_points.rds")
dat_csr_clust_polys <- readRDS("dat/dat_csr_clust_polys.rds")
dat_csr_clust_points <- readRDS("dat/dat_csr_clust_points.rds")

# For each plot
metric_gen <- function(dat_polys, dat_points) {
  metrics_df <- mclapply(seq_along(dat_polys), function(x) {
    # For each iteration of simulation
    out <- do.call(rbind, lapply(seq_along(dat_polys[[x]]), function(y) {
      message(x, "/", length(dat_polys), ", ", y, "/", length(dat_polys[[x]]))
      # For each polygon
      polyg_ex <- lapply(seq_len(nrow(dat_polys[[x]][[y]])), function(z) {

        # Extract the polygon
        polyg <- dat_polys[[x]][[y]][z,]

        # Find the vertex points of the polygon
        vertex <- st_cast(polyg, "POINT")

        # Find the corresponding point
        point <- dat_points[[x]][[y]][unlist(st_intersects(polyg, dat_points[[x]][[y]])),]

        # Distance of point to furthest vertex
        vertex_dists <- c(st_distance(point, vertex))

        # Distance between point and centre of gravity of cell
        centre_gravity_dist <- c(st_distance(st_centroid(polyg), point))

        # Cell area
        cell_area <- st_area(polyg)

        # Cell elongation (Polsby-Popper)
        pp <- (4*pi*cell_area)/(st_length(st_cast(polyg, "LINESTRING"))^2)
        
        # Cell elongation (minimum bounding rectangle ratio)      
        min_box <- unlist(minBox(polyg)[c("width", "height")])
        elong <- max(min_box) / min(min_box)

        return(list(vertex_dists, centre_gravity_dist, cell_area, pp, elong, z))
      })

      # Nearest neighbour distances
      neighb_dist <- unlist(lapply(st_nn(dat_points[[x]][[y]], dat_points[[x]][[y]], 
            k = 2, returnDist = TRUE)[[2]], 
          "[[", 2))

      # Winkelmass
      point_coords <- as.data.frame(st_coordinates(dat_points[[x]][[y]]))
      wi <- winkelmass(point_coords$X, point_coords$Y, 4)
      wi_mean <- mean(wi)
      wi_cv <- sd(wi) / mean(wi)

      # Calculate metrics 

      # Nearest neighbour distance CV measure 
      neighb_cv <- sd(neighb_dist) / mean(neighb_dist)

      # Point distribution norm (Gunzburger and Burkardt 2004, Saka et al. 2007)
      # $h = max h_{i}$
      # where $h_{i} = max |z_{i}-y|$
      # where $y$ is the point, and $z_{i}$ are the voronoi cell vertices
      # Smaller $h$ = more uniform
      point_distrib <- unlist(lapply(polyg_ex, function(z) {
        max(z[[1]])
      }))

      point_distrib_norm <- max(point_distrib)

      # Point distribution ratio (Gunzburger and Burkardt 2004)
      # $\mu = \frac{max h_{i}}{min h_{i}}$
      # where $h_{i} is same as above
      # For uniform distrib, $\mu = 1$
      point_distrib_ratio <- max(point_distrib) - min(point_distrib)

      # Regularity measure (Gunzburger and Burkardt 2004, Saka et al. 2007)
      # $\chi = max \chi_{i}$
      # where $\chi_{i} = \frac{2h_{i}}{\gamma_{i}}$
      # where $\gamma_{i}$ is the minimum distance between a point $i$ and its nearest neighbour
      # and where $h_{i}$ is same as above
      # For uniform distrib, $\chi = \chi_{i}$. Smaller is more uniform
      reg_measure <- max((2 * neighb_dist) / point_distrib)

      # Centre gravity distance CV
      centre_gravity_distrib <- unlist(lapply(polyg_ex, "[[", 2))
      centre_gravity_cv <- sd(centre_gravity_distrib) / mean(centre_gravity_distrib)

      # Cell area deviation (Gunzburger and Burkardt 2004)
      # $\upsilon = \frac{max V_{i}}{min V_{i}}$
      # where $V_{i}$ is the area of cell $i$
      # In uniform distrib, $\upsilon = 1$
      cell_area_distrib <- unlist(lapply(polyg_ex, "[[", 3))
      cell_area_dev <- max(cell_area_distrib) / min(cell_area_distrib)

      # Cell area CV
      cell_area_cv <- sd(cell_area_distrib) / mean(cell_area_distrib)

      # Polsby-Popper CV
      pp_distrib <- unlist(lapply(polyg_ex, "[[", 4))
      pp_cv <- sd(pp_distrib) / mean(pp_distrib)

      # Cell elongation CV
      elong_distrib <- unlist(lapply(polyg_ex, "[[", 5))
      elong_cv <- sd(elong_distrib) / mean(elong_distrib)

      return(data.frame(wi_mean, wi_cv, neighb_cv, point_distrib_norm, 
          point_distrib_ratio, reg_measure, centre_gravity_cv, cell_area_dev, 
          cell_area_cv, pp_cv, elong_cv, y))
    }))

    out$x <- x

    return(out)
  }, mc.cores = 4)
  return(metrics_df)
}

even_csr_metrics <- metric_gen(dat_even_csr_polys, dat_even_csr_points)
csr_clust_metrics <- metric_gen(dat_csr_clust_polys, dat_csr_clust_points)

