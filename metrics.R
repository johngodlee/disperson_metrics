# Calculate various metrics on voronoi cells
# John Godlee (johngodlee@gmail.com)
# 2021-08-23

# Packages
library(parallel)
library(RANN)

source("functions.R")

# Import data
even_csr_polys <- readRDS("dat/even_csr_polys.rds")
even_csr_points <- readRDS("dat/even_csr_points.rds")
csr_clust_polys <- readRDS("dat/csr_clust_polys.rds")
csr_clust_points <- readRDS("dat/csr_clust_points.rds")
bicuar_points <- readRDS("dat/bicuar_points.rds")
bicuar_polys <- readRDS("dat/bicuar_polys.rds")

# Define function to extract polygon metrics 
# from dataframe of polygons and points
metric_gen <- function(poly_df, point_df) {
  polyg_ex <- lapply(seq_along(poly_df), function(x) {
    # Extract polygon and point
    polyg <- poly_df[[x]]
    point <- as.matrix(point_df[x,])
    
    # Find distance of point to vertices
    vd <- nn2(point, polyg, k = 1)[[2]]

    # Find polygon area
    polygc <- rbind(polyg, polyg[1,])
    n <- nrow(polygc)
    cx <- polygc[,1]
    cy <- polygc[,2] 
    i <- 1:(n-1)
    cell_area <- sum(c(cx[i] * cy[i+1] - cx[i+1] * cy[i])) / 2

    # Find centre of gravity of polygon (centroid)
    centx <- sum((cx[i] + cx[i+1]) * 
      (cx[i] * cy[i+1] - cx[i+1] * cy[i])) / (6 * cell_area)
    centy <- sum((cy[i] + cy[i+1]) * 
      (cx[i] * cy[i+1] - cx[i+1] * cy[i])) / (6 * cell_area)

    # Distance between point and centre of gravity of cell
    cgd <- nn2(point, matrix(c(centx, centy), nrow = 1), 
      k = 1)[[2]]

    # Circumference of cell
    dist_mat <- as.matrix(dist(polyg))
    cell_perim <- sum(
      dist_mat[row(dist_mat) == col(dist_mat) + 1],
      dist_mat[nrow(dist_mat),1])

    # Cell elongation (Polsby-Popper)
    pp <- (4*pi*cell_area)/(cell_perim^2)
    
    # Cell elongation (minimum bounding rectangle ratio)      
    min_box <- unlist(minBox(polyg)[c("width", "height")])
    elong <- max(min_box) / min(min_box)

    return(list(vd, cgd, cell_area, pp, elong))
  })

  # Nearest neighbour distances
  neighb_dist <- nn2(point_df, point_df, k = 2)[[2]][,2]

  # Winkelmass
  wi <- winkelmass(point_df$x, point_df$y, 4)
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
  point_distrib <- unlist(lapply(polyg_ex, function(x) {
    max(x[[1]])
  }))

  point_distrib_norm <- max(point_distrib)

  # Point distribution ratio (Gunzburger and Burkardt 2004)
  # $\mu = \frac{max h_{i}}{min h_{i}}$
  # where $h_{i}$ is same as above
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
  centre_gravity_mean <- mean(centre_gravity_distrib)
  centre_gravity_cv <- sd(centre_gravity_distrib) / centre_gravity_mean

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
      point_distrib_ratio, reg_measure, centre_gravity_mean, centre_gravity_cv, cell_area_dev, 
      cell_area_cv, pp_cv, elong_cv))
}

# Extract metrics from simulated data
sim_metric_gen <- function(polys_list, points_list) {
  even_csr_metrics <- mclapply(seq_along(polys_list), function(x) {
    lapply(seq_along(polys_list[[x]]), function(y) {
    message(x, ":", y)
      metric_gen(
        polys_list[[x]][[y]], 
        points_list[[x]][[y]]
      )
    })
  }, mc.cores = 4)
}

even_csr_metrics <- sim_metric_gen(even_csr_polys, even_csr_points)
csr_clust_metrics <- sim_metric_gen(csr_clust_polys, csr_clust_points)

even_csr_metrics_df <- do.call(rbind, lapply(seq_along(even_csr_metrics), function(x) {
  do.call(rbind, lapply(seq_along(even_csr_metrics[[x]]), function(y) {
    even_csr_metrics[[x]][[y]]$repl <- x
    even_csr_metrics[[x]][[y]]$adj <- y-1
    return(even_csr_metrics[[x]][[y]])
  }))
}))

csr_clust_metrics_df <- do.call(rbind, lapply(seq_along(csr_clust_metrics), function(x) {
  do.call(rbind, lapply(seq_along(csr_clust_metrics[[x]]), function(y) {
    csr_clust_metrics[[x]][[y]]$repl <- x
    csr_clust_metrics[[x]][[y]]$adj <- y-1
    return(csr_clust_metrics[[x]][[y]])
  }))
}))

# Write data 
saveRDS(even_csr_metrics_df, "dat/even_csr_metrics.rds")
saveRDS(csr_clust_metrics_df, "dat/csr_clust_metrics.rds")

# Extract metrics from real data 
bicuar_metrics <- mclapply(seq_along(bicuar_polys), function(x) {
  metric_gen(bicuar_polys[[x]], bicuar_points[[x]])
}, mc.cores = 3)

bicuar_metrics_df <- do.call(rbind, bicuar_metrics)

bicuar_metrics_df$plot_id <- names(bicuar_polys)

saveRDS(bicuar_metrics_df, "dat/bicuar_metrics.rds")
