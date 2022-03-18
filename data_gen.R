# Create data for Voronoi tesselation test
# John Godlee (johngodlee@gmail.com)
# 2021-08-30

# Packages
library(parallel)
library(RANN)

# Create a regular grid of points
# 5 m buffer inside plot
# 10 m between points
# 100x100 m grid
xy_vec <- seq(5, 100, 10)
dat <- expand.grid(xy_vec, xy_vec)
names(dat) <- c("x", "y")

# Number of replicates
reps <- 50

# Replicate grid of points
even_list <- replicate(reps, dat, simplify = FALSE)

# Possible coordinates for movement
c_repls <- seq(0,100,0.1)

# For each replicate, randomly move points
# Move every point once
# Evenness -> randomness
even_csr_points <- mclapply(even_list, function(x) {
  # Randomise order in which points are moved
  ord <- sample(seq_len(nrow(x)), nrow(x))
  # Create list to fill
  x <- list(x)
  # For each point 
  for (i in seq_len(nrow(x[[1]]))) {
    x[[i + 1]] <- x[[i]]
    # Move point to a random location
    x[[i + 1]][ord[i], c("x", "y")] <- sample(c_repls, 2, replace = TRUE)
  }
  return(x)
}, mc.cores = 4)

# Write data
saveRDS(even_csr_points, "dat/even_csr_points.rds")
even_csr_points <- readRDS("dat/even_csr_points.rds")
    
# Generate random points
rand_list <- replicate(reps, 
  data.frame(
    x = sample(c_repls, 100),
    y = sample(c_repls, 100)), 
  simplify = FALSE)

# For each replicate, randomly points closer to nearest cluster point
csr_clust_points <- mclapply(rand_list, function(x) {
  # Define points to use as clusters
  clust_rows <- sample(seq_len(nrow(x)), 10)
  # Extract cluster points
  clust_pts <- x[clust_rows,]
  # Extract non-cluster points
  pts <- x[-clust_rows,]
  # Find nearest cluster point for each non-cluster point
  nn <- nn2(clust_pts[,1:2], pts[,1:2], k = 1)[[1]]
  # Randomise order in which points are moved
  ord <- sample(seq_len(nrow(pts)), nrow(pts))
  # Create list to fill
  pts <- list(pts)
  # For each point
  for (i in seq_len(nrow(pts[[1]]))) {
    pts[[i + 1]] <- pts[[i]]
    # Get midpoint between chosen point and cluster point and update coords
    pts[[i + 1]][i,1:2] <- (pts[[i + 1]][i,1:2] +
      clust_pts[nn[i],1:2]) / 2
  }
  return(pts)
}, mc.cores = 4)

# Write data
saveRDS(csr_clust_points, "dat/csr_clust_points.rds")
