# Voronoi tessellation of points
# John Godlee (johngodlee@gmail.com)
# 2022-03-17

# Packages
library(parallel)
library(deldir)

# Import data
even_csr_points <- readRDS("dat/even_csr_points.rds")
csr_clust_points <- readRDS("dat/csr_clust_points.rds")

# Get voronoi polygon vertices
vor_polys <- lapply(list(even_csr_points,csr_clust_points), function(i) {
  mclapply(seq_along(i), function(x) {
    lapply(seq_along(i[[x]]), function(y) {
      message(x, ":", y)
      vor <- deldir(i[[x]][[y]], rw = c(0,100,0,100))
      vor_tiles <- tile.list(vor)
      lapply(vor_tiles, function(z) {
        cbind(z$x, z$y)
      })
    })
  }, mc.cores = 4)
})

# Write data
saveRDS(vor_polys[[1]], "dat/even_csr_polys.rds")
saveRDS(vor_polys[[2]], "dat/csr_clust_polys.rds")
