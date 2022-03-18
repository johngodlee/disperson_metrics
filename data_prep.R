# Create an anonymised and cleaned real world dataset from Bicuar data
# John Godlee (johngodlee@gmail.com)
# 2021-11-14

# Packages
library(dplyr)
library(parallel)
library(RANN)
library(deldir)

# Import data
s <- read.csv("./dat/bicuar_stems.csv")

s_points <- s %>%
  group_by(plot_id, x_grid, y_grid) %>%
  summarise() %>% 
  filter(
    !is.na(x_grid), x_grid <= 100, x_grid >= 0,
    !is.na(y_grid), y_grid <= 100, y_grid >= 0) %>%
  rename(x = x_grid, y = y_grid) %>% 
  group_by(plot_id) %>%
  as.data.frame() %>% 
  split(., .$plot_id) %>% 
  lapply(., "[", c("x", "y"))
  
saveRDS(s_points, "dat/bicuar_points.rds")

s_polys <-  lapply(s_points, function(x) {
  vor <- deldir(x, rw = c(0,100,0,100))
  vor_tiles <- tile.list(vor)
  lapply(vor_tiles, function(z) {
    cbind(z$x, z$y)
  })
})

saveRDS(s_polys, "dat/bicuar_polys.rds")
