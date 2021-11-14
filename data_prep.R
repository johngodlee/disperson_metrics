# Create an anonymised and cleaned real world dataset from Bicuar data
# John Godlee (johngodlee@gmail.com)
# 2021-11-14

# Packages
library(dplyr)

# Import data
s <- read.csv("dat/bicuar_stems.csv")

scale01 <- function(x){(x-min(x))/(max(x)-min(x))}

s_clean <- s %>%
  filter(plot %in% 1:5) %>%
  group_by(plot, dec_longitude, dec_latitude) %>%
  slice(1) %>%
  group_by(plot) %>%
  mutate(
    tree_id = row_number(),
    dec_longitude = scale01(dec_longitude) * 100,
    dec_latitude = scale01(dec_latitude) * 100) %>% 
  dplyr::select(
    plot_id = plot, 
    tree_id, 
    x_grid = dec_longitude,
    y_grid = dec_latitude)
  
write.csv(s_clean, "dat/trees.csv")

