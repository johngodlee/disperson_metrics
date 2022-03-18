# Calculating dispersion metrics using Voronoi tessellation

## Scripts 

* `data_prep.R` - creates an anonymised and cleaned dataset of stems from Bicuar National Park permanent plots
* `data_gen.R` - creates many simulated stem location datasets with different levels of dispersion
* `voronoi.R` - performs Voronoi tessellation on simulated and real data
* `metrics.R` - calculates various metrics of dispersion from voronoi cells, using the simulated data
* `ripley_neighb.R` - explores Ripley's K and L functions, and nearest neighbour functions using the simulated data
* `vis.R` - generates plots showing the behaviour of various measures of dispersion
* `frames.R` - generates images showing evolution of simlated data as points are moved
* `functions.R` - contains miscellaneous functions
* `animate.sh` - takes the images produced by `frames.R` and creates videos
