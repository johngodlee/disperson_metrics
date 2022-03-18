# Ripley's K and L function visualisations, and nearest neighbour cumulative density distributions
# John Godlee (johngodlee@gmail.com)
# 2021-11-14

# Packages
library(ggplot2)
library(spatstat)
library(scico)
library(patchwork)

source("functions.R")

# Import data
even_csr <- readRDS("dat/even_csr_points.rds")
csr_clust <- readRDS("dat/csr_clust_points.rds")
bicuar_points <- readRDS("dat/bicuar_points.rds")

# spatstat Ripley's K and L from Bicuar data
bicuar_fil <- bicuar_points[["ABG_1"]]

s_ppp <- ppp(bicuar_fil$x, bicuar_fil$y, c(0,100), c(0,100))

s_k <- Kest(s_ppp)
s_kenv <- envelope(s_ppp,Kest)

s_l <- Lest(s_ppp)
s_lenv <- envelope(s_ppp,Lest)

s_lenv$obs_r <- s_lenv$obs - s_lenv$r
s_lenv$theo_r <- s_lenv$theo - s_lenv$r
s_lenv$lo_r <- s_lenv$lo - s_lenv$r
s_lenv$hi_r <- s_lenv$hi - s_lenv$r

ripley_k <- ggplot() + 
  geom_ribbon(data = s_kenv, aes(x = r, ymin = lo, ymax = hi), alpha = 0.5) + 
  geom_line(data = s_kenv, aes(x = r, y = theo), linetype = 2, colour = "red") + 
  geom_line(data = s_kenv, aes(x = r, y = obs)) + 
  mytheme() + 
  labs(y = expression(K(r)), x = "d (m)")
ggsave(ripley_k, width = 8, height = 5, file = "img/ripley_k.png")

ripley_l <- ggplot() + 
  geom_ribbon(data = s_lenv, aes(x = r, ymin = lo, ymax = hi), alpha = 0.5) + 
  geom_line(data = s_lenv, aes(x = r, y = theo), linetype = 2, colour = "red") + 
  geom_line(data = s_lenv, aes(x = r, y = obs)) + 
  mytheme() + 
  labs(y = expression(L(r)), x = "d (m)")
ggsave(ripley_l, width = 8, height = 5, file = "img/ripley_l.png")

ripley_l_norm <- ggplot() + 
  geom_ribbon(data = s_lenv, aes(x = r, ymin = lo_r, ymax = hi_r), alpha = 0.5) + 
  geom_line(data = s_lenv, aes(x = r, y = theo_r), linetype = 2, colour = "red") + 
  geom_line(data = s_lenv, aes(x = r, y = obs_r)) + 
  mytheme() + 
  labs(y = expression(hat("L")), x = "d (m)")
ggsave(ripley_l_norm, width = 8, height = 5, file = "img/ripley_l_norm.png")

ripley_klr <- ripley_k + ripley_l + ripley_l_norm
ggsave(ripley_klr, width = 8, height = 3.5, file = "img/ripley_klr.png")

# Evenness CSR
adj_seq <- seq(0, length(even_csr[[1]]), 10)
even_csr_fil <- even_csr[[1]][adj_seq+1]

even_csr_l <- do.call(rbind, lapply(seq_along(even_csr_fil), function(x) {
  x_ppp <- ppp(even_csr_fil[[x]]$x, even_csr_fil[[x]]$y, c(0,100), c(0,100))
  x_l <- envelope(x_ppp, Lest, verbose = FALSE)

  x_l$obs_r <- x_l$obs - x_l$r
  x_l$theo_r <- x_l$theo - x_l$r
  x_l$lo_r <- x_l$lo - x_l$r
  x_l$hi_r <- x_l$hi - x_l$r

  x_l_df <- as.data.frame(x_l)
  x_l_df$adj <- adj_seq[x]

  return(x_l_df)
}))

pdf(file = "img/even_csr_ripley_l.pdf", width = 6, height = 4)
ggplot() + 
  geom_ribbon(data = even_csr_l[even_csr_l$adj == 0,], 
    aes(x = r, ymin = lo_r, ymax = hi_r), alpha = 0.5) +
  geom_line(data = even_csr_l, 
    aes(x = r, y = obs_r, group = adj,  colour = adj)) + 
  scale_colour_scico(name = "Randomness", palette = "batlow") + 
  theme_bw() + 
  theme(
    legend.position = c(0.8,0.3), 
    legend.box.background = element_rect(colour = "black")) + 
  labs(y = expression(hat("L")), x = "d (m)")
dev.off()

csr_clust_fil <- csr_clust[[1]][adj_seq+1]
csr_clust_fil <- csr_clust_fil[lengths(csr_clust_fil) != 0]

csr_clust_l <- do.call(rbind, lapply(seq_along(csr_clust_fil), function(x) {
  x_ppp <- ppp(csr_clust_fil[[x]]$x, csr_clust_fil[[x]]$y, c(0,100), c(0,100))
  x_l <- envelope(x_ppp, Lest, verbose = FALSE)

  x_l$obs_r <- x_l$obs - x_l$r
  x_l$theo_r <- x_l$theo - x_l$r
  x_l$lo_r <- x_l$lo - x_l$r
  x_l$hi_r <- x_l$hi - x_l$r

  x_l_df <- as.data.frame(x_l)
  x_l_df$adj <- adj_seq[x]

  return(x_l_df)
}))

pdf(file = "img/csr_clust_ripley_l.pdf", width = 6, height = 4)
ggplot() + 
  geom_ribbon(data = csr_clust_l[csr_clust_l$adj == 0,], 
    aes(x = r, ymin = lo_r, ymax = hi_r), alpha = 0.5) +
  geom_line(data = csr_clust_l, 
    aes(x = r, y = obs_r, group = adj,  colour = adj)) + 
  scale_colour_scico(name = "Clustering", palette = "batlow") + 
  theme_bw() + 
  theme(
    legend.position = c(0.9,0.7), 
    legend.box.background = element_rect(colour = "black")) + 
  labs(y = expression(hat("L")), x = "d (m)")
dev.off()

# Cumulative density distributions of nearest neighbour distances - G
even_csr_g <- do.call(rbind, lapply(seq_along(even_csr_fil), function(x) {
  x_ppp <- ppp(even_csr_fil[[x]]$x, even_csr_fil[[x]]$y, c(0,100), c(0,100))
  x_g <- envelope(x_ppp, Gest, verbose = FALSE)
  x_g_df <- as.data.frame(x_g)
  x_g_df$adj <- adj_seq[x]

  return(x_g_df)
}))

pdf(file = "img/even_csr_ripley_g.pdf", width = 6, height = 4)
ggplot() + 
  geom_ribbon(data = even_csr_g[even_csr_g$adj == 0,], 
    aes(x = r, ymin = lo, ymax = hi), alpha = 0.5) +
  geom_line(data = even_csr_g[even_csr_g$adj == 0,], 
    aes(x = r, y = theo), colour = "red", linetype = 2) + 
  geom_line(data = even_csr_g, 
    aes(x = r, y = obs, group = adj,  colour = adj)) + 
  scale_colour_scico(name = "Randomness", palette = "batlow") + 
  theme_bw() + 
  theme(
    legend.position = c(0.8,0.3), 
    legend.box.background = element_rect(colour = "black")) + 
  labs(y = "G", x = "d (m)")
dev.off()

csr_clust_g <- do.call(rbind, lapply(seq_along(csr_clust_fil), function(x) {
  x_ppp <- ppp(csr_clust_fil[[x]]$x, csr_clust_fil[[x]]$y, c(0,100), c(0,100))
  x_g <- envelope(x_ppp, Gest, verbose = FALSE)
  x_g_df <- as.data.frame(x_g)
  x_g_df$adj <- adj_seq[x]

  return(x_g_df)
}))

pdf(file = "img/csr_clust_ripley_g.pdf", width = 6, height = 4)
ggplot() + 
  geom_ribbon(data = csr_clust_g[csr_clust_g$adj == 0,], 
    aes(x = r, ymin = lo, ymax = hi), alpha = 0.5) +
  geom_line(data = csr_clust_g[csr_clust_g$adj == 0,], 
    aes(x = r, y = theo), colour = "red", linetype = 2) + 
  geom_line(data = csr_clust_g, 
    aes(x = r, y = obs, group = adj,  colour = adj)) + 
  scale_colour_scico(name = "Clustering", palette = "batlow") + 
  theme_bw() + 
  theme(
    legend.position = c(0.8,0.3), 
    legend.box.background = element_rect(colour = "black")) + 
  labs(y = "G", x = "d (m)")
dev.off()

