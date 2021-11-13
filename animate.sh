convert -delay 10 -loop 0 frames/even_csr/*_even_csr.png -scale 500x500 img/even_csr.gif

convert img/even_csr.gif -coalesce +dither -colors 64 -layers optimize img/even_csr_small.gif

convert -delay 10 -loop 0 frames/csr_clust/*_csr_clust.png -scale 500x500 img/csr_clust.gif

convert img/csr_clust.gif -coalesce +dither -colors 64 -layers optimize img/csr_clust_small.gif
