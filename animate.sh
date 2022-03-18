#!/usr/bin/env sh

ffmpeg -framerate 6 -i frames/even_csr/%03d_even_csr.png -y vid/even_csr.mp4

ffmpeg -framerate 6 -i frames/csr_clust/%03d_csr_clust.png -y vid/csr_clust.mp4
