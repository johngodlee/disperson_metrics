mytheme <- function() {
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = NA),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 10),
    strip.text = element_text(size = 10),
    legend.text = element_text(size = 10),
  )
}

minBox <- function(x) {
  ## rotating calipers algorithm using the convex hull
  H <- chull(x)  ## hull indices, vertices ordered clockwise
  n <- length(H)  ## number of hull vertices
  hull <- x[H, ]  ## hull vertices
  
  ## unit basis vectors for all subspaces spanned by the hull edges
  hDir <- diff(rbind(hull, hull[1,]))  ## hull vertices are circular
  hLens <- sqrt(rowSums(hDir^2))  ## length of basis vectors
  huDir <- diag(1/hLens) %*% hDir  ## scaled to unit length
  
  ## unit basis vectors for the orthogonal subspaces
  ## rotation by 90 deg -> y' = x, x' = -y
  ouDir <- cbind(-huDir[,2], huDir[,1])
  
  ## project hull vertices on the subspaces spanned by the hull edges, and on
  ## the subspaces spanned by their orthogonal complements - in subspace coords
  projMat <- rbind(huDir, ouDir) %*% t(hull)
  
  ## range of projections and corresponding width/height of bounding rectangle
  rangeH <- matrix(numeric(n*2), ncol=2)  ## hull edge
  rangeO <- matrix(numeric(n*2), ncol=2)  ## orthogonal subspace
  widths <- numeric(n)
  heights <- numeric(n)
  
  for(i in seq(along=numeric(n))) {
    rangeH[i,] <- range(projMat[i,])
    
    ## the orthogonal subspace is in the 2nd half of the matrix
    rangeO[i,] <- range(projMat[n+i,])
    widths[i] <- abs(diff(rangeH[i,]))
    heights[i] <- abs(diff(rangeO[i,]))
  }
  
  ## extreme projections for min-area rect in subspace coordinates
  ## hull edge leading to minimum-area
  eMin <- which.min(widths*heights)
  hProj <- rbind(rangeH[eMin,], 0)
  oProj <- rbind(0, rangeO[eMin,])
  
  ## move projections to rectangle corners
  hPts <- sweep(hProj, 1, oProj[,1], "+")
  oPts <- sweep(hProj, 1, oProj[,2], "+")
  
  ## corners in standard coordinates, rows = x,y, columns = corners
  ## in combined (4x2)-matrix: reverse point order to be usable in polygon()
  ## basis formed by hull edge and orthogonal subspace
  basis <- cbind(huDir[eMin,], ouDir[eMin,])
  hCorn <- basis %*% hPts
  oCorn <- basis %*% oPts
  pts <- t(cbind(hCorn, oCorn[,c(2,1)]))
  
  ## angle of longer edge pointing up
  dPts <- diff(pts)
  e <- dPts[which.max(rowSums(dPts^2)), ]  ## one of the longer edges
  eUp <- e * sign(e[2])  ## rotate upwards 180 deg if necessary
  deg <- atan2(eUp[2], eUp[1])*180/pi  ## angle in degrees
  
  return(list(pts = pts, width = heights[eMin], 
      height = widths[eMin], angle = deg))
}

#' Winkelmass (spatial regularity of trees)
#'
#' @param x vector of individual x axis coordinates
#' @param y vector of individual y axis coordinates
#' @param k number of neighbours to consider
#'
#' @return 
#' 
#' @references von Gadow, K., Hui, G. Y. (2001). Characterising forest spatial 
#' structure and diversity. Sustainable Forestry in Temperate Regions. Proc. of 
#' an international workshop organized at the University of Lund, Sweden. 
#' Pages 20- 30.
#' 
#' @export
#' 
winkelmass <- function(x, y, k = 4, a0 = 72) {
  dat_sf <- sf::st_as_sf(data.frame(x,y), coords = c("x", "y"))

  dists <- suppressMessages(nngeo::st_nn(dat_sf, dat_sf, k = k+1, 
      progress = FALSE))

  wi <- unlist(lapply(dists, function(i) {
    focal_sfg <- sf::st_geometry(dat_sf[i[1],])[[1]]
    nb_sfg <- sf::st_geometry(dat_sf[i[-1],])
    nb_angles <- sort(unlist(lapply(nb_sfg, function(j) {
      angleCalc(focal_sfg, j)
    })))
  
    aj <- nb_angles - c(NA, head(nb_angles, -1))
    aj[1] <- nb_angles[k] - nb_angles[1]
    aj <- ifelse(aj > 180, 360 - aj, aj)
    aj <- round(aj, 1)
    sum(aj < a0)
  }))

  out <- 1 / k * wi

  return(out)
}

#' Calculate angle between two sf point objects
#'
#' @param x point feature of class 'sf'
#' @param y point feature of class 'sf'
#'
#' @return azimuthal from x to y, in degrees
#' 
#' @examples
#' 
#' @export
#' 
angleCalc <- function(x, y) {
  dst_diff <- as.numeric(x - y)
  return((atan2(dst_diff[1], dst_diff[2]) + pi) / 0.01745329)
}

