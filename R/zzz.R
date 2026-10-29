#' Project coordinates
#' 
#' @param lat latitudes 
#' @param lon longitudes
#' @param id opttional id for each location 
#' 
#' @rdname project_CRS
#' @name project_CRS
#' @export
project_CRS <- function(lat, lon, id, utm = 32632){
  # Coordinates in lat and lon 
  coords <- data.frame(lon = lon, lat = lat)
  sf_pts <- sf::st_as_sf(coords, coords = c("lon","lat"), crs = 4326) 
  # choose a local UTM 
  coords_proj <- as.data.frame(sf::st_coordinates(sf::st_transform(sf_pts, utm)))
  colnames(coords_proj) <- c("lon", "lat")
  sp::coordinates(coords_proj) = ~lon+lat
  if (missing(id)) {
    coords_proj$id <- rep("", nrow(coords))
  } else {
    coords_proj$id <- id
  }
  coords_proj
}


#' Make a matrix positive semi-defined
#'
#' The matrix is decomposed using a spectral decomposition.
#' Then the negative eigenvalues are imputed with `neg_values` 
#' and the original matrix is constructed again.
#'
#' @param x symmetric matrix
#' @param neg_values numeric
#' @noRd
#' @export
makeSemiPositive <- function(x, neg_values = 1e-5){
  mat <- x
  dec <- eigen(x)
  e <- dec$vectors
  lam <- dec$values
  if (any(lam < 0)) {
    lam[lam < 0] <- neg_values
    mat <- e %*% diag(lam) %*% t(e)
  }
  attr(mat, "index_neg_values") <- which(dec$values < 0)
  attr(mat, "original_values") <- dec$values[dec$values < 0]
  return(mat)
}


#' Haversine distance
#'
#' Compute the Haversine distance between two points.
#'
#' @param lat_1 Numeric vector, latitudes of first location.
#' @param lon_1 Numeric vector, longitudes of first location.
#' @param lat_2 Numeric vector, latitudes of second location.
#' @param lon_2 Numeric vector, longitudes of second location.
#'
#' @examples
#' havDistance(43.3, 12.1, 43.4, 12.2)
#' havDistance(c(43.35, 43.35), c(12.15, 12.1), c(43.4, 44.5), c(12.2, 13.4))
#' @name havDistance
#' @rdname havDistance
#' @return Vector of distances in kilometers.
#' @export
havDistance <- function(lat_1, lon_1, lat_2, lon_2){
  r <- 6371 # Radius of the earth in kilometres
  dist <- c()
  # Solar functions
  sf <- seasonalSolarFunctions$new()
  # Compute distances
  for(i in 1:length(lat_1)){
    # Coordinates
    coord_A <- c(lat_1[i], lon_1[i])
    coord_B <- c(lat_2[i], lon_2[i])
    # Convert in radiant
    point_A <- sf$radiant(coord_A)
    point_B <- sf$radiant(coord_B)
    # Latitudes
    phi <- c(point_A[1], point_B[1])
    delta_phi <- phi[2] - phi[1]
    # Longitudes
    lam <- c(point_A[2], point_B[2])
    delta_lam <- lam[2] - lam[1]
    # Compute distance
    dist[i] <- 2*r*asin(sqrt(0.5*(1 - cos(delta_phi) + cos(phi[1])*cos(phi[2])*(1-cos(delta_lam)))))
  }
  return(dist)
}


#' Compute the anisotropic distance 
#' 
#' @param X projected surface (see project_CRS())
#' @param theta anisotropic rotation
#' @param ratio anisotropic rotation parameter
#' @rdname dist_aniso
#' @name dist_aniso
#' @export
dist_aniso <- function(X, theta, ratio){
  R  <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)    
  S  <- diag(c(1, 1/ratio))                                    
  A  <- S %*% R
  coords_star <-  X@coords %*% t(A)   
  # plain Euclidean distances in meters
  sp::spDists(coords_star)
}
