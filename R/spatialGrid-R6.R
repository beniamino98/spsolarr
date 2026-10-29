#' Spatial Grid
#'
#' Create a grid from a range of latitudes and longitudes.
#'
#' @return a tibble with two columns `lat` and `lon`.
#'
#' @examples
#' # Initialize a spatial grid
#' sp <- spatialGrid$new()
#'
#' # Create an equally spaced grid
#' sp$make_grid(c(43.1, 44), c(9.2, 12.3), c(0.1, 0.1))
#' sp
#' # Check if a point is inside the bounds
#' sp$is_inside_bounds(44.8, 10.9)
#' sp$is_inside_bounds(43.9, 12.1)
#' 
#' # Extract a known point
#' sp$neighborhoods(43.9, 12.1)
#' 
#' # Extract its neighborhoods
#' sp$neighborhoods(43.39, 12.1)$ngb
#'
#' @rdname spatialGrid
#' @name spatialGrid
#' @export
spatialGrid <- R6::R6Class("spatialGrid",
                           public = list(
                             #' @field weights Weighting function for the distance.
                             weights = NA,
                             #' @description
                             #' Initialize a spatial grid
                             #' @param weights Weighting function for the distance.
                             initialize = function(weights = IDW(2)){
                               # Weighting function
                               self$weights <- weights
                             },
                             #' @description
                             #' Create a spatial grid
                             #' @param lat Numeric vector, from which is extracted the minimum and maximum for latitude.
                             #' @param lon Numeric vector, from which is extracted the minimum and maximum for longitude.
                             #' @param by Numeric vector, the first element is used to establish the distance between two latitudes in the grid.
                             #' The second element (if present) is used to establish the distance between two longitudes in the grid.
                             #' @param update.grid Logical, when `TRUE` the grid inside the model will be updated with the one contructed. 
                             #' @param digits Integer scalar, number of digits for latitudes and longitudes.
                             make_grid = function(lat, lon, by, update.grid = FALSE, CRS = FALSE, digits = 5){
                               # Latitudes bounds
                               lat_min <- min(lat)
                               lat_max <- max(lat)
                               lat_by  <- by[1]
                               # Longitudes bounds
                               lon_min <- min(lon)
                               lon_max <- max(lon)
                               lon_by  <- ifelse(length(by) == 1, by[1], by[2])
                               # Sequence of latitudes and longitudes 
                               lats <- round(seq(lat_min, lat_max, by = lat_by), digits = digits)
                               lons <- round(seq(lon_min, lon_max, by = lon_by), digits = digits)
                               # Matrix of latitudes and longitudes
                               grid_lat_lon <- expand.grid(lats, lons)
                               grid_lat_lon <- dplyr::as_tibble(grid_lat_lon)
                               # Standard column names
                               colnames(grid_lat_lon) <- c("lat", "lon")
                               # Assign a unique ID for each place
                               labels <- paste0("ID_", 1:nrow(grid_lat_lon))
                               grid_lat_lon <- dplyr::bind_cols(id = labels, grid_lat_lon)
                               # Projected grid 
                               grid_lat_lon_proj <- project_CRS(grid_lat_lon$lat, grid_lat_lon$lon, id = grid_lat_lon$id)
                               
                               # Store the spatial grid
                               if (update.grid) {
                                 # Latitudes bounds
                                 private[["lat_min"]] <- lat_min
                                 private[["lat_max"]] <- lat_max
                                 private[["lat_by"]]  <- lat_by
                                 # Longitudes bounds
                                 private[["lon_min"]] <- lon_min
                                 private[["lon_max"]] <- lon_max
                                 private[["lon_by"]]  <- lon_by
                                 # Sequence of latitudes and longitudes 
                                 private[["lats"]] <- lats
                                 private[["lons"]] <- lons
                                 # Grid 
                                 private[["..grid"]] <- grid_lat_lon
                                 # Projected grid 
                                 private[["..data"]] <- grid_lat_lon_proj
                               } 
                               if (CRS) {
                                 return(grid_lat_lon_proj)
                               } else {
                                 return(grid_lat_lon)
                               }
                             },
                             #' @description
                             #' Set a spatial grid
                             #' @param grid Tibble with column `id`, `lat` and `lon`.
                             set_grid = function(grid){
                               # Extract a sequence of latitudes 
                               lats <- as.numeric(unique(grid$lat))
                               lats <- lats[order(lats)]
                               # Extract a sequence of longitudes 
                               lons <- as.numeric(unique(grid$lon))
                               lons <- lons[order(lons)]
                               # Assign a unique ID for each place
                               if (is.null(grid[["id"]])) {
                                 cli::cli_alert_warning("The column `id` is missing in the grid. Assigned new ID.")
                                 labels = paste0("ID_", 1:nrow(grid))
                                 grid <- dplyr::bind_cols(id = labels, grid)
                               }
                               # ********* Update private *********
                               # Latitudes 
                               private[["lat_min"]] <- min(lats)
                               private[["lat_max"]] <- max(lats)
                               private[["lat_by"]] <- abs(diff(lats)[1])
                               private[["lats"]] <- lats
                               # Longitudes
                               private[["lon_min"]] <- min(lons)
                               private[["lon_max"]] <- max(lons)
                               private[["lon_by"]] <- abs(diff(lons)[1])
                               private[["lons"]] <- lons
                               # Spatial grid 
                               private[["..grid"]] <- grid
                               # Projected grid 
                               private[["..data"]] <- project_CRS(grid$lat, grid$lon, id = grid$id)
                             },
                             #' @description
                             #' Check if a location is inside the bounds of the grid.
                             #' @param lat Numeric vector, reference latitudes.
                             #' @param lon Numeric vector, reference longitudes.
                             #' @return `TRUE` when the point is inside the limits and `FALSE` otherwise.
                             is_inside_bounds = function(lat, lon){
                               coords <- dplyr::tibble(lat = lat, lon = lon)
                               condition <- rep(TRUE, nrow(coords))
                               for(i in 1:nrow(coords)) {
                                 condition[i] <- coords$lon[i] >= private$lon_min & coords$lon[i] <= private$lon_max
                                 condition[i] <- condition[i] & coords$lat[i] >= private$lat_min & coords$lat[i] <= private$lat_max
                               }
                               return(condition)
                             },
                             #' @description
                             #' Check if a location is a known point inside the grid.
                             #' @param lat Numeric vector, reference latitudes.
                             #' @param lon Numeric vector, reference longitudes.
                             #' @param n Integer scalar, number of neighborhoods to consider. 
                             #' @return `TRUE` when the location is known and `FALSE` otherwise.
                             neighborhoods = function(lat, lon, n = 4){
                               # Coordinates of the point 
                               coords <- dplyr::tibble(lat = lat, lon = lon)
                               coords_new <- project_CRS(coords$lat, coords$lon)
                               # Dataset of known points 
                               data <- self$data
                               res <- list()
                               isKnown <- c()
                               for(i in 1:length(lat)){
                                 coords_k <- coords_new[i,]
                                 # Distances 
                                 h <- sp::spDistsN1(data, coords_k, longlat = FALSE)
                                 # Index of the neighborhoods
                                 idx_ngb <- order(h, decreasing = FALSE)[1:n]
                                 # Add distance
                                 res[[i]] <- data[idx_ngb,]
                                 res[[i]]$h <- h[idx_ngb]
                                 # Weights 
                                 res[[i]]$wgt <- self$weights(h[idx_ngb], normalize = TRUE)
                                 # Add known point detection
                                 isKnown[i] <- any(h[idx_ngb] == 0)
                                 res[[i]]$wgt[1] <- ifelse(isKnown[i], 1, res[[i]]$wgt[1])
                               }
                               coords$isKnown <- isKnown
                               coords$ngb <- res
                               coords
                             },
                             #' @description
                             #' Method print for a `spatialGrid` object.
                             print = function(){
                               cat(paste0("--------------------- ", "spatial Grid", "--------------------- \n"))
                               cat(paste0("Latitudes:  ", "(", "\033[1;35m", round(private$lat_min, 3), "\033[0m", " - ", "\033[1;35m", round(private$lat_max, 3),
                                          "\033[0m)", " by ", round(private$lat_by, 3), "\n"))
                               cat(paste0("Longitudes: ", "(", "\033[1;35m", round(private$lon_min, 3), "\033[0m", " - ", "\033[1;35m", round(private$lon_max, 3), 
                                          "\033[0m)", " by ", round(private$lon_by, 3), "\n"))
                               cat(paste0("Number of points: ", "\033[1;35m", nrow(self$grid), "\033[0m", "\n"))
                             }
                           ),
                           private = list(
                             ..grid = NA,
                             ..data = NA,
                             lats = NA,
                             lons = NA,
                             lat_min = NA,
                             lat_max = NA,
                             lon_min = NA,
                             lon_max = NA,
                             lon_by = NA,
                             lat_by = NA
                           ),
                           active = list(
                             #' @field grid A tibble with the spatial grid.
                             grid = function(){
                               private$..grid
                             },
                             #' @field data Projected grid 
                             data = function(){
                               private$..data
                             },
                             #' @field range_lat range of latitudes
                             range_lat = function(){
                               c(min = private$lat_min, max = private$lat_max)
                             },
                             #' @field range_lon range of longitudes
                             range_lon = function(){
                               c(min = private$lon_min, max = private$lon_max)
                             }
                           )
)

#' Inverse Distance Weighting Functions
#'
#' Return a distance weighting function
#'
#' @param beta parameter used in exponential and power functions.
#' @param d0 parameter used only in exponential function.
#' @details When the parameter `d0` is not specified the function returned will be of power type otherwise of exponential type.
#'
#' @examples
#' # Power weighting
#' IDW_pow <- IDW(2)
#' IDW_pow(c(2, 3,10))
#' IDW_pow(c(2, 3,10), normalize = TRUE)
#' # Exponential weighting
#' IDW_exp <- IDW(2, d0 = 5)
#' IDW_exp(c(2, 3,10))
#' IDW_exp(c(2, 3,10), normalize = TRUE)
#' @export
IDW <- function(beta, d0 = NULL){
  # Power function
  IDW.pow <- function(beta) {
    force(beta);
    function(d, normalize = FALSE){
      w <- 1/(d^beta)
      if (normalize) {
        w <- w/sum(w)
      }
      return(w)
    }
  }
  # Exponential function
  IDW.exp <- function(beta, d0) {
    force(beta); force(d0);
    function(d, normalize = FALSE){
      w <- exp(-(d/d0)^beta)
      if (normalize) {
        w <- w/sum(w)
      }
      return(w)
    }
  }
  # Output function
  if (is.null(d0)) {
    IDW.pow(beta)
  } else {
    IDW.exp(beta, d0)
  }
}
