

plot_map <- function(data, coords = NULL, add = NULL) {
  
  lonlat <- data[, c(coords, add)]
  
  lonlat <- stats::na.omit(lonlat)
  
  rownames(lonlat) <- 1:nrow(lonlat)
  
  lonlat <- sf::st_as_sf(lonlat, coords = c("lon","lat"), crs = 4326)
  
  mapview::mapview(lonlat, alpha = 1, map.types = "OpenTopoMap")

}
 