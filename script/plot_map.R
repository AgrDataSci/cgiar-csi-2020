library(mapview)
library(sf)
library(data.table)

dt <- fread("data/wheat_pvs_data.csv")

lonlat <- dt[, c("lon","lat","project_name")]

lonlat$lat[lonlat$lat > 33] <- NA

lonlat <- na.omit(lonlat)

lonlat <- lonlat[!duplicated(lonlat$lat),]

rownames(lonlat) <- 1:nrow(lonlat)

lonlat <- st_as_sf(lonlat, coords = c("lon","lat"), crs = 4326)

mapview(lonlat, alpha = 1, map.types = "OpenTopoMap")
