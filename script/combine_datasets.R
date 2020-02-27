####################################################
## Read and clean wheat PVS data 
#Kaue de Sousa
#INN University & Bioversity International
####################################################

# Packages
library("tidyverse")
source("script/functions.R")

# Get names of csv files in data/raw
files <- as.list(list.files("data/raw",
                               pattern = "wheat", 
                               full.names = TRUE)[-1])

projnames <- lapply(files, function(x){
  
  x <- strsplit(x, "/")[[1]][[3]]
  x <- gsub("pvs_|.csv","",x)
})

projnames <- unlist(projnames)

#read cvs files and combine then
mydata <- lapply(files, function(x){
  
  file <- read_csv(x, na = c("NA","", " ","-","#DIV/0!","#VALUE!"))
  }

)

names_data <- lapply(mydata, names)

names_data <- unique(unlist(names_data))

for(i in seq_along(mydata)) {
  
  dat_i <- mydata[[i]]
  
  dat_j <- matrix(NA, 
                  nrow = nrow(dat_i),
                  ncol = length(names_data), 
                  dimnames = list(1:nrow(dat_i),
                                  names_data))
  
  dat_j <- data.frame(dat_j)
  
  index <- which(names(dat_j) %in% names(dat_i))
  
  dat_j[, index] <- dat_i
  dat_j$project_name <- projnames[i]
  mydata[[i]] <- dat_j

}

rm(dat_i, dat_j)

mydata <- do.call(rbind, mydata)

# create an id for the plots
ids <- with(mydata, paste0(farmer_name, father_name, village, district))

ids <- as.integer(as.factor(ids))

mydata$id <- ids

boxplot(mydata$wts_gm)

# yield into numeric
index <- which(grepl("yield", names_data))

mydata[, index] <- lapply(mydata[, index], as.numeric)

# convert from kg_sqrm to kg_ha
index <- which(grepl("kg_sqrm", names_data))
mydata[index,] <- mydata[, index] * 1000

# convert from q_ha to kg_ha
index <- which(grepl("q_ha", names_data))
mydata[index, ] <- mydata[, index] * 99.79

mydata$yield_kg_ha <- with(mydata, 
                           ifelse(is.na(yield_kg_sqrm) & !is.na(yield_q_ha),
                                  yield_q_ha, yield_kg_sqrm))


mydata$bhusa_yield_kg_ha <- with(mydata, 
                           ifelse(is.na(bhusa_yield_kg_sqrm) & !is.na(bhusa_yield_q_ha),
                                  bhusa_yield_q_ha, bhusa_yield_kg_sqrm))


index <- which(grepl("kg_ha", names(mydata)))

boxplot(mydata[, index])

mydata$yield_kg_ha <- with(mydata, 
                           ifelse(yield_kg_ha < 250, NA, 
                                  yield_kg_ha)) 

boxplot(mydata$yield_kg_ha)

mydata$bhusa_yield_kg_ha <- with(mydata, 
                                 ifelse(bhusa_yield_kg_ha < 250 | bhusa_yield_kg_ha > 900, 
                                        NA, bhusa_yield_kg_ha))

boxplot(mydata$bhusa_yield_kg_ha)

# .................................
# .................................
# fix lonlat coords ####
mydata[, c("lon","lat","project_name")]

unique(mydata$lon)
as.character(unique(mydata$lat))

mydata$lat <- ifelse(mydata$lat == 2467.766784, 24.67766784,
                            ifelse(mydata$lat == 2472960, 24.72960,
                                   mydata$lat))

mydata$lat <- ifelse(mydata$lat > 2400, mydata$lat / 100, mydata$lat)

mydata$lat[mydata$lat  > 33] <- NA

# put both sides as NAs
na <- is.na(mydata$lat) | is.na(mydata$lon)
mydata[na, "lon"] <- NA
mydata[na, "lat"] <- NA

sum(is.na(unlist(mydata,c("lon", "lat"))))/2

plot_map(mydata, c("lon", "lat"), "project_name")

# identify outliers in geographic coordinates per project
unique(mydata$project_name)
keep <- grepl("wheat", mydata$project_name)
mydata <- mydata[keep, ]
p <- unique(mydata$project_name)

# create a z vector which is the sum of longitude and latitude
mydata$z <- rowSums(mydata[, c("lon","lat")])

for(i in seq_along(p)){
  
  # get a logic vector for the rows within project i
  r <- mydata$project == p[i]
  
  #keep z for values in project i and NA for the other values
  z_i <- mydata$z
  z_i[!r] <- NA
  
  #if all geographical entries are missing jump to the next element in the loop
  if(sum(!is.na(mydata[r, "lat"])) == 0) next
  
  #look for outliers first with a coef of 3 them get it lower if no outliers are found
  out <- boxplot.stats(z_i[r], coef = 3)$out
  if(is_empty(out)) out <- boxplot.stats(z_i[r], coef = 2.5)$out
  if(is_empty(out)) out <- boxplot.stats(z_i[r], coef = 1.5)$out
  if(is_empty(out)) out <- boxplot.stats(z_i[r], coef = 1.1)$out
  
  #identify outliers within vector z_i
  out <- z_i %in% out
  
  # NA to the outliers
  mydata[out, "lon"] <- NA
  mydata[out, "lat"] <- NA
  mydata[out, "z"] <- NA
  
}

sum(is.na(mydata[,c("lon","lat")]))/2

# input lat and lon using closest value to the median of z within the village 
mydata$village <- tolower(mydata$village)
mydata$village <- gsub(" ","",mydata$village)
mydata$village[mydata$village=="adoli"] <- 	"addoli"
p <- unique(mydata$village)

for (i in seq_along(p)){
  
  #data frame with coordinates and z 
  xy <- mydata[mydata$village==p[i] & !is.na(mydata$village), c("lon","lat","z")]
  
  #vector with z minus the mean of z
  z <- xy$z - mean(xy$z, na.rm = TRUE)
  
  #which of these values are closest to 0 (closest to the mean)
  z <- xy[which.min(abs(z)), c("lon","lat")]
  
  #replace NA using the closest value
  mydata[,"lon"] <- ifelse(mydata[,"village"] == p[i] & is.na(mydata[,"lon"]),
                           z[,1], mydata[,"lon"] )
  mydata[,"lat"] <- ifelse(mydata[,"village"] == p[i] & is.na(mydata[,"lat"]),
                           z[,2], mydata[,"lat"] )
  
}

sum(is.na(mydata[,c("lon","lat")]))/2

plot_map(mydata, c("lon", "lat"), "project_name")

# use lonlat from CCI projects and check for lonlat in the same villages
p <- unique(mydata[is.na(mydata$lat), "village"])

lonlatcci <- read.csv("data/raw/villages_geopoints.csv", na.strings = "null")
names(lonlatcci) <- c("district","village","lon","lat")
lonlatcci$village <- tolower(lonlatcci$village)
lonlatcci <- na.omit(lonlatcci)

p <- p[p %in% unique(lonlatcci$village)]

for(i in seq_along(p)) {

  ll <- lonlatcci[lonlatcci$village == p[i], c("lon", "lat")] 
  ll <- ll[1,]

  mydata$lon[mydata$village == p[i]] <- ll$lon
  mydata$lat[mydata$village == p[i]] <- ll$lat
}

sum(is.na(mydata[,c("lon","lat")]))/2


plot_map(mydata, c("lon", "lat"), "project_name")

# try with other dataset
list.files("data/raw")
lonlatcci <- read.csv("data/raw/climmob_india27Feb2018.csv")
lonlatcci$village <- tolower(lonlatcci$village)
lonlatcci$village <- gsub(" ","",lonlatcci$village)
lonlatcci <- lonlatcci[!is.na(lonlatcci$lat) &
                         !is.na(lonlatcci$lon), ]

p <- unique(mydata[is.na(mydata$lat), "village"])
p <- as.vector(na.omit(p))

p <- p[p %in% unique(lonlatcci$village)]

for(i in seq_along(p)) {
  
  ll <- lonlatcci %>% 
    filter(village == p[i])
  
  ll <- ll[1, c("lon","lat")]
  
  mydata$lon[mydata$village == p[i]] <- ll$lon
  mydata$lat[mydata$village == p[i]] <- ll$lat
}

sum(is.na(mydata[,c("lon","lat")]))/2

plot_map(mydata, c("lon", "lat"), "project_name")

# replace NAs by district
p <- sort(unique(mydata$district))
mydata$z <- rowSums(mydata[, c("lon","lat")])

for(i in seq_along(p)) {
  xy <- 
  mydata %>% 
    filter(district == p[i]) %>% 
    select(lon, lat, z) %>% 
    filter(!is.na(lon))
  
  #vector with z minus the mean of z
  z <- xy$z - mean(xy$z, na.rm = TRUE)
  
  #which of these values are closest to 0 (closest to the mean)
  z <- xy[which.min(abs(z)), c("lon","lat")]
  
  #replace NA using the closest value
  mydata[,"lon"] <- ifelse(mydata[,"district"] == p[i] & is.na(mydata[,"lon"]),
                           z[,1], mydata[,"lon"] )
  mydata[,"lat"] <- ifelse(mydata[,"district"] == p[i] & is.na(mydata[,"lat"]),
                           z[,2], mydata[,"lat"] )
  
  
}

sum(is.na(mydata[,c("lon","lat")]))/2

plot_map(mydata, c("lon", "lat"), "project_name")






