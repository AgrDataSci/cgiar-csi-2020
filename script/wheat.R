# ..........................................................
# ..........................................................
# Analyse wheat data from station trials in India
## Kauê de Sousa
### Inland Norway University
### Bioversity International

# ..........................................................
# ..........................................................
# Packages ####
library("devtools")
# install_github("agrobioinfoservices/gosset", upgrade = "never")
library("gosset")
library("climatrends")
library("PlackettLuce")
library("qvcalc")
library("ggplot2")
library("abind")
library("foreach")
library("doParallel")
source(paste0("https://raw.githubusercontent.com/agrobioinfoservices/",
              "cgiar-csi-2020/master/script/functions.R"))

# ........................................................
# ........................................................
# Read data ####
# dt <- read.csv(paste0("https://raw.githubusercontent.com/agrobioinfoservices/",
#                       "cgiar-csi-2020/master/data/wheat_pvs.csv"),
#                stringsAsFactors = FALSE)
dt <- read.csv("data/wheat_pvs.csv", stringsAsFactors = F)
head(dt)
tail(dt)




# 255 points across northern India
length(unique(dt$id))
plot_map(dt, c("lon","lat"), map.types = "OpenTopoMap")

# with 31 bread wheat varieties
unique(dt$variety)


# here we are going to use the yield data,
# lets see how it is distributed
boxplot(dt$yield)

boxplot(dt$yield ~ dt$project_id)

boxplot(dt$yield ~ dt$id)

# this problem can be solved using rankings,
# assuming tha the permutation is correct
R <- rank_numeric(data = dt,
                  items = "variety",
                  input = "yield",
                  id = "id")

head(R)

# Let's check the network
network(R, vertex.size = 15)

# Let's check for its adjacency 
adjacency(R)

# and connectivity
connectivity(R)

# all varieties belong to the same network


# Here we can look for the dominance of each items with some 
# summary and visualization tools from gosset
# higher values indicates higher dominance of Player1 over Player2
dom <- summarise_dominance(R)
plot(dom)

# summarise_favorite shows the proportion of times the item was most or least 
# favourite. 
fav <- summarise_favorite(R)
plot(fav)

# .......................................................
# .......................................................
# PlackettLuce model ####
# Now lets fit a simple Plackett-Luce model with the rankings object
mod <- PlackettLuce(R)

summary(mod)

plot(qvcalc(mod))

# .......................................................
# .......................................................
# Add covariates ####
dt$planting_date <- as.Date(dt$planting_date, format = "%Y-%m-%d")
dt$flowering50perc_date <- as.Date(dt$flowering50perc_date, format = "%Y-%m-%d")

dt$ts <- as.integer(dt$flowering50perc_date - dt$planting_date)


geoinput <- dt[,c("id","lon","lat","planting_date","year","ts")]

geoinput <- geoinput[!duplicated(geoinput$id), ]


temp <- temperature(geoinput[,c("lon","lat")],
                    day.one = geoinput$planting_date,
                    span = geoinput$ts)

temp <- read.csv("data/temperature.csv")

# # geographic and time boundaries are too large for NASAPOWER
# s <- ifelse(geoinput$lon < 80, 1, 2)
# s <- as.integer(as.factor(paste0(s, geoinput$year)))
# 
# k <- unique(s)
# 
# temp <- as.data.frame(matrix(NA, 
#                              ncol = 8,
#                              nrow = nrow(geoinput)))
# 
# for (i in seq_along(k)) {
#   print(i)
#   g_i <- geoinput[s == k[i], ]
#   
#   tp <- temperature(g_i[, c("lon","lat")],
#                     day.one = g_i$planting_date,
#                     span = g_i$ts)
#   
#   temp[s == k[i], ] <- tp
#     
# }
# 
# names(temp) <- names(tp)
# 

geoinput <- cbind(geoinput, temp)


geoinput <- geoinput[, -match(c("id","ts","TR","CFD"), names(geoinput))]


G <- group(R, 1:length(R))


dat <- cbind(G, geoinput)

plot(density(dat$maxDT))
boxplot(dat$maxDT)

plot(density(dat$maxNT))
boxplot(dat$maxNT)

# .......................................................
# .......................................................
# Forward selection ####

# We use gosset::forward() with runs a forward stepwise selection
# with cross-validation (an internall call to gosset::crossvalidation())
# to select the best model considering the rankings 
# and explanatory variables it is exposed to

# Forward starts with an empty model (intercept-only) and add variables
# in each step 

# forward() flavours, type ?forward for details
# n rows in dat
n <- dim(dat)[[1]]
# folds based on the season
seasons <- as.integer(as.factor(dat$year))
table(seasons)
# number of folds
nk <- max(seasons)
# packages to export to parallels
pkgs <- "PlackettLuce"
# goodness-of-fit to select by
gof <- "deviance"
# number of cores
ncor <- abs(detectCores() / 2)

# PlackettLuce flavours passed to ...
# type ?PlackettLuce for details
# minimum size of each node 
mins <- round((n*0.3), -1)
# bonferroni correction
bonf <- TRUE
# the significance level for spliting the data into nodes
a <- 0.1
# pseudo rankings 
pr <- 1.5

# and then the forward stepwise selection
f <- forward(G ~ ., 
             data = dat,
             k = nk,
             folds = seasons,
             select.by = gof,
             ncores = ncor, 
             packages = pkgs,
             minsize = mins,
             bonferroni = bonf,
             alpha = a,
             npseudo = pr)
