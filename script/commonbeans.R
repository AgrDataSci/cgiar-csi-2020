# ..........................................................
# ..........................................................
# Analyse common beans data from tricot trials in Nicaragua
## Kauê de Sousa
### Inland Norway University
### Bioversity International


# ..........................................................
# ..........................................................
# Packages ####
library("devtools")
# install_github("agrobioinfoservices/gosset", upgrade = "never")
library("gosset")
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
# dt <- read.csv("data/commonbeans_tricot.csv", stringsAsFactors = FALSE)
dt <- read.csv(paste0("https://raw.githubusercontent.com/agrobioinfoservices/",
                      "cgiar-csi-2020/master/data/commonbeans_tricot.csv"),
               stringsAsFactors = FALSE)
head(dt)
tail(dt)

# 842 points across Nicaragua
plot_map(dt, c("lon","lat"), map.types = "OpenTopoMap")

# with 10 common bean varieties
unique(unlist(dt[, paste0("variety_", letters[1:3])]))

# ........................................................
# ........................................................
# Some exploration ####
# coerce into rankings using
# gosset::rank_tricot
R <- rank_tricot(data = dt,
                 items = paste0("variety_", letters[1:3]),
                 input = c("best","worst"))

head(R)
tail(R)

# We can check how these items (varieties are connected)
# even with partial rankings, we can see that all varieties connect
# together. This is due to the controlled randomisation performed 
# in ClimMobTools::randomise, when a tricot project is set up
network(R, vertex.size = 30)

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

# This fit confirms the findings in the exploratory analysis
# placing INTA Fuerte Sequia, INTA Centro Sur, BRT103-182 and INTA Rojo
# as the best varieties 


# We can plot the qasi-variance se using the package qvcalc
plot(qvcalc(mod))

# Now lets add the Local variety
R <- rank_tricot(data = dt,
                 items = paste0("variety_", letters[1:3]),
                 input = c("best","worst"),
                 additional.rank = dt[, paste0("var_", letters[1:3])])

# The Local variety is in the middle of the network since it is ranked 
# by all participants
network(R)

fav <- summarise_favorite(R)
plot(fav)

# .......................................................
# .......................................................
# PlackettLuce Trees ####

# For the PLT we use an object of class "grouped_rankings" 
G <- rank_tricot(data = dt,
                 items = paste0("variety_", letters[1:3]),
                 input = c("best","worst"),
                 additional.rank = dt[, paste0("var_", letters[1:3])],
                 group = TRUE)

# Same as a "rankings" object but "grouped_rankings" has an index that 
# allows the rankings to be linked to explanatory variables
head(G)
head(R)

# Lets use these variables
# These variables were previous computed using MODIS and CHIRPS data
# and the package climatrends
expvar <- c("lon","lat","planting_date","planting_day",
            "season","year","GDD","minNT","maxNT","MLDS",
            "MLWS","Rx5day")

# and combine the "grouped_rankings" and explanatory variables into a new 
# object
dat <- cbind(G, dt[, expvar])

head(dat)

# Too many explanatory variables
# Let's create a model that best explain these rankings and 
# that can be extrapolated across seasons (years)

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
seasons <- as.integer(as.factor(dt$season))
# number of folds
nk <- max(seasons)
# packages to export to parallels
pkgs <- "PlackettLuce"
# goodness-of-fit to select by
gof <- "deviance"
# number of cores
ncor <- 2

# PlackettLuce flavours passed to ...
# type ?PlackettLuce for details
# minimum size of each node 
mins <- round((n*0.3), -1)
# bonferroni correction
bonf <- TRUE
# the significance level for spliting the data into nodes
a <- 0.01
# pseudo rankings 
pr <- 1.5

b <- Sys.time()
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

e <- Sys.time()

print(e-b)

# these are the cross-validation estimates from the forward selection
print(f)

# the best model is 
model <- as.formula(f$raw$call)

# Lets fit a PlackettLuce tree with this model
plt <- pltree(model,
              data = dat,
              minsize = mins,
              npseudo = pr,
              bonferroni = bonf,
              alpha = a)

summary(plt)

predict(plt)

coef(plt, log = FALSE)

# we can also plot the nodes with error bars using 
# gosset::plot_nodes()
# unfortunately it does not construct the tree in a single plot
p <- plot_nodes(plt)

p[[1]]
p[[2]]

# Finally we can compute the worst regret values
# using gosset::worst_regret() which combine the estimates from 
# all nodes and compute the worst regret values
worst_regret(plt)

