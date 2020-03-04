# ..........................................................
# ..........................................................
# Investment choices and social impact investment
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
dt <- read.csv(paste0("https://raw.githubusercontent.com/agrobioinfoservices/",
                      "cgiar-csi-2020/master/data/impact_invest.csv"))

head(dt)
str(dt)

# get rankings in a different dataframe
options <- grepl("rank_", names(dt))

R <- dt[options]
names(R) <- gsub("rank_","",names(R))

head(R)

items <- names(R)

n <- dim(R)[[1]]

# convert rankings into a matrix and a PlackettLuce object of class 'rankings'
R <- as.matrix(R, dimnames = list(1:n, items))
R <- PlackettLuce::as.rankings(R, ncol = length(items), byrow = TRUE)

# ================================================
# ================================================
# Forward selection with PlackettLuce model

# tranform rankings into grouped rankings 
G <- group(R, index = 1:n)

# combine grouped rankings with explanatory variables
input <- cbind(G, dt[, !options])

head(input)

# define model parameters
a <- 0.1 
minsize <- 20 

mod <- forward(G ~ .,
               data = input,
               k = 10,
               ncores = 8,
               minsize = minsize,
               alpha = a,
               packages = "PlackettLuce", 
               seed = 1234)

mod


plt <- pltree(as.formula(mod$raw$call), 
              data = input,
              alpha = a,
              minsize = minsize)

plt

plot(plt)


worst_regret(plt)


