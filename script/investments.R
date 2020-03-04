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






# get rankings in a different dataframe
R <- mydata[items]
names(R) <- gsub("rank_","",names(R))
R[R < 0] <- 0
R[R > 7] <- 0

# also explanatory vars
exp_var <- mydata[exp_var]

# put characters as factor
exp_var[2:9] <- lapply(exp_var[2:9], function(X){
  as.factor(as.character(X))
})

# check for missing data in these two dataframes 
sum(is.na(R))
sum(is.na(exp_var))

# refresh item names 
items <- names(R)

# convert rankings into a matrix and a PlackettLuce object of class 'rankings'
R <- as.matrix(R, dimnames = list(1:n, items))
R <- PlackettLuce::as.rankings(R, ncol = length(items), byrow = TRUE)

# ================================================
# ================================================
# Forward selection with PlackettLuce model

# tranform rankings into grouped rankings 
G <- group(R, index = 1:n)

# combine grouped rankings with explanatory variables
input <- cbind(G, exp_var)

# define model parameters
a <- 0.1 # alpha
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
              minsize = 20)

plt

plot(plt, abbr = 3)

d <- summarise_dominance(R)
plot(d)

f <- summarise_favorite(R)
plot(f)

plot(plt)


worst_regret(plt)


