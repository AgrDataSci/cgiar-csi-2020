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

