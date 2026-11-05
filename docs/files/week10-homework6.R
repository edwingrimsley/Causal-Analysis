###############################################################################
# Homework 6 (Problem Set 6): Did Head Start save children's lives?
# Causal Analysis — CUNY Graduate Center, QMSS · Week 10 · due November 12
#
# Data: Ludwig & Miller (2007), "Does Head Start Improve Children's Life
# Chances? Evidence from a Regression Discontinuity Design," QJE 122(1),
# as distributed with Cattaneo, Titiunik & Vazquez-Bare (2017, JPAM).
# One row per U.S. county.
#
#   povrate60  county poverty rate in 1960 (running variable)
#              The 300 poorest counties, poverty rate >= 59.1984%, got
#              federal help writing Head Start grant applications in 1965.
#   mort_age59_related_postHS
#              deaths per 100,000 children aged 5-9, 1973-1983, from causes
#              Head Start's health services could affect (the outcome)
#   See the handout for the other variables.
#
# This script loads the data. Write the rest, using Lab 9 as your guide.
###############################################################################

# install.packages(c("rdrobust", "rddensity", "ggplot2"))
library(rdrobust)
library(rddensity)
library(ggplot2)

base <- "https://edwingrimsley.github.io/Causal-Analysis/files/"
hs <- read.csv(paste0(base, "hw6-headstart.csv"))
# Offline? Download the CSV from the course site and use:
# hs <- read.csv(file.choose())

cutoff <- 59.1984
hs$R <- hs$povrate60 - cutoff            # running variable, centered at 0
hs$Y <- hs$mort_age59_related_postHS     # main outcome

summary(hs$povrate60)
sum(hs$R >= 0, na.rm = TRUE)            # the 300 counties that got help

# Your code starts here (Q1-Q6).
