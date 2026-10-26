###############################################################################
# Homework 5 (Problem Set 5): Did COVID lockdowns cause a sourdough craze?
# Causal Analysis — CUNY Graduate Center, QMSS · Week 8 · due October 29
#
# Adapted from Nick Huntington-Klein's Chapter 18 coding homework for
# The Effect (github.com/NickCH-K/TheEffectAssignments), which builds on an
# example by Peter Nencka. Data: Google Trends search interest in the US,
# Jan 1 – Aug 1, 2020, for "sourdough" and three comparison terms.
#
# This script loads the data and sets up the variables. Write the rest.
###############################################################################

# install.packages(c("fixest", "dplyr", "ggplot2"))
library(fixest)
library(dplyr)
library(ggplot2)

# A file window opens: pick sourdough_trends.csv
sr <- read.csv(file.choose())

# Keep the variables we need and turn the date text into a real date
sr <- sr |>
  select(date, hits, keyword) |>
  mutate(date = as.Date(substr(date, 1, 10)))

# "Treatment" date: the start of COVID lockdowns, March 15, 2020
event <- as.Date("2020-03-15")

# Periods of about one month, counted from March 15:
#   -1 = Feb 14 – Mar 14 (the month just before), 0 = Mar 15 – Apr 13, ...
sr <- sr |>
  mutate(period = floor(as.numeric(date - event) / 30.44))
sr |> group_by(period) |> summarize(from = min(date), to = max(date))

# Your code starts here (Q1-Q8).
