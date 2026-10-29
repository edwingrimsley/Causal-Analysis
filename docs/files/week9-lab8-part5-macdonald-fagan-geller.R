# =====================================================================
# Causal Analysis (QMSS, CUNY GC) -- Week 9 Lab, Part 5
# Real data: the same impact zones, and their effect on CRIME and ARRESTS
#
# Replicates MacDonald, Fagan & Geller (2016), "The Effects of Local
# Police Surges on Crime and Arrests in New York City," PLOS ONE 11(6).
# Data: github.com/macdonaldjohn/Impact-Zone-Data (released under the
# Floyd v. City of New York court order). Extracted here to the columns
# we need: lab8-mfg-crime.rds and lab8-mfg-arrests.rds on the course site.
#
# Unit: census block group (anonymized "fid") within a precinct and
# neighborhood (NTA), by month, 2004-2012. The paper reports about
# -12% crime and +53% arrests in impact zones.
#
# Packages: install.packages("fixest")
# =====================================================================

library(fixest)

base <- "https://edwingrimsley.github.io/Causal-Analysis/files/"
get_rds <- function(f) readRDS(gzcon(url(paste0(base, f))))
cr <- get_rds("lab8-mfg-crime.rds")
ar <- get_rds("lab8-mfg-arrests.rds")
# Offline? Download both .rds files and use readRDS(file.choose()) for each.

# treatment  = block group is in an impact zone this month
# neighbor   = block group borders a zone this month (not itself zoned)
# eventneg1, eventneg2 = one and two months before the zone starts
mean(cr$treatment)                               # about 0.11 of rows
table(cr$year, cr$treatment)


# ---------------------------------------------------------------------
# Model 1, as published: Poisson with precinct-by-month fixed effects
# ---------------------------------------------------------------------
m1_crime  <- fepois(crimes  ~ treatment | pct^year^month, cr, cluster = ~pct^year^month)
m1_arrest <- fepois(arrests ~ treatment | pct^year^month, ar, cluster = ~pct^year^month)
etable(m1_crime, m1_arrest)
exp(coef(m1_crime)) - 1      # about -0.12: 12% fewer crimes
exp(coef(m1_arrest)) - 1     # about +0.52: 52% more arrests

# Q1. With only precinct-by-month fixed effects, what comparison produces
#     this estimate? (Hint: which block groups are compared with which,
#     and in which month?) Is that a difference-in-differences?


# ---------------------------------------------------------------------
# Add block-group fixed effects: compare each area with ITSELF over time
# ---------------------------------------------------------------------
m5_crime  <- fepois(crimes  ~ treatment | fid^pct + pct^year^month, cr, cluster = ~pct^year^month)
m5_arrest <- fepois(arrests ~ treatment | fid^pct + pct^year^month, ar, cluster = ~pct^year^month)
etable(m1_crime, m5_crime, m1_arrest, m5_arrest)
# Expect: crime about -6%, arrests about +37%.

# Q2. Adding block-group fixed effects cuts the crime effect in half.
#     What was Model 1 picking up that the fixed effects remove? Which
#     estimate would you report, and why? (The paper's Model 5 uses
#     block-group fixed effects with time splines as a robustness check.)


# ---------------------------------------------------------------------
# Model 2: spillover to neighboring block groups
# ---------------------------------------------------------------------
m2 <- fepois(crimes ~ treatment + neighbor | pct^year^month, cr, cluster = ~pct^year^month)
exp(coef(m2)) - 1            # zones about -12%; neighbors about -3%
# Q3. Why does a spillover to neighbors threaten the comparison in
#     Model 1? Which direction would it bias the zone estimate?


# ---------------------------------------------------------------------
# Model 3: do zoned areas change BEFORE the zone arrives?
# ---------------------------------------------------------------------
m3 <- fepois(crimes ~ eventneg2 + eventneg1 + treatment | pct^year^month, cr,
             cluster = ~pct^year^month)
exp(coef(m3)) - 1            # about -5% and -2% in the two months before
# Q4. Crime is already lower in the two months before designation.
#     Two stories: police surged before the official start, or zones
#     were chosen after a crime spike that was already fading (mean
#     reversion). Which one biases the estimate, and in which direction?


# ---------------------------------------------------------------------
# Connecting back to Legewie & Fagan
# ---------------------------------------------------------------------
# Q5. These are the same zones. Monthly, they switch on and off as
#     phases begin and end. Legewie & Fagan treat education as affected
#     only while a student's area is zoned. MacDonald et al. treat crime
#     the same way. For which outcome is "only while zoned" more
#     plausible, and what would carryover do to each estimate?
