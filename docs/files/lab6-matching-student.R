# ============================================================
# Week 6 Lab: Matching and Propensity Scores
# Question: Did job training raise earnings?
# Data: LaLonde job-training data (comes with the MatchIt package)
# ============================================================

# ---- 0. Setup (run once) -----------------------------------
# install.packages(c("MatchIt", "cobalt", "chk"), dependencies = TRUE)
# (if you see "no package called chk", run install.packages("chk"))
library(MatchIt)   # does the matching
library(cobalt)    # makes balance plots

data("lalonde")
head(lalonde)
# treat  = 1 if the person got job training, 0 if not
# re78   = earnings in 1978 (the outcome we care about)
# age, educ, race, married, nodegree, re74, re75 = background traits
#   (re74 and re75 = earnings BEFORE the program)

# ---- Step 1. Are the two groups alike? ---------------------
# Average of each trait, for trained (1) vs not trained (0)
aggregate(cbind(age, educ, married, nodegree, re74, re75) ~ treat,
          data = lalonde, FUN = mean)
# Look at re74 and re75. Who earned more BEFORE the program?

# ---- Step 2. The naive comparison --------------------------
mean(lalonde$re78[lalonde$treat == 1]) -
  mean(lalonde$re78[lalonde$treat == 0])
# The real experiment found training RAISED earnings by about $1,794.
# Is the naive number even pointing the right way?

# ---- Step 3. Balance BEFORE matching -----------------------
# method = NULL means "don't match yet, just show me the balance"
before <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
                  data = lalonde, method = NULL, distance = "glm")
summary(before)
# DELIVERABLE 1: save the balance table above
# Std. Mean Diff. = how far apart the groups are (0 = identical).
# Rule of thumb: bigger than 0.1 (in absolute value) is a problem.

# ---- Step 4. Match -----------------------------------------
# For each trained person, find the untrained person with the
# most similar propensity score (= chance of getting training,
# given their traits).
m1 <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
              data = lalonde, method = "nearest", distance = "glm")
m1
summary(m1)   # look at "Summary of Balance for Matched Data"

# A picture of balance before vs after (dots should move toward 0)
love.plot(m1, binary = "std", thresholds = c(m = 0.1))
# DELIVERABLE 2: save this plot.  DELIVERABLE 3: save the Sample Sizes table from summary(m1)

# ---- Step 5. Estimate the effect on the matched data -------
md <- match.data(m1)
fit <- lm(re78 ~ treat, data = md, weights = weights)
summary(fit)
# The coefficient on treat is our matched estimate.
# DELIVERABLE 4: fill in the estimates table (naive, matched, caliper)
# Compare to: naive estimate (Step 2) and the $1,794 benchmark.

# ---- Step 6. Try one change (your turn) --------------------
# Add a caliper: refuse matches that are too far apart.
m2 <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
              data = lalonde, method = "nearest", distance = "glm",
              caliper = 0.2)
m2                      # how many treated people were dropped?
summary(m2)
fit2 <- lm(re78 ~ treat, data = match.data(m2), weights = weights)
coef(fit2)["treat"]
# DELIVERABLE 5: save the Sample Sizes table from summary(m2)

# ---- Step 7. Extension: more than one twin per person ------
# replace = TRUE lets an untrained person be used more than once
m3 <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
              data = lalonde, method = "nearest", distance = "glm",
              ratio = 2, replace = TRUE, caliper = 0.2)
m3
summary(m3)
coef(lm(re78 ~ treat, data = match.data(m3), weights = weights))["treat"]

# Now try 3 twins
m3b <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
               data = lalonde, method = "nearest", distance = "glm",
               ratio = 3, replace = TRUE, caliper = 0.2)
m3b
summary(m3b)
coef(lm(re78 ~ treat, data = match.data(m3b), weights = weights))["treat"]

# ---- Step 8. Extension: force exact balance on key traits --
# Trained people are only matched to untrained people with the
# SAME race and marital status.
m4 <- matchit(treat ~ age + educ + race + married + nodegree + re74 + re75,
              data = lalonde, method = "nearest", distance = "glm",
              exact = ~ race + married, caliper = 0.2)
m4
summary(m4)
coef(lm(re78 ~ treat, data = match.data(m4), weights = weights))["treat"]
love.plot(m4, binary = "std", thresholds = c(m = 0.1))
# DELIVERABLE 6: save the Sample Sizes tables from summary(m3) and summary(m4)
# DELIVERABLE 7: save the love plot for m4
# DELIVERABLE 4 now also includes the m3, m3b and m4 estimates

# ---- Step 9. Regression with and without matching ----------
covs <- "age + educ + race + married + nodegree + re74 + re75"
f0 <- re78 ~ treat
f1 <- as.formula(paste("re78 ~ treat +", covs))

models <- list(
  "No matching, no controls" = lm(f0, data = lalonde),
  "No matching + controls"   = lm(f1, data = lalonde),
  "Caliper (m2)"             = lm(f0, data = match.data(m2),  weights = weights),
  "Caliper (m2) + controls"  = lm(f1, data = match.data(m2),  weights = weights),
  "3 twins (m3b)"            = lm(f0, data = match.data(m3b), weights = weights),
  "3 twins (m3b) + controls" = lm(f1, data = match.data(m3b), weights = weights),
  "Exact (m4)"               = lm(f0, data = match.data(m4),  weights = weights),
  "Exact (m4) + controls"    = lm(f1, data = match.data(m4),  weights = weights)
)

results <- data.frame(
  estimate  = sapply(models, function(m) coef(m)["treat"]),
  std_error = sapply(models, function(m) summary(m)$coefficients["treat", "Std. Error"])
)
round(results, 0)
# The real experiment found about +1,794.
# (Standard errors ignore the matching step, so they are rough.)
# DELIVERABLE 8: save this table

# TAKEAWAYS
# - A well-balanced match and a regression with controls can agree. That
#   agreement is reassuring: the answer doesn't hinge on one method.
# - 1:1 matching is not always optimal. Compare several specifications
#   (caliper, more twins, exact matching) and check balance for each.
