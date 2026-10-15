###############################################################################
# Lab 6: Fixed Effects, Growth Curves, and Matching with Panel Data
# Does a marijuana possession charge change earnings trajectories?
#
# Causal Analysis — CUNY Graduate Center, QMSS · Week 7
#
# Data: National Longitudinal Survey of Youth 1997 (NLSY97), extracts
# prepared for Professor Grimsley's paper on low-level criminal justice
# contact and labor-market trajectories.
#   panel   person-year panel, 2001–2019 (every 2 years)
#   person  one row per respondent (background traits)
# Both come in one file, lab6-nlsy97.rda, which loads straight from the
# course site (no download or folder setup needed).
#
#   Part 1  The panel and who gets "treated"
#   Part 2  Pooled OLS
#   Part 3  Person fixed effects, then two-way fixed effects
#   Part 4  Who identifies a fixed-effects estimate?
#   Part 5  Event study: earnings before and after a charge
#   Part 6  Growth curves (the paper's approach)
#   Part 7  Propensity score matching + fixed effects
#   Part 8  Does the penalty differ by race?
#   Part 9  Robustness: how the outcome is measured
###############################################################################


## ---- setup ----
# install.packages(c("dplyr", "ggplot2", "plm", "lmtest", "sandwich", "lme4", "MatchIt", "estimatr"))
library(dplyr)
library(ggplot2)
library(plm)        # fixed-effects (within) models for panel data
library(lmtest)     # coeftest() for clustered standard errors
library(sandwich)   # cluster-robust variance
library(lme4)       # mixed-effects growth models
library(MatchIt)    # propensity score matching
library(estimatr)   # lm_robust(): weighted FE with clustered SEs


# Helper: coefficient, standard error clustered by person, and p-value
# for a plm model.
cl <- function(model, terms = "charge") {
  V  <- vcovHC(model, cluster = "group", type = "HC1")
  b  <- coef(model)[terms]
  se <- sqrt(diag(V))[terms]
  data.frame(term = terms, estimate = round(b, 3), std.error = round(se, 3),
             p.value = round(2 * pnorm(-abs(b / se)), 3), row.names = NULL)
}


###############################################################################
# PART 1. THE PANEL AND WHO GETS "TREATED"
###############################################################################

## ---- load ----
# One line loads both data frames (panel and person) from the course site.
load(url("https://edwingrimsley.github.io/Causal-Analysis/files/lab6-nlsy97.rda"))
# Offline? Download lab6-nlsy97.rda from the course site and run instead:
# load(file.choose())

panel <- panel |> mutate(charge = marijuana_charge_primary)

dim(panel)    # 89,840 rows = 8,984 people x 10 survey years
dim(person)   # 8,984 people
table(panel$year)

# Look at one person over time
panel |> filter(ID == 46) |> select(ID, year, age, charge, income, log_income)

## ---- treatment-timing ----
# charge is time-varying: 0 before someone's first marijuana charge, 1 from
# then on. Once charged, always "charged" (an absorbing treatment).
timing <- panel |>
  group_by(ID) |>
  summarise(first_charge = if (any(charge == 1)) min(year[charge == 1]) else NA_real_,
            group = case_when(all(charge == 0) ~ "Never charged",
                              all(charge == 1) ~ "Charged before 2001",
                              TRUE             ~ "Charged during panel"))
table(timing$group)
table(timing$first_charge)   # first-charge year among those charged during the panel

## ---- trajectories ----
# Average log income by age for people ever vs. never charged
panel <- panel |> left_join(timing |> select(ID, group, first_charge), by = "ID")

panel |>
  filter(!is.na(log_income)) |>
  mutate(ever = ifelse(group == "Never charged", "Never charged", "Ever charged")) |>
  group_by(ever, age) |>
  summarise(mean_log_income = mean(log_income), n = n(), .groups = "drop") |>
  filter(n >= 50) |>
  ggplot(aes(age, mean_log_income, colour = ever)) +
  geom_line(linewidth = 1) + geom_point() +
  labs(title = "Earnings rise steeply with age for everyone",
       x = "Age", y = "Mean log income", colour = NULL) +
  theme_minimal()

# Q1. The NLSY97 follows people from their late teens into their thirties.
#     Why does that make AGE a confounder for any variable that switches on
#     partway through the panel?


###############################################################################
# PART 2. POOLED OLS
###############################################################################

## ---- pooled ----
est <- panel |> filter(!is.na(log_income))

pooled <- lm(log_income ~ charge + age_c + I(age_c^2) + race_rc + sex_rc + factor(year) +
               youth_delinquency + substance_use_index + mental_health_problem + chronic_condition,
             data = est)
coeftest(pooled, vcov = vcovCL(pooled, cluster = ~ID))["charge", ]

# Q2. Pooled OLS compares charged person-years with uncharged person-years,
#     across different people. List two unmeasured differences between
#     people who do and don't get charged that could bias this estimate.


###############################################################################
# PART 3. PERSON FIXED EFFECTS, THEN TWO-WAY FIXED EFFECTS
###############################################################################

## ---- fe ----
pdat <- pdata.frame(as.data.frame(est), index = c("ID", "year"))

# (a) Person fixed effects only: compare each person with themselves
fe_person <- plm(log_income ~ charge, data = pdat, model = "within", effect = "individual")
cl(fe_person)

# (b) Person + year fixed effects (two-way fixed effects)
fe_twoway <- plm(log_income ~ charge, data = pdat, model = "within", effect = "twoways")
cl(fe_twoway)

# Q3. The person-FE estimate is POSITIVE. Look back at the trajectory plot.
#     Charges switch on as people age into higher earnings. What is (a)
#     comparing, and why does adding year fixed effects fix it?
# Q4. Why can't a person-FE model estimate the effect of race or sex?


###############################################################################
# PART 4. WHO IDENTIFIES A FIXED-EFFECTS ESTIMATE?
###############################################################################

## ---- identification ----
est |>
  group_by(ID) |>
  summarise(switch = n_distinct(charge) > 1) |>
  count(switch)

# Q5. Only people whose charge status CHANGES during the panel move the
#     two-way FE estimate. How many people is that? What happened to the
#     241 people charged before 2001? Who does the estimate describe?


###############################################################################
# PART 5. EVENT STUDY: EARNINGS BEFORE AND AFTER A CHARGE
###############################################################################

## ---- event-study ----
# Event time = years since first charge (surveys are every 2 years).
# Reference period: the survey 2 years BEFORE the charge.
# People charged before 2001 are dropped: we never see them before treatment.
es <- est |>
  filter(group != "Charged before 2001") |>
  mutate(et = year - first_charge,
         e_m6 = as.numeric(!is.na(et) & et <= -6),
         e_m4 = as.numeric(!is.na(et) & et == -4),
         e_p0 = as.numeric(!is.na(et) & et == 0),
         e_p2 = as.numeric(!is.na(et) & et == 2),
         e_p4 = as.numeric(!is.na(et) & et == 4),
         e_p6 = as.numeric(!is.na(et) & et == 6),
         e_p8 = as.numeric(!is.na(et) & et == 8),
         e_p10 = as.numeric(!is.na(et) & et >= 10))

es_p <- pdata.frame(as.data.frame(es), index = c("ID", "year"))
es_mod <- plm(log_income ~ e_m6 + e_m4 + e_p0 + e_p2 + e_p4 + e_p6 + e_p8 + e_p10,
              data = es_p, model = "within", effect = "twoways")
es_tab <- cl(es_mod, c("e_m6", "e_m4", "e_p0", "e_p2", "e_p4", "e_p6", "e_p8", "e_p10"))
es_tab$event_time <- c(-6, -4, 0, 2, 4, 6, 8, 10)
es_tab <- rbind(es_tab, data.frame(term = "ref", estimate = 0, std.error = 0, p.value = NA, event_time = -2))
es_tab

ggplot(es_tab, aes(event_time, estimate)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_vline(xintercept = -1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = estimate - 1.96 * std.error, ymax = estimate + 1.96 * std.error)) +
  labs(title = "Log income relative to 2 years before a first marijuana charge",
       subtitle = "Two-way FE event study; -6 bin = 6+ years before; 10 bin = 10+ years after",
       x = "Years since first charge", y = "Difference in log income") +
  theme_minimal()

# Q6. Are the pre-charge coefficients (-6, -4) near zero? What would a
#     DECLINE in earnings before the charge (an "Ashenfelter dip") mean for
#     the two-way FE estimate in Part 3?
# Q7. These estimates are noisy. How many switchers sit in each bin?
colSums(es[, grep("^e_", names(es))])


###############################################################################
# PART 6. GROWTH CURVES (THE PAPER'S APPROACH)
###############################################################################

## ---- growth ----
# The paper (Grimsley, "NLSY ASR submission draft") estimates linear
# mixed-effects growth models: person-years nested in people, a random
# intercept for each person, linear age, and age x race x charge
# interactions. Table 2 adds covariates in six steps.

# (a) A simple growth model with no race interactions, for comparison
#     with pooled OLS and FE
growth_simple <- lmer(log_income ~ age_c + race_rc + sex_rc + charge + (1 | ID),
                      data = panel, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
round(summary(growth_simple)$coefficients, 3)

# (b) REPLICATE TABLE 2, MODEL 1 of the paper
growth_m1 <- lmer(log_income ~ age_c * race_rc * charge + sex_rc + (1 | ID),
                  data = panel, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
m1 <- summary(growth_m1)$coefficients
published_m1 <- c("(Intercept)" = 7.410, "age_c" = 0.151, "race_rcBlack" = -1.278,
                  "race_rcHispanic" = -0.660, "charge" = -0.528, "sex_rcMale" = 0.799,
                  "race_rcBlack:charge" = -1.524, "race_rcHispanic:charge" = -0.429,
                  "age_c:race_rcBlack:charge" = -0.022, "age_c:race_rcHispanic:charge" = 0.070)
data.frame(ours = round(m1[names(published_m1), "Estimate"], 3), published = published_m1)
nobs(growth_m1)   # 64,849 person-years

# (c) REPLICATE ALL SIX MODELS OF TABLE 2. Each adds a block of controls.
panel <- panel |>
  mutate(youth_delinquency = factor(youth_delinquency,
                                    c("No Deliquency Incidents", "Low Deliquency Incidents",
                                      "Medium Deliquency Incidents", "High Deliquency Incidents")))
blocks <- c(M1 = "",
            M2 = "+ highschool_dropout",
            M3 = "+ youth_delinquency + substance_use_index",
            M4 = "+ mental_health_problem + chronic_condition",
            M5 = "+ assault_charge + robbery_charge + burglary_charge",
            M6 = "+ incarceration")
rhs <- "log_income ~ age_c * race_rc * charge + sex_rc"
table2 <- list()
for (m in names(blocks)) {
  rhs <- paste(rhs, blocks[m])
  table2[[m]] <- lmer(as.formula(paste(rhs, "+ (1 | ID)")), data = panel, REML = FALSE,
                      control = lmerControl(optimizer = "bobyqa"))
}
published_bxc <- c(-1.524, -1.610, -1.680, -1.681, -1.524, -1.379)
data.frame(model = names(table2),
           N = sapply(table2, nobs),
           charge = round(sapply(table2, function(m) fixef(m)["charge"]), 3),
           black_x_charge = round(sapply(table2, function(m) fixef(m)["race_rcBlack:charge"]), 3),
           published_black_x_charge = published_bxc,
           hispanic_x_charge = round(sapply(table2, function(m) fixef(m)["race_rcHispanic:charge"]), 3),
           row.names = NULL)

# Q8. The simple growth model's charge coefficient is close to POOLED OLS,
#     not to two-way FE. Random-intercept models use both between-person and
#     within-person variation. What do they assume about the random
#     intercepts and charge that FE does not? Why might that fail here?
# Q9. The paper's key result is the Black x charge interaction (-1.524 in
#     Model 1, -1.379 in Model 6). Which block of controls moves it most? Compare it with the FE race interactions
#     in Part 8. What does each model buy you, and what does it cost?

# Optional: Hausman test comparing FE and RE
re <- plm(log_income ~ charge + factor(year), data = pdat, model = "random")
fe <- plm(log_income ~ charge + factor(year), data = pdat, model = "within")
phtest(fe, re)


###############################################################################
# PART 7. PROPENSITY SCORE MATCHING + FIXED EFFECTS
###############################################################################

## ---- psm ----
# Build a person-level file of background traits measured in adolescence
ever <- timing |> mutate(ever = as.numeric(group != "Never charged")) |> select(ID, ever)

pp <- person |>
  left_join(ever, by = "ID") |>
  mutate(race_rc = factor(race, 0:2, c("White", "Black", "Hispanic")),
         male    = as.numeric(sex == "Male"),
         parent_income = ifelse(is.na(parents_household_income), "Missing", parents_household_income),
         delinquency = factor(youth_delinquency,
                              c("No Deliquency Incidents", "Low Deliquency Incidents",
                                "Medium Deliquency Incidents", "High Deliquency Incidents")),
         chronic = as.numeric(chronic_condition == "Chronic Condition")) |>
  filter(!is.na(race_rc))

# Covariates: adolescent background only. We deliberately leave out high
# school dropout and marriage. (Q10 asks why.)
# 2:1 nearest-neighbor matching (two controls per charged respondent), as in
# Professor Grimsley's dissertation, with a caliper of 0.2 SD of the score.
m_out <- matchit(ever ~ race_rc + male + parent_income + delinquency +
                   substance_use_index + mental_health_problem + chronic,
                 data = pp, method = "nearest", distance = "glm", caliper = 0.2, ratio = 2)
summary(m_out)$nn

# Balance: standardized mean differences before and after
bal <- summary(m_out, standardize = TRUE)
round(cbind(before = bal$sum.all[, "Std. Mean Diff."],
            after  = bal$sum.matched[, "Std. Mean Diff."]), 3)
plot(summary(m_out), abs = TRUE)   # Love plot

# How many controls did each charged respondent get within the caliper?
md <- match.data(m_out)
table(table(md$subclass[md$ever == 0]))

# Two-way FE on the matched sample, using the matching weights.
# (Controls matched to a treated person who found only one partner get
# more weight.) lm_robust handles weights with clustered SEs.
matched_panel <- est |> inner_join(md |> select(ID, weights), by = "ID")
fe_matched <- lm_robust(log_income ~ charge, data = matched_panel, weights = weights,
                        fixed_effects = ~ ID + year, clusters = ID, se_type = "stata")
summary(fe_matched)$coefficients["charge", ]

# Compare: 1:1 matching (the more common default)
m_11 <- matchit(ever ~ race_rc + male + parent_income + delinquency +
                  substance_use_index + mental_health_problem + chronic,
                data = pp, method = "nearest", distance = "glm", caliper = 0.2, ratio = 1)
p11 <- est |> inner_join(match.data(m_11) |> select(ID, weights), by = "ID")
summary(lm_robust(log_income ~ charge, data = p11, weights = weights,
                  fixed_effects = ~ ID + year, clusters = ID, se_type = "stata"))$coefficients["charge", ]

# Q10. Why might high school dropout and marriage be bad controls here?
#      (Week 2: post-treatment variables.)
# Q11. Matching balances ADOLESCENT traits; fixed effects remove STABLE
#      traits. What kind of confounder would survive both?
# Q12. Two-way FE now compares charged people with matched never-charged
#      people of similar background, instead of with everyone. Why might
#      that make the "parallel trends" assumption more plausible?
# Q12b. Compare 2:1 with 1:1 matching. What do you gain from a second
#      control per charged person, and what do you risk? (Bias vs. variance.)


###############################################################################
# PART 8. DOES THE PENALTY DIFFER BY RACE?
###############################################################################

## ---- race ----
fe_race <- plm(log_income ~ charge + charge:race_rc, data = pdat, model = "within", effect = "twoways")
cl(fe_race, c("charge", "charge:race_rcBlack", "charge:race_rcHispanic"))

# Q13. Race itself drops out of a FE model, but charge x race does not. Why?
#      Compare these interactions with the growth-model results in Part 6.


###############################################################################
# PART 9. ROBUSTNESS: HOW THE OUTCOME IS MEASURED
###############################################################################

## ---- robustness ----
# log_income is log(earnings + 1), so it is 0 when earnings are 0 and the
# outcome mixes "any income" with "how much." (The paper separates the two
# with a two-part model in its Table 3.) And in 2007 no one has missing income but 44% report zero.
panel |> group_by(year) |>
  summarise(share_zero = mean(income == 0, na.rm = TRUE), share_missing = mean(is.na(income)))

no07 <- pdata.frame(as.data.frame(filter(est, year != 2007)), index = c("ID", "year"))
rob <- rbind(
  cbind(model = "Log income, all years",      cl(plm(log_income ~ charge, pdat, model = "within", effect = "twoways"))),
  cbind(model = "Log income, drop 2007",      cl(plm(log_income ~ charge, no07, model = "within", effect = "twoways"))),
  cbind(model = "Any income (0/1), drop 2007", cl(plm(I(as.numeric(income > 0)) ~ charge, no07, model = "within", effect = "twoways"))),
  cbind(model = "Income in $1,000s, drop 2007", cl(plm(I(income / 1000) ~ charge, no07, model = "within", effect = "twoways")))
)
rob

# Q14. Which of these outcomes would you report as the main result, and why?
#      What does the "any income" result suggest about where the penalty
#      comes from?

## ---- end ----
