## =====================================================================
## Week 10 Lab: Regression Discontinuity in Time (RDiT)
## Did use-of-force policies reduce police use of force?
## Five cities: Chicago, Minneapolis, Indianapolis, Austin, Cincinnati
##
## Data: uof_city_month.csv  (city x month x force type x subject race)
##       built from each department's open use-of-force data.
##
## The design: time is the running variable. The cutoff is the month a
## city's use-of-force policy took effect. If nothing else changed at
## that exact moment, a jump in force right at the cutoff is the
## policy's effect.
##
## Key reading: Hausman & Rapson (2018), "Regression Discontinuity in
## Time: Considerations for Empirical Applications," Annual Review of
## Resource Economics 10: 533-552.
## =====================================================================

## ---- 0. Setup --------------------------------------------------------
# install.packages(c("rdrobust", "dplyr", "tidyr", "ggplot2",
#                    "data.table", "lubridate", "sandwich", "lmtest"))
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(lubridate)
library(sandwich)
library(lmtest)
has_rdrobust <- requireNamespace("rdrobust", quietly = TRUE)
if (has_rdrobust) library(rdrobust)

## Download uof_city_month.csv from the course site first:
##   https://edwingrimsley.github.io/Causal-Analysis/files/uof_city_month.csv
## A window opens: pick uof_city_month.csv wherever you saved it.
## (Falls back to the working directory when run non-interactively.)
DATA_FILE    <- if (interactive()) file.choose() else "uof_city_month.csv"
WINDOW       <- 36      # months on each side of the cutoff kept in the sample
H            <- 12      # default bandwidth (months) for the hand-built RD
INCLUDE_2020 <- FALSE   # COVID + post-Floyd 2020 contaminates every city

## Policy effective months (the cutoffs).
## "verified" = effective date confirmed in a public source.
## The others are best available and are part of the lesson (Part 7).
policy_dates <- tibble::tribble(
  ~city,          ~policy_date,  ~date_status,
  "Chicago",      "2017-10-01",  "verified",   # G03-02 suite, effective mid-Oct 2017
  "Minneapolis",  "2016-08-01",  "verified",   # sanctity of life / de-escalation / duty to intervene
  "Indianapolis", "2016-08-01",  "proposal",   # ACLU proposal Table 1
  "Austin",       "2016-10-01",  "uncertain",  # year known (2016), month not
  "Cincinnati",   "2015-10-01",  "uncertain"   # year known (2015), month not
) %>% mutate(policy_date = as.Date(policy_date))

## ---- 1. Load data and build the analysis file -----------------------
raw <- fread(DATA_FILE, data.table = FALSE) %>%
  mutate(month = as.Date(month)) %>%
  select(city, month, force_type, race, n)       # drop old cutoff columns

## Re-attach cutoffs and build the running variable:
##   x = months since the policy took effect (0 = first treated month)
raw <- raw %>%
  left_join(policy_dates, by = "city") %>%
  mutate(x    = (year(month) - year(policy_date)) * 12 +
                (month(month) - month(policy_date)),
         post = as.integer(x >= 0))

## Collapse to one row per city-month, with total, Black, and White counts
cm <- raw %>%
  group_by(city, month, x, post, date_status) %>%
  summarise(total = sum(n),
            black = sum(n[race == "Black"]),
            white = sum(n[race == "White"]),
            .groups = "drop") %>%
  mutate(black_share = black / pmax(black + white, 1),
         cal_month   = month(month))

## Drop each city's "phase-in" months: leading months where reporting
## had not started yet (e.g., Indianapolis, Jan-Jul 2014, ~0 reports).
cm <- cm %>%
  group_by(city) %>% arrange(month) %>%
  filter(cumsum(total >= 0.2 * median(total)) > 0) %>%
  ungroup()

## Outcome: log(count + 1). Within a city, population barely moves over
## a few years, so a jump in log counts is approximately a jump in rates.
## A coefficient of -0.10 means roughly a 10% drop in monthly force.
cm <- cm %>%
  mutate(y_total = log(total + 1),
         y_black = log(black + 1),
         y_white = log(white + 1))

## Analysis sample: +/- WINDOW months around each city's cutoff
an <- cm %>%
  filter(abs(x) <= WINDOW) %>%
  filter(INCLUDE_2020 | month < as.Date("2020-03-01"))

cat("\nMonths in the analysis sample, by city and side of the cutoff:\n")
print(table(an$city, ifelse(an$post == 1, "post", "pre")))

## (Force-type categories are NOT comparable across cities -- e.g.,
## Chicago's export misaligns its weapon flags -- so this lab uses total
## force and subject race, not force type.)


## =====================================================================
## PART 0. Descriptive tables, every city
## =====================================================================
## Always describe the data before estimating anything. Four tables:
##   D1 design and coverage, D2 monthly force by subject race,
##   D3 population and force rates, D4 force-type mix.
## All use the analysis window (+/- 36 months, before March 2020).
## Tables are printed and saved as CSV files (desc_*.csv).

policy_elements <- c(
  Austin       = "Force continuum; warning before shooting; duty to intervene",
  Chicago      = "7 of 8 Campaign Zero elements (all but exhaust alternatives)",
  Cincinnati   = "All 8 Campaign Zero elements",
  Indianapolis = "De-escalation; force continuum; deadly force, vehicles, chokeholds",
  Minneapolis  = "De-escalation; force continuum; duty to intervene")

## D1. Design and coverage
d1 <- cm %>%
  group_by(city) %>%
  summarise(data_from = format(min(month), "%b %Y"),
            data_to   = format(max(month), "%b %Y"), .groups = "drop") %>%
  left_join(policy_dates, by = "city") %>%
  left_join(an %>% group_by(city) %>%
              summarise(months_pre = sum(post == 0), months_post = sum(post == 1),
                        .groups = "drop"), by = "city") %>%
  mutate(policy_month = format(policy_date, "%b %Y"),
         elements = policy_elements[city]) %>%
  select(city, data_from, data_to, policy_month, date_status,
         months_pre, months_post, elements)
cat("\nTable D1. Design and data coverage\n"); print(as.data.frame(d1))

## D2. Monthly use of force by subject race, pre vs. post
d2 <- an %>%
  group_by(city, period = ifelse(post == 1, "post", "pre")) %>%
  summarise(reports    = sum(total),
            mean_total = mean(total), sd_total = sd(total),
            mean_black = mean(black), mean_white = mean(white),
            pct_black  = 100 * sum(black) / sum(total),
            pct_white  = 100 * sum(white) / sum(total),
            .groups = "drop") %>%
  arrange(city, desc(period)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))
cat("\nTable D2. Monthly use-of-force reports by subject race (+/- 36 months)\n")
print(as.data.frame(d2))

## D3. Population and force rates per 100,000 residents (monthly)
## CAUTION: the population columns in this file are COUNTY totals
## (e.g., Chicago = Cook County, ~5.2M; Cincinnati = Hamilton County),
## not city populations. Levels of these rates are therefore too low,
## by different amounts in each city, and should not be compared across
## cities. Within a city, the Black/White RATIO is less affected, but
## still reflects county, not city, demographics.
## Second caution: Census "White" includes Hispanic Whites, but some
## departments (Austin, Indianapolis) record Hispanic subjects separately,
## so the White rate is understated there and the B/W ratio overstated.
pop <- fread(DATA_FILE, data.table = FALSE) %>%
  distinct(city, year, Total_Population, Black_Pop, White_Pop) %>%
  group_by(city) %>% arrange(year) %>%
  fill(Total_Population, Black_Pop, White_Pop, .direction = "downup") %>%
  ungroup()
d3 <- an %>%
  mutate(year = year(month)) %>%
  left_join(pop, by = c("city", "year")) %>%
  group_by(city, period = ifelse(post == 1, "post", "pre")) %>%
  summarise(county_pop   = mean(Total_Population),
            pct_pop_black = 100 * mean(Black_Pop / Total_Population),
            pct_pop_white = 100 * mean(White_Pop / Total_Population),
            rate_black   = mean(1e5 * black / Black_Pop),
            rate_white   = mean(1e5 * white / White_Pop),
            .groups = "drop") %>%
  mutate(bw_rate_ratio = rate_black / rate_white) %>%
  arrange(city, desc(period)) %>%
  mutate(county_pop = round(county_pop),
         across(c(pct_pop_black, pct_pop_white), ~ round(.x, 1)),
         across(c(rate_black, rate_white, bw_rate_ratio), ~ round(.x, 2)))
cat("\nTable D3. Monthly force per 100,000 county residents, by race\n")
print(as.data.frame(d3))

## D4. Force-type mix (% of reports in the window)
## Chicago is excluded: its export misaligns the weapon-flag columns.
## Cincinnati is excluded: its INCIDENT_DESCRIPTION codes are not
## captured by the harmonized force categories (most fall in "Other").
d4 <- raw %>%
  semi_join(an, by = c("city", "month")) %>%
  filter(!city %in% c("Chicago", "Cincinnati")) %>%
  count(city, force_type, wt = n) %>%
  group_by(city) %>% mutate(pct = round(100 * n / sum(n), 1)) %>% ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = force_type, values_from = pct, values_fill = 0)
cat("\nTable D4. Force-type mix, % of reports (Chicago and Cincinnati excluded)\n")
print(as.data.frame(d4))

write.csv(d1, "desc_D1_design.csv", row.names = FALSE)
write.csv(d2, "desc_D2_force_by_race.csv", row.names = FALSE)
write.csv(d3, "desc_D3_rates.csv", row.names = FALSE)
write.csv(d4, "desc_D4_force_type.csv", row.names = FALSE)

## Q0. Before any RD: in which cities is force against Black subjects
##     most disproportionate (D2 shares vs. D3 population shares)? Did
##     the Black share of force change from pre to post in any city?


## =====================================================================
## PART 1. The naive pre/post comparison (what the original study did)
## =====================================================================
## Q1. Compare mean monthly force before vs. after the policy, using the
##     full window. What does this comparison assume about trends?

naive <- an %>%
  group_by(city) %>%
  summarise(pre_mean  = mean(total[post == 0]),
            post_mean = mean(total[post == 1]),
            pct_change = 100 * (post_mean / pre_mean - 1),
            .groups = "drop")
print(naive)
## It assumes the only reason force differs before vs. after is the
## policy -- no trend, no seasonality, no other change in 3 years.


## =====================================================================
## PART 2. Look at the data
## =====================================================================
## Q2. Plot monthly force for each city with a vertical line at the
##     cutoff. Where do you see a jump? Where do you see a trend that
##     would fool a naive pre/post comparison?

p_raw <- ggplot(an, aes(x = month, y = total)) +
  geom_line(colour = "grey40") +
  geom_vline(data = policy_dates, aes(xintercept = policy_date),
             colour = "firebrick", linetype = "dashed") +
  facet_wrap(~ city, scales = "free_y", ncol = 2) +
  labs(x = NULL, y = "Use-of-force reports per month",
       title = "Monthly use of force, +/- 3 years around policy adoption") +
  theme_minimal()
print(p_raw)
ggsave("fig1_raw_series.png", p_raw, width = 9, height = 8, dpi = 150)

## Pooled RD plot. Cities differ hugely in size, so first remove each
## city's average level (a city fixed effect), then plot the monthly
## mean residual against months-since-adoption with a line on each side.
an <- an %>%
  group_by(city) %>%
  mutate(r_total = y_total - mean(y_total),
         r_black = y_black - mean(y_black),
         r_white = y_white - mean(y_white)) %>%
  ungroup()

bins <- an %>% group_by(x) %>% summarise(r = mean(r_total), .groups = "drop")
p_rd <- ggplot(bins, aes(x, r)) +
  geom_point(colour = "grey30") +
  geom_smooth(data = filter(bins, x < 0),  method = "lm", se = TRUE, colour = "steelblue") +
  geom_smooth(data = filter(bins, x >= 0), method = "lm", se = TRUE, colour = "steelblue") +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "firebrick") +
  labs(x = "Months since policy took effect",
       y = "log(force + 1), city-demeaned",
       title = "Pooled RD plot, five cities") +
  theme_minimal()
print(p_rd)
ggsave("fig2_pooled_rdplot.png", p_rd, width = 7, height = 4.5, dpi = 150)


## =====================================================================
## PART 3. Local linear RD -- by hand, then with rdrobust
## =====================================================================
## An RD estimate is just a weighted regression near the cutoff:
##   y = a + tau*post + b1*x + b2*post*x + city FE,
## using only months with |x| <= h, weighting months closer to the
## cutoff more (triangular kernel). tau is the jump at the cutoff.
##
## Standard errors: monthly counts are serially correlated, and all
## cities share calendar shocks. We use Driscoll-Kraay (panel
## Newey-West) errors via sandwich::vcovPL. With only five cities,
## treat any p-value with humility.

rd_ll <- function(d, y, h = H, cutoff = 0, donut = FALSE, pooled = TRUE,
                  extra = NULL) {
  d <- d %>% mutate(xc = x - cutoff, D = as.integer(xc >= 0)) %>%
    filter(abs(xc) <= h)
  if (donut) d <- d %>% filter(!xc %in% c(-1, 0))
  d$w <- 1 - abs(d$xc) / (h + 1)
  rhs <- "D * xc"
  if (pooled && length(unique(d$city)) > 1) rhs <- paste(rhs, "+ factor(city)")
  if (!is.null(extra)) rhs <- paste(rhs, "+", extra)
  fit <- lm(as.formula(paste(y, "~", rhs)), data = d, weights = w)
  V <- if (length(unique(d$city)) > 1)
         vcovPL(fit, cluster = ~ city, order.by = ~ month, lag = 3) else
         NeweyWest(fit, lag = 3, prewhite = FALSE)
  ct <- coeftest(fit, vcov. = V)["D", ]
  tibble(est = ct[1], se = ct[2], p = ct[4], n_months = nrow(d), h = h)
}

## Q3. Estimate the pooled jump in total force at H = 12 months.
##     Interpret the coefficient as a percent change.

pooled_total <- rd_ll(an, "y_total")
print(pooled_total)
cat(sprintf("Pooled jump: %.1f%% (SE %.1f)\n",
            100 * (exp(pooled_total$est) - 1), 100 * pooled_total$se))

## The same thing with rdrobust (data-driven bandwidth, bias-corrected
## robust confidence intervals). rdrobust has no city fixed effects
## option, so we feed it the city-demeaned outcome.
if (has_rdrobust) {
  rd1 <- rdrobust(y = an$r_total, x = an$x, c = 0, kernel = "triangular",
                  p = 1, bwselect = "mserd")
  summary(rd1)
  rdplot(y = an$r_total, x = an$x, c = 0, p = 1,
         x.label = "Months since policy", y.label = "log force (demeaned)",
         title = "rdplot: pooled")
}


## =====================================================================
## PART 4. City by city
## =====================================================================
## Q4. Estimate the jump separately for each city. Which results are
##     you most willing to believe, and why? (Hint: date_status.)

by_city <- an %>% group_split(city) %>%
  lapply(function(d) rd_ll(d, "y_total", pooled = FALSE) %>%
           mutate(city = d$city[1], date_status = d$date_status[1])) %>%
  bind_rows() %>% select(city, date_status, everything())
print(by_city)

p_city <- ggplot(by_city, aes(x = est, y = reorder(city, est))) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_errorbarh(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se), height = 0.2) +
  geom_point(aes(shape = date_status), size = 3) +
  labs(x = "Jump in log(force) at cutoff (95% CI)", y = NULL,
       title = "City-specific RD estimates, h = 12 months") +
  theme_minimal()
print(p_city)
ggsave("fig3_city_estimates.png", p_city, width = 7, height = 4, dpi = 150)


## =====================================================================
## PART 5. Black vs. White subjects
## =====================================================================
## Three ways to ask whether the policy changed force differently for
## Black and White subjects:
##   (a) separate RDs for force against Black and against White subjects
##   (b) a "difference in discontinuities": stack the two series and
##       test whether the Black jump differs from the White jump
##   (c) RD on the Black share of force (Black / (Black + White))
##
## Note on denominators: within one city over a few years, the Black
## and White populations are close to constant, so jumps in log counts
## track jumps in population rates. Levels (the Black-White gap in
## rates) need census denominators; jumps at the cutoff do not.

## Q5a. Separate pooled RDs for Black and White subjects.
race_sep <- bind_rows(
  rd_ll(an, "y_black") %>% mutate(group = "Black subjects"),
  rd_ll(an, "y_white") %>% mutate(group = "White subjects"))
print(race_sep)

## Q5b. Difference in discontinuities. Stack the two series; the
##      coefficient on D:black is (Black jump) - (White jump).
stacked <- an %>%
  select(city, month, x, y_black, y_white) %>%
  pivot_longer(c(y_black, y_white), names_to = "grp", values_to = "y") %>%
  mutate(black = as.integer(grp == "y_black"))

did_rd <- function(d, h = H) {
  d <- d %>% mutate(D = as.integer(x >= 0)) %>% filter(abs(x) <= h)
  d$w <- 1 - abs(d$x) / (h + 1)
  fit <- lm(y ~ D * x * black + factor(city) * black, data = d, weights = d$w)
  V <- vcovPL(fit, cluster = ~ city, order.by = ~ month, lag = 3)
  coeftest(fit, vcov. = V)[c("D", "D:black"), ]
}
print(did_rd(stacked))
## Row "D"       = jump for White subjects
## Row "D:black" = how much larger (more negative) the Black jump is

## Q5c. RD on the Black share of force. Did the racial composition of
##      force change at the cutoff?
share_rd <- rd_ll(an, "black_share")
print(share_rd)

## Q5d. City-by-city Black vs. White plot
race_city <- an %>% group_split(city) %>%
  lapply(function(d) bind_rows(
    rd_ll(d, "y_black", pooled = FALSE) %>% mutate(group = "Black"),
    rd_ll(d, "y_white", pooled = FALSE) %>% mutate(group = "White")) %>%
      mutate(city = d$city[1])) %>% bind_rows()
print(race_city)

p_race <- ggplot(race_city, aes(x = est, y = city, colour = group)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_pointrange(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se),
                  position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c(Black = "#1b6ca8", White = "#d1495b")) +
  labs(x = "Jump in log(force) at cutoff (95% CI)", y = NULL, colour = "Subject race",
       title = "RD estimates by subject race, h = 12 months") +
  theme_minimal()
print(p_race)
ggsave("fig4_race_estimates.png", p_race, width = 7, height = 4.5, dpi = 150)


## =====================================================================
## PART 6. Is it real? Robustness for RD in time
## =====================================================================

## Q6a. Bandwidth sensitivity: re-estimate the pooled jump for
##      h = 6, 9, 12, 18, 24, 36. Why does a wide bandwidth worry us
##      more in RDiT than in a cross-sectional RD?
bw <- lapply(c(6, 9, 12, 18, 24, 36), function(h)
  bind_rows(rd_ll(an, "y_total", h = h) %>% mutate(outcome = "Total"),
            rd_ll(an, "y_black", h = h) %>% mutate(outcome = "Black"),
            rd_ll(an, "y_white", h = h) %>% mutate(outcome = "White"))) %>%
  bind_rows()
print(bw)
p_bw <- ggplot(bw, aes(h, est, colour = outcome)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_pointrange(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se),
                  position = position_dodge(width = 1.5)) +
  labs(x = "Bandwidth (months)", y = "Jump in log(force)",
       title = "Bandwidth sensitivity") + theme_minimal()
print(p_bw)
ggsave("fig5_bandwidth.png", p_bw, width = 7, height = 4, dpi = 150)

## Q6b. Placebo cutoffs: pretend the policy started 12 or 18 months
##      earlier (using only pre-period data) or later (post-period data
##      only). What should you find if the design is sound?
placebo <- bind_rows(
  rd_ll(filter(an, x < 0),  "y_total", cutoff = -18) %>% mutate(fake = -18),
  rd_ll(filter(an, x < 0),  "y_total", cutoff = -12) %>% mutate(fake = -12),
  rd_ll(filter(an, x >= 0), "y_total", cutoff =  12) %>% mutate(fake =  12),
  rd_ll(filter(an, x >= 0), "y_total", cutoff =  18) %>% mutate(fake =  18))
print(placebo)

## Q6c. Donut RD: drop the adoption month and the month before it
##      (training, rollout, mid-month effective dates, anticipation).
print(bind_rows(rd_ll(an, "y_total") %>% mutate(spec = "baseline"),
                rd_ll(an, "y_total", donut = TRUE) %>% mutate(spec = "donut")))

## Q6d. Seasonality. Use of force has a calendar cycle: in Chicago,
##      monthly reports average ~515 in July but ~355 in December and
##      January. A cutoff in October therefore sits right before the
##      seasonal decline. Add city-by-calendar-month fixed effects and
##      re-run the pooled, race-specific, and city estimates. With only
##      12 months on each side, every calendar month appears about once
##      per side, so use h = 24 to separate season from policy.
##      Which conclusions survive?
seas_rd <- function(d, y, h = 24, pooled = TRUE)
  rd_ll(d, y, h = h, pooled = pooled,
        extra = if (pooled) "factor(city):factor(cal_month)" else "factor(cal_month)")
cm <- cm %>% mutate(cal_month = month(month))
cm_s <- cm %>% filter(INCLUDE_2020 | month < as.Date("2020-03-01"))

seas_pooled <- bind_rows(lapply(c(y_total = "Total", y_black = "Black", y_white = "White") %>%
  names(), function(y) bind_rows(
    rd_ll(cm_s, y, h = 24) %>% mutate(spec = "no seasonal control"),
    seas_rd(cm_s, y)       %>% mutate(spec = "city x calendar-month FE")) %>%
    mutate(outcome = y)))
print(seas_pooled)

seas_city <- lapply(unique(cm_s$city), function(cty)
  lapply(c("y_total", "y_black", "y_white"), function(y)
    seas_rd(filter(cm_s, city == cty), y, pooled = FALSE) %>%
      mutate(city = cty, outcome = y)) %>% bind_rows()) %>% bind_rows()
print(seas_city)

## Q6e. Why no McCrary density test? In a cross-sectional RD we test
##      whether units bunch just above the cutoff. Here the running
##      variable is calendar time: every month appears exactly once per
##      city. There is nothing to manipulate, so the test is meaningless.
##      The threat in RDiT is different: OTHER THINGS that changed at the
##      same moment. Name two for these cities.


## =====================================================================
## PART 7. Uncertain cutoffs and recording changes
## =====================================================================
## Q7a. For Austin and Cincinnati we know the year but not the month.
##      Shift the cutoff from -4 to +4 months and re-estimate. What
##      happens to the estimate when the cutoff is wrong?
## (uses the full city series, not the +/-36 window, so shifted
##  cutoffs still have 12 months on each side)
shift_cities <- intersect(c("Chicago", "Austin", "Cincinnati", "Indianapolis"), unique(cm$city))
shift <- lapply(shift_cities, function(cty)
  lapply(-4:4, function(s)
    rd_ll(filter(cm, city == cty, month < as.Date("2020-03-01")), "y_total",
          cutoff = s, pooled = FALSE) %>% mutate(city = cty, s = s)) %>%
    bind_rows()) %>% bind_rows()
print(shift)
p_shift <- ggplot(shift, aes(s, est)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_pointrange(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se)) +
  facet_wrap(~ city, nrow = 1) +
  labs(x = "Cutoff shifted by (months)", y = "Jump in log(force)",
       title = "Sensitivity to the assumed adoption month") + theme_minimal()
print(p_shift)
ggsave("fig6_cutoff_shift.png", p_shift, width = 9, height = 3.5, dpi = 150)
## Chicago (verified date) is the benchmark: a real effect should peak
## near s = 0 and fade as the cutoff moves away from the true date.

## Q7b. Austin changed its data format in 2018 (from "Weapon Used"
##      fields to "Used ..." flags). Plot Austin's physical force by
##      month. What would a 36-month bandwidth attribute to the policy?
aus_phys <- raw %>%
  filter(city == "Austin", force_type == "Physical") %>%
  count(month, wt = n)
p_aus <- ggplot(aus_phys, aes(month, n)) +
  geom_line() +
  geom_vline(xintercept = as.Date("2016-10-01"), colour = "firebrick", linetype = "dashed") +
  geom_vline(xintercept = as.Date("2018-01-01"), colour = "darkorange", linetype = "dotted") +
  annotate("text", x = as.Date("2016-10-01"), y = Inf, label = "policy", vjust = 1.5, hjust = 1.1) +
  annotate("text", x = as.Date("2018-01-01"), y = Inf, label = "new data format", vjust = 1.5, hjust = -0.1) +
  labs(x = NULL, y = "Physical force reports per month",
       title = "Austin: a recording change inside the post window") + theme_minimal()
print(p_aus)
ggsave("fig7_austin_recording.png", p_aus, width = 8, height = 3.5, dpi = 150)


## =====================================================================
## PART 8. Serious vs. minor force
## =====================================================================
## Theory: use-of-force policies (de-escalation, force continuums,
## warnings, chokehold bans) target the most dangerous force. They
## should reduce SERIOUS force, while MINOR hands-on force (grabs,
## takedowns, handcuffing resistance) may not change -- or may rise if
## officers substitute down the force continuum.
##
## Definition used here (harmonized categories):
##   serious = Firearm + Less-lethal weapon (Taser, baton, impact) + Chemical
##   minor   = Physical (weaponless)
##   "Other" is excluded from both.
## Only cities whose force types are reliably coded are used: Chicago's
## export misaligns its weapon flags, and Cincinnati's description codes
## are not mapped to these categories yet.
SEVERITY_CITIES <- c("Austin", "Indianapolis", "Minneapolis")

sev <- raw %>%
  filter(city %in% SEVERITY_CITIES) %>%
  mutate(severity = case_when(
    force_type %in% c("Firearm", "Less-lethal weapon", "Chemical") ~ "serious",
    force_type == "Physical" ~ "minor",
    TRUE ~ NA_character_)) %>%
  filter(!is.na(severity)) %>%
  group_by(city, month, severity) %>%
  summarise(all = sum(n), black = sum(n[race == "Black"]),
            white = sum(n[race == "White"]), .groups = "drop") %>%
  pivot_wider(names_from = severity, values_from = c(all, black, white),
              values_fill = 0) %>%
  semi_join(cm, by = c("city", "month")) %>%        # same phase-in trimming
  left_join(policy_dates %>% select(city, policy_date), by = "city") %>%
  mutate(x = (year(month) - year(policy_date)) * 12 +
             (month(month) - month(policy_date)),
         post = as.integer(x >= 0),
         cal_month = month(month),
         across(c(all_serious, all_minor, black_serious, black_minor,
                  white_serious, white_minor), ~ log(.x + 1), .names = "y_{.col}")) %>%
  filter(INCLUDE_2020 | month < as.Date("2020-03-01"))

## Table D5. Monthly serious and minor force, pre vs. post (+/- 36 months)
d5 <- sev %>%
  filter(abs(x) <= WINDOW) %>%
  group_by(city, period = ifelse(post == 1, "post", "pre")) %>%
  summarise(serious = mean(all_serious), minor = mean(all_minor),
            pct_serious = 100 * sum(all_serious) / sum(all_serious + all_minor),
            serious_black = mean(black_serious), serious_white = mean(white_serious),
            minor_black = mean(black_minor), minor_white = mean(white_minor),
            .groups = "drop") %>%
  arrange(city, desc(period)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))
cat("\nTable D5. Mean monthly serious vs. minor force\n"); print(as.data.frame(d5))
write.csv(d5, "desc_D5_severity.csv", row.names = FALSE)

## Q8a. Estimate the jump for serious and for minor force, overall and
##      by subject race, with and without seasonal controls.
outcomes <- c(Serious_all   = "y_all_serious",   Minor_all   = "y_all_minor",
              Serious_Black = "y_black_serious", Serious_White = "y_white_serious",
              Minor_Black   = "y_black_minor",   Minor_White   = "y_white_minor")
sev_res <- lapply(names(outcomes), function(o) bind_rows(
  rd_ll(sev, outcomes[[o]], h = 12) %>% mutate(spec = "h=12, no seasonal"),
  rd_ll(sev, outcomes[[o]], h = 24, extra = "factor(city):factor(cal_month)") %>%
    mutate(spec = "h=24, seasonal FE")) %>%
  mutate(outcome = o)) %>% bind_rows() %>%
  select(outcome, spec, est, se, p, n_months)
print(sev_res, n = 30)

p_sev <- sev_res %>%
  tidyr::separate(outcome, c("severity", "group"), sep = "_") %>%
  ggplot(aes(x = est, y = group, colour = severity)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_pointrange(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se),
                  position = position_dodge(width = 0.5)) +
  facet_wrap(~ spec) +
  scale_colour_manual(values = c(Serious = "#b2182b", Minor = "#7f7f7f")) +
  labs(x = "Jump in log(force) at cutoff (95% CI)", y = NULL, colour = NULL,
       title = "Serious vs. minor force (Austin, Indianapolis, Minneapolis)") +
  theme_minimal()
print(p_sev)
ggsave("fig8_severity.png", p_sev, width = 9, height = 4, dpi = 150)

## Q8b. Did serious force fall RELATIVE to minor force? Stack the two
##      series; D:serious is (serious jump) - (minor jump). This also
##      nets out anything that shifts all reporting at the cutoff.
sev_stack <- sev %>%
  select(city, month, x, cal_month, y_all_serious, y_all_minor) %>%
  pivot_longer(c(y_all_serious, y_all_minor), names_to = "grp", values_to = "y") %>%
  mutate(serious = as.integer(grp == "y_all_serious"))
sev_did <- function(d, h, seasonal) {
  d <- d %>% mutate(D = as.integer(x >= 0)) %>% filter(abs(x) <= h)
  d$w <- 1 - abs(d$x) / (h + 1)
  f <- if (seasonal) y ~ D * x * serious + factor(city):factor(cal_month) * serious
       else          y ~ D * x * serious + factor(city) * serious
  fit <- lm(f, data = d, weights = d$w)
  coeftest(fit, vcov. = vcovPL(fit, cluster = ~ city, order.by = ~ month, lag = 3))[c("D", "D:serious"), ]
}
cat("\nh = 12, no seasonal control:\n"); print(sev_did(sev_stack, 12, FALSE))
cat("\nh = 24, seasonal FE:\n");          print(sev_did(sev_stack, 24, TRUE))
## Row "D"         = jump in minor force
## Row "D:serious" = extra jump in serious force beyond minor force

## Q8c. Robustness: chemical agents (OC spray) are low on many force
##      continuums. Re-define serious = Firearm + Less-lethal weapon only.
sev_alt <- raw %>%
  filter(city %in% SEVERITY_CITIES,
         force_type %in% c("Firearm", "Less-lethal weapon")) %>%
  group_by(city, month) %>% summarise(n = sum(n), .groups = "drop") %>%
  right_join(sev %>% select(city, month, x, cal_month), by = c("city", "month")) %>%
  mutate(n = coalesce(n, 0L), y_alt = log(n + 1))
print(bind_rows(
  rd_ll(sev_alt, "y_alt", h = 12) %>% mutate(spec = "h=12, no seasonal"),
  rd_ll(sev_alt, "y_alt", h = 24, extra = "factor(city):factor(cal_month)") %>%
    mutate(spec = "h=24, seasonal FE")))


## =====================================================================
## Discussion (write 1-2 paragraphs)
## =====================================================================
## 1. RDiT identifies the effect at the moment of adoption. Is that the
##    effect a city council cares about? What about effects that build
##    over months as officers are retrained?
## 2. A policy that requires officers to REPORT more types of force can
##    raise recorded force while lowering actual force. Which of the
##    eight policy elements could do this?
## 3. Part 8: the theory says policies cut serious force, not minor
##    force. What pattern in Q8a-Q8b would support it? What would it mean
##    if minor force ROSE (substitution down the continuum)? And if
##    serious force fell more for White than Black subjects?
## 3b. Compare this design with the DiD from Weeks 8-9. What would a
##    comparison city buy you here?
## 4. The Black-White difference in discontinuities: what would have to
##    be true for a policy to reduce force more for Black subjects?
##    What else could produce that pattern at the cutoff? (Hint: if
##    Black residents' encounters with police are more concentrated in
##    summer, seasonality alone can manufacture a racial gap.)
## 5. Chicago's drop survives seasonal controls. What else happened in
##    Chicago in late 2017 (the DOJ findings and the move toward a
##    consent decree, new reporting forms, training)? Can RDiT tell the
##    policy apart from those?
