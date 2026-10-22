###############################################################################
# Lab 7: Difference-in-Differences and Parallel Trends
# Causal Analysis — CUNY Graduate Center, QMSS · Week 8 · October 22
#
# Two studies:
#   Card & Krueger (1994): New Jersey raised its minimum wage on April 1, 1992;
#     Pennsylvania didn't. Fast-food employment, Feb 1992 vs. Nov 1992.
#   Kessler & Roth (2014): California switched to "active choice" organ-donor
#     sign-up in July 2011; 26 other states didn't. (The Effect, ch. 18.)
#
# WHAT TO SUBMIT (paste a screenshot or copy of each into one document):
#   Output 1  Your 2 x 2 table and difference-in-differences (Part 1)
#   Output 2  Your regression table (Part 2)
#   Output 3  Your event-study plot for organ donation (Part 3)
# Under each output, answer its one question in 2-3 sentences.
#
# Data files: card-krueger-1994.csv and organ-donations.csv
###############################################################################

# install.packages(c("fixest", "dplyr", "ggplot2"))
library(fixest)
library(dplyr)
library(ggplot2)


###############################################################################
# PART 1. CARD & KRUEGER: THE 2 x 2 BY HAND                       -> Output 1
###############################################################################

# A file window opens: pick card-krueger-1994.csv
ck <- read.csv(file.choose())
ck <- ck |> mutate(nj = as.numeric(state == "NJ"),        # 1 = New Jersey (treated)
                   post = as.numeric(wave == 2))         # 1 = November (after)

# The four averages: employment per store, by state and wave
cells <- ck |> group_by(state, period) |>
  summarize(mean_fte = round(mean(fte, na.rm = TRUE), 2), .groups = "drop")
as.data.frame(cells)

# Difference-in-differences by hand
nj_change <- mean(ck$fte[ck$nj == 1 & ck$post == 1], na.rm = TRUE) -
             mean(ck$fte[ck$nj == 1 & ck$post == 0], na.rm = TRUE)
pa_change <- mean(ck$fte[ck$nj == 0 & ck$post == 1], na.rm = TRUE) -
             mean(ck$fte[ck$nj == 0 & ck$post == 0], na.rm = TRUE)
round(c(NJ_change = nj_change, PA_change = pa_change, DiD = nj_change - pa_change), 2)

# Picture: the dashed line is New Jersey IF it had followed Pennsylvania's trend
nj_before <- mean(ck$fte[ck$nj == 1 & ck$post == 0], na.rm = TRUE)
pic <- rbind(
  data.frame(group = "New Jersey", period = c("Feb", "Nov"),
             fte = c(nj_before, nj_before + nj_change)),
  data.frame(group = "Pennsylvania", period = c("Feb", "Nov"),
             fte = c(mean(ck$fte[ck$nj == 0 & ck$post == 0], na.rm = TRUE),
                     mean(ck$fte[ck$nj == 0 & ck$post == 1], na.rm = TRUE))),
  data.frame(group = "NJ if parallel to PA", period = c("Feb", "Nov"),
             fte = c(nj_before, nj_before + pa_change)))
ggplot(pic, aes(period, fte, colour = group, group = group, linetype = group == "NJ if parallel to PA")) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  guides(linetype = "none") + labs(x = NULL, y = "Employment per store", colour = NULL) + theme_bw()

# OUTPUT 1: the four averages and the DiD line above.
# Question 1: Card & Krueger report +2.76. Do you match it? In one sentence,
#             what does the parallel-trends assumption say about NJ and PA?


###############################################################################
# PART 2. THE SAME NUMBER FROM A REGRESSION                      -> Output 2
###############################################################################

# (a) The 2 x 2 as a regression: the coefficient on nj:post is the DiD
m1 <- feols(fte ~ nj * post, data = ck, cluster = ~ store)

# (b) Two-way fixed effects: a dummy for every store and for each wave
m2 <- feols(fte ~ nj:post | store + wave, data = ck, cluster = ~ store)

etable(m1, m2, headers = c("2 x 2 regression", "Two-way fixed effects"))

# OUTPUT 2: the etable() output.
# Question 2: Find the DiD in each column. Why do the two models give
#             (almost) the same answer? What do the store and wave fixed
#             effects each remove?


###############################################################################
# PART 3. ORGAN DONATION: CHECKING PARALLEL TRENDS                -> Output 3
###############################################################################

# A file window opens: pick organ-donations.csv
od <- read.csv(file.choose())
od <- od |> mutate(california = as.numeric(State == "California"),
                   after = as.numeric(Quarter_Num >= 4),     # Q3 2011 and later
                   quarter = factor(Quarter, levels = c("Q42010", "Q12011", "Q22011",
                                                        "Q32011", "Q42011", "Q12012")))

# (a) Look at the trends: California vs. the average of the other states
od |> group_by(group = ifelse(california == 1, "California", "Other states"), quarter) |>
  summarize(rate = mean(Rate), .groups = "drop") |>
  ggplot(aes(quarter, rate, colour = group, group = group)) +
  geom_line(linewidth = 1) + geom_point() + geom_vline(xintercept = 3.5, linetype = "dashed") +
  labs(x = NULL, y = "Organ-donor registration rate", colour = NULL) + theme_bw()

# (b) The difference-in-differences
feols(Rate ~ california:after | State + Quarter, data = od, cluster = ~ State)

# (c) Event study: one coefficient per quarter, relative to Q2 2011
es <- feols(Rate ~ i(quarter, california, ref = "Q22011") | State + Quarter,
            data = od, cluster = ~ State)
iplot(es, main = "California vs. other states, relative to Q2 2011")

# (d) Placebo: pretend the policy started in Q1 2011 (before-quarters only)
od_pre <- od |> filter(Quarter_Num <= 3) |> mutate(fake_after = as.numeric(Quarter_Num >= 2))
feols(Rate ~ california:fake_after | State + Quarter, data = od_pre, cluster = ~ State)

# OUTPUT 3: the event-study plot from (c).
# Question 3: Look at the two quarters BEFORE the policy. Do they support
#             parallel trends? What could Card & Krueger NOT check that this
#             study can?


###############################################################################
# OPTIONAL (not graded): more ways to probe parallel trends in Card & Krueger
###############################################################################

# Did the law actually raise wages? (Should be yes in NJ, no in PA)
ck |> group_by(state, period) |> summarize(mean_wage = round(mean(wage_st, na.rm = TRUE), 2), .groups = "drop")

# A second comparison group inside NJ: stores already paying $5.00+ in February
nj <- ck |> filter(nj == 1) |> group_by(store) |> mutate(wage_feb = wage_st[wave == 1]) |> ungroup() |>
  filter(wage_feb == 4.25 | wage_feb >= 5.00) |> mutate(low = as.numeric(wage_feb == 4.25))
feols(fte ~ low * post, data = nj, cluster = ~ store)
