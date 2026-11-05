# =====================================================================
# Causal Analysis (QMSS, CUNY GC) -- Week 10 Lab 9
# Regression Discontinuity
#
# Two classic sharp RD designs, both with public data:
#   Part 1. The minimum legal drinking age (Carpenter & Dobkin 2009;
#           Angrist & Pischke, Mastering 'Metrics, ch. 4). Running
#           variable: age. Cutoff: 21. Outcome: deaths per 100,000.
#   Part 2. Close U.S. Senate elections, 1914-2010 (Cattaneo, Frandsen &
#           Titiunik 2015, building on Lee 2008). Running variable: the
#           Democratic margin of victory. Cutoff: 0. Outcome: the
#           Democratic vote share in the next election for that seat.
#   Part 3. Is the design believable? Density, covariates, placebos.
#
# Packages: install.packages(c("rdrobust", "rddensity", "ggplot2"))
#
# ---------------------------------------------------------------------
# WHAT TO SUBMIT ON BRIGHTSPACE: THREE OUTPUTS
#   Upload one PDF or Word file with three sections, labeled Output 1-3.
#   Each section has the file the code below saves, plus a short answer
#   (about 150 words) to the question listed with it.
#
#   OUTPUT 1  Figure  output1_mlda.png          Part 1, answer Q1
#                     (report the jumps for all deaths and motor vehicle
#                     deaths, and the ages 20-22 estimate)
#   OUTPUT 2  Figure  output2_senate.png        Part 2, answer Q2
#                     (report the main estimate, its robust 95% CI, and
#                     the bandwidth)
#   OUTPUT 3  Table   output3_validity.txt      Part 3, answer Q3
#
#   Q4 and Q5 are for discussion in class.
#   Look for the lines marked  >>> OUTPUT  below.
# ---------------------------------------------------------------------
# =====================================================================

library(rdrobust)
library(rddensity)
library(ggplot2)

base <- "https://edwingrimsley.github.io/Causal-Analysis/files/"
mlda   <- read.csv(paste0(base, "lab9-mlda.csv"))
senate <- read.csv(paste0(base, "lab9-senate.csv"))
# Offline? Download both CSVs from the course site and use read.csv(file.choose()).


# ---------------------------------------------------------------------
# PART 1. Does turning 21 raise the death rate?
#   One row per month of age, ages 19-23. Rates are deaths per 100,000
#   person-years. Columns: all, mva (motor vehicle), internal (disease),
#   suicide, homicide, alcohol, drugs.
# ---------------------------------------------------------------------
mlda <- subset(mlda, !is.na(all))
mlda$over21 <- as.integer(mlda$agecell >= 21)   # treatment: legal to drink
mlda$age_c  <- mlda$agecell - 21                # running variable, centered

# Sharp RD by regression: a jump at the cutoff, separate slopes on each side
rd_all <- lm(all ~ over21 * age_c, data = mlda)
rd_mva <- lm(mva ~ over21 * age_c, data = mlda)
coef(rd_all)["over21"]          # about 7.7 more deaths per 100,000
coef(rd_mva)["over21"]          # about 4.5
mean(mlda$all[mlda$over21 == 0])   # baseline just under 21: about 93

# Same model, closer to the cutoff (ages 20-22 only)
near <- subset(mlda, agecell >= 20 & agecell < 22)
coef(lm(all ~ over21 * age_c, data = near))["over21"]   # about 9.8

# A placebo outcome: deaths from disease should not jump at 21
coef(lm(internal ~ over21 * age_c, data = mlda))["over21"]  # about 0.4

# >>> OUTPUT 1: all deaths and motor vehicle deaths by age, with the
#     fitted lines on each side of 21.
long <- rbind(data.frame(cause = "All causes",      age = mlda$agecell, rate = mlda$all, over21 = mlda$over21),
              data.frame(cause = "Motor vehicle",   age = mlda$agecell, rate = mlda$mva, over21 = mlda$over21))
p1 <- ggplot(long, aes(age, rate)) +
  geom_point(colour = "grey35") +
  geom_smooth(aes(group = over21), method = "lm", formula = y ~ x, se = FALSE, colour = "#C8553D") +
  geom_vline(xintercept = 21, linetype = "dashed") +
  facet_wrap(~cause, scales = "free_y") +
  labs(x = "Age (months, in years)", y = "Deaths per 100,000",
       title = "Mortality and the minimum legal drinking age") +
  theme_minimal()
ggsave("output1_mlda.png", p1, width = 9, height = 4.2, dpi = 150)
p1

# Q1 (Output 1). How large is the jump in all-cause deaths, relative to
#     the rate just under 21, and which cause drives it? The design
#     assumes nothing else that affects mortality changes sharply at 21.
#     Name one thing that does change at 21 and say whether it threatens
#     the estimate. What does the "internal causes" result add?


# ---------------------------------------------------------------------
# PART 2. Close Senate elections: is there an incumbency advantage?
#   margin = Democratic margin at election t (running variable)
#   vote   = Democratic vote share at the next election for the seat (t+2)
#   Democrats who barely won vs. barely lost: as good as a coin flip?
# ---------------------------------------------------------------------
summary(senate$margin)
mean(senate$margin > 0)                         # Democrats won about 54%

# The RD plot: binned means, with local linear fits within 30 points
rdp <- rdplot(senate$vote, senate$margin, c = 0, p = 1,
              subset = abs(senate$margin) <= 30,
              x.label = "Democratic margin at election t",
              y.label = "Democratic vote share at t+2",
              title = "Close Senate elections")

# The main estimate: local linear, triangular kernel, MSE-optimal bandwidth
est <- rdrobust(senate$vote, senate$margin, c = 0)
summary(est)
# Expect: conventional estimate about 7.4 (SE 1.5); robust 95% CI about
# [4.1, 10.9]; bandwidth h about 17.8 points; about 680 elections used.

# How much does the answer depend on the bandwidth and the polynomial?
grid <- expand.grid(h = c(5, 10, 15, 20, 25, 30), p = 1:2)
grid[c("est", "lo", "hi")] <- NA
for (i in seq_len(nrow(grid))) {
  r <- rdrobust(senate$vote, senate$margin, c = 0, h = grid$h[i], p = grid$p[i])
  grid$est[i] <- r$coef[1]
  grid$lo[i]  <- r$ci[3, 1]              # robust CI
  grid$hi[i]  <- r$ci[3, 2]
}
grid
# Expect, local linear: about 12.3, 8.0, 7.5, 7.3, 7.1, 7.3 as h widens.

# >>> OUTPUT 2: the RD plot (left) and the bandwidth check (right).
grid$poly <- ifelse(grid$p == 1, "Local linear", "Local quadratic")
p2b <- ggplot(grid, aes(h, est, colour = poly)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 2)) +
  scale_colour_manual(values = c("#2E6F73", "#C8553D")) +
  labs(x = "Bandwidth (points of margin)", y = "Estimate and robust 95% CI",
       colour = NULL, title = "Sensitivity to bandwidth") +
  theme_minimal() + theme(legend.position = "bottom")
png("output2_senate.png", width = 1700, height = 750, res = 150)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
print(rdp$rdplot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2b,        vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

# Q2 (Output 2). Interpret the main estimate in one sentence. For which
#     seats and candidates does it apply, and for which does it say
#     nothing? Why do estimates get larger and noisier as the bandwidth
#     shrinks, and which bandwidth would you report?


# ---------------------------------------------------------------------
# PART 3. Is the design believable?
#   (a) No sorting: is the number of elections smooth across 0?
#   (b) Balance: are characteristics fixed BEFORE election t smooth at 0?
#   (c) Placebo cutoffs: no jumps where there is no treatment.
# ---------------------------------------------------------------------
# (a) Density (manipulation) test
dens <- rddensity(X = senate$margin, c = 0)
summary(dens)
dens_p <- dens$test$p_jk                        # expect about 0.38
rdplotdensity(dens, senate$margin)

# (b) Covariate balance: run the same RD on pre-election characteristics
covs <- c(demvoteshlag1 = "Dem. vote share, election t-1",
          demvoteshlag2 = "Dem. vote share, election t-2",
          presdemvoteshlag1 = "Dem. presidential vote share",
          demwinprv1 = "Dem. won election t-1",
          dopen = "Open seat",
          termshouse = "House terms served",
          termssenate = "Senate terms served",
          population = "State population")
bal <- do.call(rbind, lapply(names(covs), function(v) {
  r <- rdrobust(senate[[v]], senate$margin, c = 0)
  data.frame(check = covs[[v]], estimate = signif(r$coef[1], 3),
             robust_p = round(r$pv[3], 3), bandwidth = round(r$bws[1, 1], 1))
}))
bal
# Expect robust p-values above 0.2 for all but "Open seat" (about 0.02).

# (c) Placebo cutoffs, using only one side of 0 and a fixed bandwidth of 10
plac <- do.call(rbind, lapply(c(-20, -10, 10, 20), function(cut) {
  side <- if (cut < 0) subset(senate, margin < 0) else subset(senate, margin >= 0)
  r <- rdrobust(side$vote, side$margin, c = cut, h = 10)
  data.frame(check = paste("Placebo cutoff at", cut), estimate = round(r$coef[1], 2),
             robust_p = round(r$pv[3], 3), bandwidth = 10)
}))
plac
# Expect no robust p-value below 0.1.

# >>> OUTPUT 3: one table with all three checks.
out3 <- rbind(data.frame(check = "Density test (no sorting at 0)", estimate = NA,
                         robust_p = round(dens_p, 3), bandwidth = NA),
              bal, plac)
writeLines(capture.output(print(out3, row.names = FALSE)), "output3_validity.txt")
out3

# Q3 (Output 3). What would each check catch if the design were broken?
#     One covariate of eight is "significant" at 0.05. Is that a problem?
#     What would you look at next?


# ---------------------------------------------------------------------
# For discussion
# ---------------------------------------------------------------------
# Q4. Try the placebo cutoff at -10 WITHOUT fixing the bandwidth:
side <- subset(senate, margin < 0)
summary(rdrobust(side$vote, side$margin, c = -10))
#     The automatic bandwidth collapses to under 2 points and the
#     "effect" is about -23. What went wrong, and what does it teach you
#     about reading any single RD estimate?
# Q5. Last week's DiD rested on parallel trends. What plays that role in
#     RD? Why can't anyone choose to land just above or just below the
#     cutoff here, and in what RD settings could they?
