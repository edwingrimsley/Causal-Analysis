# =====================================================================
# Causal Analysis (QMSS, CUNY GC) -- Week 9 Lab
# Staggered, switching treatments: Operation Impact
#
# Reading: Legewie & Fagan (2019), ASR.
# Data: the ZONES are real (MacDonald, Fagan & Geller 2016 replication
# data: which areas were impact zones, month by month, 2004-2012).
# The STUDENTS and SCORES are simulated, so we know the truth:
#   Black boys aged 13-15 lose about 0.13 SD per full year of exposure.
#   No one else is affected. (Column true_effect holds each student's.)
#
# The data load from the course site below. week9-lab8-simulate-data.R shows
# how they were simulated, and week9-lab8-build-zones.R how the real zone
# histories were built (both optional).
# Packages: install.packages(c("fixest", "did", "ggplot2"))
#
# Three outcomes, one lesson each:
#   score          the design works: parallel trends hold
#   score_trend    zoned areas were already declining; never-zoned were not
#   score_persist  effects linger after a zone is removed
#
# ---------------------------------------------------------------------
# WHAT TO SUBMIT ON BRIGHTSPACE: THREE OUTPUTS
#   Upload one PDF or Word file with three sections, labeled Output 1-3.
#   Each section has the table or figure the code below saves, plus a
#   short answer (about 150 words) to the questions listed with it.
#
#   OUTPUT 1  Table   output1_lf_models.txt        Part 2, answer Q2
#   OUTPUT 2  Figure  output2_event_studies.png    Part 3, answer Q3 + Q4
#                     (include the four pre-trend p-values from 3c)
#   OUTPUT 3  Figure  output3_cs_dynamic.png       Part 4, answer Q5
#                     (include the TWFE and C&S estimates for A and C)
#
#   The other questions (Q1, Q6, Q7) are for discussion in class.
#   Look for the lines marked  >>> OUTPUT  below.
# ---------------------------------------------------------------------
# =====================================================================

library(fixest)
library(did)
library(ggplot2)

base <- "https://edwingrimsley.github.io/Causal-Analysis/files/"
get_rds <- function(f) readRDS(gzcon(url(paste0(base, f))))
d <- get_rds("lab8-operation-impact-sim.rds")
# Offline? Download the .rds from the course site and use: d <- readRDS(file.choose())
d$group <- ifelse(d$race == "Black"    & d$female == 0, "Black boys",
           ifelse(d$race == "Black"    & d$female == 1, "Black girls",
           ifelse(d$race == "Hispanic" & d$female == 0, "Hispanic boys",
           ifelse(d$race == "Hispanic" & d$female == 1, "Hispanic girls", "White"))))
d$agegrp <- cut(d$age, c(8, 10, 11, 12, 15), labels = c("9-10", "11", "12", "13-15"))

# Legewie & Fagan's sample: students who were EVER exposed to a zone.
s <- subset(d, student_id %in% unique(d$student_id[d$exposure > 0]))
nrow(s); length(unique(s$student_id))   # about 211,000 student-years, 44,000 students
mean(s$exposure[s$exposure > 0])        # about 0.70 (paper: 0.77)
mean(s$exposure)                        # about 0.41 (paper: 0.26; these zones last longer)


# ---------------------------------------------------------------------
# PART 1. What does the treatment look like? (real zones)
# ---------------------------------------------------------------------
nb <- get_rds("lab8-operation-impact-nbhd.rds")
table(nb$first_zone_year[nb$year == 2004])      # first designation, by school year

set.seed(1)
show <- sample(unique(nb$nbhd[nb$first_zone_year > 0]), 30)
ggplot(subset(nb, nbhd %in% show), aes(year, factor(nbhd), fill = exposure)) +
  geom_tile() +
  scale_fill_gradient(low = "grey95", high = "firebrick") +
  labs(y = "Neighborhood", x = "School year (spring)",
       title = "Real impact zones: staggered, partial, and reversible") +
  theme_minimal()

# Q1. Name two ways this treatment differs from Card & Krueger's (Week 8).


# ---------------------------------------------------------------------
# PART 2. Legewie & Fagan's two main models (ever-exposed sample)
#   Eq. 1: neighborhood FE + age-by-year FE
#   Eq. 2: student FE + age-by-year FE
# ---------------------------------------------------------------------
eq1 <- feols(score ~ exposure:i(agegrp) | nbhd + age^year,
             data = s, split = ~group, cluster = ~nbhd)
eq2 <- feols(score ~ exposure:i(agegrp) | student_id + age^year,
             data = s, split = ~group, cluster = ~nbhd)
etable(eq1)
etable(eq2)
# Expect, Black boys 13-15: about -0.11 (Eq. 1) and -0.12 (Eq. 2).
# Truth: -0.13. Other groups: near zero.

pooled <- feols(score ~ exposure | nbhd + age^year, data = s, cluster = ~nbhd)
pooled                                    # about -0.02

# >>> OUTPUT 1: Eq. 1 and Eq. 2 for Black boys and Black girls, plus the
#     pooled estimate, in one table.
out1 <- etable(eq1[sample = "Black boys"], eq2[sample = "Black boys"],
               eq1[sample = "Black girls"], eq2[sample = "Black girls"], pooled,
               headers = c("Black boys, Eq. 1", "Black boys, Eq. 2",
                           "Black girls, Eq. 1", "Black girls, Eq. 2", "Pooled"))
writeLines(capture.output(print(out1)), "output1_lf_models.txt")
out1

# Q2 (Output 1). Why is the pooled estimate so small? What does that imply
#     about reporting a single average effect for a policy like this?


# ---------------------------------------------------------------------
# PART 3. Event studies: which comparison group?
#   Sample: Black boys aged 13-15. Reference: the year before designation.
# ---------------------------------------------------------------------
bb <- subset(d, group == "Black boys" & age >= 13)
bb$rel <- ifelse(bb$first_zone_year > 0,
                 pmin(pmax(bb$year - bb$first_zone_year, -4), 4), -1000)

# 3a. Compare zoned areas to areas NEVER zoned
es_a <- feols(score       ~ i(rel, ref = c(-1, -1000)) | nbhd + age^year, bb, cluster = ~nbhd)
es_b <- feols(score_trend ~ i(rel, ref = c(-1, -1000)) | nbhd + age^year, bb, cluster = ~nbhd)
iplot(list(es_a, es_b), main = "Never-zoned controls: A (black) vs. B (red)")
# Expect: A flat before 0. B leads of about +0.13, +0.05, +0.01, and
# "effects" that keep growing, to about -0.22 by year 4.
# Legewie & Fagan report that alternative control groups failed exactly
# this check. So they compared zoned areas only with each other, by TIMING.

# 3b. Their approach: only areas zoned at some point
bz <- subset(bb, first_zone_year > 0)
es_a2 <- feols(score       ~ i(rel, ref = -1) | nbhd + age^year, bz, cluster = ~nbhd)
es_b2 <- feols(score_trend ~ i(rel, ref = -1) | nbhd + age^year, bz, cluster = ~nbhd)
iplot(list(es_a2, es_b2), main = "Ever-zoned only: A (black) vs. B (red)")
# Expect: A and B now look alike, with flat leads.

# 3c. A formal pre-trend test: are the leads (-4, -3, -2) jointly zero?
wald(es_a,  keep = "rel::-[234]$")   # never-zoned controls, A: p about 0.43
wald(es_b,  keep = "rel::-[234]$")   # never-zoned controls, B: p < 0.001
wald(es_a2, keep = "rel::-[234]$")   # ever-zoned only, A:      p about 0.40
wald(es_b2, keep = "rel::-[234]$")   # ever-zoned only, B:      p about 0.32
# (Check that 'keep' printed three coefficients; adjust the pattern if not.)

# >>> OUTPUT 2: both event-study panels side by side, saved as one figure.
png("output2_event_studies.png", width = 1600, height = 650, res = 130)
par(mfrow = c(1, 2))
iplot(list(es_a, es_b), main = "Never-zoned controls", xlab = "Years since designation")
legend("bottomleft", c("A: score", "B: score_trend"), col = 1:2, pch = 20, bty = "n")
iplot(list(es_a2, es_b2), main = "Ever-zoned only", xlab = "Years since designation")
legend("bottomleft", c("A: score", "B: score_trend"), col = 1:2, pch = 20, bty = "n")
dev.off()
# Add the four p-values from 3c to your write-up.

# Q3 (Output 2). Why does restricting to ever-zoned areas rescue version B?
#     What must be true about areas zoned in 2006 vs. 2009 for this to
#     work? What kind of selection would it NOT fix?
# Q4 (Output 2). A p-value of 0.32 for B means "we can't reject parallel pre-trends."
#     Does it mean the trends ARE parallel? What would a small or noisy
#     sample do to this test? And if we only keep designs that pass it,
#     what happens to the estimates we end up reporting?


# ---------------------------------------------------------------------
# PART 4. Staggered timing: TWFE vs. Callaway & Sant'Anna
#   Ever-zoned areas only; compare each cohort to areas NOT YET zoned.
#   Treatment = "has been designated" (on, and stays on).
#   Students age in and out, so this is a repeated cross-section.
# ---------------------------------------------------------------------
bz <- subset(bz, first_zone_year != 2004)        # no pre-period for 2004
bz$designated <- as.integer(bz$year >= bz$first_zone_year)

twfe <- lapply(c("score", "score_persist"), function(y)
  feols(as.formula(paste(y, "~ designated | nbhd + age^year")), bz, cluster = ~nbhd))
etable(twfe)

cs <- lapply(c("score", "score_persist"), function(y)
  att_gt(yname = y, tname = "year", gname = "first_zone_year", data = bz,
         panel = FALSE, control_group = "notyettreated", clustervars = "nbhd"))
lapply(cs, aggte, type = "simple")
summary(cs[[1]])   # bottom line: "P-value for pre-test of parallel trends
                   # assumption" -- the C&S version of the test in 3c
ggdid(aggte(cs[[1]], type = "dynamic", min_e = -3, max_e = 4))
ggdid(aggte(cs[[2]], type = "dynamic", min_e = -3, max_e = 4))

# Expect (simple ATT):    A       C
#   TWFE designated     -0.08   -0.09
#   Callaway-Sant'Anna  -0.08   -0.14
# True average effect among designated students, version A: about -0.06.

# >>> OUTPUT 3: C&S dynamic effects for versions A and C, side by side.
#     Report the four simple ATTs (TWFE and C&S, A and C) in your write-up.
dyn_df <- function(x, label) {
  a <- aggte(x, type = "dynamic", min_e = -3, max_e = 4)
  data.frame(version = label, e = a$egt, att = a$att.egt, se = a$se.egt)
}
dyn <- rbind(dyn_df(cs[[1]], "A: score"), dyn_df(cs[[2]], "C: score_persist"))
p3 <- ggplot(dyn, aes(e, att)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(ymin = att - 1.96 * se, ymax = att + 1.96 * se)) +
  facet_wrap(~version) +
  labs(x = "Years since designation", y = "ATT (SD)",
       title = "Callaway & Sant'Anna, not-yet-treated controls") +
  theme_minimal()
ggsave("output3_cs_dynamic.png", p3, width = 9, height = 4.5, dpi = 150)
p3

# Q5 (Output 3). In version A, TWFE and C&S agree. In version C they don't. What is
#     different about C, and why does that hurt TWFE? (Who serves as a
#     "control" for a 2009 zone in TWFE?)
# Q6. In version A, C&S effects are largest in years 0-1 and then fade,
#     but the true per-year effect never fades. Why? Is "has been
#     designated" the estimand we want? Compare with Part 2.
# Q7. Legewie & Fagan code areas as UNEXPOSED once a zone is removed.
#     If effects linger (version C), what does that do to Part 2's
#     exposure estimates? Try it:
etable(feols(score_persist ~ exposure:i(agegrp) | nbhd + age^year,
             data = subset(s, group == "Black boys"), cluster = ~nbhd))
# Expect about -0.14 for 13-15, vs. -0.11 in version A: lingering effects
# get credited to whoever is currently zoned.
