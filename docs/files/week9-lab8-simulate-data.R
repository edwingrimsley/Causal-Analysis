# =====================================================================
# Causal Analysis (QMSS, CUNY GC) -- Week 9
# Simulated students on REAL Operation Impact zone histories.
#
# Modeled on Legewie & Fagan (2019, ASR). Their student records are
# restricted, but the zones are not: MacDonald, Fagan & Geller (2016,
# PLOS ONE) posted block-group-by-month zone data for 2004-2012.
#
#   REAL (week9-lab8-build-zones.R -> lab8-impact-zone-neighborhoods.csv):
#     - which areas were zoned, in which months, phases III-XVII,
#       including zones that end and areas re-designated later
#     - precincts, block-group counts, crimes and stops per area
#     - racial mix of stops by precinct (NYPD files; Stanton 2024)
#   SIMULATED:
#     - students, their race, and their test scores
#     - the effect: Black boys only, growing with age, about -0.15 at 15
#
# Neighborhoods = block groups in the same precinct with the same
# month-by-month zone history (Legewie & Fagan split precincts along
# zone boundaries). Students per neighborhood scale with its number of
# block groups; never-zoned areas are sampled at a lower rate.
#
# Three outcome versions, one lesson each:
#   score          parallel trends hold; effect only while zoned
#   score_trend    zoned areas were already declining, never-zoned were not
#   score_persist  effects linger after a zone is removed (carryover)
#
# Inputs : lab8-impact-zone-neighborhoods.csv, lab8-nypd-precinct-stops-2003-2012.csv
#          (read from the course site)
# Outputs: impact_students.csv, impact_nbhd_year.csv in your working folder.
#          The lab itself loads the same data from the site as .rds files.
# Base R only.
# =====================================================================

set.seed(4242)

# ---- 1. Real neighborhoods and zone histories ---------------------------
base <- "https://edwingrimsley.github.io/Causal-Analysis/files/"
nb <- read.csv(paste0(base, "lab8-impact-zone-neighborhoods.csv"), colClasses = c(history = "character"))
N  <- nrow(nb)
Tm <- t(sapply(strsplit(nb$history, ""), as.integer))  # N x 108, month 1 = Jan 2004

stops <- read.csv(paste0(base, "lab8-nypd-precinct-stops-2003-2012.csv"))
s03   <- subset(stops, year == 2003)
bsh_p <- (s03$stops_black / s03$stops)[match(nb$precinct, s03$precinct)]
bsh_p[is.na(bsh_p)] <- mean(bsh_p, na.rm = TRUE)

crime_rate <- as.numeric(scale(log((nb$crimes + 1) / nb$n_blockgroups)))
black_sh <- plogis(qlogis(pmin(pmax(bsh_p, .02), .98)) * 0.9 + 0.5 * crime_rate +
                   rnorm(N, 0, 0.5))
hisp_sh  <- 0.8                                       # Hispanic share of non-Black

mu_n <- -0.25 * crime_rate + rnorm(N, 0, 0.15)       # neighborhood score level

# ---- 2. Exposure by school year -----------------------------------------
# School year 2004 = Sep 2003-Jun 2004. Zones before Jan 2004 (phases I-II)
# are not in the data and count as unexposed, as in Legewie & Fagan.
years <- 2004:2012
sy <- lapply(years, function(y) {
  m <- c((y - 2005) * 12 + 9:12, (y - 2004) * 12 + 1:6)   # Sep..Dec, Jan..Jun
  m[m >= 1]
})
exposure <- sapply(sy, function(m) rowSums(Tm[, m, drop = FALSE]) / 10)
# "different shocks": the early phases (2004-2006) were more intense
month_year <- 2004 + (0:107) %/% 12
Dm  <- sweep(Tm, 2, ifelse(month_year <= 2006, 1.25, 0.88), `*`)
dose <- sapply(sy, function(m) rowSums(Dm[, m, drop = FALSE]) / 10)
dose <- dose / (sum(dose) / sum(exposure))          # one full year averages 1.0
colnames(exposure) <- colnames(dose) <- years

carry <- matrix(0, N, length(years))                  # carryover, decays 40%/yr
for (k in 2:length(years)) carry[, k] <- 0.6 * (carry[, k - 1] + dose[, k - 1])

first_zone_year <- apply(exposure, 1, function(x) if (any(x > 0)) years[which(x > 0)[1]] else 0)

# Version B: areas that were ever zoned were already declining (about
# 0.03 SD a year), whenever their zone arrived; never-zoned areas were not
lambda_n <- ifelse(first_zone_year > 0, -0.03, 0) + rnorm(N, 0, 0.008)

nbhd_year <- data.frame(
  nbhd = rep(nb$nbhd, times = length(years)), precinct = rep(nb$precinct, times = length(years)),
  year = rep(years, each = N), exposure = round(as.vector(exposure), 3),
  dose = round(as.vector(dose), 3), first_zone_year = rep(first_zone_year, times = length(years)),
  n_blockgroups = rep(nb$n_blockgroups, times = length(years)))

# ---- 3. Students ----------------------------------------------------------
# Per block group per cohort: 2 students in zoned areas, 0.5 elsewhere.
cohorts <- 1998:2012                                   # year turning 9
rate <- ifelse(nb$ever_zone == 1, 2, 0.5) * nb$n_blockgroups
n_nc <- rpois(N * length(cohorts), rep(rate, each = length(cohorts)))
stu  <- data.frame(nbhd   = rep(rep(seq_len(N), each = length(cohorts)), n_nc),
                   cohort = rep(rep(cohorts, times = N), n_nc))
stu$student_id <- seq_len(nrow(stu))
stu$race <- ifelse(runif(nrow(stu)) < black_sh[stu$nbhd], "Black",
            ifelse(runif(nrow(stu)) < hisp_sh, "Hispanic", "White"))
stu$female <- rbinom(nrow(stu), 1, 0.5)
stu$alpha  <- rnorm(nrow(stu), 0, 0.6)

panel <- merge(stu, data.frame(year = years))
panel$age <- panel$year - panel$cohort + 9
panel <- panel[panel$age >= 9 & panel$age <= 15, ]
panel <- panel[order(panel$student_id, panel$year), ]
idx <- cbind(panel$nbhd, match(panel$year, years))
panel$exposure <- exposure[idx]; panel$dose <- dose[idx]; panel$carry <- carry[idx]

# ---- 4. True effects --------------------------------------------------------
beta_age <- c(`9` = 0, `10` = 0, `11` = -0.02, `12` = -0.06,
              `13` = -0.10, `14` = -0.13, `15` = -0.15)
panel$black_boy <- as.integer(panel$race == "Black" & panel$female == 0)
panel$true_beta <- panel$black_boy * beta_age[as.character(panel$age)]
panel$true_effect         <- panel$true_beta * panel$dose
panel$true_effect_persist <- panel$true_beta * (panel$dose + panel$carry)

# ---- 5. Outcomes -------------------------------------------------------------
ay_fe <- rnorm(length(years) * 7, 0, 0.05)
ay    <- ay_fe[(match(panel$year, years) - 1) * 7 + (panel$age - 8)]
race_gap <- ifelse(panel$race == "White", 0.35, 0)
eps <- rnorm(nrow(panel), 0, 0.45)
panel$score <- panel$alpha + mu_n[panel$nbhd] + race_gap + ay + panel$true_effect + eps
panel$score_trend   <- panel$score + lambda_n[panel$nbhd] * (panel$year - 2008)
panel$score_persist <- panel$score - panel$true_effect + panel$true_effect_persist

panel$nbhd <- nb$nbhd[panel$nbhd]
panel$precinct <- nb$precinct[match(panel$nbhd, nb$nbhd)]
panel$first_zone_year <- first_zone_year[match(panel$nbhd, nb$nbhd)]
panel$ever_zone <- as.integer(panel$first_zone_year > 0)

out <- panel[, c("student_id", "nbhd", "precinct", "year", "age", "race", "female",
                 "black_boy", "exposure", "dose", "first_zone_year", "ever_zone",
                 "score", "score_trend", "score_persist", "true_effect")]
for (v in c("score", "score_trend", "score_persist", "true_effect")) out[[v]] <- round(out[[v]], 4)
out$exposure <- round(out$exposure, 3); out$dose <- round(out$dose, 3)
write.csv(out, "impact_students.csv", row.names = FALSE)
write.csv(nbhd_year, "impact_nbhd_year.csv", row.names = FALSE)

# ---- 6. Simulation vs. paper -------------------------------------------------
ez <- out[out$student_id %in% unique(out$student_id[out$exposure > 0]), ]
cat("Neighborhoods:", N, " zoned:", sum(nb$ever_zone),
    " precincts with zones:", length(unique(nb$precinct[nb$ever_zone == 1])), "of 76\n")
cat("Mean exposure, students in zones:", round(mean(ez$exposure[ez$exposure > 0]), 2), "(paper 0.77)\n")
cat("Mean exposure, students ever exposed:", round(mean(ez$exposure), 2), "(paper 0.26)\n")
cat("Students ever exposed:", length(unique(ez$student_id)), " student-years:", nrow(ez), "\n")
cat("Neighborhoods first zoned, by school year:\n"); print(table(first_zone_year[first_zone_year > 0]))
