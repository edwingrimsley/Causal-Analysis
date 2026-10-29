# =====================================================================
# Build real Operation Impact zone histories from MacDonald, Fagan &
# Geller (2016, PLOS ONE) public replication data:
#   github.com/macdonaldjohn/Impact-Zone-Data
#   file: blocks2004_2012_crime_fid_impactzones.csv (download and unzip
#   blocks2004_2012_crime_fid_impactzones.zip from that repository first)
#
# Unit in their data: census block group x precinct x month, 2004-2012.
# Block groups are anonymized ("fid"). Each unit carries flags for
# membership in impact-zone phases III-XVII (impact3 ... impact17).
#
# Treatment timing follows the authors' Stata code: a phase's areas are
# zoned for that calendar year; mid-year phases add areas from June,
# July or August through December. This rule reproduces their
# 'treatment' variable in all but 109 of 840,297 rows.
#
# Neighborhoods follow Legewie & Fagan's idea of splitting precincts
# along zone boundaries: block groups in the same precinct with the
# same month-by-month zone history form one neighborhood.
#
# Output: lab8-impact-zone-neighborhoods.csv (one row per neighborhood)
# =====================================================================

f <- "blocks2004_2012_crime_fid_impactzones.csv"
hdr <- names(read.csv(f, nrows = 1))
phase <- c(paste0("impact", 3:15), "impact16_a", "impact16_h", "impact17")
keep <- c("fid", "pct", "borough", "crimes1", "stops", phase)
cc <- ifelse(hdr %in% keep, NA, "NULL")
d <- read.csv(f, colClasses = cc)

# one row per block-group x precinct unit
for (p in phase) d[[p]] <- ifelse(is.na(d[[p]]), 0, d[[p]] > 0)
d$stops[is.na(d$stops)] <- 0
key <- paste(d$fid, d$pct)
U <- data.frame(key = unique(key))
U$fid <- d$fid[match(U$key, key)]; U$pct <- d$pct[match(U$key, key)]
U$borough <- d$borough[match(U$key, key)]
U$crimes <- tapply(d$crimes1, key, sum)[U$key]
U$stops  <- tapply(d$stops,   key, sum)[U$key]
for (p in phase) U[[p]] <- tapply(d[[p]], key, max)[U$key] > 0

# phase start (year, month); each runs through December of that year
start <- list(impact3 = c(2004, 1), impact4 = c(2005, 1), impact5 = c(2005, 7),
              impact6 = c(2006, 1), impact7 = c(2006, 6), impact8 = c(2007, 1),
              impact9 = c(2007, 7), impact10 = c(2008, 1), impact11 = c(2008, 7),
              impact12 = c(2009, 1), impact13 = c(2009, 7), impact14 = c(2010, 1),
              impact15 = c(2010, 8), impact16_a = c(2011, 1), impact16_h = c(2011, 8),
              impact17 = c(2012, 1))
Tm <- matrix(0L, nrow(U), 108)                       # month 1 = Jan 2004
for (p in phase) {
  s <- (start[[p]][1] - 2004) * 12 + start[[p]][2]
  e <- (start[[p]][1] - 2004) * 12 + 12
  Tm[U[[p]], s:e] <- 1L
}
U$history <- apply(Tm, 1, paste, collapse = "")

# neighborhoods: precinct x identical zone history
U$nbhd_key <- paste(U$pct, U$history)
nb <- aggregate(cbind(n_blockgroups = 1, crimes, stops) ~ nbhd_key + pct + borough + history,
                data = U, FUN = sum)
nb <- nb[order(nb$pct, nb$history), ]
nb$nbhd <- seq_len(nrow(nb))
nb$ever_zone <- as.integer(grepl("1", nb$history))
out <- nb[, c("nbhd", "pct", "borough", "n_blockgroups", "crimes", "stops", "ever_zone", "history")]
names(out)[2] <- "precinct"
write.csv(out, "lab8-impact-zone-neighborhoods.csv", row.names = FALSE)

cat("Block-group units:", nrow(U), " ever zoned:", sum(grepl("1", U$history)), "\n")
cat("Neighborhoods:", nrow(out), " zoned:", sum(out$ever_zone),
    " precincts with zones:", length(unique(out$precinct[out$ever_zone == 1])), "\n")
