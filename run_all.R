# =====================================================================
# One command for the whole R part of the pilot.
# Before this, run once in a terminal:  python make_hourly.py
# Then in RStudio (working directory = this folder):  source("run_all.R")
# Takes about 10 minutes. Writes every results file and session_info.txt.
# =====================================================================

set.seed(2026)                       # fixed seed, in case any step draws random numbers
source("build_A.R")                  # builds and checks A, saves A.rds
source("tests/test_pilot.R")         # stops here if the data or building blocks are wrong
source("pilot_core.R")               # run 1: 28-day window, 2 base models, 95% coverage
source("pilot_training_length.R")    # run 2: training windows of 2, 4, 7 and 28 days
source("pilot_upgrade.R")            # run 3: 80/90/95%, interval score, per level, plot

writeLines(capture.output(sessionInfo()), "session_info.txt")
cat("\nDone. Exact R and package versions saved in session_info.txt\n")
