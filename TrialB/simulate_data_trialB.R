# =====================================================================
# simulate_data_trialB.R  (publication version)
# Synthetic data generation for the competing-risks worked example
# ("Trial B") of the "Median Not Reached" article.
#
# One randomized two-arm trial, ONE data cut at month 48, one record
# per participant. Output variables: USUBJID, TRT, AVAL, CNSR,
# EVNTDESC. Nothing else is computed.
#
# Design intent:
#   - two competing first events: Progression (cause 1) and Death
#     without progression (cause 2);
#   - treatment barely changes the progression RATE (csHR ~ 0.90) but
#     strongly reduces the competing death rate (csHR ~ 0.35), so more
#     treated participants remain at risk of progressing; the
#     cumulative incidence of progression is slightly LOWER under
#     treatment at early horizons, the curves cross (around month 13),
#     and the treated incidence is higher thereafter - the
#     direction-discordance example (subdistribution HR > 1);
#   - naive 1-KM for progression (deaths censored) overestimates the
#     absolute probability, more strongly in the Control arm.
#
# Requires: base R only.
# =====================================================================

## Reproducibility: fix the seed AND the generator kinds explicitly.
set.seed(20260917, kind = "Mersenne-Twister",
         normal.kind = "Inversion", sample.kind = "Rejection")

## 1. Design parameters (all times in months) --------------------------
n_per_arm   <- 500
accrual_max <- 12               # uniform staggered enrolment window
dco         <- 48               # single data cut, months from study start
l1_ctrl     <- 0.020            # cause-specific hazard, Progression, Control
l2_ctrl     <- 0.030            # cause-specific hazard, Death,       Control
cshr1_trt   <- 0.90             # treatment multiplier on cause-1 hazard
cshr2_trt   <- 0.35             # treatment multiplier on cause-2 hazard

stopifnot(accrual_max < dco)    # no negative potential follow-up

## 2. Latent subject-level quantities ----------------------------------
## RNG contract: exactly THREE vectorized draws, in this order.
N   <- 2 * n_per_arm
trt <- rep(c("Control", "Treatment"), each = n_per_arm)

t_prog  <- rexp(N, l1_ctrl * ifelse(trt == "Treatment", cshr1_trt, 1))
t_death <- rexp(N, l2_ctrl * ifelse(trt == "Treatment", cshr2_trt, 1))
enrl    <- runif(N, 0, accrual_max)

## 3. Observed data at the data cut ------------------------------------
potential_fu <- dco - enrl
stopifnot(all(potential_fu > 0))
t_first  <- pmin(t_prog, t_death)
is_event <- t_first <= potential_fu            # tie -> event
first_is_prog <- t_prog <= t_death             # tie -> Progression

trialB <- data.frame(
  USUBJID  = sprintf("TRB-%04d", 1:N),
  TRT      = trt,
  AVAL     = round(pmin(t_first, potential_fu), 6),
  CNSR     = as.integer(!is_event),
  EVNTDESC = ifelse(!is_event, "Administrative censoring",
             ifelse(first_is_prog, "Progression",
                                   "Death without progression"))
)
write.csv(trialB, "trialB_example.csv", row.names = FALSE)
cat("Written trialB_example.csv:", nrow(trialB), "records\n")

## 4. Dataset fingerprint check ----------------------------------------
fingerprint <- sum(trialB$AVAL)
reference   <- 21167.150705
cat(sprintf("Sum of AVAL: %.6f (reference %.6f)\n", fingerprint, reference))
if (abs(fingerprint - reference) > 1e-6)
  warning("Dataset fingerprint does NOT match the frozen reference.")

sink("sessionInfo_simulation_trialB.txt"); print(sessionInfo()); sink()
