# =====================================================================
# simulate_data.R  (publication version)
# Synthetic data generation for the "Median Not Reached" worked example
# (Trial A). One randomized two-arm trial, TWO data cuts (Interim /
# Final), long format. Output variables: USUBJID, TRT, AVAL, CNSR,
# EVNTDESC, DCOID. No other variables are computed.
#
# Design intent:
#   Interim : median NR in BOTH arms (immature follow-up);
#   Final   : Control median finite; Treatment median NR with a
#             sustained plateau supported by a large risk set.
#
# Mechanism: mixture cure (Bernoulli-by-arm incidence; Weibull latency
# with an AFT treatment effect). Requires: base R only.
# =====================================================================

## Reproducibility: fix the seed AND the generator kinds explicitly.
## These are the current R defaults; naming them guards against any
## future change of defaults.
set.seed(20261011, kind = "Mersenne-Twister",
         normal.kind = "Inversion", sample.kind = "Rejection")

## 1. Design parameters (all times in months) --------------------------
n_per_arm    <- 140
accrual_max  <- 12              # uniform staggered enrolment window
dco_interim  <- 16
dco_final    <- 48
dropout_rate <- 0.005           # exponential dropout hazard per month
p_ltef_ctrl  <- 0.38            # P(long-term event-free), Control
p_ltef_trt   <- 0.55            # P(long-term event-free), Treatment
w_shape      <- 1.3
lat_med_ctrl <- 11              # median latency among Control susceptibles
b_trt_aft    <- 0.25            # AFT log time ratio, Treatment vs Control

stopifnot(accrual_max < dco_interim)   # every participant enrolled before
                                       # the earliest data cut (no negative
                                       # potential follow-up)

## 2. Latent subject-level quantities ----------------------------------
## RNG contract: exactly four vectorized draws, in this order.
N    <- 2 * n_per_arm
trt  <- rep(c("Control", "Treatment"), each = n_per_arm)

longterm   <- rbinom(N, 1, ifelse(trt == "Treatment",
                                  p_ltef_trt, p_ltef_ctrl))
scale0     <- lat_med_ctrl / log(2)^(1 / w_shape)
wb_draws   <- rweibull(N, shape = w_shape,
                       scale = scale0 * exp(b_trt_aft * (trt == "Treatment")))
enrl_time  <- runif(N, 0, accrual_max)
t_dropout  <- rexp(N, dropout_rate)

t_event <- wb_draws
t_event[longterm == 1] <- Inf   # long-term event-free: event never occurs

## 3. One record per participant per data cut --------------------------
make_cut <- function(dco, dcoid) {
  potential_fu <- dco - enrl_time
  stopifnot(all(potential_fu > 0))
  t_obs        <- pmin(t_event, t_dropout, potential_fu)
  is_event     <- t_event <= pmin(t_dropout, potential_fu)  # tie -> event
  data.frame(
    USUBJID  = sprintf("TTE-%03d", 1:N),
    TRT      = trt,
    AVAL     = round(t_obs, 6),
    CNSR     = as.integer(!is_event),
    EVNTDESC = ifelse(is_event, "Event",
               ifelse(potential_fu <= t_dropout,
                      "Administrative censoring", "Dropout")),
    DCOID    = dcoid
  )
}

tte <- rbind(make_cut(dco_interim, "Interim"),
             make_cut(dco_final,   "Final"))
write.csv(tte, "tte_example.csv", row.names = FALSE)
cat("Written tte_example.csv:", nrow(tte), "records\n")

## 4. Dataset fingerprint check ----------------------------------------
## The reference dataset is frozen; any deviation indicates an RNG or
## derivation change and must be investigated before the data are used.
fingerprint <- sum(tte$AVAL)
reference   <- 9474.695804
cat(sprintf("Sum of AVAL: %.6f (reference %.6f)\n", fingerprint, reference))
if (abs(fingerprint - reference) > 1e-6)
  warning("Dataset fingerprint does NOT match the frozen reference.")

sink("sessionInfo_simulation.txt"); print(sessionInfo()); sink()
