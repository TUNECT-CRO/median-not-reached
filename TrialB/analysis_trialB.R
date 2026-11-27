# ----------------------------------------------------------
# Sponsor: TUNECT
# Project Name: Median Not Reached. Supplementary Materials. Trial B.
#
# Program Name: analysis_trialB
# Description : R implementation of the methods presented in the
#               article for survival analysis with competing risks.
#
# Programmer : Elizaveta Tikhonova
#
# R Version : 4.3.3
# Operating System: Windows 10 x64
# ----------------------------------------------------------

# =====================================================================
# analysis_trialB.R  (publication version)
# Reference R analysis for the competing-risks worked example (Trial B).
# Reads ONLY trialB_example.csv. Writes display-ready tables (all
# values rounded to 3 decimal places) and the figure to outputs_trialB/.
#
# Analyses:
#   B1  First-event counts by arm.
#   B2  Probability of progression at 6/12/24/36 months: Aalen-Johansen
#       CIF versus the naive 1-KM with deaths censored; CIF differences
#       (Treatment - Control) by horizon, and the crossing time.
#   B3  Cause-specific Cox models (progression; death).
#   B4  Fine-Gray subdistribution model for progression.
#   Figure: CIF (solid) vs 1-KM (dashed) for progression, by arm,
#       with numbers at risk.
#
# Methods conventions:
#   - the event type is derived from EVNTDESC; CNSR = 1 means censored;
#     inputs are validated and any inconsistency stops the script;
#   - CIF: Aalen-Johansen estimator (multistate survfit);
#   - naive comparator: KM with competing deaths treated as censored,
#     reported as 1 - KM (intentionally incorrect, for illustration);
#   - cause-specific Cox: competing events censored; ties = Efron;
#     interpretation is a rate ratio among participants still free of
#     any event;
#   - Fine-Gray: survival::finegray weighted dataset + Cox model with
#     participant-level clustering (robust sandwich variance by
#     USUBJID); the SHR conveys the direction of the covariate effect
#     on the CIF, not its magnitude;
#   - group contrasts are matched by stratum NAME, never by position;
#   - all displayed numbers are rounded to 3 decimal places.
# Requires: R >= 4.0, package 'survival' (>= 3.0 for finegray()).
# =====================================================================

library(survival)
dir.create("outputs_trialB", showWarnings = FALSE)

## ---- Input and derived variables ------------------------------------
d <- read.csv("trialB_example.csv", stringsAsFactors = FALSE)

allowed_desc <- c("Progression", "Death without progression",
                  "Administrative censoring")
stopifnot(nrow(d) == 1000,
          all(d$EVNTDESC %in% allowed_desc),
          all((d$CNSR == 1) == (d$EVNTDESC == "Administrative censoring")))

d$TRT    <- factor(d$TRT, c("Control", "Treatment"))
d$TRTN   <- as.integer(d$TRT == "Treatment")        # derived, not stored
d$STATUS <- ifelse(d$CNSR == 1, 0L,
            ifelse(d$EVNTDESC == "Progression", 1L, 2L))
d$fstatus <- factor(d$STATUS, 0:2,
                    c("censored", "progression", "death"))
arms <- levels(d$TRT)
fmt3 <- function(x) sprintf("%.3f", x)
fmtp <- function(p) ifelse(p < 0.0001, "<.0001", sprintf("%.4f", p))

## ---- B1. First-event counts -----------------------------------------
npct <- function(n, N) sprintf("%d (%.1f%%)", n, 100 * n / N)  # "180 (36.0%)"
b1 <- do.call(rbind, lapply(arms, function(a) {
  x <- subset(d, TRT == a)
  data.frame(TRT = a,
             Progression = npct(sum(x$EVNTDESC == "Progression"), nrow(x)),
             Death_without_progression =
               npct(sum(x$EVNTDESC == "Death without progression"), nrow(x)),
             Administrative_censoring =
               npct(sum(x$EVNTDESC == "Administrative censoring"), nrow(x)))
}))
names(b1) <- c("Treatment", "Progression", "Death without progression",
               "Administrative censoring")
write.csv(b1, "outputs_trialB/tB1_counts.csv", row.names = FALSE)
cat("=== B1 First-event counts ===\n"); print(b1, row.names = FALSE)

## ---- B2. CIF (Aalen-Johansen) vs naive 1-KM -------------------------
b2_times <- c(6, 12, 24, 36)
aj  <- survfit(Surv(AVAL, fstatus) ~ TRT, data = d,
               conf.type = "log-log")   # multistate (Aalen-Johansen)
kmn <- survfit(Surv(AVAL, STATUS == 1) ~ TRT, data = d,
               conf.type = "log-log")  # deaths censored

cif_at <- function(a, tt) {    # name-based lookup, one arm and time
  s <- summary(aj[paste0("TRT=", a), ], times = tt)
  s$pstate[, which(aj$states == "progression")]
}
km_at <- function(a, tt) 1 - summary(kmn[paste0("TRT=", a)], times = tt)$surv

## Interval conventions: log(-log) transformation for BOTH estimators,
## matching the project convention for KM and the SAS PROC LIFETEST
## default for CIFs. CIF limits: log-log on the Aalen-Johansen
## estimate; 1-KM limits: the log-log KM interval for S(t) with both
## limits transformed (1 - upper; 1 - lower), so the limits cross over.
fmt_ci <- function(e, l, u)
  sprintf("%s [%s, %s]", fmt3(e), fmt3(l), fmt3(u))
cif_cells <- function(a, tt) {
  s  <- summary(aj[paste0("TRT=", a), ], times = tt)
  ix <- which(aj$states == "progression")
  fmt_ci(s$pstate[, ix], s$lower[, ix], s$upper[, ix])
}
km_cells <- function(a, tt) {
  s <- summary(kmn[paste0("TRT=", a)], times = tt)
  fmt_ci(1 - s$surv, 1 - s$upper, 1 - s$lower)   # limits cross over
}
b2 <- do.call(rbind, lapply(arms, function(a)
  data.frame(Treatment = a, Timepoints = b2_times,
             CIF = cif_cells(a, b2_times),
             KM  = km_cells(a, b2_times))))
names(b2) <- c("Treatment", "Timepoints", "CIF [95%]", "1 - KM [95%]")
write.csv(b2, "outputs_trialB/tB2_cif_vs_km.csv", row.names = FALSE)
cat("\n=== B2 Probability of progression: CIF vs naive 1-KM ===\n")
print(b2, row.names = FALSE)

b2d <- data.frame(Parameter = "Cumulative Incidence Function",
                  Timepoints = b2_times,
                  Trt = fmt3(cif_at("Treatment", b2_times)),
                  Ctl = fmt3(cif_at("Control",   b2_times)),
                  Diff = fmt3(cif_at("Treatment", b2_times) -
                              cif_at("Control",   b2_times)))
names(b2d) <- c("Parameter", "Timepoints", "Treatment", "Control",
                "Difference, Treatment - Control")
write.csv(b2d, "outputs_trialB/tB2_cif_diff.csv", row.names = FALSE)
cat("\n=== B2 CIF difference (Treatment - Control) by horizon ===\n")
print(b2d, row.names = FALSE)

## Crossing time of the two CIF curves: the last sign change on the
## merged grid of observed event times (single early events make the
## difference oscillate around zero, so the first sign change is noise;
## the curves separate for good only after the last change)
evt_grid <- sort(unique(d$AVAL[d$STATUS > 0]))
dd <- cif_at("Treatment", evt_grid) - cif_at("Control", evt_grid)
neg <- which(dd < 0)
cross_t <- NA
if (length(neg) > 0 && max(neg) < length(dd))
  cross_t <- evt_grid[max(neg) + 1]

## ---- B3. Cause-specific Cox models ----------------------------------
cs1 <- coxph(Surv(AVAL, STATUS == 1) ~ TRTN, data = d, ties = "efron")
cs2 <- coxph(Surv(AVAL, STATUS == 2) ~ TRTN, data = d, ties = "efron")
ci_of <- function(m) summary(m)$conf.int   # HR, lower, upper
## One display row per model: log-HR and its SE (robust SE and robust
## p-value when the model carries a cluster term), HR with 95% CI, p.
hr_row <- function(m, label, hr_name = "HR") {
  sm <- summary(m); cf <- sm$coefficients; ci <- sm$conf.int
  se_col <- if ("robust se" %in% colnames(cf)) "robust se" else "se(coef)"
  out <- data.frame(Model = label,
                    coef = fmt3(cf[1, "coef"]),
                    SE   = fmt3(cf[1, se_col]),
                    HR   = fmt3(ci[1]),
                    LCL  = fmt3(ci[3]), UCL = fmt3(ci[4]),
                    p_value = fmtp(cf[1, "Pr(>|z|)"]))
  names(out)[names(out) == "HR"] <- hr_name
  out
}
b3 <- rbind(
  hr_row(cs1, "Cause-specific Cox: progression (deaths censored)"),
  hr_row(cs2, "Cause-specific Cox: death without progression"))
names(b3) <- c("Model", "Coefficient", "Standard Error", "Hazard Ratio",
               "Hazard Ratio, 95% CI Lower Limit",
               "Hazard Ratio, 95% CI Upper Limit", "p-value")
write.csv(b3, "outputs_trialB/tB3_cs_cox.csv", row.names = FALSE)
cat("\n=== B3 Cause-specific hazard ratios (Treatment vs Control) ===\n")
print(b3, row.names = FALSE)

## ---- B4. Fine-Gray subdistribution model for progression ------------
fgd <- finegray(Surv(AVAL, fstatus) ~ ., data = d, etype = "progression")
fg  <- coxph(Surv(fgstart, fgstop, fgstatus) ~ TRTN + cluster(USUBJID),
             weights = fgwt, data = fgd, ties = "efron")
cif_fg <- ci_of(fg)
b4 <- hr_row(fg, "Fine-Gray subdistribution model: progression",
             hr_name = "SHR")
names(b4) <- c("Model", "Coefficient", "Standard Error", "Hazard Ratio",
               "Hazard Ratio, 95% CI Lower Limit",
               "Hazard Ratio, 95% CI Upper Limit", "p-value")
write.csv(b4, "outputs_trialB/tB4_finegray.csv", row.names = FALSE)
cat("\n=== B4 Fine-Gray subdistribution HR (Treatment vs Control) ===\n")
print(b4, row.names = FALSE)
cat("   (participant-clustered robust variance, cluster = USUBJID)\n")
cat(sprintf("log-likelihoods: csCox progression %.4f; csCox death %.4f; Fine-Gray %.4f\n",
            logLik(cs1), logLik(cs2), logLik(fg)))

## ---- Interpretation (generated from the estimates, not hard-coded) --
desc_hr <- function(ci) {
  if (ci[3] < 1 && ci[4] > 1) "no clear difference (the 95% CI includes 1)"
  else if (ci[4] < 1) "a lower rate"
  else "a higher rate"
}
cat("\nINTERPRETATION (direction, not magnitude):\n")
cat(sprintf(" csHR progression = %s (%s-%s): %s in the progression rate\n",
            fmt3(ci_of(cs1)[1]), fmt3(ci_of(cs1)[3]), fmt3(ci_of(cs1)[4]),
            desc_hr(ci_of(cs1))))
cat("   among participants still free of any event;\n")
cat(sprintf(" csHR death = %s (%s-%s): treatment strongly lowers the\n",
            fmt3(ci_of(cs2)[1]), fmt3(ci_of(cs2)[3]), fmt3(ci_of(cs2)[4])))
cat("   competing death rate, keeping more participants at risk;\n")
cat(sprintf(" SHR progression = %s (%s-%s): the cumulative incidence of\n",
            fmt3(cif_fg[1]), fmt3(cif_fg[3]), fmt3(cif_fg[4])))
cat("   progression is higher under treatment over the follow-up as a whole.\n")
if (!is.na(cross_t))
  cat(sprintf(" The CIF curves cross at about month %.2f: the difference is %s\n",
              cross_t, paste(b2d$CIF_diff_TminusC, collapse = " / ")))
cat("   at months", paste(b2_times, collapse = "/"),
    "- the subdistribution hazards are not proportional, so the single\n")
cat("   SHR is a follow-up-weighted average. Neither hazard ratio\n")
cat("   measures the change in absolute risk: read that from the CIFs\n")
cat("   at prespecified time points (tB2).\n")

## ---- Figure: CIF vs 1-KM for progression, with numbers at risk ------
cols <- c(Control = "grey45", Treatment = "black")
fig_grid <- seq(0, 48, 12)
draw_figB <- function() {
  par(mar = c(8, 4, 3, 1))
  plot(aj[, "progression"], col = cols, lwd = 2, lty = 1,
       xlim = c(0, 48), ylim = c(0, 0.8), xaxt = "n",
       xlab = "Months since randomization",
       ylab = "Probability of progression",
       main = "Progression: Aalen-Johansen CIF vs naive 1-KM")
  axis(1, at = fig_grid)
  lines(kmn, fun = "event", col = cols, lwd = 2, lty = 2)
  legend("topleft",
         c("Control CIF", "Treatment CIF",
           "Control 1-KM (naive)", "Treatment 1-KM (naive)"),
         col = rep(cols, 2), lty = c(1, 1, 2, 2), lwd = 2,
         bty = "n", cex = 0.9)
  s <- summary(aj, times = fig_grid, extend = TRUE)  # at risk of ANY event
  ## multistate n.risk is a matrix (time x state): take the initial
  ## state's column = participants still free of any event
  s0 <- which(aj$states == "(s0)"); if (length(s0) == 0) s0 <- 1
  lab_at <- -0.065 * 48
  mtext("At risk", side = 1, line = 3.9, at = lab_at, adj = 1,
        cex = 0.75, font = 2)
  for (a in arms) {
    ln <- 3.9 + which(arms == a)
    mtext(a, side = 1, line = ln, at = lab_at, adj = 1, cex = 0.75)
    nr <- s$n.risk[s$strata == paste0("TRT=", a), s0]
    mtext(sprintf("%d", nr), side = 1, line = ln, at = fig_grid, cex = 0.75)
  }
}
png("outputs_trialB/figB_cif_vs_km.png", width = 7.5, height = 6,
    units = "in", res = 300); draw_figB(); dev.off()
pdf("outputs_trialB/figB_cif_vs_km.pdf", width = 7.5, height = 6)
draw_figB(); dev.off()

sink("outputs_trialB/sessionInfo_analysis_trialB.txt")
print(sessionInfo()); sink()
cat("\nAll tables and the figure written to outputs_trialB/\n")
