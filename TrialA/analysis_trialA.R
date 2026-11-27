# ----------------------------------------------------------
# Sponsor: TUNECT
# Project Name: Median Not Reached. Supplementary Materials. Trial A.
#
# Program Name: analysis_R
# Description : R implementation of the methods presented in the
#               article for survival analysis when the median is
#               not reached.
#
# Programmer : Elizaveta Tikhonova
#
# R Version : 4.3.3
# Operating System: Windows 10 x64
# ----------------------------------------------------------

# =====================================================================
# analysis_R.R  (publication version)
# Reference R analysis for the "Median Not Reached" worked example
# (Trial A). Reads ONLY tte_example.csv. Writes display-ready tables
# (all values rounded to 3 decimal places; NR/NE conventions applied)
# and the figures to outputs/.
#
# Analyses, per data cut unless stated otherwise:
#   5.1  Descriptive summary and numbers at risk on frozen time grids.
#   5.2  Kaplan-Meier medians with NR / NE display rules.
#   5.3  Follow-up by the reverse Kaplan-Meier method.
#   5.4  Fixed-time event-free probabilities and differences (Final).
#   5.5  RMST; the horizon tau is defined by a prespecified RULE:
#        the largest time supported by the observed data in BOTH arms
#        (minimum over arms of the largest observed time), so that no
#        segment of either curve is extrapolated.
#   5.6  Logistic-Weibull mixture cure model (Final): incidence
#        (probability of being susceptible) and latency components.
#
# Conventions:
#   - one survfit() call per data cut, conf.type = "log-log"
#     (Brookmeyer-Crowley intervals for the median) reused everywhere;
#   - a non-estimable point estimate is displayed as NR, a
#     non-estimable confidence limit as NE;
#   - group contrasts are matched by stratum NAME, never by position;
#   - all displayed numbers are rounded to 3 decimal places.
# Requires: R >= 4.0, package 'survival'.
# =====================================================================

library(survival)
dir.create("outputs", showWarnings = FALSE)

## ---- Input and derived variables ------------------------------------
tte <- read.csv("tte_example.csv", stringsAsFactors = FALSE)

allowed_desc <- c("Event", "Administrative censoring", "Dropout")
stopifnot(nrow(tte) == 560,
          all(tte$EVNTDESC %in% allowed_desc),
          all((tte$CNSR == 0) == (tte$EVNTDESC == "Event")))

tte$EVENT <- 1 - tte$CNSR
tte$TRTN  <- as.integer(tte$TRT == "Treatment")  # derived, not stored
tte$TRT   <- factor(tte$TRT, c("Control", "Treatment"))
cuts <- c("Interim", "Final")
arms <- levels(tte$TRT)

## ---- Display helpers (3 decimal places everywhere) ------------------
fmt3  <- function(x) sprintf("%.3f", x)
fmt1  <- function(x) sprintf("%.1f", x)    # journal convention for times
fmt2  <- function(x) sprintf("%.2f", x)
dnum  <- function(x, miss) ifelse(is.na(x), miss, fmt3(x))  # NR or NE
dnum1 <- function(x, miss) ifelse(is.na(x), miss, fmt1(x))  # NR or NE, 1 dp
dnum2 <- function(x, miss) ifelse(is.na(x), miss, fmt2(x))  # NR or NE, 2 dp

## ---- One KM fit per data cut, log-log CI, reused everywhere ---------
km_fit <- lapply(cuts, function(cc)
  survfit(Surv(AVAL, EVENT) ~ TRT, data = subset(tte, DCOID == cc),
          conf.type = "log-log"))
names(km_fit) <- cuts
stratum <- function(km, arm) paste0("TRT=", arm)  # name-based row lookup

## =====================================================================
## 5.1 Descriptive summary and numbers at risk (frozen grids)
## =====================================================================
npct <- function(n, N) sprintf("%d (%.1f%%)", n, 100 * n / N)  # "41 (29.3%)"
desc <- do.call(rbind, lapply(cuts, function(cc) {
  do.call(rbind, lapply(arms, function(a) {
    d <- subset(tte, DCOID == cc & TRT == a)
    data.frame(DCOID = cc, TRT = a, N = nrow(d),
               EVENT    = npct(sum(d$EVENT), nrow(d)),
               CENSORED = npct(sum(d$CNSR),  nrow(d)),
               Admin   = sum(d$EVNTDESC == "Administrative censoring"),
               Dropout = sum(d$EVNTDESC == "Dropout"))
  }))
}))
names(desc) <- c("Analysis Type", "Treatment", "Number of Subjects",
                 "Number of Events (%)", "Number of Censored (%)",
                 "Number of Administrative Censoring",
                 "Number of Dropout Censoring")
write.csv(desc, "outputs/tA1_descriptive.csv", row.names = FALSE)
cat("=== 5.1 Descriptive summary ===\n"); print(desc, row.names = FALSE)

grids <- list(Interim = c(0, 6, 12, 15), Final = c(0, 12, 24, 36, 42, 47))
at_risk <- do.call(rbind, lapply(cuts, function(cc) {
  s <- summary(km_fit[[cc]], times = grids[[cc]])
  data.frame(DCOID = cc, TRT = sub("TRT=", "", s$strata),
             time = s$time, n_risk = s$n.risk)
}))
names(at_risk) <- c("Analysis Type", "Treatment", "Timepoint",
                    "Number at Risk")
write.csv(at_risk, "outputs/tA1_at_risk.csv", row.names = FALSE)
cat("\n=== 5.1 Numbers at risk (frozen grids) ===\n")
print(at_risk, row.names = FALSE)

## =====================================================================
## 5.2 Kaplan-Meier medians with NR / NE display rules
## =====================================================================
med <- do.call(rbind, lapply(cuts, function(cc) {
  tb <- summary(km_fit[[cc]])$table
  data.frame(DCOID = cc, TRT = sub("TRT=", "", rownames(tb)),
             Median     = dnum1(tb[, "median"],  "NR"),
             Median_LCL = dnum1(tb[, "0.95LCL"], "NE"),
             Median_UCL = dnum1(tb[, "0.95UCL"], "NE"))
}))
names(med) <- c("Analysis Type", "Treatment", "Median",
                "95% Lower Limit", "95% Upper Limit")
write.csv(med, "outputs/tA2_median.csv", row.names = FALSE)
cat("\n=== 5.2 KM median (log-log / Brookmeyer-Crowley CI) ===\n")
print(med, row.names = FALSE)

## =====================================================================
## 5.3 Follow-up: reverse Kaplan-Meier (median, quartiles)
## =====================================================================
fu <- do.call(rbind, lapply(cuts, function(cc) {
  d  <- subset(tte, DCOID == cc)
  rk <- survfit(Surv(AVAL, CNSR) ~ TRT, data = d, conf.type = "log-log")
  tb <- summary(rk)$table
  q  <- quantile(rk, probs = c(0.25, 0.75))
  data.frame(DCOID = cc, TRT = sub("TRT=", "", rownames(tb)),
             FU_median     = dnum2(tb[, "median"],  "NR"),
             FU_median_LCL = dnum2(tb[, "0.95LCL"], "NE"),
             FU_median_UCL = dnum2(tb[, "0.95UCL"], "NE"),
             FU_Q1 = dnum2(q$quantile[, "25"], "NE"),
             FU_Q3 = dnum2(q$quantile[, "75"], "NE"))
}))
names(fu) <- c("Analysis Type", "Treatment", "Follow-Up Median",
               "Follow-Up Median 95% CI Lower Limit",
               "Follow-Up Median 95% CI Upper Limit",
               "Follow-Up Q1", "Follow-Up Q3")
write.csv(fu, "outputs/tA3_followup.csv", row.names = FALSE)
cat("\n=== 5.3 Potential follow-up (reverse Kaplan-Meier) ===\n")
print(fu, row.names = FALSE)

## =====================================================================
## 5.4 Fixed-time event-free probabilities (Final cut)
## =====================================================================
fix_times <- c(12, 24, 36)
dF <- subset(tte, DCOID == "Final")
## Guard: every requested time point must lie within the observed range
## of BOTH arms; otherwise the estimate would be unsupported.
max_obs <- tapply(dF$AVAL, dF$TRT, max)
if (any(fix_times > min(max_obs)))
  stop("A fixed time point exceeds the observed follow-up in an arm.")

sF <- summary(km_fit[["Final"]], times = fix_times)
fx <- data.frame(TRT = sub("TRT=", "", sF$strata), t = sF$time,
                 S   = fmt3(sF$surv),
                 LCL = dnum(sF$lower, "NE"),
                 UCL = dnum(sF$upper, "NE"))
names(fx) <- c("Treatment", "Timepoint", "Event-Free Probability",
               "95% CI Lower Limit", "95% CI Upper Limit")
write.csv(fx, "outputs/tA4_fixedtime.csv", row.names = FALSE)
cat("\n=== 5.4 Fixed-time event-free probabilities (Final) ===\n")
print(fx, row.names = FALSE)

## Difference Treatment - Control, matched by stratum name
fxd <- do.call(rbind, lapply(fix_times, function(tt) {
  i_t <- sF$strata == stratum(km_fit[["Final"]], "Treatment") & sF$time == tt
  i_c <- sF$strata == stratum(km_fit[["Final"]], "Control")   & sF$time == tt
  dr  <- sF$surv[i_t] - sF$surv[i_c]
  se  <- sqrt(sF$std.err[i_t]^2 + sF$std.err[i_c]^2)  # Greenwood, plain scale
  data.frame(Parameter = "Difference of Probabilities (Treatment - Control)",
             t = tt, diff = fmt3(dr),
             LCL = fmt3(dr - 1.96 * se), UCL = fmt3(dr + 1.96 * se),
             SE = fmt3(se))
}))
names(fxd) <- c("Parameter", "Timepoint", "Difference",
                "95% CI Lower Limit", "95% CI Upper Limit",
                "Standard Error")
write.csv(fxd, "outputs/tA4_fixedtime_diff.csv", row.names = FALSE)
cat("\n=== 5.4 Difference in event-free probability (T - C) ===\n")
print(fxd, row.names = FALSE)

## =====================================================================
## 5.5 RMST: tau = last commonly observed time (prespecified rule)
## =====================================================================
taus <- sapply(cuts, function(cc)
  min(tapply(tte$AVAL[tte$DCOID == cc], tte$TRT[tte$DCOID == cc], max)))
names(taus) <- cuts
cat(sprintf("\nRMST horizons (last commonly observed time): %s\n",
            paste(sprintf("%s %.3f", cuts, taus), collapse = ", ")))

rm_tab <- do.call(rbind, lapply(cuts, function(cc) {
  tb <- summary(km_fit[[cc]], rmean = taus[cc])$table
  r_t <- tb[stratum(km_fit[[cc]], "Treatment"), ]
  r_c <- tb[stratum(km_fit[[cc]], "Control"),   ]
  dr  <- r_t["rmean"] - r_c["rmean"]
  se  <- sqrt(r_t["se(rmean)"]^2 + r_c["se(rmean)"]^2)
  data.frame(DCOID = cc, tau = fmt3(taus[cc]),
             RMST_Control   = fmt3(r_c["rmean"]),
             RMST_Treatment = fmt3(r_t["rmean"]),
             Diff = fmt3(dr),
             LCL = fmt3(dr - 1.96 * se), UCL = fmt3(dr + 1.96 * se))
}))
names(rm_tab) <- c("Analysis Type", "Tau", "RMST, Control",
                   "RMST, Treatment",
                   "RMST Difference, Treatment - Control",
                   "RMST Difference, 95% CI Lower Limit",
                   "RMST Difference, 95% CI Upper Limit")
write.csv(rm_tab, "outputs/tA5_rmst.csv", row.names = FALSE)
cat("\n=== 5.5 RMST (area under KM up to tau; no extrapolation) ===\n")
print(rm_tab, row.names = FALSE)

## =====================================================================
## 5.6 Mixture cure model (Final cut): logistic incidence + Weibull AFT
## =====================================================================
## Shared model functions: used by BOTH the likelihood and Figure 2,
## so the fitted curves always reflect the estimated model.
cure_p  <- function(par, trtn) plogis(par[1] + par[2] * trtn)
cure_su <- function(par, t, trtn) {
  shp <- exp(-par[5]); sc <- exp(par[3] + par[4] * trtn)
  exp(-(t / sc)^shp)
}
cure_fu <- function(par, t, trtn) {
  shp <- exp(-par[5]); sc <- exp(par[3] + par[4] * trtn)
  (shp / sc) * (t / sc)^(shp - 1) * cure_su(par, t, trtn)
}
cure_spop <- function(par, t, trtn) {         # population survival
  p <- cure_p(par, trtn)
  1 - p + p * cure_su(par, t, trtn)
}

negll <- function(par, d) {
  p  <- cure_p(par, d$TRTN)
  Su <- cure_su(par, d$AVAL, d$TRTN)
  fu <- cure_fu(par, d$AVAL, d$TRTN)
  -sum(d$EVENT * log(pmax(p * fu, 1e-300)) +
       (1 - d$EVENT) * log(pmax(1 - p + p * Su, 1e-300)))
}

opt <- optim(c(0, 0, log(12), 0, log(1 / 1.3)), negll, d = dF,
             method = "BFGS", hessian = TRUE,
             control = list(reltol = 1e-12, maxit = 1000))
if (opt$convergence != 0)
  stop("Cure model did not converge (code ", opt$convergence, ").")
ch <- tryCatch(chol(opt$hessian), error = function(e) NULL)
if (is.null(ch))
  stop("Cure model Hessian is not positive definite at the optimum.")
V   <- chol2inv(ch)
est <- opt$par
se  <- sqrt(diag(V))

## Reported quantities. Confidence intervals are Wald intervals on the
## ESTIMATION scale, with limits transformed monotonically to the
## reporting scale (never the delta method on the reporting scale):
##   LTEF  = 1 - plogis(eta): decreasing, limits cross over;
##   OR    = exp(-g1): decreasing, limits cross over;
##   shape = exp(-log_sigma): decreasing, limits cross over.
eta_c <- est[1];           se_c <- se[1]
eta_t <- est[1] + est[2];  se_t <- sqrt(V[1, 1] + V[2, 2] + 2 * V[1, 2])
ltef  <- function(eta, se) c(1 - plogis(eta),
                             1 - plogis(eta + 1.96 * se),
                             1 - plogis(eta - 1.96 * se))
lc <- ltef(eta_c, se_c); lt <- ltef(eta_t, se_t)

cure_tab <- data.frame(
  Parameter = c("Long-term event-free fraction, Control",
                "Long-term event-free fraction, Treatment",
                "Odds ratio for long-term event-free status (Treatment vs Control)",
                "Latency log time ratio (Treatment vs Control)",
                "Weibull shape"),
  Estimate = fmt3(c(lc[1], lt[1], exp(-est[2]), est[4], exp(-est[5]))),
  LCL = fmt3(c(lc[2], lt[2], exp(-(est[2] + 1.96 * se[2])),
               est[4] - 1.96 * se[4], exp(-(est[5] + 1.96 * se[5])))),
  UCL = fmt3(c(lc[3], lt[3], exp(-(est[2] - 1.96 * se[2])),
               est[4] + 1.96 * se[4], exp(-(est[5] - 1.96 * se[5])))))
names(cure_tab) <- c("Parameter", "Estimate",
                     "95% CI, Lower Limit", "95% CI, Upper Limit")
write.csv(cure_tab, "outputs/tA6_curemodel.csv", row.names = FALSE)

cat("\n=== 5.6 Mixture cure model (Final) ===\n")
print(cure_tab, row.names = FALSE)
cat(sprintf("log-likelihood: %.6f; convergence code: %d\n",
            -opt$value, opt$convergence))

## =====================================================================
## Figures
## =====================================================================
cols <- c(Control = "grey45", Treatment = "black")

km_panel <- function(cc, xmax, xticks) {
  km <- km_fit[[cc]]; grid <- grids[[cc]]
  plot(km, col = cols, lwd = 2, lty = 1, conf.int = FALSE,
       xlim = c(0, xmax), ylim = c(0, 1), xaxt = "n",
       xlab = "Months since randomization",
       ylab = "Event-free probability", main = paste(cc, "data cut"))
  axis(1, at = xticks)
  lines(km, conf.int = "only", col = cols, lty = 3, lwd = 1)
  abline(h = 0.5, col = "grey70", lty = 2)
  legend("bottomleft",
         c("Control", "Treatment", "95% CI (log-log)", "0.50 level"),
         col = c(cols, "grey45", "grey70"), lty = c(1, 1, 3, 2),
         lwd = c(2, 2, 1, 1), bty = "n", cex = 0.85)
  s <- summary(km, times = grid)
  lab_at <- -0.065 * xmax            # row labels in the left margin,
  mtext("At risk", side = 1, line = 3.9, at = lab_at, adj = 1,
        cex = 0.75, font = 2)        # clear of the t = 0 counts
  for (a in arms) {
    ln <- 3.9 + which(arms == a)
    mtext(a, side = 1, line = ln, at = lab_at, adj = 1, cex = 0.75)
    nr <- s$n.risk[s$strata == paste0("TRT=", a)]
    mtext(sprintf("%d", nr), side = 1, line = ln, at = grid, cex = 0.75)
  }
}
draw_fig1 <- function() {
  par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))
  km_panel("Interim", 16, c(0, 6, 12, 16))
  km_panel("Final",   48, seq(0, 48, 12))
}
png("outputs/fig1_km_interim_final.png", width = 11, height = 5.8,
    units = "in", res = 300); draw_fig1(); dev.off()
pdf("outputs/fig1_km_interim_final.pdf", width = 11, height = 5.8)
draw_fig1(); dev.off()

draw_fig2 <- function() {
  par(mar = c(5, 4, 3, 1))
  km <- km_fit[["Final"]]
  plot(km, col = cols, lwd = 2, lty = 1, conf.int = FALSE,
       xlim = c(0, 48), ylim = c(0, 1), xaxt = "n",
       xlab = "Months since randomization",
       ylab = "Event-free probability",
       main = "Final cut: Kaplan-Meier and fitted mixture cure model")
  axis(1, at = seq(0, 48, 12))
  tt <- seq(0, 48, by = 0.1)
  lines(tt, cure_spop(est, tt, 0), col = cols["Control"],   lty = 2, lwd = 1.6)
  lines(tt, cure_spop(est, tt, 1), col = cols["Treatment"], lty = 2, lwd = 1.6)
  text(48, lc[1] - 0.045, sprintf("LTEF %s", fmt3(lc[1])),
       col = cols["Control"],   adj = 1, cex = 0.85)
  text(48, lt[1] + 0.045, sprintf("LTEF %s", fmt3(lt[1])),
       col = cols["Treatment"], adj = 1, cex = 0.85)
  legend("bottomleft",
         c("Control, KM", "Treatment, KM",
           "Control, model", "Treatment, model"),
         col = cols[c(1, 2, 1, 2)], lty = c(1, 1, 2, 2),
         lwd = c(2, 2, 1.6, 1.6), bty = "n", cex = 0.85)
}
png("outputs/fig2_cure_fit.png", width = 7, height = 5.5,
    units = "in", res = 300); draw_fig2(); dev.off()
pdf("outputs/fig2_cure_fit.pdf", width = 7, height = 5.5)
draw_fig2(); dev.off()

sink("outputs/sessionInfo_analysis.txt"); print(sessionInfo()); sink()
cat("\nAll tables and figures written to outputs/\n")
