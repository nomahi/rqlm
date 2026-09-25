# rqlm: practical binary-outcome tutorial
# Author: Hisashi Noma

# Install once if necessary:
# install.packages("rqlm")

# 3.2. Load the synthetic target trial emulation datasets
library(rqlm)
packageVersion("rqlm")

data(SGLT2i01, package = "rqlm")
data(SGLT2i02, package = "rqlm")
rr_data <- SGLT2i01
rd_data <- SGLT2i02

dim(rr_data)
head(rr_data)
table(rr_data$treatment)

# Basic input checks; no data are imputed, recoded, or dropped here.
required <- c(
  "id", "A", "treatment", "death36", "age", "hba1c", "egfr",
  "proteinuria", "prior_heart_failure", "prior_stroke",
  "recent_hospitalization", "risk_score", "followup_months"
)
for (d in list(rr_data, rd_data)) {
  if (!all(required %in% names(d)))
    stop("A required tutorial variable is missing.", call. = FALSE)
  if (anyNA(d[required]))
    stop("The supplied tutorial data should have no missing required values.",
         call. = FALSE)
  if (anyDuplicated(d$id) || !all(d$followup_months == 36))
    stop("The tutorial requires one row per person and complete 36-month outcomes.",
         call. = FALSE)
  if (!all(d$A %in% c(0, 1)) || !all(d$death36 %in% c(0, 1)))
    stop("Treatment and outcome must be numeric 0/1.", call. = FALSE)
  if (!is.numeric(d$risk_score))
    stop("risk_score must be numeric.", call. = FALSE)
}
if (!identical(rr_data$id, rd_data$id) || !identical(rr_data$A, rd_data$A))
  stop("The two teaching versions should have identical IDs and assignments.",
       call. = FALSE)
if (nrow(rr_data) != 700L || nrow(rd_data) != 700L)
  warning("These data differ in size from the 700-person manuscript examples.")
rm(d)

# 4.1. Baseline differences and crude risks
aggregate(
  cbind(prior_heart_failure, recent_hospitalization) ~ treatment,
  data = rr_data, FUN = mean
)
with(rr_data, table(treatment, death36))
with(rr_data, tapply(death36, treatment, mean))

# 4.2. Modified Poisson: adjusted RR
f_rr <- death36 ~ A + age + hba1c + egfr + proteinuria +
  prior_heart_failure + prior_stroke + recent_hospitalization

fit_p <- rqlm(
  f_rr, data = rr_data, family = poisson,
  eform = TRUE, var.method = "standard"
)
fit_p

if (!isTRUE(fit_p$model$converged) || anyNA(coef(fit_p)))
  stop("fit_p did not yield a converged full-rank fit; inspect the output.",
       call. = FALSE)
if (!is.finite(fit_p$se[["A"]]))
  stop("The treatment standard error is unavailable in fit_p.", call. = FALSE)

# 4.3. Modified logistic: adjusted RR
fit_q <- qlogist(
  f_rr, data = rr_data,
  eform = TRUE, var.method = "standard"
)
fit_q

if (!isTRUE(fit_q$model$converged) || anyNA(coef(fit_q)))
  stop("fit_q did not yield a converged full-rank fit; inspect the output.",
       call. = FALSE)
if (!is.finite(fit_q$se[["A"]]))
  stop("The treatment standard error is unavailable in fit_q.", call. = FALSE)

# 4.3. Treatment-only RR table
rr_results <- rbind(
  Modified_Poisson = exp(c(
    coef(fit_p)[["A"]], fit_p$cl[["A"]], fit_p$cu[["A"]]
  )),
  Modified_logistic = exp(c(
    coef(fit_q)[["A"]], fit_q$cl[["A"]], fit_q$cu[["A"]]
  ))
)
colnames(rr_results) <- c("RR", "Lower95", "Upper95")
round(rr_results, 3)

# 5.1. Crude RD in the alternative teaching dataset
with(rd_data, table(treatment, death36))
rd_risks <- with(rd_data, tapply(death36, treatment, mean))
100 * (rd_risks[["SGLT2i"]] - rd_risks[["DPP-4i"]])

# 5.1. Modified least squares: adjusted RD
rd_data$baseline_risk_index <- exp(2.4 * rd_data$risk_score)
f_rd <- death36 ~ A + baseline_risk_index

fit_rd <- rqlm(
  f_rd, data = rd_data, family = gaussian,
  eform = FALSE, var.method = "standard"
)
fit_rd

100 * c(
  RD = coef(fit_rd)[["A"]],
  Lower95 = fit_rd$cl[["A"]],
  Upper95 = fit_rd$cu[["A"]]
)

if (!isTRUE(fit_rd$model$converged) || anyNA(coef(fit_rd)))
  stop("fit_rd did not yield a converged full-rank fit; inspect the output.",
       call. = FALSE)
if (!is.finite(fit_rd$se[["A"]]))
  stop("The treatment standard error is unavailable in fit_rd.", call. = FALSE)

# 5.2. HAD uncertainty estimate
fit_had <- rqlm(
  f_rd, data = rd_data, family = gaussian,
  eform = FALSE, var.method = "HAD"
)
fit_had

if (!isTRUE(fit_had$model$converged) || anyNA(coef(fit_had)))
  stop("fit_had did not yield a converged full-rank fit; inspect the output.",
       call. = FALSE)
if (!is.finite(fit_had$se[["A"]]))
  stop("The treatment standard error is unavailable in fit_had.", call. = FALSE)

# Check the fitted risks for the least-squares example; do not clip them.
rd_fitted <- stats::fitted(fit_rd$model)
if (any(rd_fitted < 0 | rd_fitted > 1))
  warning("Some least-squares fitted risks are outside [0, 1]; review the model.")

# 5.2. Treatment-only RD table, in percentage points
rd_results <- 100 * rbind(
  Standard = c(
    coef(fit_rd)[["A"]], fit_rd$cl[["A"]], fit_rd$cu[["A"]]
  ),
  HAD = c(
    coef(fit_had)[["A"]], fit_had$cl[["A"]], fit_had$cu[["A"]]
  )
)
colnames(rd_results) <- c("RD_pp", "Lower95_pp", "Upper95_pp")
round(rd_results, 2)

# 6. Optional independent-person Poisson MBN comparison
fit_p_mbn <- rqlm(
  f_rr, data = rr_data, family = poisson,
  eform = TRUE, var.method = "MBN"
)
fit_p_mbn

if (!isTRUE(fit_p_mbn$model$converged) || anyNA(coef(fit_p_mbn)))
  stop("fit_p_mbn did not yield a converged full-rank fit; inspect the output.",
       call. = FALSE)
if (!is.finite(fit_p_mbn$se[["A"]]))
  stop("The treatment standard error is unavailable in fit_p_mbn.", call. = FALSE)

# 7.1. Extracting stored model results
exp(coef(fit_p)[["A"]])
fit_p$se[["A"]]
vcov(fit_p)
family(fit_p)

# 7.1. Save results and the R session information
out_dir <- "rqlm_tutorial_results"
dir.create(out_dir, showWarnings = FALSE)
write.csv(rr_results, file.path(out_dir, "risk_ratios.csv"))
write.csv(rd_results, file.path(out_dir, "risk_differences_pp.csv"))
saveRDS(
  list(poisson = fit_p, logistic = fit_q,
       least_squares = fit_rd, had = fit_had, mbn = fit_p_mbn),
  file.path(out_dir, "fitted_models.rds")
)
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo.txt"))

# Changing only the variance estimator must not change point estimates.
stopifnot(
  isTRUE(all.equal(coef(fit_p), coef(fit_p_mbn), tolerance = 1e-8)),
  isTRUE(all.equal(coef(fit_rd), coef(fit_had), tolerance = 1e-8))
)
