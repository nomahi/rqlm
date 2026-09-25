# rqlm package

## Estimation of risk ratios and risk differences for binary outcomes

`rqlm` helps epidemiologists and clinical researchers fit **modified Poisson**, **modified least-squares**, and **modified logistic (augmented logistic)** regression models. It returns effect estimates, robust standard errors, confidence intervals, and a full covariance matrix through familiar R formula syntax.

This repository accompanies a practical tutorial built around two **fully synthetic datasets based on a target trial emulation** comparing SGLT2-inhibitor and DPP-4-inhibitor initiation. The outcome is death within 36 months. The data contain **no actual patient records** and do not provide clinical evidence about either drug.

The worked examples describe **rqlm 4.5-1**. The example scripts and this documentation do not modify the package implementation.

- [CRAN package](https://CRAN.R-project.org/package=rqlm)
- [Practical tutorial](https://doi.org/10.51094/jxiv.1054)
- [Full worked-example R script](https://github.com/nomahi/rqlm/blob/main/Examplecode.r)

## Installation

Install the available CRAN release, then check its version:

```r
install.packages("rqlm")
library(rqlm)
packageVersion("rqlm")
```

The examples require **4.5-1 or later** and the `SGLT2i01` and `SGLT2i02` datasets.


## Which method should I use?

| Goal | Command | Treatment coefficient to report |
| --- | --- | --- |
| Adjusted RR | `rqlm(..., family = poisson, eform = TRUE)` | `exp(coef)` is the RR |
| Adjusted RR | `qlogist(..., eform = TRUE)` | `exp(coef)` is an RR, not an ordinary logistic odds ratio |
| Adjusted RD | `rqlm(..., family = gaussian, eform = FALSE)` | `Estimate` is the RD; multiply by 100 for percentage points |

The interpretations assume an appropriate model without treatment interactions. A regression coefficient is not automatically a causal effect.


## Quick start: adjusted risk ratios

```r
library(rqlm)
data(SGLT2i01, package = "rqlm")
rr_data <- SGLT2i01

f_rr <- death36 ~ A + age + hba1c + egfr + proteinuria +
  prior_heart_failure + prior_stroke + recent_hospitalization

fit_p <- rqlm(
  f_rr, data = rr_data, family = poisson,
  eform = TRUE, var.method = "standard"
)
fit_p

fit_q <- qlogist(
  f_rr, data = rr_data,
  eform = TRUE, var.method = "standard"
)
fit_q
```

Read the `A` row. `A = 1` means SGLT2i and `A = 0` means DPP-4i. Give `qlogist()` the **original data**: it constructs pseudo-observations internally.

The supplied manuscript describes a crude RR near 1.59 and adjusted RRs near 0.91 and 0.90. The direction changes because baseline prognosis differs substantially between treatment groups. The adjusted intervals remain wide and include one. These rounded model values are inherited from an independent reconstruction of the supplied implementation, not yet verified by running the new R script.


## Quick start: adjusted risk differences

```r
data(SGLT2i02, package = "rqlm")
rd_data <- SGLT2i02
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

fit_had <- rqlm(
  f_rd, data = rd_data, family = gaussian,
  eform = FALSE, var.method = "HAD"
)
fit_had
```


## Variance options and version-specific cautions

`rqlm()` and `qlogist()` default to `var.method = "MBN"`. The main tutorial **explicitly uses `"standard"`**. An optional MBN comparison is limited to the independent-person modified Poisson fit.

- `standard`: standard robust covariance; used in the main examples.
- `MBN`, `GST`, `WL`: small-sample correction options; they change uncertainty, not the point estimate.
- `HAD`: finite-sample unbiased covariance under its assumptions, for independent, unweighted Gaussian-identity fits only.
- `HAD0`, `HAD025`: truncated HAD variants; they do **not** preserve exact unbiasedness.

