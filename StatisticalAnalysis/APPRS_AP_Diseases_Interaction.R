
suppressPackageStartupMessages({
	library(data.table)
	library(dplyr)
	library(survival)
})

# ------------------------------
# Helpers
# ------------------------------

ensure_columns_exist <- function(dataset, required_columns, label) {
	missing_columns <- setdiff(required_columns, colnames(dataset))
	if (length(missing_columns) > 0) {
		stop(sprintf("Missing required %s columns: %s", label, paste(missing_columns, collapse = ", ")))
	}
}

quantile_group <- function(x, n, na.rm = TRUE) {
	x <- as.numeric(x)
	if (na.rm) x_na <- is.na(x) else x_na <- rep(FALSE, length(x))
	probs <- seq(0, 1, length.out = n + 1)
	qs <- unique(quantile(x[!is.na(x)], probs = probs, na.rm = TRUE, type = 7))
	# If too many ties collapse, fallback to rank-based grouping
	if (length(qs) <= 2 && n > 2) {
		r <- rank(x, na.last = "keep", ties.method = "average")
		breaks <- quantile(r[!is.na(r)], probs = probs, na.rm = TRUE, type = 7)
		grp <- cut(r, breaks = breaks, include.lowest = TRUE, labels = FALSE)
	} else {
		grp <- cut(x, breaks = qs, include.lowest = TRUE, labels = FALSE)
	}
	grp[is.na(grp)] <- NA_integer_
	grp
}

make_combined_group <- function(apprs_grp, exposure_grp, apprs_labels, exposure_labels) {
	apprs_fac <- factor(apprs_grp, levels = sort(unique(apprs_grp)), labels = apprs_labels[seq_len(nlevels(factor(apprs_grp)))])
	exp_fac <- factor(exposure_grp, levels = sort(unique(exposure_grp)), labels = exposure_labels[seq_len(nlevels(factor(exposure_grp)))])
	interaction(apprs_fac, exp_fac, sep = " & ")
}

fit_cox_groups <- function(df, outcome, duration, group_col, covariates, ref_level = NULL) {
	ensure_columns_exist(df, c(outcome, duration, group_col, covariates), "cox vars")
	dd <- df[, c(outcome, duration, group_col, covariates), drop = FALSE]
	dd <- dd[complete.cases(dd), , drop = FALSE]
	if (!nrow(dd)) return(NULL)
	dd[[outcome]] <- as.numeric(dd[[outcome]])
	dd[[duration]] <- as.numeric(dd[[duration]])
	dd[[group_col]] <- as.factor(dd[[group_col]])
	if (!is.null(ref_level) && ref_level %in% levels(dd[[group_col]])) {
		dd[[group_col]] <- relevel(dd[[group_col]], ref = ref_level)
	}
	mySurv <- Surv(time = dd[[duration]], event = dd[[outcome]])
	form <- as.formula(paste0("mySurv ~ ", group_col, if (length(covariates)) paste0(" + ", paste(covariates, collapse = " + ")) else ""))
	fit <- coxph(form, data = dd)
	ci <- exp(cbind(coef(fit), confint.default(fit))) %>% as.data.frame()
	sm <- summary(fit)
	pvals <- as.data.frame(sm$coefficients)[, 5, drop = FALSE]
	res <- cbind(ci, pvals)
	colnames(res) <- c("HR", "Lower", "Upper", "P_value")
	res$Group <- rownames(res)
	res
}

# ------------------------------
# Public API
# ------------------------------

#' Group APPRS and pollutant, construct combined groups, run Cox for a single outcome
#'
#' @param data data.frame or path; must include columns: exposure, apprs, outcome, duration, covariates
#' @param exposure single pollutant column
#' @param apprs single APPRS column
#' @param outcome binary status (0/1)
#' @param duration time-to-event (in years)
#' @param covariates vector of covariate names
#' @param exposure_quantiles number of quantile groups for exposure (default 3)
#' @param apprs_quantiles number of quantile groups for APPRS (default 2)
#' @param write_path optional CSV output
#' @return data.frame of HR, CI, P and group labels
run_cox_apprs_exposure_single_public <- function(
	data,
	exposure,
	apprs,
	outcome,
	duration,
	covariates = c("sex", "age", "smoking", "alcohol", "BMI"),
	exposure_quantiles = 3,
	apprs_quantiles = 2,
	write_path = NULL
) {
	df <- if (is.character(data)) fread(data) %>% as.data.frame() else data
	ensure_columns_exist(df, c(exposure, apprs, outcome, duration, covariates), "input")

	df$apprs_grp <- quantile_group(df[[apprs]], apprs_quantiles)
	df$exp_grp <- quantile_group(df[[exposure]], exposure_quantiles)

	apprs_labels <- if (apprs_quantiles == 2) c("Low APPRS", "High APPRS") else paste("APPRS Q", seq_len(apprs_quantiles))
	exposure_labels <- c("Low", if (exposure_quantiles >= 3) "Intermediate" else NULL, "High")[seq_len(exposure_quantiles)]

	df$group <- make_combined_group(df$apprs_grp, df$exp_grp, apprs_labels, exposure_labels)
	ref <- paste(apprs_labels[1], exposure_labels[1], sep = " & ")

	res <- fit_cox_groups(df, outcome = outcome, duration = duration, group_col = "group", covariates = covariates, ref_level = ref)
	if (!is.null(res)) {
		res$Outcome <- outcome
		res$Exposure <- exposure
		res$APPRS <- apprs
	}
	if (!is.null(write_path) && !is.null(res)) fwrite(res, file = write_path)
	res
}

#' Batch runner over multiple exposures and outcomes
run_cox_apprs_exposure_batch_public <- function(
	data,
	exposures,
	apprs,
	outcomes,
	durations,
	covariates = c("sex", "age", "smoking", "alcohol", "BMI"),
	exposure_quantiles = 3,
	apprs_quantiles = 2,
	write_dir = NULL,
	basename = "APPRS_Exposure_Cox"
) {
	df <- if (is.character(data)) fread(data) %>% as.data.frame() else data
	stopifnot(length(outcomes) == length(durations))
	all_res <- list()
	for (ex in exposures) {
		for (i in seq_along(outcomes)) {
			outc <- outcomes[i]
			dur <- durations[i]
			r <- run_cox_apprs_exposure_single_public(
				data = df,
				exposure = ex,
				apprs = apprs,
				outcome = outc,
				duration = dur,
				covariates = covariates,
				exposure_quantiles = exposure_quantiles,
				apprs_quantiles = apprs_quantiles
			)
			if (!is.null(r)) {
				all_res[[length(all_res) + 1]] <- r
				if (!is.null(write_dir)) {
					if (!dir.exists(write_dir)) dir.create(write_dir, recursive = TRUE, showWarnings = FALSE)
					fwrite(r, file = file.path(write_dir, paste0(basename, "_", ex, "_", outc, ".csv")))
				}
			}
		}
	}
	if (length(all_res)) dplyr::bind_rows(all_res) else data.frame()
}

# ------------------------------
# Example (commented)
# ------------------------------

# res <- run_cox_apprs_exposure_batch_public(
# 	data = "./inputs/ukb_with_apprs_exposures_outcomes.csv",
# 	exposures = c("PM25_2010", "NO_2010", "NO2_average", "PM10_average", "ap_score_healthy_stat"),
# 	apprs = "ApPRS_HeartFailure",
# 	outcomes = c("HeartFailure", "CKD"),
# 	durations = c("duration_HeartFailure", "duration_CKD"),
# 	covariates = c("sex", "age", "smoking", "alcohol", "BMI"),
# 	exposure_quantiles = 3,
# 	apprs_quantiles = 2,
# 	write_dir = "./outputs",
# 	basename = "APPRS_Exposure_Cox"
# )


