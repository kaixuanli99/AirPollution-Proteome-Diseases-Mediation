
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

coerce_numeric <- function(df, cols) {
	for (nm in intersect(cols, colnames(df))) df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
	df
}

get_outcome_pairs <- function(outcome_cols, duration_cols) {
	stopifnot(length(outcome_cols) == length(duration_cols))
	data.frame(outcome = outcome_cols, duration = duration_cols, stringsAsFactors = FALSE)
}

fit_one_cox <- function(df, feature, outcome, duration, covariates, do_zph = FALSE) {
	# Build formula: Surv(duration, outcome) ~ feature + covariates
	# Ensure required columns
	needed <- c(feature, outcome, duration, covariates)
	ensure_columns_exist(df, needed, "cox variables")

	# Clean rows: non-missing in key fields
	rows_ok <- complete.cases(df[, needed, drop = FALSE])
	if (!any(rows_ok)) return(NULL)
	dd <- df[rows_ok, , drop = FALSE]

	# Coerce numeric where appropriate
	dd <- coerce_numeric(dd, c(feature, outcome, duration, covariates))

	# Cox fit
	form <- as.formula(paste0("Surv(", duration, ", ", outcome, ") ~ ", paste(c(feature, covariates), collapse = "+")))
	fit <- try(coxph(form, data = dd), silent = TRUE)
	if (inherits(fit, "try-error")) return(NULL)

	# Extract results for the feature term (first coefficient)
	ci <- try(exp(cbind(coef(fit), confint.default(fit))), silent = TRUE)
	if (inherits(ci, "try-error")) return(NULL)
	ci <- as.data.frame(ci)
	row1 <- ci[1, , drop = FALSE]
	HR <- as.numeric(row1[1, 1])
	LCI <- as.numeric(row1[1, 2])
	UCI <- as.numeric(row1[1, 3])
	pval <- try(summary(fit)$coefficients[1, 5], silent = TRUE)
	if (inherits(pval, "try-error")) pval <- NA_real_
	nevent <- try(fit$nevent[1], silent = TRUE)
	if (inherits(nevent, "try-error")) nevent <- NA_real_
	n_total <- try(fit$n[1], silent = TRUE)
	if (inherits(n_total, "try-error")) n_total <- sum(rows_ok)
	n_controls <- n_total - nevent

	p_zph_feature <- NA_real_
	p_zph_global <- NA_real_
	if (isTRUE(do_zph)) {
		z <- try(cox.zph(fit), silent = TRUE)
		if (!inherits(z, "try-error")) {
			pvals <- try(z$table[, "p"], silent = TRUE)
			if (!inherits(pvals, "try-error")) {
				p_zph_feature <- suppressWarnings(as.numeric(pvals[1]))
				p_zph_global <- suppressWarnings(as.numeric(pvals["GLOBAL"]))
			}
		}
	}

	data.frame(
		Predictor = feature,
		Outcome = outcome,
		Hazard.Ratio = round(HR, 2),
		LCI = round(LCI, 2),
		UCI = round(UCI, 2),
		P.Value = pval,
		Cases = nevent,
		Controls = n_controls,
		cox.zph_feature = p_zph_feature,
		cox.zph_global = p_zph_global,
		stringsAsFactors = FALSE
	)
}

run_cox_for_outcome <- function(df, features, outcome, duration, covariates, do_zph = FALSE) {
	res_list <- vector("list", length(features))
	for (j in seq_along(features)) {
		ft <- features[j]
		r <- fit_one_cox(df, feature = ft, outcome = outcome, duration = duration, covariates = covariates, do_zph = do_zph)
		res_list[[j]] <- r
	}
	res <- dplyr::bind_rows(res_list)
	res <- res[order(res$P.Value), , drop = FALSE]
	res
}

# ------------------------------
# Public API
# ------------------------------

#' Run Cox regressions for many proteins against multiple disease outcomes
#'
#' @param data Either a data.frame with all columns or a path to CSV/TSV
#' @param feature_cols Character vector of protein column names to test as predictors
#' @param outcome_cols Character vector of binary outcome columns (same length as duration_cols)
#' @param duration_cols Character vector of time-to-event columns (same length as outcome_cols)
#' @param covariates Character vector of covariate column names
#' @param do_zph Logical; whether to compute cox.zph tests (default FALSE)
#' @param write_dir Optional directory to write per-outcome CSVs; if NULL, nothing is written
#' @param basename Base name for output files when write_dir is provided
#' @return A named list of data.frames: one result per outcome
run_cox_proteome_disease_public <- function(
	data,
	feature_cols,
	outcome_cols,
	duration_cols,
	covariates = c("age", "sex", "smoking", "alcohol", "BMI", "ethnic"),
	do_zph = FALSE,
	write_dir = NULL,
	basename = "Proteome_Disease_Cox"
) {
	df <- if (is.character(data)) fread(data) %>% as.data.frame() else data
	ensure_columns_exist(df, c(feature_cols, outcome_cols, duration_cols, covariates), "input data")

	pairs <- get_outcome_pairs(outcome_cols, duration_cols)
	results <- list()
	if (!is.null(write_dir) && !dir.exists(write_dir)) dir.create(write_dir, recursive = TRUE, showWarnings = FALSE)

	for (i in seq_len(nrow(pairs))) {
		outc <- pairs$outcome[i]
		dur <- pairs$duration[i]
		cat("Running:", outc, "with", length(feature_cols), "proteins\n")
		res <- run_cox_for_outcome(df, features = feature_cols, outcome = outc, duration = dur, covariates = covariates, do_zph = do_zph)
		results[[outc]] <- res
		if (!is.null(write_dir)) {
			fwrite(res, file = file.path(write_dir, paste0(basename, "_", outc, ".csv")))
		}
	}
	results
}

# Example (commented):
# res_list <- run_cox_proteome_disease_public(
# 	data = "./inputs/UKB_all_covariates_proteome.csv",
# 	feature_cols = paste0("prot_", 1:1463),
# 	outcome_cols = c("DiseaseA_2010_2020", "DiseaseB_2010_2020"),
# 	duration_cols = c("date_DiseaseA_2010_2020", "date_DiseaseB_2010_2020"),
# 	covariates = c("age", "sex", "smoking", "alcohol", "BMI", "ethnic"),
# 	write_dir = "./outputs",
# 	basename = "Proteome_Disease_Cox"
# )


