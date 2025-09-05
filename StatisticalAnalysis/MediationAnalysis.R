
suppressPackageStartupMessages({
	library(data.table)
	library(dplyr)
	library(mediation)
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

coerce_binary <- function(x) {
	if (is.factor(x)) x <- as.character(x)
	x <- suppressWarnings(as.numeric(x))
	# If values are not 0/1, try to map unique values to 0/1
	ux <- unique(x[!is.na(x)])
	if (!all(ux %in% c(0, 1))) {
		if (length(ux) == 2) {
			m <- setNames(c(0, 1), sort(ux))
			x <- unname(m[as.character(x)])
		}
	}
	x
}

# Significance filtering helper
filter_significant <- function(df, p_col = "P.Value", adjust_method = NULL, alpha = 0.05, out_col_name = NULL) {
	stopifnot(p_col %in% colnames(df))
	dd <- df
	if (!is.null(adjust_method)) {
		padj <- p.adjust(dd[[p_col]], method = adjust_method)
		keep <- padj < alpha
	} else {
		keep <- dd[[p_col]] < alpha
	}
	dd[keep, , drop = FALSE]
}

# Compute overlapping mediators between PWAS (AP~Proteome) and Cox (Proteome~Disease)
compute_overlap_mediators <- function(
	pwas_results,
	cox_results,
	exposure,
	disease,
	pwas_exposure_col = "exposure",
	pwas_protein_col = "feature",
	pwas_p_col = "P.Value",
	pwas_adjust = "fdr",
	pwas_alpha = 0.05,
	cox_outcome_col = "Outcome",
	cox_protein_col = "Predictor",
	cox_p_col = "P.Value",
	cox_adjust = "bonferroni",
	cox_alpha = 0.05
) {
	stopifnot(pwas_exposure_col %in% colnames(pwas_results))
	stopifnot(cox_outcome_col %in% colnames(cox_results))
	# Filter PWAS for specific exposure and significance
	pwas_sub <- pwas_results[pwas_results[[pwas_exposure_col]] == exposure, , drop = FALSE]
	pwas_sig <- filter_significant(pwas_sub, p_col = pwas_p_col, adjust_method = pwas_adjust, alpha = pwas_alpha)
	proteins_pwas <- unique(as.character(pwas_sig[[pwas_protein_col]]))
	# Filter Cox for specific disease and significance
	cox_sub <- cox_results[cox_results[[cox_outcome_col]] == disease, , drop = FALSE]
	cox_sig <- filter_significant(cox_sub, p_col = cox_p_col, adjust_method = cox_adjust, alpha = cox_alpha)
	proteins_cox <- unique(as.character(cox_sig[[cox_protein_col]]))
	intersect(proteins_pwas, proteins_cox)
}

# ------------------------------
# Single mediation run and tidy summary
# ------------------------------

fit_mediation_single <- function(
	data,
	mediator,
	exposure,
	outcome,
	covariates = c("age", "sex", "smoking", "BMI", "alcohol"),
	sims = 100,
	outcome_family = c("poisson", "binomial"),
	offset_col = NULL
) {
	outcome_family <- match.arg(outcome_family)
	needed <- c(mediator, exposure, outcome, covariates)
	if (!is.null(offset_col)) needed <- c(needed, offset_col)
	ensure_columns_exist(data, needed, "mediation variables")

	dd <- data[, needed, drop = FALSE]
	dd <- na.omit(dd)
	if (!nrow(dd)) return(NULL)

	# Coerce types
	dd[[mediator]] <- suppressWarnings(as.numeric(dd[[mediator]]))
	dd[[exposure]] <- suppressWarnings(as.numeric(dd[[exposure]]))
	if (outcome_family == "binomial") dd[[outcome]] <- coerce_binary(dd[[outcome]])

	# Build models
	med_fml <- stats::as.formula(paste0(mediator, " ~ ", exposure, if (length(covariates)) paste0(" + ", paste(covariates, collapse = " + ")) else ""))
	out_fml <- stats::as.formula(paste0(outcome, " ~ ", mediator, " + ", exposure, if (length(covariates)) paste0(" + ", paste(covariates, collapse = " + ")) else ""))

	med.fit <- lm(med_fml, data = dd)
	if (outcome_family == "poisson") {
		if (!is.null(offset_col)) {
			dd[[offset_col]] <- suppressWarnings(as.numeric(dd[[offset_col]]))
			dd <- dd[is.finite(dd[[offset_col]]) & dd[[offset_col]] > 0, , drop = FALSE]
			if (!nrow(dd)) return(NULL)
			out.fit <- glm(out_fml, data = dd, family = "poisson", offset = log(dd[[offset_col]]))
		} else {
			out.fit <- glm(out_fml, data = dd, family = "poisson")
		}
	} else {
		out.fit <- glm(out_fml, data = dd, family = binomial())
	}

	md <- mediate(med.fit, out.fit, treat = exposure, mediator = mediator, sims = sims, boot = TRUE, boot.ci.type = "bca")
	sm <- summary(md)

	# Tidy like original structure
	res <- data.frame(
		Term = c("ACME1", "ACME2", "ADE1", "ADE2", "Prop_Mediated1", "Prop_Mediated2", "ACME_average", "ADE_average", "Total_Effect", "Prop_Mediated_average"),
		Estimate = c(sm$d0, sm$d1, sm$z0, sm$z1, sm$n0, sm$n1, (sm$d0 + sm$d1) / 2, (sm$z0 + sm$z1) / 2, sm$tau.coef, (sm$n0 + sm$n1) / 2),
		P_value = suppressWarnings(as.numeric(c(sm$d0.p, sm$d1.p, sm$z0.p, sm$z1.p, sm$n0.p, sm$n1.p, (as.numeric(sm$d0.p) + as.numeric(sm$d1.p)) / 2, (as.numeric(sm$z0.p) + as.numeric(sm$z1.p)) / 2, sm$tau.p, (as.numeric(sm$n0.p) + as.numeric(sm$n1.p)) / 2))),
		stringsAsFactors = FALSE
	)
	attr(res, "meta") <- list(N = nrow(dd), mediator = mediator, exposure = exposure, outcome = outcome)
	res
}

# ------------------------------
# Batch runner and aggregation
# ------------------------------

run_mediation_public <- function(
	data,
	exposures,
	outcomes,
	mediators,
	covariates = c("age", "sex", "smoking", "BMI", "alcohol"),
	sims = 100,
	outcome_family = c("poisson", "binomial"),
	offset_col = NULL,
	# Optional: provide prior results to derive mediator overlaps per (exposure, outcome)
	pwas_results = NULL,
	cox_results = NULL,
	use_overlap = FALSE,
	pwas_exposure_col = "exposure",
	pwas_protein_col = "feature",
	pwas_p_col = "P.Value",
	pwas_adjust = "fdr",
	pwas_alpha = 0.05,
	cox_outcome_col = "Outcome",
	cox_protein_col = "Predictor",
	cox_p_col = "P.Value",
	cox_adjust = "bonferroni",
	cox_alpha = 0.05,
	write_dir = NULL,
	basename = "Mediation"
) {
	dd <- if (is.character(data)) fread(data) %>% as.data.frame() else data
	outcome_family <- match.arg(outcome_family)

	all_rows <- list()
	for (y in outcomes) {
		for (x in exposures) {
			# Determine mediator set: either user-specified or overlap-derived
			mediators_set <- mediators
			if (isTRUE(use_overlap)) {
				stopifnot(!is.null(pwas_results), !is.null(cox_results))
				mediators_overlap <- compute_overlap_mediators(
					pwas_results = if (is.character(pwas_results)) fread(pwas_results) %>% as.data.frame() else pwas_results,
					cox_results = if (is.character(cox_results)) fread(cox_results) %>% as.data.frame() else cox_results,
					exposure = x,
					disease = y,
					pwas_exposure_col = pwas_exposure_col,
					pwas_protein_col = pwas_protein_col,
					pwas_p_col = pwas_p_col,
					pwas_adjust = pwas_adjust,
					pwas_alpha = pwas_alpha,
					cox_outcome_col = cox_outcome_col,
					cox_protein_col = cox_protein_col,
					cox_p_col = cox_p_col,
					cox_adjust = cox_adjust,
					cox_alpha = cox_alpha
				)
				if (!is.null(mediators)) {
					mediators_set <- intersect(mediators, mediators_overlap)
				} else {
					mediators_set <- mediators_overlap
				}
			}
			if (length(mediators_set) == 0) next
			for (m in mediators_set) {
				res <- try(fit_mediation_single(dd, mediator = m, exposure = x, outcome = y, covariates = covariates, sims = sims, outcome_family = outcome_family, offset_col = offset_col), silent = TRUE)
				if (inherits(res, "try-error") || is.null(res)) next
				meta <- attr(res, "meta")
				res$Mediator <- meta$mediator
				res$Exposure <- meta$exposure
				res$Outcome <- meta$outcome
				res$N <- meta$N
				all_rows[[length(all_rows) + 1]] <- res
			}
		}
	}

	combined <- if (length(all_rows)) bind_rows(all_rows) else data.frame()
	if (!is.null(write_dir)) {
		if (!dir.exists(write_dir)) dir.create(write_dir, recursive = TRUE, showWarnings = FALSE)
		fwrite(combined, file = file.path(write_dir, paste0(basename, "_long.csv")))
	}
	combined
}

summarise_mediation_public <- function(
	results_long,
	adjust_method = "fdr",
	alpha = 0.05,
	min_prop = 0
) {
	if (!nrow(results_long)) return(results_long)
	# Pivot to a compact wide-like summary per (Mediator, Exposure, Outcome)
	key_rows <- results_long %>%
		filter(Term %in% c("ACME_average", "ADE_average", "Total_Effect", "Prop_Mediated_average")) %>%
		select(Mediator, Exposure, Outcome, Term, Estimate, P_value, N)
	wide <- tidyr::pivot_wider(key_rows, names_from = Term, values_from = c(Estimate, P_value))
	colnames(wide) <- gsub("Estimate_", "", colnames(wide))
	colnames(wide) <- gsub("P_value_", "P.", colnames(wide))
	# FDR adjustments
	wide$ACME_FDR <- p.adjust(wide$`P.ACME_average`, method = adjust_method)
	wide$ADE_FDR <- p.adjust(wide$`P.ADE_average`, method = adjust_method)
	wide$Total_FDR <- p.adjust(wide$`P.Total_Effect`, method = adjust_method)
	wide$Prop_FDR <- p.adjust(wide$`P.Prop_Mediated_average`, method = adjust_method)
	# Basic filtering
	filtered <- wide %>% filter(ACME_FDR < alpha, Total_FDR < alpha, `Prop_Mediated_average` > min_prop)
	filtered
}

# ------------------------------
# Example (commented):
# ------------------------------

# res_long <- run_mediation_public(
# 	data = "./inputs/ukb_ap_proteome_disease.csv",  # contains columns for mediators, exposures, outcomes, covariates
# 	exposures = c("Nitrogen_oxides_2010", "NO2_mean", "PM25_2010", "PM10_mean"),
# 	outcomes = c("Asthma_2010_2020", "COPD_2010_2020"),
# 	mediators = c("GDF15", "ITGAV"),
# 	covariates = c("age", "sex", "smoking", "BMI", "alcohol"),
# 	sims = 100, outcome_family = "poisson",
# 	write_dir = "./outputs",
# 	basename = "Mediation_AP_Proteome_Disease"
# )
#
# res_summary <- summarise_mediation_public(res_long, adjust_method = "fdr", alpha = 0.05, min_prop = 0.04)
# data.table::fwrite(res_summary, file = "./outputs/mediation_summary.csv")


