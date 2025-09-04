## Public, cleaned Cox analysis: Air pollution ~ Disease outcomes
## - Removes AP index construction and any correlation plots
## - No hard-coded paths; inputs/outputs are parameterized

suppressPackageStartupMessages({
	library(data.table)
	library(dplyr)
	library(survival)
})

# ------------------------------
# Helpers
# ------------------------------

IQR_transfer <- function(dataset, x) {
	for (z in x) {
		if (!z %in% colnames(dataset)) next
		new_x_IQR <- paste0(z, "_", "IQR")
		dataset[[new_x_IQR]] <- as.numeric(dataset[[z]]) / stats::IQR(as.numeric(dataset[[z]]), na.rm = TRUE)
	}
	dataset
}

cox_IQR_lkx <- function(dataset, outcome, duration, covariates, exposure) {
	vars <- c(exposure, covariates)
	mySurv <- Surv(time = dataset[[duration]], event = dataset[[outcome]])
	FML <- as.formula(paste0("mySurv~", paste(vars, collapse = "+")))
	fit_cox <- coxph(FML, data = dataset)
	res_cox <- exp(cbind(coef(fit_cox), confint.default(fit_cox))) %>% as.data.frame()
	res_cox[, 4] <- paste0(round(res_cox[, 1], 2), '(', round(res_cox[, 2], 2), '~', round(res_cox[, 3], 2), ')')
	sm <- summary(fit_cox)
	b <- as.data.frame(sm$coefficients)
	res_cox <- cbind(res_cox, b[, 5])
	colnames(res_cox) <- c("HR", "Lowerlimit", "Upperlimit", "HR_95%CI", "P_value")
	res_cox$exposure <- exposure
	res_cox[1, , drop = FALSE]
}

derive_pollutant_means <- function(pheno) {
	pheno %>%
		mutate(
			NO2_2005 = as.numeric(NO2_2005),
			NO2_2006 = as.numeric(NO2_2006),
			NO2_2007 = as.numeric(NO2_2007),
			NO2_2010 = as.numeric(NO2_2010),
			PM10_2007 = as.numeric(PM10_2007),
			PM10_2010 = as.numeric(PM10_2010)
		) %>%
		mutate(
			NO2_mean = (NO2_2005 + NO2_2006 + NO2_2007 + NO2_2010) / 4,
			PM10_mean = (PM10_2007 + PM10_2010) / 2
		)
}

get_disease_pairs <- function(disease_table, baseline_date_col = "53-0.0") {
	cols <- colnames(disease_table)
	status_cols <- setdiff(cols[!startsWith(cols, "date_")], c("eid", baseline_date_col))
	status_cols <- status_cols[status_cols != ""]
	status_cols <- status_cols[order(status_cols)]
	data.frame(
		status = status_cols,
		date = paste0("date_", status_cols),
		stringsAsFactors = FALSE
	)
}

# ------------------------------
# Public API
# ------------------------------

#' Run Cox regression: Air pollution exposures ~ incident diseases
#'
#' @param pheno Either a data.frame of phenotypes (with exposures and covariates) or a file path
#' @param disease_table Either a data.frame of disease status/date columns or a file path
#' @param start_date String date for left-truncation of disease onset (e.g., '2010-01-01')
#' @param exposures Character vector of exposure column names (without _IQR suffix)
#' @param covariates Character vector of covariate column names
#' @param baseline_date_col Column name for baseline date present in disease_table (default '53-0.0')
#' @param write_path Optional CSV output file; if NULL, no output is written
#' @return data.frame of Cox results with HR, CI, p-value, disease, exposure
run_cox_airpollution_disease_public <- function(
	pheno,
	disease_table,
	start_date = "2010-01-01",
	exposures = c("PM25_2010", "Nitrogen_oxides_2010", "NO2_mean", "PM10_mean"),
	covariates = c("Age", "Sex", "Smoking_status", "Drinking_freq", "BMI", "Ethnic"),
	baseline_date_col = "53-0.0",
	write_path = NULL
) {
	# Load if paths provided
	pheno_df <- if (is.character(pheno)) fread(pheno) %>% as.data.frame() else pheno
	disease_df <- if (is.character(disease_table)) fread(disease_table) %>% as.data.frame() else disease_table
	stopifnot("eid" %in% colnames(pheno_df), "eid" %in% colnames(disease_df))

	# Prepare exposures (add NO2_mean, PM10_mean if components exist)
	pheno_df <- derive_pollutant_means(pheno_df)

	# Merge
	ph <- merge(pheno_df, disease_df, by = "eid", all.x = TRUE) %>% distinct()

	# IQR scale exposures
	ph <- IQR_transfer(ph, exposures)
	exposures_iqr <- paste0(exposures, "_IQR")

	# Disease pairs (status/date)
	pairs <- get_disease_pairs(disease_df, baseline_date_col = baseline_date_col)

	# Prepare baseline date
	if (baseline_date_col %in% colnames(ph)) ph[[baseline_date_col]] <- as.Date(ph[[baseline_date_col]])

	results <- list()
	for (k in seq_len(nrow(pairs))) {
		status_col <- pairs$status[k]
		date_col <- pairs$date[k]
		if (!all(c(status_col, date_col) %in% colnames(ph))) next
		# Duration since start_date
		ph[[date_col]] <- as.Date(ph[[date_col]]) - as.Date(start_date)
		tmp <- ph[ph[[date_col]] > 0, , drop = FALSE]
		if (!nrow(tmp)) next
		tmp[[status_col]] <- as.numeric(tmp[[status_col]])
		tmp[[date_col]] <- tmp[[date_col]] / 365.25

		for (ex in exposures_iqr) {
			if (!ex %in% colnames(tmp)) next
			# Drop rows with missing exposure
			tmp2 <- tmp[!is.na(tmp[[ex]]), , drop = FALSE]
			if (!nrow(tmp2)) next
			res <- try(cox_IQR_lkx(tmp2, outcome = status_col, duration = date_col, covariates = covariates, exposure = ex), silent = TRUE)
			if (inherits(res, "try-error")) next
			res$Disease <- status_col
			results[[length(results) + 1]] <- res
		}
	}

	res_table <- if (length(results)) bind_rows(results) else data.frame()
	if (nrow(res_table)) {
		res_table$Disease <- gsub('^', '', gsub('date_', '', res_table$Disease))
	}

	if (!is.null(write_path)) fwrite(res_table, file = write_path)
	res_table
}

# Example (commented):
# res <- run_cox_airpollution_disease_public(
# 	pheno = "./inputs/UKB_pheno.csv",
# 	disease_table = "./inputs/UKB_DiseaseTable.csv",
# 	start_date = "2010-01-01",
# 	write_path = "./outputs/cox_airpollution_disease.csv"
# )


