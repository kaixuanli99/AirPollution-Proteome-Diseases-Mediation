
suppressPackageStartupMessages({
	library(data.table)
	library(dplyr)
	library(limma)
})

# ------------------------------
# Helpers
# ------------------------------

# Ensure columns exist
ensure_columns_exist <- function(dataset, required_columns, label) {
	missing_columns <- setdiff(required_columns, colnames(dataset))
	if (length(missing_columns) > 0) {
		stop(sprintf("Missing required %s columns: %s", label, paste(missing_columns, collapse = ", ")))
	}
}

# Prepare pollutant summary variables if components exist
derive_pollutant_means <- function(pheno) {
	pheno <- pheno %>% mutate(
		NO2_2005 = suppressWarnings(as.numeric(NO2_2005)),
		NO2_2006 = suppressWarnings(as.numeric(NO2_2006)),
		NO2_2007 = suppressWarnings(as.numeric(NO2_2007)),
		NO2_2010 = suppressWarnings(as.numeric(NO2_2010)),
		PM10_2007 = suppressWarnings(as.numeric(PM10_2007)),
		PM10_2010 = suppressWarnings(as.numeric(PM10_2010))
	)
	if (all(c("NO2_2005","NO2_2006","NO2_2007","NO2_2010") %in% colnames(pheno))) {
		pheno$NO2_mean <- rowMeans(pheno[, c("NO2_2005","NO2_2006","NO2_2007","NO2_2010")], na.rm = TRUE)
	}
	if (all(c("PM10_2007","PM10_2010") %in% colnames(pheno))) {
		pheno$PM10_mean <- rowMeans(pheno[, c("PM10_2007","PM10_2010")], na.rm = TRUE)
	}
	pheno
}

# Convert proteome to ExpressionSet-like matrix: features x samples
# Supports input where rows are samples and a column 'eid' identifies samples
proteome_to_matrix <- function(proteome_df, sample_id_col = "eid") {
	ensure_columns_exist(proteome_df, c(sample_id_col), "proteome")
	sample_ids <- as.character(proteome_df[[sample_id_col]])
	E <- as.matrix(proteome_df[, setdiff(colnames(proteome_df), sample_id_col), drop = FALSE])
	E <- t(E)  # features x samples
	colnames(E) <- sample_ids
	E
}

# Align proteome matrix (features x samples) and trait dataframe by common sample IDs
align_by_samples <- function(E, trait_df, sample_id_col = "eid") {
	ids_trait <- if (sample_id_col %in% colnames(trait_df)) as.character(trait_df[[sample_id_col]]) else rownames(trait_df)
	stopifnot(!is.null(ids_trait))
	common <- intersect(colnames(E), ids_trait)
	if (length(common) == 0) stop("No overlapping sample IDs between proteome and phenotype.")
	E2 <- E[, common, drop = FALSE]
	trait2 <- trait_df[match(common, ids_trait), , drop = FALSE]
	list(E = E2, trait = trait2)
}

# Build design matrix
build_design <- function(trait_df, exposure, covariates = character(), adjust_exposures = character()) {
	vars <- c(exposure, covariates, adjust_exposures)
	ensure_columns_exist(trait_df, vars, "trait/design")
	stats::model.matrix(
		stats::as.formula(paste0("~ ", paste(c(exposure, covariates, adjust_exposures), collapse = " + "))),
		data = trait_df
	)
}

# Run limma for a single exposure
run_limma_single <- function(E, trait_df, exposure, covariates = character(), adjust_exposures = character()) {
	design <- build_design(trait_df, exposure, covariates, adjust_exposures)
	fit <- lmFit(E, design)
	fit2 <- eBayes(fit)
	res <- topTable(fit2, num = Inf, coef = 2) # coef 2 corresponds to the first predictor (exposure)
	res$feature <- rownames(res)
	res$exposure <- exposure
	res
}

# ------------------------------
# Public API
# ------------------------------

#' Run PWAS (limma) for multiple air pollution exposures vs proteome
#'
#' @param proteome Either a data.frame with samples in rows (must contain 'eid') and proteins in columns, or a file path (CSV/TSV)
#' @param pheno Either a data.frame (must contain 'eid' and covariates/exposures) or a file path (CSV/TSV)
#' @param exposures Character vector of exposure variable names (in pheno) to test
#' @param covariates Character vector of covariate names (in pheno)
#' @param adjust_map Optional named list mapping exposure -> character vector of additional exposures to adjust for
#' @param write_dir Optional directory to write per-exposure results and a combined CSV; if NULL, no files are written
#' @param basename Base name for output files when write_dir is provided
#' @return data.frame with per-feature results across exposures
run_pwas_airpollution_public <- function(
	proteome,
	pheno,
	exposures = c("Nitrogen_oxides_2010", "NO2_mean", "PM25_2010", "PM10_mean"),
	covariates = c("Age", "Sex", "BMI", "Smoking_status", "Drinking_status", "Ethnic"),
	adjust_map = NULL,
	write_dir = NULL,
	basename = "PWAS_AirPollution"
) {
	prot_df <- if (is.character(proteome)) fread(proteome) %>% as.data.frame() else proteome
	pheno_df <- if (is.character(pheno)) fread(pheno) %>% as.data.frame() else pheno
	ensure_columns_exist(prot_df, c("eid"), "proteome")
	ensure_columns_exist(pheno_df, c("eid"), "pheno")

	# Derive pollutant means if possible
	pheno_df <- derive_pollutant_means(pheno_df)

	# Convert proteome and align
	E <- proteome_to_matrix(prot_df, sample_id_col = "eid")
	aligned <- align_by_samples(E, pheno_df, sample_id_col = "eid")
	E2 <- aligned$E
	trait <- aligned$trait

	# Iterate exposures
	all_res <- list()
	for (ex in exposures) {
		if (!ex %in% colnames(trait)) next
		adj <- if (!is.null(adjust_map) && ex %in% names(adjust_map)) adjust_map[[ex]] else character()
		adj <- intersect(adj, colnames(trait))
		res <- run_limma_single(E2, trait, exposure = ex, covariates = covariates, adjust_exposures = adj)
		all_res[[length(all_res) + 1]] <- res
		if (!is.null(write_dir)) {
			if (!dir.exists(write_dir)) dir.create(write_dir, recursive = TRUE, showWarnings = FALSE)
			fwrite(res, file = file.path(write_dir, paste0(basename, "_", ex, ".csv")))
		}
	}

	combined <- if (length(all_res)) bind_rows(all_res) else data.frame()
	if (!is.null(write_dir)) fwrite(combined, file = file.path(write_dir, paste0(basename, "_combined.csv")))
	combined
}

# Example (commented):
# res <- run_pwas_airpollution_public(
# 	proteome = "./inputs/olink_instance_0.csv",   # rows=samples with 'eid', cols=proteins
# 	pheno = "./inputs/pheno_for_limma.csv",       # rows=samples with 'eid' and exposures/covariates
# 	exposures = c("Nitrogen_oxides_2010", "NO2_mean", "PM25_2010", "PM10_mean"),
# 	covariates = c("Age", "Sex", "BMI", "Smoking_status", "Drinking_status", "Ethnic"),
# 	adjust_map = list(
# 		Nitrogen_oxides_2010 = c("PM25_2010"),
# 		NO2_mean = c("PM10_mean"),
# 		PM25_2010 = c("Nitrogen_oxides_2010"),
# 		PM10_mean = c("NO2_mean")
# 	),
# 	write_dir = "./outputs",
# 	basename = "PWAS_AirPollution"
# )


