
suppressPackageStartupMessages({
	library(data.table)
	library(dplyr)
})

# ------------------------------
# Helpers
# ------------------------------

# Ensure a character vector of column names exists in a data.frame
ensure_columns_exist <- function(dataset, required_columns, label) {
	missing_columns <- setdiff(required_columns, colnames(dataset))
	if (length(missing_columns) > 0) {
		stop(sprintf("Missing required %s columns: %s", label, paste(missing_columns, collapse = ", ")))
	}
}

# Extract incident disease events from UKB ICD-10 longitudinal columns (41270 codes, 41280 dates)
# dataset: data.frame with 'eid', code columns (41270-*) and date columns (41280-*)
# icd_code_cols: character vector of code-column names (e.g., grep('41270', colnames(dataset), value=TRUE))
# icd_date_cols: character vector of date-column names (aligned with icd_code_cols; same length and order)
# icd_pattern: regex to match ICD-10 code(s), e.g., 'I50' or 'J12|J13|...'
# disease_name: base name for output columns
# date_origin: date origin for numeric UKB dates ('1970-01-01' in original code)
extract_incident_by_icd <- function(dataset,
		icd_code_cols,
		icd_date_cols,
		icd_pattern,
		disease_name,
		date_origin = "1970-01-01") {

	stopifnot(length(icd_code_cols) == length(icd_date_cols))
	ensure_columns_exist(dataset, c("eid", icd_code_cols, icd_date_cols), "ICD")

	# Subset for speed
	data_sub <- dataset[, c("eid", icd_code_cols, icd_date_cols)]

	matched_eids <- vector(mode = "list", length = length(icd_code_cols))
	matched_dates <- vector(mode = "list", length = length(icd_code_cols))

	for (i in seq_along(icd_code_cols)) {
		codes <- data_sub[[icd_code_cols[i]]]
		rows <- grep(icd_pattern, codes)
		if (length(rows) == 0) next
		matched_eids[[i]] <- data_sub$eid[rows]
		matched_dates[[i]] <- data_sub[[icd_date_cols[i]]][rows]
	}

	if (length(Filter(Negate(is.null), matched_eids)) == 0) {
		# No matches; return empty structure with expected columns
		out <- data.frame(
			eid = integer(0),
			status = integer(0),
			date = as.Date(character(0))
		)
		colnames(out) <- c("eid", disease_name, paste0("date_", disease_name))
		return(out)
	}

	eids <- unlist(matched_eids, use.names = FALSE)
	dates_raw <- unlist(matched_dates, use.names = FALSE)
	dates <- as.Date(dates_raw, origin = date_origin)

	events <- data.frame(
		eid = eids,
		status = 1L,
		date = dates,
		stringsAsFactors = FALSE
	)

	# For multiple occurrences per eid, keep earliest date
	events <- events %>%
		group_by(eid) %>%
		summarise(status = first(status), date = min(date, na.rm = TRUE), .groups = "drop")

	colnames(events) <- c("eid", disease_name, paste0("date_", disease_name))
	return(events)
}

# Merge list of disease event tables into a single wide table keyed by eid
merge_disease_tables <- function(base_df, disease_tables) {
	stopifnot("eid" %in% colnames(base_df))
	merged <- base_df
	for (tbl in disease_tables) {
		if (nrow(tbl) == 0) next
		merged <- merge(merged, tbl, by = "eid", all.x = TRUE)
	}
	merged
}

# Build a base frame with eid and baseline date if available
build_base_frame <- function(dataset, baseline_date_col = NULL) {
	ensure_columns_exist(dataset, c("eid"), "base")
	base <- dataset[, "eid", drop = FALSE]
	if (!is.null(baseline_date_col) && baseline_date_col %in% colnames(dataset)) {
		base[[baseline_date_col]] <- dataset[[baseline_date_col]]
	}
	base
}

# UKB disease definitions (ICD-10 patterns)
ukb_disease_definitions <- function() {
	list(
		# Unreported list in the original
		LiverFibrosis = "K70|K71|K74",
		Schizophrenia = "F20|F21|F22|F25",
		PD = "F023|G20",
		VD = "F01",
		BrainCancer = "C70|C71|C72|C751|C752|C753|C754|C755",
		IschaemicHeartDisease = "I250|I251|I253|I254|I255|I256|I258|I259",
		Hypertension = "I10|I11|I12|I13|I14|I15",
		Arrhythmias = "I46|I47|I48|I49",
		Anemia = "D50|D51|D52|D53|D54|D55|D56|D57|D58|D59|D60|D61|D62|D63|D64",
		SkinDisorder = "^L", # '\nL' -> anchored to start to avoid over-match
		MultipleSclerosis = "G35",
		PeripheralArteryDisease = "I70|I71|I72|I73|I74|I75|I76|I77|I78|I79",
		LupusErythematosus = "L93|M32",
		InflammatoryBowelDisease = "K50|K51",
		Endometriosis = "N80",
		SleepDisorder = "G47",
		Obesity = "E66",
		GynaecologicalCancer = "C54|C55|C56|D06|N87|C53|C51|C52|C57",
		ColorectalCancer = "C18|C19|C20|C21",
		ProstateCancer = "C61",

		# Reported 16 in the original
		Depression = "F32|F33",
		BipolarDisorder = "F31",
		AD = "F00|G30",
		RA = "M05|M06",
		Osteoporosis = "M80|M81|M82",
		AllergicRhinitis = "J301|J302|J303|J304",
		HeartFailure = "I50",
		IschamicStroke = "I630|I631|I632|I633|I634|I635|I638|I639|I693",
		Asthma = "J45|J46",
		COPD = "J41|J42|J43|J44",
		LungCancer = "C33|C34",
		Pneumonia = "J12|J13|J14|J15|J16|J17|J18",
		T2D = "E11|E12",
		NAFLD = "K758|K759|K760",
		CKD = "I120|I131|I132|N180|N181|N182|N183|N184|N185|N188|N189",
		BreastCancer = "C50|D05"
	)
}

# Derive disease tables for a dataset given the standard UKB 41270/41280 columns
derive_disease_table <- function(dataset,
		baseline_date_col = "53-0.0",
		code_column_selector = function(cols) grep("41270", cols, value = TRUE),
		date_column_selector = function(cols) grep("41280", cols, value = TRUE),
		disease_patterns = ukb_disease_definitions(),
		date_origin = "1970-01-01") {

	all_cols <- colnames(dataset)
	icd_code_cols <- code_column_selector(all_cols)
	icd_date_cols <- date_column_selector(all_cols)

	if (length(icd_code_cols) == 0 || length(icd_date_cols) == 0) {
		stop("No ICD code/date columns detected. Please adjust selectors.")
	}

	if (length(icd_code_cols) != length(icd_date_cols)) {
		stop("ICD code and date columns must have the same length and order.")
	}

	base <- build_base_frame(dataset, baseline_date_col)

	disease_tables <- lapply(names(disease_patterns), function(name) {
		pattern <- disease_patterns[[name]]
		extract_incident_by_icd(
			dataset = dataset,
			icd_code_cols = icd_code_cols,
			icd_date_cols = icd_date_cols,
			icd_pattern = pattern,
			disease_name = name,
			date_origin = date_origin
		)
	})
	names(disease_tables) <- names(disease_patterns)

	DiseaseTable <- merge_disease_tables(base, disease_tables)
	DiseaseTable
}

# Create simple summary (cases and non-cases) for each disease
summarise_disease_counts <- function(disease_table, exclude_columns = c("eid", "53-0.0")) {
	status_cols <- setdiff(colnames(disease_table), c(exclude_columns, grep("^date_", colnames(disease_table), value = TRUE)))
	status_cols <- status_cols[status_cols != ""]
	N <- nrow(disease_table)
	res <- lapply(status_cols, function(col) {
		cases <- sum(as.numeric(disease_table[[col]] %in% 1), na.rm = TRUE)
		data.frame(Disease = col, Cases = cases, Normal = N - cases, stringsAsFactors = FALSE)
	})
	bind_rows(res)
}

# Optional: Compute AP index (air pollution score) as in the original, provided all required columns exist
compute_ap_index <- function(pheno_df) {
	required <- c(
		"Overall_health_rating", # '2178-0.0'
		"Age", "Sex", "BMI", "Smoking_status", "Drinking_freq",
		"PM25_2010", "Nitrogen_oxides_2010", "NO2_2005", "NO2_2006", "NO2_2007", "NO2_2010", "PM10_2007", "PM10_2010"
	)
	ensure_columns_exist(pheno_df, required, "AP index")

	pheno_df <- pheno_df %>%
		mutate(
			NO2_mean = (as.numeric(NO2_2005) + as.numeric(NO2_2006) + as.numeric(NO2_2007) + as.numeric(NO2_2010)) / 4,
			PM10_mean = (as.numeric(PM10_2007) + as.numeric(PM10_2010)) / 2
		)

	ap_predictors <- c("PM25_2010", "Nitrogen_oxides_2010", "NO2_mean", "PM10_mean")
	weights <- numeric(length(ap_predictors))
	coeffs <- numeric(length(ap_predictors))
	for (i in seq_along(ap_predictors)) {
		fml <- stats::as.formula(paste0("Overall_health_rating ~ ", ap_predictors[i], " + Age + Sex + BMI + Drinking_freq + Smoking_status"))
		fit <- lm(fml, data = pheno_df)
		sm <- summary(fit)
		weights[i] <- abs(sm$coefficients[2, 3])
		coeffs[i] <- sm$coefficients[2, 1]
	}

	score <- 0
	for (i in seq_along(ap_predictors)) {
		term <- as.numeric(pheno_df[[ap_predictors[i]]]) * coeffs[i]
		term[is.na(term)] <- 0
		score <- score + term
	}
	score <- score / sum(coeffs) * length(ap_predictors)
	pheno_df$ap_score <- score
	pheno_df
}

# Build phenotype frame with standard columns if present
build_pheno_for_cox <- function(dataset) {
	cols <- colnames(dataset)
	get <- function(id) {
		col <- cols[cols == id]
		if (length(col) == 1) dataset[[col]] else NA
	}
	data.frame(
		eid = dataset$eid,
		NO2_2010 = get("24003-0.0"),
		NO2_2007 = get("24018-0.0"),
		NO2_2006 = get("24017-0.0"),
		NO2_2005 = get("24016-0.0"),
		Nitrogen_oxides_2010 = get("24004-0.0"),
		PM10_2007 = get("24019-0.0"),
		PM10_2010 = get("24005-0.0"),
		PM25_2010 = get("24006-0.0"),
		PM25_absorbance_2010 = get("24007-0.0"),
		PM25_10um_2010 = get("24008-0.0"),
		Traffic_intensity1 = get("24009-0.0"),
		Traffic_intensity2 = get("24011-0.0"),
		Inverse_distance1 = get("24010-0.0"),
		Inverse_distance2 = get("24012-0.0"),
		Total_traffic_load = get("24013-0.0"),
		Close_to_major_load = get("24014-0.0"),
		Overall_health_rating = get("2178-0.0"),
		Sex = get("31-0.0"),
		Age = get("21003-0.0"),
		BMI = get("23104-0.0"),
		Ethnic = get("21000-0.0"),
		Date_of_attending_assessment_centre = get("53-0.0"),
		Smoking_status = get("20116-0.0"),
		Drinking_status = get("20117-0.0"),
		Drinking_freq = get("1558-0.0"),
		stringsAsFactors = FALSE
	)
}

# ------------------------------
# Public API
# ------------------------------

#' Derive DiseaseTable and optional summary from a UKB-like dataset
#'
#' @param input Either a data.frame or a path to a CSV/TSV file readable by data.table::fread
#' @param baseline_date_col Baseline date column name (default '53-0.0')
#' @param write_dir Optional directory to write outputs; if NULL, nothing is written
#' @param write_basename Base name for output files (without extension)
#' @param return_pheno Whether to return a phenotype frame and computed ap_score
#' @return A list with DiseaseTable, Summary (and optionally Pheno)
derive_ukb_diseases_public <- function(
	input,
	baseline_date_col = "53-0.0",
	write_dir = NULL,
	write_basename = "UKB-Disease-Outputs",
	return_pheno = FALSE
) {
	if (is.character(input)) {
		dataset <- fread(input) %>% as.data.frame()
	} else if (is.data.frame(input)) {
		dataset <- input
	} else {
		stop("input must be a data.frame or a path to a delimited file")
	}

	ensure_columns_exist(dataset, c("eid"), "input")

	DiseaseTable <- derive_disease_table(dataset, baseline_date_col = baseline_date_col)
	Summary <- summarise_disease_counts(DiseaseTable, exclude_columns = c("eid", baseline_date_col))

	result <- list(DiseaseTable = DiseaseTable, Summary = Summary)

	if (isTRUE(return_pheno)) {
		pheno <- build_pheno_for_cox(dataset)
		# ap_score is optional; compute only if all required columns exist
		required_ap_cols <- c("2178-0.0", "21003-0.0", "31-0.0", "23104-0.0", "20116-0.0", "1558-0.0",
			"24006-0.0", "24004-0.0", "24016-0.0", "24017-0.0", "24018-0.0", "24003-0.0", "24019-0.0", "24005-0.0")
		if (all(required_ap_cols %in% colnames(dataset))) {
			# Rename to the expected names used by compute_ap_index
			colnames(pheno)[colnames(pheno) == "Overall_health_rating"] <- "Overall_health_rating"
			pheno <- compute_ap_index(pheno)
		}
		result$Pheno <- pheno
	}

	if (!is.null(write_dir)) {
		if (!dir.exists(write_dir)) {
			dir.create(write_dir, recursive = TRUE, showWarnings = FALSE)
		}
		fwrite(result$DiseaseTable, file = file.path(write_dir, paste0(write_basename, "_DiseaseTable.csv")))
		fwrite(result$Summary, file = file.path(write_dir, paste0(write_basename, "_Summary.csv")))
		if (!is.null(result$Pheno)) {
			fwrite(result$Pheno, file = file.path(write_dir, paste0(write_basename, "_Pheno.csv")))
		}
	}

	result
}

# ------------------------------
# Example usage (commented)
# ------------------------------

# res <- derive_ukb_diseases_public(
# 	input = "path/to/HH_data_220812_0512.csv",
# 	baseline_date_col = "53-0.0",
# 	write_dir = "./outputs",
# 	write_basename = "UKB_Disease_50w",
# 	return_pheno = TRUE
# )


