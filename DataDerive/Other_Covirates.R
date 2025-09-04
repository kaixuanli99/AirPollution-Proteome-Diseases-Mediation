
suppressPackageStartupMessages({
	library(data.table)
	library(dplyr)
	library(stringr)
})

ensure_columns_exist <- function(dataset, required_columns, label) {
	missing_columns <- setdiff(required_columns, colnames(dataset))
	if (length(missing_columns) > 0) {
		stop(sprintf("Missing required %s columns: %s", label, paste(missing_columns, collapse = ", ")))
	}
}

# Build core covariates per original logic (excluding PRS)
derive_covariates_core <- function(dataset) {
	data <- dataset

	# assessment_center (54-0.0)
	if ("X54.0.0" %in% colnames(data)) {
		data$assessment_center <- data[["X54.0.0"]]
	}

	# education (6138)
	edu_cols <- colnames(data)[grep("^X6138\\.0\\.", colnames(data))]
	if (length(edu_cols) > 0) {
		cond <- Reduce(`|`, lapply(edu_cols, function(col) data[[col]] %in% c(1, 5, 6)))
		data$education <- ifelse(cond, 1, 0)
	}

	# healthy_diet_score removed per request

	# income (738-0.0)
	if ("X738.0.0" %in% colnames(data)) {
		data$income <- ifelse(data[["X738.0.0"]] %in% c(1, 2, 3, 4, 5), data[["X738.0.0"]], 0)
	}

	# physical_activity removed per request

	# townsend_index (189-0.0)
	if ("X189.0.0" %in% colnames(data)) data$townsend_index <- data[["X189.0.0"]]

	# FEV1 (take row-wise max across arrays)
	fev_cols <- c("X3063.0.0", "X3063.0.1", "X3063.0.2")
	fev_cols <- fev_cols[fev_cols %in% colnames(data)]
	if (length(fev_cols) > 0) data$FEV1 <- do.call(pmax, c(data[fev_cols], na.rm = TRUE))

	# FVC (take row-wise max across arrays)
	fvc_cols <- c("X3062.0.0", "X3062.0.1", "X3062.0.2")
	fvc_cols <- fvc_cols[ fvc_cols %in% colnames(data) ]
	if (length(fvc_cols) > 0) data$FVC <- do.call(pmax, c(data[fvc_cols], na.rm = TRUE))

	# race (21000-0.0) recode
	if ("X21000.0.0" %in% colnames(data)) {
		val <- data[["X21000.0.0"]]
		data$race <- ifelse(val %in% c(1, 1001, 2001, 3001, 4001), 1,
			ifelse(val %in% c(2, 1002, 2002, 3002, 4002), 2,
				ifelse(val %in% c(3, 1003, 2003, 3003, 4003, 5), 3,
					ifelse(val %in% c(4, 2004, 3004), 4, NA))))
	}

	# Drop chronic disease baselines and baseline NO2 per public spec

	# self-report health rating (2178-0.0), recode 1<->4, 2<->3; set -1/-3 to 0
	if ("X2178.0.0" %in% colnames(data)) {
		data$health_rating <- data$X2178.0.0
		idx1 <- which(data$health_rating == 1)
		idx2 <- which(data$health_rating == 2)
		idx3 <- which(data$health_rating == 3)
		idx4 <- which(data$health_rating == 4)
		data$health_rating[idx1] <- 4
		data$health_rating[idx4] <- 1
		data$health_rating[idx2] <- 3
		data$health_rating[idx3] <- 2
		data$health_rating <- ifelse(data$health_rating %in% c(-1, -3), 0, data$health_rating)
	}

	# Select output columns that exist
	keep <- c("eid", "assessment_center", "education", "income",
		"townsend_index", "FEV1", "FVC", "race", "health_rating")
	keep <- intersect(keep, colnames(data))
	data[, keep, drop = FALSE]
}

# Public API: derive and optionally write to CSV
derive_covariates_public <- function(input, write_path = NULL) {
	if (is.character(input)) {
		dataset <- fread(input) %>% as.data.frame()
	} else if (is.data.frame(input)) {
		dataset <- input
	} else stop("input must be a data.frame or a path to a delimited file")
	stopifnot("eid" %in% colnames(dataset))
	covars <- derive_covariates_core(dataset)
	if (!is.null(write_path)) fwrite(covars, file = write_path)
	covars
}

# Example usage (commented):
# covars <- derive_covariates_public(
# 	input = "path/to/HH_data_220812_0512.csv",
# 	write_path = "./outputs/UKB_covariates_public.csv"
# )


