suppressPackageStartupMessages({
	library(data.table)
	library(dplyr)
	library(tidyr)
	library(glmnet)
	library(caret)
	library(pROC)
})

# ------------------------------
# APPRS construction
# ------------------------------

# Build APPRS score matrix per disease using ACME weights
build_apprs_scores <- function(
	weights_df,          # data.frame with columns: mediation_factor, Disease, ACME
	proteome_df,         # data.frame with sample rows, protein columns (lowercase names), includes 'eid'
	disease_list = NULL, # optional subset of Disease to build
	mediator_col = "mediation_factor",
	disease_col = "Disease",
	weight_col = "ACME",
	sample_id_col = "eid",
	write_dir = NULL,
	basename = "APPRS"
) {
	wd <- if (is.character(weights_df)) fread(weights_df) %>% as.data.frame() else weights_df
	prot <- if (is.character(proteome_df)) fread(proteome_df) %>% as.data.frame() else proteome_df
	stopifnot(sample_id_col %in% colnames(prot))
	if (!is.null(disease_list)) wd <- wd[wd[[disease_col]] %in% disease_list, , drop = FALSE]

	# Ensure lowercase protein names to match proteome columns
	wd[[mediator_col]] <- toupper(wd[[mediator_col]])
	prot_cols <- setdiff(colnames(prot), sample_id_col)
	colnames(prot)[match(prot_cols, colnames(prot))] <- tolower(prot_cols)

	# Prepare output matrix: samples x diseases
	diseases <- unique(wd[[disease_col]])
	scores <- matrix(0, nrow = nrow(prot), ncol = length(diseases))
	colnames(scores) <- paste0("ApPRS_", diseases)

	for (d in diseases) {
		wd_d <- wd[wd[[disease_col]] == d, c(mediator_col, weight_col), drop = FALSE]
		wd_d[[mediator_col]] <- tolower(wd_d[[mediator_col]])
		wd_d <- wd_d[wd_d[[mediator_col]] %in% colnames(prot), , drop = FALSE]
		if (!nrow(wd_d)) next
		W <- as.matrix(wd_d[[weight_col]])
		rownames(W) <- wd_d[[mediator_col]]
		X <- as.matrix(prot[, rownames(W), drop = FALSE])
		X[is.na(X)] <- 0
		# Multiply sample-by-protein (X) with protein-by-1 weights (W) => sample scores
		s <- X %*% W
		# Normalize by sum of weights as in original
		denom <- sum(wd_d[[weight_col]], na.rm = TRUE)
		if (denom == 0) denom <- 1
		s <- as.numeric(s) / denom
		scores[, paste0("ApPRS_", d)] <- s
	}

	res <- cbind(prot[, sample_id_col, drop = FALSE], as.data.frame(scores))
	if (!is.null(write_dir)) {
		if (!dir.exists(write_dir)) dir.create(write_dir, recursive = TRUE, showWarnings = FALSE)
		fwrite(res, file = file.path(write_dir, paste0(basename, "_scores.csv")))
	}
	res
}

# ------------------------------
# LASSO prediction (three model variants)
# ------------------------------

run_lasso_three_models <- function(
	data,
	outcome,
	basic_covariates = c("age", "sex", "smoking", "alcohol", "BMI"),
	clinical_covariates = c("ApoB_ApoA", "Albumin", "Creatinine", "HbA1c", "Glucose"),
	apprs_var = NULL,   # e.g., "ApPRS_HeartFailure"
	sample_id_col = "eid",
	seed = 123,
	kfold = 10
) {
	df <- if (is.character(data)) fread(data) %>% as.data.frame() else data
	stopifnot(all(c(outcome, basic_covariates, clinical_covariates) %in% colnames(df)))
	if (!is.null(apprs_var)) stopifnot(apprs_var %in% colnames(df))
	set.seed(seed)

	# Prepare response as factor for ROC
	y <- factor(df[[outcome]], labels = c("Class0", "Class1"))

	# Helper to fit one glmnet via caret with ROC metric
	fit_glmnet <- function(X, y) {
		ctrl <- trainControl(method = "cv", number = kfold, summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = TRUE)
		train(x = as.matrix(scale(X)), y = y, method = "glmnet",
			trControl = ctrl,
			tuneGrid = expand.grid(alpha = 1, lambda = 10^seq(-5, -2, length = 100)),
			metric = "ROC",
			family = "binomial")
	}

	# Model 1: Basic only
	X1 <- df[, basic_covariates, drop = FALSE]
	m1 <- fit_glmnet(X1, y)
	roc1 <- pROC::roc(y, predict(m1, as.matrix(scale(X1)), type = "prob")[, 2])

	# Model 2: Basic + Clinical
	X2 <- df[, c(basic_covariates, clinical_covariates), drop = FALSE]
	m2 <- fit_glmnet(X2, y)
	roc2 <- pROC::roc(y, predict(m2, as.matrix(scale(X2)), type = "prob")[, 2])

	# Model 3: Basic + Clinical + ApPRS (if provided)
	m3 <- NULL; roc3 <- NULL
	if (!is.null(apprs_var)) {
		X3 <- df[, c(basic_covariates, clinical_covariates, apprs_var), drop = FALSE]
		m3 <- fit_glmnet(X3, y)
		roc3 <- pROC::roc(y, predict(m3, as.matrix(scale(X3)), type = "prob")[, 2])
	}

	list(
		model_basic = m1,
		auc_basic = as.numeric(pROC::auc(roc1)),
		model_basic_clinical = m2,
		auc_basic_clinical = as.numeric(pROC::auc(roc2)),
		model_basic_clinical_apprs = m3,
		auc_basic_clinical_apprs = if (!is.null(roc3)) as.numeric(pROC::auc(roc3)) else NA_real_
	)
}

# ------------------------------
# Example (commented)
# ------------------------------

# 1) Build APPRS
# apprs <- build_apprs_scores(
# 	weights_df = "./inputs/mediation_acme.csv",           # columns: mediation_factor, Disease, ACME
# 	proteome_df = "./inputs/proteome_matrix.csv",          # rows=samples, cols=proteins (lowercase), plus 'eid'
# 	disease_list = c("HeartFailure","CKD","IHD"),
# 	write_dir = "./outputs",
# 	basename = "APPRS"
# )

# 2) Merge APPRS with clinical covariates and run LASSO
# merged <- merge(apprs, read.csv("./inputs/clinical_outcomes.csv"), by = "eid")
# res <- run_lasso_three_models(
# 	data = merged,
# 	outcome = "HeartFailure",
# 	basic_covariates = c("ApoB_ApoA","Albumin","Creatinine","HbA1c","Glucose"),
# 	clinical_covariates = c("age","sex","smoking","alcohol","BMI"),
# 	apprs_var = "ApPRS_HeartFailure",
# 	kfold = 10
# )
# print(res$auc_basic); print(res$auc_basic_clinical); print(res$auc_basic_clinical_apprs)


