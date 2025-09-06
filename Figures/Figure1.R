
rm(list = ls())
suppressPackageStartupMessages({
	library(openxlsx)
	library(dplyr)
	library(stringr)
	library(ggplot2)
	library(data.table)
	library(ggpubr)
	library(GGally)
	library(tidyr)
	library(scales)
})

# ------------------------------
# I/O (edit these)
# ------------------------------

reported_xlsx <- "./inputs/Reported AP-Disease HR in UKB.xlsx"   # sheet indices same as original
unreported_csv <- "./inputs/500k-case-IQR-Unreported.csv"
apscore_csv <- "./inputs/500k-case-IQR-APScore.csv"

out_reported_pdf <- "./Figure1_reported.pdf"
out_reported_png <- "./Figure1_reported.png"

out_unreported_pdf <- "./Figure1_unreported.pdf"
out_unreported_png <- "./Figure1_unreported.png"

out_apscore_pdf <- "./Figure1B.pdf"
out_apscore_png <- "./Figure1B.png"
out_heatmap_pdf <- "./Figure1_heatmap.pdf"
out_heatmap_png <- "./Figure1_heatmap.png"

# ------------------------------
# Reported forest (left part of fig1)
# ------------------------------

data <- read.xlsx(reported_xlsx, sheet = 4) %>% as.data.frame()

unreported <- fread(unreported_csv)
apscore <- fread(apscore_csv)

unreported$Disease <- factor(unreported$Disease, levels = data$Incident.Disease[17:36])
unreported <- arrange(unreported, Disease)
apscore$Disease <- factor(apscore$Disease, levels = data$Incident.Disease)
apscore <- arrange(apscore, Disease)

for (i in 3:7) {
	data[, i] <- gsub(' ', '', data[, i])
	data[, i] <- gsub('Unreported', '1.00[1.00,1.00]', data[, i])
}

extracted_numbers <- str_extract_all(data[, 3], "[0-9]+\\.[0-9]+")
extrat_df <- do.call(rbind, lapply(extracted_numbers, function(x) as.numeric(x)))
result <- data.frame(extrat_df)
colnames(result) <- c("PM25_Value", "PM25_Lower", "PM25_Upper")
data <- cbind(data, result)
colnames(data)[3:9] <- c('HR_PM2_5', 'HR_PM2_5_10', 'HR_PM10', 'HR_NO2', 'HR_NOx', 'HR_O3', 'HR_AP_Score')

data_long <- data.frame(
	disease = rep(data$Incident.Disease, 5),
	HR = c(data$HR_PM2_5, data$HR_PM2_5_10, data$HR_PM10, data$HR_NO2, data$HR_NOx),
	Pollutants = c(rep('PM2.5', 36), rep('PM2.5-10', 36), rep('PM10', 36), rep('NO2', 36), rep('NOx', 36))
)

extracted_numbers <- str_extract_all(data_long$HR, "[0-9]+\\.[0-9]+")
extrat_df <- do.call(rbind, lapply(extracted_numbers, function(x) as.numeric(x)))
result <- data.frame(extrat_df)
colnames(result) <- c("value", "lower", "high")
data_long <- cbind(data_long, result)

diseaselist <- unique(data_long$disease)
data_long$disease <- data_long$disease %>% factor(levels = diseaselist)
data_long$Pollutants <- data_long$Pollutants %>% factor(levels = c('NOx', 'NO2', 'PM10', 'PM2.5-10', 'PM2.5'))
data_long$sigificance <- ifelse((data_long$lower - 1) * (data_long$high - 1) > 0, 1, 0)
data_long <- data_long %>% arrange(disease)
data_long <- data_long[data_long$Pollutants != 'PM2.5-10', ]
data_long$Pollutants <- data_long$Pollutants %>% factor(levels = c('NOx', 'NO2', 'PM10', 'PM2.5'))

diseaselist <- c('Depression', 'Bipolar Disorder', 'AD', 'RA', 'Osteoporosis', 'Allergic Rhinitis',
				'Heart Failure', 'Ischamic Stroke',
				'Asthma', 'COPD', 'Lung Cancer', 'Pneumonia', 'T2D', 'NAFLD', 'Chronic Kidney Disease', 'Breast Cancer')

data_long$disease <- gsub('AllergicRhinitis', 'Allergic Rhinitis', data_long$disease)
data_long$disease <- gsub('BipolarDisorder', 'Bipolar Disorder', data_long$disease)
data_long$disease <- gsub('BreastCancer', 'Breast Cancer', data_long$disease)
data_long$disease <- gsub('HeartFailure', 'Heart Failure', data_long$disease)
data_long$disease <- gsub('IschamicStroke', 'Ischamic Stroke', data_long$disease)
data_long$disease <- gsub('LungCancer', 'Lung Cancer', data_long$disease)
data_long$disease <- gsub('BrainCancer', 'Brain Cancer', data_long$disease)
data_long$disease <- gsub('IschaemicHeartDisease', 'Ischaemic Heart Disease', data_long$disease)
data_long$disease <- gsub('LiverFibrosis', 'Liver Fibrosis', data_long$disease)
data_long$disease <- gsub('MultipleSclerosis', 'Multiple Sclerosis', data_long$disease)
data_long$disease <- gsub('SkinDisorder', 'Skin Disorder', data_long$disease)
data_long$disease <- gsub('PeripheralArteryDisease', 'Peripheral Artery Disease', data_long$disease)
data_long$disease <- gsub('SleepDisorder', 'Sleep Disorder', data_long$disease)
data_long$disease <- gsub('ColorectalCancer', 'Colorectal Cancer', data_long$disease)
data_long$disease <- gsub('LupusErythematosus', 'Lupus Erythematosus', data_long$disease)
data_long$disease <- gsub('GynaecologicalCancer', 'Gynaecological Cancer', data_long$disease)
data_long$disease <- gsub('ProstateCancer', 'Prostate Cancer', data_long$disease)
data_long$disease <- gsub('InflammatoryBowelDisease', 'Inflammatory Bowel Disease', data_long$disease)
data_long$disease <- gsub('ChronicKidneyDisease', 'Chronic Kidney Disease', data_long$disease)

data_long1 <- data_long[data_long$disease %in% diseaselist[1:8], ]
data_long1$disease <- factor(data_long1$disease, levels = rev(diseaselist[1:8]))
data_long1$text_position <- data_long1$high + 0.015

p1 <- ggplot(data_long1, aes(x = disease, y = value, group = as.factor(Pollutants), colour = as.factor(Pollutants))) +
	geom_hline(yintercept = 1, linetype = 2, color = "grey4", size = 1.15) +
	geom_point(position = position_dodge(width = 0.9), size = 4) +
	geom_errorbar(aes(ymax = lower, ymin = high), width = 0, alpha = 0.8, position = position_dodge(width = 0.9), size = 2) +
	ylim(0, 3) +
	theme_minimal() +
	geom_text(aes(y = text_position, label = ifelse(sigificance, '***', '')), position = position_dodge(width = 0.9), hjust = -0.2, vjust = 0.76, size = 6) +
	theme_bw() +
	coord_flip() +
	theme(axis.title.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank(), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14), axis.title.x = element_blank(), plot.title = element_text(size = 20.5), legend.position = 'none') +
	labs(colour = "Pollutants", y = "Hazard Ratio") +
	scale_colour_manual(values = c("#54A966", '#C55054', '#4D71AF', '#4C4C4C')) +
	guides(colour = guide_legend(reverse = TRUE)) +
	geom_stripped_cols(even = "#CCCCCC", alpha = .1)

data_long2 <- data_long[data_long$disease %in% diseaselist[9:16], ]
data_long2$disease <- factor(data_long2$disease, levels = rev(diseaselist[9:16]))
data_long2$text_position <- data_long2$high + 0.015

p2 <- ggplot(data_long2, aes(x = disease, y = value, group = as.factor(Pollutants), colour = as.factor(Pollutants))) +
	geom_hline(yintercept = 1, linetype = 2, color = "grey4", size = 1.15) +
	geom_point(position = position_dodge(width = 0.9), size = 4) +
	geom_errorbar(aes(ymax = lower, ymin = high), width = 0, alpha = 0.8, position = position_dodge(width = 0.9), size = 2) +
	ylim(0, 3) +
	theme_minimal() +
	geom_text(aes(y = text_position, label = ifelse(sigificance, '***', '')), position = position_dodge(width = 0.9), hjust = -0.2, vjust = 0.76, size = 6) +
	theme_bw() +
	coord_flip() +
	theme(axis.title.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank(), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14), axis.title.x = element_blank(), plot.title = element_text(size = 20.5)) +
	labs(colour = "Pollutants", y = "Hazard Ratio") +
	scale_colour_manual(values = c("#54A966", '#C55054', '#4D71AF', '#4C4C4C')) +
	guides(colour = guide_legend(reverse = TRUE)) +
	geom_stripped_cols(even = "#CCCCCC", alpha = .1)

p_reported <- ggarrange(p1, p2, widths = c(1, 1.46))
p_reported <- annotate_figure(p_reported, top = text_grob("Disease Reported", face = "bold", size = 14), bottom = text_grob("Hazard Ratio", face = "bold", vjust = 0.5, hjust = 0.5, size = 14))
ggsave(out_reported_png, plot = p_reported, width = 10, height = 8, dpi = 350)
ggsave(out_reported_pdf, plot = p_reported, device = "pdf", width = 10, height = 8)

# (O3 panel removed as requested)

# ------------------------------
# AP Score forest (Figure 1B)
# ------------------------------

APScore_disease <- fread(apscore_csv) %>% as.data.frame()
colnames(APScore_disease)[1:3] <- c("value", "lower", "high")
APScore_disease$exposure <- gsub('ap_score_IQR', 'AP Score', APScore_disease$exposure)

APScore_disease$Disease <- gsub('AllergicRhinitis', 'Allergic Rhinitis', APScore_disease$Disease)
APScore_disease$Disease <- gsub('BipolarDisorder', 'Bipolar Disorder', APScore_disease$Disease)
APScore_disease$Disease <- gsub('BreastCancer', 'Breast Cancer', APScore_disease$Disease)
APScore_disease$Disease <- gsub('HeartFailure', 'Heart Failure', APScore_disease$Disease)
APScore_disease$Disease <- gsub('IschamicStroke', 'Ischamic Stroke', APScore_disease$Disease)
APScore_disease$Disease <- gsub('LungCancer', 'Lung Cancer', APScore_disease$Disease)
APScore_disease$Disease <- gsub('BrainCancer', 'Brain Cancer', APScore_disease$Disease)
APScore_disease$Disease <- gsub('IschaemicHeartDisease', 'Ischaemic Heart Disease', APScore_disease$Disease)
APScore_disease$Disease <- gsub('LiverFibrosis', 'Liver Fibrosis', APScore_disease$Disease)
APScore_disease$Disease <- gsub('MultipleSclerosis', 'Multiple Sclerosis', APScore_disease$Disease)
APScore_disease$Disease <- gsub('SkinDisorder', 'Skin Disorder', APScore_disease$Disease)
APScore_disease$Disease <- gsub('PeripheralArteryDisease', 'Peripheral Artery Disease', APScore_disease$Disease)
APScore_disease$Disease <- gsub('SleepDisorder', 'Sleep Disorder', APScore_disease$Disease)
APScore_disease$Disease <- gsub('ColorectalCancer', 'Colorectal Cancer', APScore_disease$Disease)
APScore_disease$Disease <- gsub('LupusErythematosus', 'Lupus Erythematosus', APScore_disease$Disease)
APScore_disease$Disease <- gsub('GynaecologicalCancer', 'Gynaecological Cancer', APScore_disease$Disease)
APScore_disease$Disease <- gsub('ProstateCancer', 'Prostate Cancer', APScore_disease$Disease)
APScore_disease$Disease <- gsub('InflammatoryBowelDisease', 'Inflammatory Bowel Disease', APScore_disease$Disease)
APScore_disease$Disease <- gsub('ChronicKidneyDisease', 'Chronic Kidney Disease', APScore_disease$Disease)

APScore_disease$sigificance <- ifelse((APScore_disease$lower - 1) * (APScore_disease$high - 1) > 0, 1, 0)
APScore_disease$text_position <- APScore_disease$high + 0.01

diseaselist <- c('Depression', 'Schizophrenia', 'Bipolar Disorder', 'Sleep Disorder', 'PD', 'VD', 'AD', 'Brain Cancer',
				'RA', 'Multiple Sclerosis', 'Peripheral Artery Disease', 'Lupus Erythematosus', 'Inflammatory Bowel Disease',
				'Osteoporosis', 'Allergic Rhinitis', 'Endometriosis',
				'Heart Failure', 'Ischamic Stroke', 'Ischaemic Heart Disease', 'Hypertension', 'Arrhythmias', 'Anemia',
				'Asthma', 'COPD', 'Lung Cancer', 'Pneumonia', 'Skin Disorder',
				'T2D', 'NAFLD', 'Chronic Kidney Disease', 'Liver Fibrosis', 'Obesity',
				'Breast Cancer', 'Gynaecological Cancer', 'Colorectal Cancer', 'Prostate Cancer')
APScore_disease$Disease <- factor(APScore_disease$Disease, levels = rev(diseaselist))

p_apscore <- ggplot(APScore_disease, aes(x = Disease, y = value)) +
	geom_hline(yintercept = 1, linetype = 2, color = "black", size = 1.15) +
	geom_point(position = position_dodge(width = 0.9), size = 4, color = '#F5B879') +
	geom_errorbar(aes(ymax = lower, ymin = high), width = 0, alpha = 0.8, position = position_dodge(width = 0.9), size = 2, color = '#F5B879') +
	ylim(0.8, 1.4) +
	theme_minimal() +
	geom_text(aes(y = text_position, label = ifelse(sigificance, '***', '')), position = position_dodge(width = 0.9), hjust = -0.2, vjust = 0.76, size = 6, color = '#F5B879') +
	theme_bw() +
	coord_flip() +
	theme(axis.title.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank(), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), legend.text = element_text(size = 14), axis.title.x = element_text(size = 14, color = "black"), plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
	labs(title = 'AP index and Disease', colour = "Pollutants", y = "Hazard Ratio") +
	guides(colour = guide_legend(reverse = TRUE)) +
	geom_stripped_cols(even = "#CCCCCC", alpha = .1)

ggsave(out_apscore_png, plot = p_apscore, width = 6.4, height = 9.6, dpi = 350)
ggsave(out_apscore_pdf, plot = p_apscore, device = "pdf", width = 6.4, height = 9.6)

# ------------------------------
# Unreported pollutants forest (two panels)
# ------------------------------

data2 <- read.csv(unreported_csv) %>% as.data.frame()
data2$exposure <- gsub('PM25_2010_IQR', 'PM2.5', data2$exposure)
data2$exposure <- gsub('PM10_mean_IQR', 'PM10', data2$exposure)
data2$exposure <- gsub('NO2_mean_IQR', 'NO2', data2$exposure)
data2$exposure <- gsub('Nitrogen_oxides_2010_IQR', 'NOx', data2$exposure)
colnames(data2)[1:3] <- c('value', 'lower', 'high')
data2$exposure <- data2$exposure %>% factor(levels = rev(c('PM2.5', 'PM2.5-10', 'PM10', 'NO2', 'NOx')))

diseaselist_unreported <- c('Schizophrenia', 'Sleep Disorder', 'PD', 'VD', 'Brain Cancer',
							 'Ischaemic Heart Disease', 'Hypertension', 'Arrhythmias', 'Anemia', 'Skin Disorder',
							 'Multiple Sclerosis', 'Peripheral Artery Disease', 'Lupus Erythematosus', 'Inflammatory Bowel Disease', 'Endometriosis',
							 'Liver Fibrosis', 'Obesity', 'Gynaecological Cancer', 'Colorectal Cancer', 'Prostate Cancer')

data2$Disease <- gsub('AllergicRhinitis', 'Allergic Rhinitis', data2$Disease)
data2$Disease <- gsub('BipolarDisorder', 'Bipolar Disorder', data2$Disease)
data2$Disease <- gsub('BreastCancer', 'Breast Cancer', data2$Disease)
data2$Disease <- gsub('HeartFailure', 'Heart Failure', data2$Disease)
data2$Disease <- gsub('IschamicStroke', 'Ischamic Stroke', data2$Disease)
data2$Disease <- gsub('LungCancer', 'Lung Cancer', data2$Disease)
data2$Disease <- gsub('BrainCancer', 'Brain Cancer', data2$Disease)
data2$Disease <- gsub('IschaemicHeartDisease', 'Ischaemic Heart Disease', data2$Disease)
data2$Disease <- gsub('LiverFibrosis', 'Liver Fibrosis', data2$Disease)
data2$Disease <- gsub('MultipleSclerosis', 'Multiple Sclerosis', data2$Disease)
data2$Disease <- gsub('SkinDisorder', 'Skin Disorder', data2$Disease)
data2$Disease <- gsub('PeripheralArteryDisease', 'Peripheral Artery Disease', data2$Disease)
data2$Disease <- gsub('SleepDisorder', 'Sleep Disorder', data2$Disease)
data2$Disease <- gsub('ColorectalCancer', 'Colorectal Cancer', data2$Disease)
data2$Disease <- gsub('LupusErythematosus', 'Lupus Erythematosus', data2$Disease)
data2$Disease <- gsub('GynaecologicalCancer', 'Gynaecological Cancer', data2$Disease)
data2$Disease <- gsub('ProstateCancer', 'Prostate Cancer', data2$Disease)
data2$Disease <- gsub('InflammatoryBowelDisease', 'Inflammatory Bowel Disease', data2$Disease)
data2$Disease <- gsub('ChronicKidneyDisease', 'Chronic Kidney Disease', data2$Disease)

APall_disease1 <- data2[data2$Disease %in% diseaselist_unreported[1:10], ]
APall_disease1$Disease <- factor(APall_disease1$Disease, levels = rev(diseaselist_unreported[1:10]))
APall_disease1$text_position <- APall_disease1$high + 0.015

p_u1 <- ggplot(APall_disease1, aes(x = Disease, y = value, group = as.factor(exposure), colour = as.factor(exposure))) +
	geom_hline(yintercept = 1, linetype = 2, color = "grey4", size = 1.15) +
	geom_point(position = position_dodge(width = 0.9), size = 3.2) +
	geom_errorbar(aes(ymax = lower, ymin = high), width = 0, alpha = 0.8, position = position_dodge(width = 0.9), size = 2) +
	ylim(0, 2) +
	theme_minimal() +
	geom_text(aes(y = text_position, label = ifelse((lower - 1) * (high - 1) > 0, '***', '')), position = position_dodge(width = 0.9), hjust = -0.2, vjust = 0.76, size = 6) +
	theme_bw() +
	coord_flip() +
	theme(axis.title.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank(), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14), axis.title.x = element_blank(), plot.title = element_text(size = 14), legend.position = 'none') +
	labs(colour = "Pollutants") +
	scale_colour_manual(values = c("#54A966", '#C55054', '#4D71AF', '#4C4C4C')) +
	guides(colour = guide_legend(reverse = TRUE)) +
	geom_stripped_cols(even = "#CCCCCC", alpha = .1)

APall_disease2 <- data2[data2$Disease %in% diseaselist_unreported[11:20], ]
APall_disease2$Disease <- factor(APall_disease2$Disease, levels = rev(diseaselist_unreported[11:20]))
APall_disease2$text_position <- APall_disease2$high + 0.015

p_u2 <- ggplot(APall_disease2, aes(x = Disease, y = value, group = as.factor(exposure), colour = as.factor(exposure))) +
	geom_hline(yintercept = 1, linetype = 2, color = "grey4", size = 1.15) +
	geom_point(position = position_dodge(width = 0.9), size = 3.2) +
	geom_errorbar(aes(ymax = lower, ymin = high), width = 0, alpha = 0.8, position = position_dodge(width = 0.9), size = 2) +
	ylim(0, 2) +
	theme_minimal() +
	geom_text(aes(y = text_position, label = ifelse((lower - 1) * (high - 1) > 0, '***', '')), position = position_dodge(width = 0.9), hjust = -0.2, vjust = 0.76, size = 6) +
	theme_bw() +
	coord_flip() +
	theme(axis.title.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank(), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14), axis.title.x = element_blank(), plot.title = element_text(size = 14)) +
	labs(colour = "Pollutants") +
	scale_colour_manual(values = c("#54A966", '#C55054', '#4D71AF', '#4C4C4C')) +
	guides(colour = guide_legend(reverse = TRUE)) +
	geom_stripped_cols(even = "#CCCCCC", alpha = .1)

p_unreported <- ggarrange(p_u1, p_u2, widths = c(1, 1.32))
p_unreported <- annotate_figure(p_unreported, top = text_grob("Disease Unreported", face = "bold", size = 14), bottom = text_grob("Hazard Ratio", face = "bold", vjust = 0.5, hjust = 0.5, size = 14))
ggsave(out_unreported_png, plot = p_unreported, width = 10, height = 8, dpi = 350)
ggsave(out_unreported_pdf, plot = p_unreported, device = "pdf", width = 10, height = 8)

# ------------------------------
# Heatmap (as in original Fig1 heatmap part)
# ------------------------------

data_hm <- read.xlsx(reported_xlsx, sheet = 5)
data_hm <- pivot_longer(data_hm[1:16, c(1, 3:6)], cols = !Incident.Disease, names_to = 'polltants', values_from = 'HR')
data_hm$polltants <- gsub('Hazard\\.Ratio\\[95\\%\\.CI\\]\\.', '', data_hm$polltants)
data_hm <- data_hm[data_hm$polltants != 'PM2.5,10', ]
data_hm$HR <- gsub('~1)', '~1.00)', data_hm$HR)
data_hm$HR <- gsub('\\(1~', '\\(1.00~', data_hm$HR)
extracted_numbers <- str_extract_all(data_hm$HR, "[0-9]+\\.[0-9]+")
extrat_df <- do.call(rbind, lapply(extracted_numbers, function(x) as.numeric(x)))
result <- data.frame(extrat_df)
colnames(result) <- c("value", "lower", "high")
data_hm$label <- ifelse((result$lower - 1) * (result$high - 1) > 0, '*', '')
data_hm$HR <- result$value %>% as.numeric()

unrep_hm <- fread(unreported_csv) %>% as.data.frame()
if (!all(c('value','lower','high') %in% colnames(unrep_hm)) && all(c('HR','Lowerlimit','Upperlimit') %in% colnames(unrep_hm))) {
	unrep_hm$value <- as.numeric(unrep_hm$HR)
	unrep_hm$lower <- as.numeric(unrep_hm$Lowerlimit)
	unrep_hm$high <- as.numeric(unrep_hm$Upperlimit)
}
if (!'label' %in% colnames(unrep_hm)) {
	if ('P_value' %in% colnames(unrep_hm)) unrep_hm$label <- ifelse(unrep_hm$P_value < 0.05, '*', '') else unrep_hm$label <- ifelse((unrep_hm$lower - 1) * (unrep_hm$high - 1) > 0, '*', '')
}
unrep_hm <- unrep_hm[, c('Disease', 'exposure', 'value', 'label')]
colnames(unrep_hm) <- c('Incident.Disease', 'polltants', 'HR', 'label')

ap_hm <- fread(apscore_csv) %>% as.data.frame()
colnames(ap_hm)[1:3] <- c('HR', 'lower', 'high')
ap_hm <- ap_hm[, c('Disease', 'exposure', 'HR', 'label')]
colnames(ap_hm) <- c('Incident.Disease', 'polltants', 'HR', 'label')

data_hm <- rbind(data_hm[, c('Incident.Disease', 'polltants', 'HR', 'label')], unrep_hm, ap_hm)

data_hm$polltants <- gsub('PM25_2010_IQR', 'PM2.5', data_hm$polltants)
data_hm$polltants <- gsub('PM10_mean_IQR', 'PM10', data_hm$polltants)
data_hm$polltants <- gsub('NO2_mean_IQR', 'NO2', data_hm$polltants)
data_hm$polltants <- gsub('Nitrogen_oxides_2010_IQR', 'NOx', data_hm$polltants)
data_hm$polltants <- gsub('ap_score_IQR', 'AP.Score', data_hm$polltants)

diseaselist_hm <- c('Depression','Schizophrenia','BipolarDisorder','SleepDisorder','PD','VD','AD','BrainCancer','Asthma',
					 'RA','MultipleSclerosis','PeripheralArteryDisease','LupusErythematosus','InflammatoryBowelDisease','Osteoporosis','AllergicRhinitis','Endometriosis',
					 'HeartFailure','IschamicStroke','IschaemicHeartDisease','Hypertension','Arrhythmias','Anemia',
					 'COPD','LungCancer','Pneumonia','SkinDisorder',
					 'T2D','NAFLD','ChronicKidneyDisease','LiverFibrosis','Obesity',
					 'BreastCancer','GynaecologicalCancer','ColorectalCancer','ProstateCancer')

data_hm$Incident.Disease <- factor(data_hm$Incident.Disease, levels = diseaselist_hm)
data_hm$polltants <- factor(data_hm$polltants, levels = rev(c('PM2.5','PM10','NO2','NOx','AP.Score')))
data_hm <- na.omit(data_hm)

p_hm <- ggplot(data = data_hm, aes(x = Incident.Disease, y = polltants)) +
	geom_tile(aes(fill = HR), color = 'white', size = 0.25) +
	scale_fill_gradientn(colours = c("#053264", "white",  "#820823"), values = rescale(c(0.83, 1, 2)), guide = "colorbar", breaks = c(0.5, 1, 1.5)) +
	geom_text(aes(label = label), size = 6, vjust = 0.8) +
	scale_y_discrete(position = "left") +
	theme(axis.text.x = element_text(angle = 30, hjust=1, size = 11, face = 'bold'), axis.text.y = element_text(size = 11, face = 'bold'), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(0.25, 0.25, 0.25, 1, "cm"), legend.key.size = unit(0.7, "cm"), legend.text = element_text(size = 10))

ggsave(out_heatmap_png, plot = p_hm, width = 13, height = 2.9, dpi = 350)
ggsave(out_heatmap_pdf, plot = p_hm, device = "pdf", width = 13, height = 2.9)


