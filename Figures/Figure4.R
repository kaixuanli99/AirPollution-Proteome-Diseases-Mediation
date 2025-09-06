############################################################

## Dependencies ---------------------------------------------------
library(ggplot2)
library(dplyr)
library(data.table)
library(ggpubr)     # for ggarrange, annotate_figure, stat_cor
library(ggsignif)   # for geom_signif

## Paths (parameterized) -------------------------------------------
input_dir  <- "<YOUR_INPUT_DIR>"      # e.g., "/path/to/data"
output_dir <- "<YOUR_OUTPUT_DIR>"     # e.g., "/path/to/output"

## Load required RData (expects UKB_RA_df_clean, result_df) --------
load(file.path(input_dir, "fig4.RData"))


############################################################
# Block A: QQ-style percentile plots per APPRS
############################################################

# Disease/APPRS lists
disease_list_name <- c('COPD','Pneumonia',
                       'RA','LupusErythematosus','Osteoporosis',
                       'T2D','NAFLD','CKD',
                       'IHD','Arrhythmias','HeartFailure','PeripheralArteryDisease','Anemia')

# APPRS columns order derived from result_df columns
apprs <- colnames(result_df)[c(2,12,5,4,13,6,11,10,3,1,7,8,9)]

color_list <- c('#B55D60','#B55D60','#5875A4','#5875A4','#9C7BA5',
                '#5F9E6E','#5F9E6E','#5F9E6E',
                '#9C7BA5','#9C7BA5','#9C7BA5','#9C7BA5','#9C7BA5')

# Build plots for k in c(1:4,6:13)
plot_indices <- c(1:4,6:13)
plots <- vector("list", length(plot_indices))
names(plots) <- paste0("p", plot_indices)

for(idx in seq_along(plot_indices)) {
  k <- plot_indices[idx]
  disease <- disease_list_name[k]
  disease_prs <- apprs[k]

  plotdf_female <- UKB_RA_df_clean[UKB_RA_df_clean$sex==0, c('eid','sex','age',disease,disease_prs)] %>% na.omit()
  plotdf_male   <- UKB_RA_df_clean[UKB_RA_df_clean$sex==1, c('eid','sex','age',disease,disease_prs)] %>% na.omit()
  colnames(plotdf_female)[5] <- 'ProRS'
  colnames(plotdf_male)[5]   <- 'ProRS'
  colnames(plotdf_female)[4] <- 'Y'
  colnames(plotdf_male)[4]   <- 'Y'
  colnames(plotdf_female)[3] <- 'Age'
  colnames(plotdf_male)[3]   <- 'Age'

  ob_female_lst <- ob_male_lst <- ob_female_age_lst <- ob_male_age_lst <- numeric(100)
  for (i in seq(0, 99)) {
    l_cut_female <- quantile(plotdf_female$ProRS, i/100)
    u_cut_female <- quantile(plotdf_female$ProRS, (i+1)/100)
    tmp_female <- plotdf_female %>% dplyr::filter(ProRS > l_cut_female & ProRS <= u_cut_female)
    ob_female_lst[i+1] <- sum(tmp_female$Y) / nrow(tmp_female)
    ob_female_age_lst[i+1] <- mean(tmp_female$Age)

    l_cut_male <- quantile(plotdf_male$ProRS, i/100)
    u_cut_male <- quantile(plotdf_male$ProRS, (i+1)/100)
    tmp_male <- plotdf_male %>% dplyr::filter(ProRS > l_cut_male & ProRS <= u_cut_male)
    ob_male_lst[i+1] <- sum(tmp_male$Y) / nrow(tmp_male)
    ob_male_age_lst[i+1] <- mean(tmp_male$Age)
  }

  scatter_df <- data.frame(ob_prop = c(ob_female_lst, ob_male_lst),
                           pro_q_prop = rep(seq(0, 99), 2),
                           gender = factor(rep(c("Female", "Male"), each = 100)),
                           age = c(ob_female_age_lst, ob_male_age_lst))

  p_title <- disease_list_name[k]
  p_title <- gsub('PeripheralArteryDisease','PAD', p_title)
  p_title <- gsub('LupusErythematosus','LE', p_title)
  label.y <- ifelse(k %in% c(3,4,7), 0.04, 0.15)

  p <- ggplot(scatter_df, aes(x = pro_q_prop, y = ob_prop, size = age, alpha = age)) +
    geom_point(color = color_list[k]) +
    labs(title = p_title) +
    scale_alpha_continuous(range = c(0.1, 1)) +
    scale_size(range = c(0.6, 1.8)) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          plot.title = element_text(size = 14),
          legend.position = "none",
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm")) +
    stat_cor(method = "pearson", label.x = 3, label.y = label.y)

  plots[[idx]] <- p
}

p_all <- ggarrange(plotlist = plots, ncol = 4, nrow = 3, widths = c(1,1,1,1))
p_all <- annotate_figure(p_all,
                         left = text_grob("Observed event rate", face = "bold", rot = 90, vjust = 1, size = 14),
                         bottom = text_grob("APPRS Percentile (%)", face = "bold", vjust = 0.5, hjust = 0.5, size = 14))

ggsave(file.path(output_dir, 'qq_percentile_all.pdf'), plot = p_all, device = 'pdf', width = 12, height = 8)
ggsave(file.path(output_dir, 'qq_percentile_all.tiff'), plot = p_all, device = 'tiff', width = 12, height = 8, dpi = 350)


############################################################
# Block B: Violin+box plots (APPRS in Cases vs Normal)
############################################################

tmp_df <- NULL
combine_result <- function(tmp_df, outcome, apprs) {
  tmp <- UKB_RA_df_clean[, c(outcome, apprs)]
  colnames(tmp) <- c('Disease','ApPRS')
  tmp$Disease <- gsub('1','Cases', tmp$Disease)
  tmp$Disease <- gsub('0','Normal', tmp$Disease)
  tmp$Disease <- factor(tmp$Disease)
  tmp <- na.omit(tmp)
  tmp$Group <- outcome
  tmp$Group <- gsub('_2010_2020','', tmp$Group)
  tmp$Group <- gsub('LupusErythematosus','LE', tmp$Group)
  tmp$Group <- gsub('PeripheralArteryDisease','PAD', tmp$Group)
  tmp_df <- rbind(tmp_df, tmp)
  return(tmp_df)
}

diseaselist <- c('HeartFailure','T2D','IHD','CKD','Anemia','RA','Arrhythmias',
                 'COPD','LupusErythematosus','PeripheralArteryDisease','NAFLD','Pneumonia')
apprslist <- c('ApPRS_HeartFailure','ApPRS_T2D','ApPRS_IHD','ApPRS_CKD','ApPRS_Anemia','ApPRS_RA','ApPRS_Arrhythmias',
               'ApPRS_COPD','ApPRS_LupusErythematosus','ApPRS_PeripheralArteryDisease','ApPRS_NAFLD','ApPRS_Pneumonia')

for(k in 1:12) {
  tmp_df <- combine_result(tmp_df, diseaselist[k], apprslist[k])
}

tmp_df$Group <- factor(tmp_df$Group)
tmp_df$Disease <- factor(tmp_df$Disease, levels = c('Normal','Cases'))

# significance bars positions
xmin <- seq(0,11) + 0.8
xmax <- seq(1,12) + 0.2
max_apprs_per_group <- tmp_df %>% group_by(Group) %>% summarise(max_apprs = max(ApPRS))
y_position <- max_apprs_per_group$max_apprs + 0.3

p_cases <- ggplot(tmp_df, aes(x = Group, y = ApPRS, fill = Disease, color = Disease)) +
  geom_violin(trim = FALSE, size = 0.8, alpha = 0.5) +
  geom_boxplot(width = 0.2, alpha = 1, outlier.size = 0.00001, position = position_dodge(0.88), color = "black") +
  labs(y = 'APPRS', title = 'UKB') +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 12, face = "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 12, angle = 0, vjust = 1, hjust = 1),
        axis.text.x = element_text(size = 12, angle = 30, vjust = 1, hjust = 1),
        axis.title.y = element_text(size = 12, face = 'bold')) +
  ylim(min(tmp_df$ApPRS) - 0.05, max(tmp_df$ApPRS) + 0.4) +
  scale_color_manual(name = "Disease", values = c("Normal" = "#5875A4", "Cases" = "#B55D60")) +
  scale_fill_manual(name = "Disease", values = c("Normal" = "#5875A4", "Cases" = "#B55D60")) +
  geom_signif(y_position = y_position, xmin = xmin, xmax = xmax, annotation = rep('***', 12), tip_length = .01, color = 'black')

ggsave(file.path(output_dir, 'boxplot-UKB-Each.png'), plot = p_cases, width = 11, height = 4, dpi = 350)
ggsave(file.path(output_dir, 'boxplot-UKB-Global.png'), plot = p_cases, width = 11, height = 4, dpi = 350)

############################################################
# Done
############################################################




############################################################
# Block C (from Replication.R L265–523): Replication panels
# p_all = ggarrange(p_pan, p_vo2_arr, p_vo2_pad)
############################################################

library(tidyr)

# Prepare data_wide_acme (coefficients) from mediation_factor_df_FDR.csv ---------
med_df <- fread(file.path(input_dir, 'mediation_factor_df_FDR.csv')) %>% as.data.frame()
med_df <- med_df[med_df$ACME_FDR < 0.05 & med_df$Total_FDR < 0.05 & med_df$Prop_mediated_ave > 0.04, ]
ap_df <- med_df[med_df$exposure == 'ap_score_healthy_stat', ]

ap_freq <- ap_df$mediation_factor %>% table() %>% as.data.frame()
ap_freq <- ap_freq[ap_freq$Freq >= 2, ]
ap_df <- ap_df[ap_df$mediation_factor %in% ap_freq$. , ]
ap_df <- ap_df %>% group_by(mediation_factor) %>% mutate(sum_proportion = sum(Prop_mediated_ave))
ord_vec <- ap_df %>% dplyr::select(mediation_factor, sum_proportion) %>% distinct() %>% arrange(desc(sum_proportion)) %>% pull(mediation_factor)
disease_ord <- ap_df$Disease %>% unique()

build_df <- med_df %>% dplyr::filter(mediation_factor %in% ord_vec, Disease %in% disease_ord, exposure == 'ap_score_healthy_stat') %>%
  dplyr::select(mediation_factor, Disease, ACME)

data_wide_acme <- build_df %>%
  pivot_wider(names_from = mediation_factor, values_from = ACME) %>% as.data.frame()
data_wide_acme[is.na(data_wide_acme)] <- 0
row.names(data_wide_acme) <- data_wide_acme$Disease
data_wide_acme <- data_wide_acme[,-1]
data_wide_acme <- data_wide_acme[, ord_vec, drop = FALSE]
data_wide_acme <- t(data_wide_acme) %>% as.data.frame()   # rows: proteins, cols: diseases


# C1. Yale replication: build p_pan --------------------------------------------
yale_raw <- fread(file.path(input_dir, 'proteome-yale.csv')) %>% as.data.frame()
yale_grp <- yale_raw %>%
  group_by(EntrezGeneSymbol) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = 'drop') %>%
  as.data.frame()
row.names(yale_grp) <- yale_grp$EntrezGeneSymbol
yale_mat <- yale_grp[, -1, drop = FALSE] %>% t() %>% as.data.frame()
yale_mat <- scale(yale_mat) %>% as.data.frame()

co_pro_yale <- intersect(row.names(data_wide_acme), colnames(yale_mat))
data_wide_acme_yale <- as.matrix(data_wide_acme[co_pro_yale, , drop = FALSE])
yale_mat_use <- as.matrix(yale_mat[, co_pro_yale, drop = FALSE])

yale_ApPRS <- yale_mat_use %*% data_wide_acme_yale   # samples × diseases
for(i in seq_len(ncol(yale_ApPRS))) {
  denom <- sum(data_wide_acme_yale[, i], na.rm = TRUE)
  if (denom != 0) yale_ApPRS[, i] <- yale_ApPRS[, i] / denom
}

# Aggregate sample groups (requires sample columns with these labels)
AutoInflammatoryDiseases <- c('PAPA','TRAPS','DADA2','FMF','FCAS','HIDS','Muckle-Wells')
Telo <- c('TERT','TERC')
PrimaryImmunodeficiencies <- c('p47-CGD','CARD14 DN','CTLA4','GATA2','STAT3 DN','LAD1','NEMO','NEMO carrier',
                               'PGM3','APDS1','STAT1 GOF','TERC','X-CGD')

df_with_rownames <- colnames(yale_ApPRS)
yale_agg <- t(yale_ApPRS) %>% as.data.frame() %>%
  rowwise() %>%
  mutate(
    AID_ApPRS_mean = mean(c_across(all_of(AutoInflammatoryDiseases)), na.rm = TRUE),
    Telo_ApPRS_mean = mean(c_across(all_of(Telo)), na.rm = TRUE),
    PID_ApPRS_mean = mean(c_across(all_of(PrimaryImmunodeficiencies)), na.rm = TRUE)
  ) %>%
  ungroup() %>% as.data.frame()
row.names(yale_agg) <- df_with_rownames

# Expect a 'Healthy' column in the sample set
yale_show <- yale_agg[, c('Healthy','AID_ApPRS_mean','Telo_ApPRS_mean','PID_ApPRS_mean')]
yale_show$Disease <- row.names(yale_show)
yale_long <- yale_show %>% pivot_longer(cols = c('Healthy','AID_ApPRS_mean','Telo_ApPRS_mean','PID_ApPRS_mean'), names_to = 'Group')
yale_long$Group <- gsub('AID_ApPRS_mean','AID',yale_long$Group)
yale_long$Group <- gsub('Telo_ApPRS_mean','Telo',yale_long$Group)
yale_long$Group <- gsub('PID_ApPRS_mean','PID',yale_long$Group)
yale_long$Group <- factor(yale_long$Group, levels = c('Healthy','AID','PID','Telo'))

p_pan <- ggplot(yale_long, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot(width = 0.55, alpha = 1, outlier.shape = NA, position = position_dodge(0.88), color = 'black') +
  labs(title = 'Panmonogenic', y = 'Global APPRS') +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    plot.title = element_text(size = 12, face = 'bold'),
    text = element_text(size = 8),
    axis.title = element_text(face = 'bold'),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 10, face = 'bold')
  ) +
  geom_signif(y_position = c(0.01, 0.095, 0.34), xmin = c(1, 1, 1), xmax = c(2, 3, 4), annotation = rep('***', 3), tip_length = .01, color = 'black') +
  scale_fill_manual(values = c('#B55D60', '#5875A4', '#9C7BA5', '#5F9E6E'))


# C2. HERITAGE replication: build p_vo2_arr, p_vo2_pad --------------------------
heritage_proteome <- fread(file.path(input_dir, 'heritage_proteome.csv')) %>% as.data.frame()
co_pro_h <- intersect(colnames(heritage_proteome)[10:ncol(heritage_proteome)], row.names(data_wide_acme))

heritage_proteome[,10:ncol(heritage_proteome)] <- scale(heritage_proteome[,10:ncol(heritage_proteome)])
colnames(heritage_proteome)[6] <- 'Baseline_VO2max'
hp_df <- heritage_proteome[, c('Baseline_VO2max','BMI', co_pro_h)]

# Quantile groups for VO2max (q3: tertiles; q2: median split)
q <- quantile(hp_df$Baseline_VO2max, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
hp_df$Baseline_VO2max_q3 <- cut(hp_df$Baseline_VO2max, breaks = unique(q), include.lowest = TRUE, labels = FALSE)
med_q <- median(hp_df$Baseline_VO2max, na.rm = TRUE)
hp_df$Baseline_VO2max_q2 <- ifelse(hp_df$Baseline_VO2max <= med_q, 1, 2)

mat_expr <- as.matrix(hp_df[, co_pro_h, drop = FALSE])
mat_coef <- as.matrix(data_wide_acme[co_pro_h, , drop = FALSE])
Heritage_ApPRS <- mat_expr %*% mat_coef %>% as.data.frame()
for(i in seq_len(ncol(Heritage_ApPRS))) {
  denom <- sum(mat_coef[, i], na.rm = TRUE)
  if (denom != 0) Heritage_ApPRS[, i] <- Heritage_ApPRS[, i] / denom
}

Heritage_ApPRS <- Heritage_ApPRS[, 1:min(12, ncol(Heritage_ApPRS)), drop = FALSE]
Heritage_ApPRS$Global_ApPRS <- rowMeans(Heritage_ApPRS, na.rm = TRUE)
Heritage_ApPRS$BMI <- hp_df$BMI
Heritage_ApPRS$Baseline_VO2max_q3 <- hp_df$Baseline_VO2max_q3
Heritage_ApPRS$Baseline_VO2max_q2 <- hp_df$Baseline_VO2max_q2

tmp_df_h <- Heritage_ApPRS[, c('Global_ApPRS','Baseline_VO2max_q3')]
tmp_df_h$VO2max <- ifelse(tmp_df_h$Baseline_VO2max_q3 == 1, 'Low CRF', ifelse(tmp_df_h$Baseline_VO2max_q3 == 2, 'Mid CRF', ifelse(tmp_df_h$Baseline_VO2max_q3 == 3, 'High CRF', NA)))
tmp_df_h$VO2max <- factor(tmp_df_h$VO2max, levels = c('High CRF','Mid CRF','Low CRF'))

# Extract Arrhythmias and PAD columns if present
arr_col <- grep('Arrhythmias', colnames(Heritage_ApPRS), value = TRUE)[1]
pad_col <- grep('PeripheralArteryDisease|PAD', colnames(Heritage_ApPRS), value = TRUE)[1]

plot_df <- data.frame(
  VO2max = tmp_df_h$VO2max,
  Arrhythmias = if (!is.na(arr_col)) Heritage_ApPRS[[arr_col]] else NA,
  PAD = if (!is.na(pad_col)) Heritage_ApPRS[[pad_col]] else NA
)

p_vo2_arr <- ggplot(plot_df, aes(x = VO2max, y = Arrhythmias, fill = VO2max)) +
  geom_boxplot(width = 0.55, alpha = 1, outlier.shape = NA, position = position_dodge(0.88), color = 'black') +
  labs(y = 'Arrhythmias APPRS', title = 'HERITAGE') +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    plot.title = element_text(size = 12, face = 'bold'),
    text = element_text(size = 8),
    axis.title = element_text(face = 'bold'),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 10, face = 'bold'),
    legend.position = 'none'
  ) +
  ylim(-0.7, 1) +
  scale_fill_manual(values = c('#5875A4','#B55D60','#5F9E6E')) +
  geom_signif(y_position = c(0.65, 0.85), xmin = c(1, 1), xmax = c(2, 3), annotation = rep('***', 2), tip_length = .01, color = 'black')

p_vo2_pad <- ggplot(plot_df, aes(x = VO2max, y = PAD, fill = VO2max)) +
  geom_boxplot(width = 0.55, alpha = 1, outlier.shape = NA, position = position_dodge(0.88), color = 'black') +
  labs(y = 'PAD APPRS', title = 'HERITAGE') +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    plot.title = element_text(size = 12, face = 'bold'),
    text = element_text(size = 8),
    axis.title = element_text(face = 'bold'),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title.y = element_text(size = 10, face = 'bold')
  ) +
  ylim(-0.7, 1) +
  scale_fill_manual(values = c('#5875A4','#B55D60','#5F9E6E')) +
  geom_signif(y_position = c(0.65, 0.75), xmin = c(1, 1), xmax = c(2, 3), annotation = rep('***', 2), tip_length = .01, color = 'black')

# Compose p_all
p_all <- ggarrange(p_pan, p_vo2_arr, p_vo2_pad, nrow = 1, ncol = 3, widths = c(1.25, 1, 1.2), labels = c('A','B','C'), align = 'h')

ggsave(file.path(output_dir, 'Replication.png'), plot = p_all, width = 11, height = 3.5, dpi = 350)
ggsave(file.path(output_dir, 'Figure4.pdf'), plot = p_all, device = 'pdf', width = 11, height = 3.5)
ggsave(file.path(output_dir, 'Figure4.tiff'), plot = p_all, device = 'tiff', width = 11, height = 3.5, dpi = 350)
