
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(ggpubr)


input_dir  <- "."
output_dir <- "."

# Helper to draw AUC compare scatter
draw_auc_scatter <- function(plot_list_data, model1_label, model2_label, disease_order, disease_colors) {
  auc_results <- lapply(plot_list_data, function(rl) sapply(rl, auc))
  auc_df <- as.data.frame(do.call(cbind, auc_results))
  colnames(auc_df) <- disease_order
  auc_df$Models <- c('Model1','Model2')
  auc_df_long <- t(auc_df[, disease_order]) %>% as.data.frame()
  auc_df_long$Disease <- disease_order
  auc_df_long$Model1 <- as.numeric(auc_df_long$Model1)
  auc_df_long$Model2 <- as.numeric(auc_df_long$Model2)
  auc_df_long$Disease <- factor(auc_df_long$Disease, levels = rev(disease_order))

  p <- ggplot(auc_df_long, aes(x = Model2, y = Model1, color = Disease)) +
      geom_point(size = 5) +
      geom_abline(intercept = 0, slope = 1, linetype = 'solid', color = 'black') +
      geom_hline(yintercept = 0.6, linetype = 'dashed', color = 'darkred') +
      geom_vline(xintercept = 0.6, linetype = 'dashed', color = 'darkred') +
      labs(x = model2_label, y = model1_label, title = '') +
      theme_classic() +
      theme(
        plot.title = element_text(size = 18, face = 'bold'),
        axis.title = element_text(face = 'bold'),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16, face = 'bold'),
        axis.title.x = element_text(size = 16, face = 'bold'),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.9, 'cm')
      ) +
      scale_color_manual(values = disease_colors) +
      labs(tag = '')
  list(plot = p, data = auc_df_long)
}

############################################################
# Inputs required for B–D
# Expect 3 RData files each defines plot_list_data
############################################################

load(file.path(input_dir, 'LASSO_basic.RData'))            # defines plot_list_data
plot_list_data_basic <- plot_list_data
load(file.path(input_dir, 'LASSO_clinical.RData'))         # defines plot_list_data
plot_list_data_clinical <- plot_list_data
load(file.path(input_dir, 'LASSO_basic_clinical.RData'))   # defines plot_list_data
plot_list_data_basic_clinical <- plot_list_data

disease_order <- c('HeartFailure','CKD','IHD','T2D','Anemia','RA','Arrhythmias','COPD','Pneumonia','NAFLD','PAD','LE')
disease_colors <- c(
  'Anemia' = '#9C7BA5','HeartFailure' = '#A88BB0','PAD' = '#B49CBB','IHD' = '#C1ACC6','Arrhythmias' = '#CDBDD2',
  'CKD' = '#5F9E6E','T2D' = '#73AA80','NAFLD' = '#87B692','RA' = '#5875A4','LE' = '#6C86AF','Pneumonia' = '#B55D60','COPD' = '#BE7173'
)

# B: Basic vs Basic+APPRS
res_B <- draw_auc_scatter(plot_list_data_basic, 'Basic', 'Basic + APPRS', disease_order, disease_colors)
pB <- res_B$plot + labs(tag = 'B') + theme(plot.tag = element_text(size = 22, face = 'bold'), plot.tag.position = c(0,1))
ggsave(file.path(output_dir, 'Figure5_B.pdf'), pB, width = 7, height = 6)

# C: Clinical vs Clinical+APPRS
res_C <- draw_auc_scatter(plot_list_data_clinical, 'Clinical', 'Clinical + APPRS', disease_order, disease_colors)
pC <- res_C$plot + labs(tag = 'C') + theme(plot.tag = element_text(size = 22, face = 'bold'), plot.tag.position = c(0,1))
ggsave(file.path(output_dir, 'Figure5_C.pdf'), pC, width = 7, height = 6)

# D: Basic+Clinical vs +APPRS
res_D <- draw_auc_scatter(plot_list_data_basic_clinical, 'Basic + Clinical', 'Basic + Clinical + APPRS', disease_order, disease_colors)
pD <- res_D$plot + labs(tag = 'D') + theme(plot.tag = element_text(size = 22, face = 'bold'), plot.tag.position = c(0,1))
ggsave(file.path(output_dir, 'Figure5_D.pdf'), pD, width = 7, height = 6)

############################################################
# E: Contribution stacked bar (Model3)
############################################################

load(file.path(input_dir, 'LASSO_coef_df_result.RData'))  # expects coef_df
coef_df <- coef_df %>% filter(Model == 'model3', Variable != '(Intercept)')
coef_df <- coef_df[coef_df$Disease != 'Osteoporosis',]
coef_df$Variable <- gsub('ApPRS','APPRS', coef_df$Variable)
coef_df$Variable <- gsub('ApoB_ApoA','ApoB/ApoA', coef_df$Variable)
coef_df$Variable <- gsub('sex','Sex', coef_df$Variable)
coef_df$Variable <- gsub('age','Age', coef_df$Variable)
coef_df$Variable <- gsub('smoking','Smoking', coef_df$Variable)
coef_df$Variable <- gsub('alcohol','Alcohol', coef_df$Variable)
coef_df$Disease <- gsub('_2010_2020','', coef_df$Disease)
coef_df$Disease <- gsub('PeripheralArteryDisease','PAD', coef_df$Disease)
coef_df$Disease <- gsub('MultipleSclerosis','MS', coef_df$Disease)
coef_df$Disease <- gsub('LupusErythematosus','LE', coef_df$Disease)
coef_df$Variable <- factor(coef_df$Variable, levels = c('APPRS','ApoB/ApoA','Sex','Age','Smoking','Alcohol','BMI','Albumin','Creatinine','HbA1c','Glucose'))

palette_E <- colorRampPalette(c('#B55D60','#DAC060','#5F9E6E','#5875A4','#9C7BA5'))(11)

pE <- ggplot(coef_df, aes(y = Disease, x = abs(Coefficient))) +
  geom_col(aes(fill = Variable), position = 'fill') +
  theme_classic() +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 16, face = 'bold'),
        legend.title = element_text(size = 16),
        plot.margin = unit(c(0.3,0.3,0.3,0.3), 'cm')) +
  labs(y = '', x = 'Contribution in LASSO Prediction model', fill = 'Variables', tag = 'E') +
  scale_fill_manual(values = palette_E) +
  theme(plot.tag = element_text(size = 22, face = 'bold'), plot.tag.position = c(0,1))
ggsave(file.path(output_dir, 'Figure5_E.pdf'), pE, width = 10, height = 6)

############################################################
# F: Replication AUC grouped bars (Arrhythmias, PAD)
############################################################

load(file.path(input_dir, 'LASSO_Replication_result.RData'))  # expects plot_list_data
load(file.path(input_dir, 'Heritage_ApPRS.RData'))            # expects Heritage_ApPRS

auc_results <- lapply(plot_list_data, function(rl) sapply(rl, auc))
auc_df <- as.data.frame(do.call(cbind, auc_results))
colnames(auc_df) <- colnames(Heritage_ApPRS)[1:13]
auc_df <- auc_df[, c(1,8)]
auc_df$Models <- c('Model1','Model2')
colnames(auc_df) <- c('Arrhythmias','PAD','Models')
result_auc_longer <- auc_df %>% pivot_longer(cols = c('Arrhythmias','PAD'), names_to = 'Disease')
result_auc_longer$value <- as.numeric(result_auc_longer$value)

pF <- ggplot(result_auc_longer, aes(x = Disease, y = value-0.5, fill = Models)) +
  geom_col(position = position_dodge(width = 0.5), width = 0.5) +
  scale_fill_manual(values = c('#5F9E6E','#5875A4'), labels = c('Basic','Basic + APPRS')) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.position = 'top',
        plot.margin = unit(c(0.3,0.3,0.3,0.3), 'cm')) +
  labs(y = 'AUC', title = '', tag = 'F') +
  scale_y_continuous(limits = c(0,0.2), breaks = seq(0,0.2,0.1), labels = c('0.5','0.6','0.7')) +
  theme(plot.tag = element_text(size = 22, face = 'bold'), plot.tag.position = c(0,1))
ggsave(file.path(output_dir, 'Figure5_F.pdf'), pF, width = 6, height = 6)

############################################################
# G: Replication contribution stacked bars (Arrhythmias, PAD)
############################################################

load(file.path(input_dir, 'LASSO_coef_df_Replication.RData'))  # expects coef_df
coef_df <- coef_df %>% filter(Model == 'model2', Variable != '(Intercept)')
coef_df <- coef_df[c(1:5,36:40),]
coef_df$Variable[c(5,10)] <- 'APPRS'
coef_df$Disease <- rep(c('Arrhythmias','PAD'), each = 5)
coef_df$Variable <- factor(coef_df$Variable, levels = c('APPRS','Age','Sex','Race','BMI'))

palette_G <- c('#9E0142','#FEE08B','#ABDDA4','#66C2A5','#3288BD')

pG <- ggplot(coef_df, aes(y = abs(Coefficient), x = Disease)) +
  geom_col(aes(fill = Variable), position = 'fill', width = 0.6) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title.x = element_text(size = 18),
        legend.title = element_text(size = 16),
        axis.title.y = element_text(size = 16, face = 'bold'),
        plot.margin = unit(c(0.3,0.3,0.3,0.3), 'cm')) +
  labs(y = 'Contribution in LASSO Prediction model', x = '', fill = 'Disease', tag = 'G') +
  scale_fill_manual(values = palette_G) +
  theme(plot.tag = element_text(size = 22, face = 'bold'), plot.tag.position = c(0,1))
ggsave(file.path(output_dir, 'Figure5_G.pdf'), pG, width = 5, height = 6)


