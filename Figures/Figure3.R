
## Dependencies ---------------------------------------------------
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(circlize)
library(grid)
library(gridExtra)
library(openxlsx)

## Paths (parameterized) -------------------------------------------
input_dir  <- "."      # e.g., "/path/to/data"
output_dir <- "."     # e.g., "/path/to/output"


############################################################
# C. Chord diagram (AP index mediators ↔ diseases)
############################################################

# Read and filter data 
data <- fread(file.path(input_dir, "mediation_factor_df_FDR.csv")) %>% as.data.frame()
data <- data[data$Prop_mediated_ave >= 0.04 & data$ACME_FDR < 0.05 & data$Total_FDR < 0.05 & data$Prop_mediated_ave_FDR < 0.05, ]
colnames(data)[1:2] <- c("mediation_factor", "exposure")

# Keep AP index only
dat <- data[data$exposure == 'ap_score_healthy_stat', c("mediation_factor", "Disease", "Prop_mediated_ave")]

# Harmonize disease names
dat$Disease <- gsub('PeripheralArteryDisease','PAD', dat$Disease)
dat$Disease <- gsub('LupusErythematosus','LE', dat$Disease)
dat$Disease <- gsub('HeartFailure','Heart Failure', dat$Disease)

# Disease order 
disease_order <- c('COPD','Pneumonia', 'IHD','Osteoporosis','PAD','Arrhythmias','Heart Failure','Anemia',
                   'LE','RA','T2D','CKD','NAFLD')
dat$Disease <- factor(dat$Disease, levels = disease_order)
dat <- dat[order(dat$Disease),]
dat$Disease <- as.vector(dat$Disease)

# Disease system groups
lung <- c('COPD','LungCancer','Pneumonia')
cardio <- c('Arrhythmias','Hypertension','IschamicStroke','Heart Failure','Anemia','Osteoporosis','PAD')
immune <- c('LE','RA','Asthma')
neuron <- c('DP','Schizophrenia','AD')
metabolic <- c('T2D','Liverfibrosis','CKD','Obesity','NAFLD')
skin <- c('SkinDisorder')

# Grouping and colors (same logic)
name <- unique(c(dat$mediation_factor, dat$Disease))
group <- rep(3, length(name))
group[which(name %in% unique(dat$mediation_factor))] <- 1
group[which(name %in% unique(dat$Disease))] <- 2
names(group) <- name

grid.col <- rep("#5F9E6E", length(name))
grid.col[which(name %in% unique(dat$mediation_factor))] <- "grey4"
grid.col[which(name %in% lung)] <- "#B55D60"
grid.col[which(name %in% cardio)] <- "#9C7BA5"
grid.col[which(name %in% immune)] <- "#5875A4"
grid.col[which(name %in% neuron)] <- "#DAC060"
grid.col[which(name %in% metabolic)] <- "#5F9E6E"
grid.col[which(name %in% skin)] <- "#595959"
names(grid.col) <- name

# Link colors
dat$col <- "#CCCCCC"
dat[dat$Disease %in% lung,]$col <- '#B55D60'
dat[dat$Disease %in% cardio,]$col <- '#9C7BA5'
dat[dat$Disease %in% immune,]$col <- '#5875A4'
dat[dat$Disease %in% neuron,]$col <- '#DAC060'
dat[dat$Disease %in% metabolic,]$col <- '#5F9E6E'
dat[dat$Disease %in% skin,]$col <- '#595959'

# Plot and export
pdf(file.path(output_dir, 'Fig3_C.pdf'), width = 11, height = 11)
par(mar = c(0.5,0.5,0.5,0.5))
chordDiagram(dat, transparency = 0, group = group, grid.col = grid.col,
             big.gap = 5, small.gap = 1.5, annotationTrack = "grid", col = dat$col,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(dat)))))
)
circos.track(
  track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.85)
  }, bg.border = NA
)
mtext("C", side = 3, adj = 0, line = -2, cex = 1.8, font = 2)
circos.clear()
dev.off()


############################################################
# D. Percentage stacked bar (blue-yellow; specific/shared)
############################################################

data <- fread(file.path(input_dir, "mediation_factor_df_FDR.csv")) %>% as.data.frame()
data <- data[data$ACME_FDR < 0.05 & data$Total_FDR < 0.05 & data$Prop_mediated_ave > 0.04, ]
data <- data[data$Disease != 'Osteoporosis', ]
data <- data[data$exposure != 'ap_score_healthy_stat', ]
data$Disease <- gsub('PeripheralArteryDisease','PAD', data$Disease)
data$Disease <- gsub('LupusErythematosus','LE', data$Disease)

# Count per disease-protein across exposures
protein_exposure_count <- data %>%
  group_by(mediation_factor, Disease) %>%
  summarise(exposures = list(unique(exposure)),
            exposure_count = n_distinct(exposure), .groups = 'drop')

# Specific proteins (appear in one exposure only)
specific_proteins <- protein_exposure_count %>%
  filter(exposure_count == 1) %>%
  mutate(exposure = unlist(exposures))

# Specific protein count per disease-exposure
specific_counts <- specific_proteins %>%
  group_by(Disease, exposure) %>%
  summarise(specific_protein_count = n(), .groups = 'drop') %>%
  spread(exposure, specific_protein_count, fill = 0)

# Shared proteins (appear in ≥2 exposures)
shared_proteins <- protein_exposure_count %>%
  filter(exposure_count >= 2) %>%
  group_by(Disease) %>%
  summarise(shared_protein_count = n())

# Merge
final_result <- specific_counts %>%
  left_join(shared_proteins, by = "Disease") %>%
  replace_na(list(shared_protein_count = 0))

# Standardize column names (for mapping/legend)
colnames(final_result) <- c('Disease', "nox", "no2", "pm10", "pm2.5", "shared_protein_count")

# Build plot data
plot_data <- final_result %>%
  mutate(total = pm2.5 + pm10 + no2 + nox + shared_protein_count) %>%
  pivot_longer(cols = c("no2", "nox", "pm2.5", "pm10", "shared_protein_count"),
               names_to = "exposure", values_to = "count") %>%
  group_by(Disease) %>%
  mutate(freq = count / sum(count))

# Disease order: descending by total
disease_order <- plot_data %>%
  group_by(Disease) %>% summarise(total = sum(count)) %>%
  arrange(desc(total)) %>% pull(Disease)
plot_data$Disease <- factor(plot_data$Disease, levels = disease_order)
plot_data$exposure <- factor(plot_data$exposure, levels = c("pm2.5", "pm10","no2", "nox", "shared_protein_count"))

# Plot
fig_D <- ggplot(plot_data, aes(x = Disease, y = freq, fill = exposure)) +
  geom_bar(stat = "identity", position = position_stack(vjust = 0)) +
  scale_fill_manual(values = c(
    "no2" = "#E5C548",
    "nox" = "#F3E581",
    "pm2.5" = "#5875A4",
    "pm10" = "#97AACB",
    "shared_protein_count" = "#D9D9D9"
  ), labels = c(
    "PM2.5 Specific",
    "PM10 Specific",
    "NO2 Specific",
    "NOx Specific",
    "Shared"
  )) +
  ylim(0, 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "top",
    legend.justification = c(0.95, 0.5),
    legend.key.size = unit(0.4, "cm"),
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = c(0, 1)
  ) +
  labs(tag = "D")
ggsave(file.path(output_dir, 'Fig3_D.pdf'), fig_D, device = 'pdf', width = 9, height = 5, dpi = 350)


############################################################
# E. Bubble plot (Disease × Exposure; size=n mediators; fill=global mediation)
############################################################

data <- fread(file.path(input_dir, "mediation_factor_df_FDR.csv")) %>% as.data.frame()
data <- data[data$ACME_FDR < 0.05 & data$Total_FDR < 0.05 & data$Prop_mediated_ave > 0.04, ]
data <- data[data$Disease != 'Osteoporosis', ]
data <- data[data$exposure == 'ap_score_healthy_stat', ]
data$ACME_P_Value[data$ACME_P_Value == 0] <- 2e-16

fisher_combine_p <- function(p_values) {
  p_values <- p_values[!is.na(p_values)]
  X2 <- -2 * sum(log(p_values))
  df <- 2 * length(p_values)
  pchisq(X2, df, lower.tail = FALSE)
}

result <- data %>%
  group_by(exposure, Disease) %>%
  summarise(
    fisher_p = -log(fisher_combine_p(ACME_P_Value)),
    median_proportion = median(Prop_mediated_ave) * 100,
    n = n(), .groups = 'drop'
  ) %>% arrange(exposure, Disease)

result$fisher_p[result$fisher_p == Inf] <- 744

global_mediation <- fread(file.path(input_dir, "PCMA_Global_PC_indirect_effect_10pc.csv")) %>% as.data.frame()
global_mediation <- global_mediation %>%
  pivot_longer(cols = c('ap_score_healthy_stat','PM25_2010','PM10_average','NO2_average','NO_2010'),
               names_to = 'exposure', values_to  = 'Global_Mediation')
result <- inner_join(result, global_mediation, by = c("Disease", "exposure"))
result$Global_Mediation <- result$Global_Mediation * 100

result$exposure <- gsub('NO_2010','NOx',result$exposure)
result$exposure <- gsub('NO2_average','NO2',result$exposure)
result$exposure <- gsub('PM10_average','PM10',result$exposure)
result$exposure <- gsub('PM25_2010','PM2.5',result$exposure)
result$exposure <- gsub('ap_score_healthy_stat','AP index',result$exposure)
result$Disease <- gsub('PeripheralArteryDisease','PAD',result$Disease)
result$Disease <- gsub('LupusErythematosus','LE',result$Disease)
result$Disease <- gsub('HeartFailure','Heart Failure',result$Disease)

result <- result[result$Disease %in% c('Anemia','Heart Failure','PAD','IHD','Arrhythmias',
                                       'CKD','T2D','NAFLD','RA','LE','Pneumonia','COPD'),]
result$Disease <- factor(result$Disease, levels = rev(c('Anemia','Heart Failure','PAD','IHD','Arrhythmias',
                                                        'CKD','T2D','NAFLD','RA','LE','Pneumonia','COPD')))

fig_E <- ggplot(result, aes(x = Disease, y = exposure)) +
  geom_point(aes(size = n, fill = Global_Mediation), shape = 21, color = "black", alpha = 0.9) +
  scale_size_continuous(range = c(3, 23), name = "Mediator Count") +
  scale_fill_gradientn(colors = c('white', "#820823"),
                       values = scales::rescale(c(0, max(result$Global_Mediation))),
                       limits = c(0, max(result$Global_Mediation)),
                       name = "Global Mediated Proportion") +
  labs(y = "Exposure", x = "Disease", tag = "E") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 19),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 19),
    axis.title = element_text(size = 19, face = "bold"),
    plot.margin = margin(0.25, 0.25, 0.25, 0.25),
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
    panel.grid.minor = element_line(color = "gray95", linetype = "dotted"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", size = 1),
    legend.position = "top",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.key = element_rect(fill = "white"),
    legend.background = element_rect(fill = "white"),
    legend.box = "horizontal",
    legend.box.just = "center",
    legend.margin = margin(0, 10, 0, 10),
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = c(0, 1)
  )
ggsave(file.path(output_dir, 'Fig3_E.pdf'), fig_E, width = 11, height = 6.2)


############################################################
# F. Tile heatmap 
############################################################

specific_mediators_df <- read.xlsx(file.path(input_dir, "go_plot_specific_mediators.xlsx"), sheet = 1) %>% as.data.frame()
specific_mediators_df$Description <- factor(specific_mediators_df$Description, levels = rev(unique(specific_mediators_df$Description)))
specific_mediators_df$DiseaseType <- factor(specific_mediators_df$DiseaseType, levels = c('Cardio','Immune','Metabolic','Lung'))


go_all <- specific_mediators_df %>%
  group_by(Description, DiseaseType) %>% summarise(min_LogP = min(LogP), .groups = 'drop') %>% distinct()
go_all$Description <- factor(go_all$Description, levels = rev(unique(specific_mediators_df$Description)))
go_all$DiseaseType <- factor(go_all$DiseaseType, levels = c('Cardio','Immune','Metabolic','Lung'))

# Wide to long, fill NA with 0, back to long
go_all <- go_all %>% spread(DiseaseType, min_LogP)
go_all[is.na(go_all)] <- 0
go_all <- go_all %>% gather(DiseaseType, min_LogP, Cardio:Lung)
go_all$DiseaseType <- factor(go_all$DiseaseType, levels = c('Cardio','Immune','Metabolic','Lung'))

fig_F <- ggplot(data = go_all, aes(x = DiseaseType, y = Description)) +
  geom_tile(aes(fill = -min_LogP), color = 'white', size = 0.1) +
  scale_fill_gradientn(name = '-log10(P.Value)', colours = c("white", "#9C7BA5", "#9C7BA5"), breaks = c(2, 5, 7.5)) +
  theme_bw() +
  scale_y_discrete(position = "left") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(0.25, 0.25, 0.25, 1, "cm"),
    legend.key.size = unit(0.9, "cm"),
    legend.text = element_text(size = 16),
    legend.position = 'bottom',
    legend.direction = "horizontal",
    legend.title = element_text(size = 16),
    panel.grid = element_blank(),
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = c(0, 1)
  ) +
  labs(title = '', tag = 'F') +
  scale_y_discrete(position = "right")
ggsave(file.path(output_dir, 'Fig3_F.pdf'), fig_F, width = 8, height = 11)


############################################################
# G. Stacked bar (x=protein, y=mediated proportion; stacked by disease)
############################################################

data <- fread(file.path(input_dir, "mediation_factor_df_FDR.csv")) %>% as.data.frame()
data <- data[data$ACME_FDR < 0.05 & data$Total_FDR < 0.05 & data$Prop_mediated_ave > 0.04, ]
data <- data[data$Disease != 'Osteoporosis', ]

data_apscore <- data[data$exposure == 'ap_score_healthy_stat', ]
data_apscore_freq <- data_apscore$mediation_factor %>% table() %>% as.data.frame()
data_apscore_freq <- data_apscore_freq[data_apscore_freq$Freq >= 3, ]  # 出现≥3疾病的蛋白

data_apscore <- data_apscore[data_apscore$mediation_factor %in% data_apscore_freq$. , ]
data_apscore <- data_apscore %>% group_by(mediation_factor) %>% mutate(sum_proportion = sum(Prop_mediated_ave))
tmp <- data_apscore %>% dplyr::select(mediation_factor, sum_proportion) %>% distinct() %>% arrange(desc(sum_proportion))
order_vec <- tmp$mediation_factor

data_apscore$mediation_factor <- factor(data_apscore$mediation_factor, levels = order_vec)

data_apscore$Disease <- gsub('PeripheralArteryDisease','PAD', data_apscore$Disease)
data_apscore$Disease <- gsub('LupusErythematosus','LE', data_apscore$Disease)
data_apscore$Disease <- factor(data_apscore$Disease, levels = rev(c('Anemia','HeartFailure','PAD','IHD','Arrhythmias',
                                                                   'CKD','T2D','NAFLD','RA','LE','Pneumonia','COPD')))
data_apscore <- na.omit(data_apscore)

# Color palette 
color_palette <- c()
tmp <- colorRampPalette(c("#9C7BA5", "white"))(9); color_palette <- c(color_palette, tmp[1:5])
tmp <- colorRampPalette(c("#5F9E6E", "white"))(9); color_palette <- c(color_palette, tmp[1:3])
tmp <- colorRampPalette(c("#5875A4", "white"))(9); color_palette <- c(color_palette, tmp[1:2])
tmp <- colorRampPalette(c("#B55D60", "white"))(9); color_palette <- c(color_palette, tmp[1:2])

fig_G <- ggplot(data = data_apscore, aes(x = mediation_factor, fill = Disease)) +
  geom_bar(aes(y = Prop_mediated_ave), stat = "identity") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title.y = element_text(size = 18, face = 'bold'),
    axis.title.x = element_blank(),
    legend.position = c(0.9, 0.95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    plot.tag = element_text(size = 22, face = "bold"),
    plot.tag.position = c(0, 1)
  ) +
  labs(x = "Protein", y = "Stacked Mediated Proportion", fill = "Disease", tag = 'G') +
  scale_fill_manual(values = rev(color_palette))
ggsave(file.path(output_dir, 'Fig3_G.pdf'), fig_G, width = 13, height = 6.315789)



