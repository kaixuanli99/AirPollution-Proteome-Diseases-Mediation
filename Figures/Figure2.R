
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
data_vol <- fread(ap_score_pwas_csv) %>% as.data.frame()
data_vol$eid <- toupper(data_vol$eid)
data_vol$log2FC <- data_vol$logFC
data_vol$log10fdr <- -log10(data_vol$adj.P.Val)
cut_off_FDR <- 0.05
cut_off_log2FC <- quantile(abs(data_vol$log2FC), 0.8, na.rm = TRUE)
data_vol$Sig <- ifelse(data_vol$adj.P.Val < cut_off_FDR & abs(data_vol$log2FC) >= cut_off_log2FC, ifelse(data_vol$log2FC > cut_off_log2FC, 'Up', 'Down'), 'no')

p_vol <- ggplot(data_vol, aes(x = log2FC, y = log10fdr, colour = Sig)) +
	geom_point(alpha = 1, size = 2) +
	scale_color_manual(values = c("#5875A4", "#d2dae2", "#B55D60")) +
	geom_vline(xintercept = c(-cut_off_log2FC, cut_off_log2FC), col = "black", lwd = 0.8, linetype = "dashed") +
	geom_hline(yintercept = -log10(cut_off_FDR), col = "black", lwd = 0.8, linetype = "dashed") +
	labs(x = "Effect size") + ylab(expression(-log[10](FDR))) +
	geom_text_repel(data = subset(data_vol, adj.P.Val < 0.05 & abs(log2FC) > cut_off_log2FC), aes(label = eid), size = 6) +
	ggtitle("AP index") + theme_bw() +
	theme(plot.title = element_text(size = 16, face = "bold"), legend.position = 'none', legend.title = element_blank(), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16), panel.grid = element_blank()) +
	theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

# -----