
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

go_df <- fread(go_enrichment_csv) %>% as.data.frame()
go_df$exposure <- gsub('Nitrogen_oxides_2010', 'NOx', go_df$exposure)
go_df$exposure <- gsub('NO2_mean', 'NO2', go_df$exposure)
go_df$exposure <- gsub('PM25_2010', 'PM2.5', go_df$exposure)
go_df$exposure <- gsub('PM10_mean', 'PM10', go_df$exposure)
go_df$exposure <- gsub('ap_score_healthy_stat', 'AP index', go_df$exposure)

top10_df <- go_df %>% group_by(exposure) %>% slice_min(order_by = p.adjust, n = 10)
top10_df$logP <- -log10(top10_df$p.adjust)
top10_df$Description <- stringr::str_wrap(top10_df$Description, 42)

p_go <- ggplot(data = top10_df, aes(x = exposure, y = Description)) +
	geom_tile(aes(fill = logP), color = 'white', size = 0.1) +
	scale_fill_gradientn(name = '-log10(P.adjust)', colours = c("white",  "#B55D60")) +
	theme_bw() +
	scale_y_discrete(position = "left") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16.8), plot.title = element_text(size = 16, face = "bold"), axis.text.y = element_text(size = 16.8), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(0.25, 0.25, 0.25, 1, "cm"), legend.key.size = unit(0.4, "cm"), legend.title = element_text(size = 16), legend.text = element_text(size = 16), legend.position = "bottom", panel.grid = element_blank()) +
	labs(title = 'GO') +
	scale_y_discrete(position = "right")
