
rm(list = ls())
suppressPackageStartupMessages({
	library(ggplot2)
	library(ggrepel)
	library(dplyr)
	library(data.table)
	library(ggpubr)
	library(scales)
})



gene_annotations_csv <- "./inputs/Proteome_Annotation.csv"                 # eid -> chr/pos
pwas_pm25_csv <- "./inputs/PM25_2010_BasicCovar+NO.csv"                    # one pollutant PWAS (PM2.5)
ap_score_pwas_csv <- "./inputs/ap_score_healthy_stat_BasicCovar.csv"       # AP index PWAS
go_enrichment_csv <- "./inputs/go_ap_enrichment_df.csv"                    # GO enrichment prepared
kegg_enrichment_csv <- "./inputs/kegg_ap_enrichment_df.csv"                # KEGG enrichment prepared
reactome_enrichment_csv <- "./inputs/reactome_ap_enrichment_df.csv"        # Reactome enrichment prepared
fig2c_counts_csv <- "./inputs/fig2c_counts.csv"                            # counts table for Fig2C (pollutants, Common, Specific, All)

out_pdf <- "./Figure2.pdf"
out_png <- "./Figure2.png"

# ------------------------------


gene_annotations <- fread(gene_annotations_csv)

data <- fread(pwas_pm25_csv) %>% as.data.frame()
data$eid <- toupper(data$eid)
data <- merge(data, gene_annotations[, c(1, 2, 3)], by.x = 'eid', by.y = 'hgnc_symbol', all.y = TRUE)
colnames(data)[c(8, 9)] <- c('chr', 'pos')

rpinf <- function(x) { x[is.infinite(x)] <- 323.3062; x }
color_pm25 <- '#DAC060'

t <- as.data.table(data)
t[, mlog10p := rpinf(-log10(P.Value))]
t[chr == "chrX"]$chr <- 23
t[chr == "chrY"]$chr <- 24
t$chr <- gsub("chr", "", t$chr) %>% as.numeric()
t$chr <- factor(t$chr, levels = 1:24)
t <- t[order(chr, pos)]
t$fdr <- p.adjust(t$P.Value, method = "fdr")
thres_fdr <- t[fdr >= 0.05, max(mlog10p)]

# Build cumulative position
t$pos1 <- NA_real_; t$index <- NA_integer_; ind <- 0
setDF(t)
for (i in unique(t$chr)) { ind <- ind + 1; t[t$chr == i, "index"] <- ind }
lastbase <- 0; ticks <- NULL; whole_col <- NULL
cols_chr <- c("gray", "gray")
for (i in unique(t$index)) {
	if (i == 1) {
		t[t$index == i, "pos1"] <- t[t$index == i, "pos"]
	} else {
		lastbase <- lastbase + tail(t[t$index == i - 1, "pos"], 1)
		t[t$index == i, "pos1"] <- t[t$index == i, "pos"] + lastbase
	}
	ticks <- c(ticks, (min(t[t$index == i, "pos1"]) + max(t[t$index == i, "pos1"])) / 2 + 1)
	whole_col <- c(whole_col, rep(cols_chr[i %% 2 + 1], length(t[t$index == i, "pos"])))
}
xlabel <- "Chromosome"
labs_chr <- unique(t$chr)
xmax <- ceiling(max(t$pos1) * 1.03)
xmin <- floor(max(t$pos1) * -0.03)

sig <- which((t$fdr <= 0.05) & (!is.na(t$P.Value)))
tmp <- t %>% as.data.frame() %>% arrange(fdr)
cutoff <- tmp$fdr[20]
ticks_break <- c(); for (i in 1:(length(ticks) - 1)) { ticks_break <- c(ticks_break, c(ticks[i] + ticks[i + 1]) / 2) }

p_mht <- ggplot(t, aes(x = pos1, y = mlog10p)) +
	geom_point(color = whole_col, size = 1.5) +
	geom_point(data = t[sig, ], color = color_pm25, size = 2.1) +
	geom_hline(yintercept = thres_fdr, color = color_pm25) +
	geom_hline(yintercept = -log10(0.05), linetype = "twodash", color = 'grey') +
	xlab(xlabel) + ylab(expression(-log[10](italic(p)))) + ggtitle("PM2.5") +
	scale_x_continuous(breaks = ticks, labels = labs_chr, limits = c(xmin, xmax), expand = c(0.02, 4)) +
	scale_y_continuous(breaks = seq(2, ceiling(max(-log10(t$P.Value), na.rm = TRUE)), 5)) +
	theme_bw() +
	theme(panel.grid = element_blank(), plot.title = element_text(size = 16, face = "bold"), axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 14), axis.title.y = element_text(size = 16), axis.title.x = element_text(size = 16)) +
	geom_text_repel(data = t[t$fdr <= cutoff, ], aes(label = eid), color = color_pm25, size = 6) +
	theme(plot.margin = unit(c(0.6, 0.6, 0.6, 0.6), "cm"))
for (tick_i in ticks_break) { p_mht <- p_mht + geom_vline(xintercept = tick_i, linetype = "dashed", color = "grey70", linewidth = 0.2) }

# ------------------------------


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

# ------------------------------


kegg_df <- fread(kegg_enrichment_csv) %>% as.data.frame()
kegg_df$exposure <- gsub('Nitrogen_oxides_2010', 'NOx', kegg_df$exposure)
kegg_df$exposure <- gsub('NO2_mean', 'NO2', kegg_df$exposure)
kegg_df$exposure <- gsub('PM25_2010', 'PM2.5', kegg_df$exposure)
kegg_df$exposure <- gsub('PM10_mean', 'PM10', kegg_df$exposure)
kegg_df$exposure <- gsub('ap_score_healthy_stat', 'AP index', kegg_df$exposure)

# If not provided, create Neg_logP
if (!'Neg_logP' %in% colnames(kegg_df) && 'pvalue' %in% tolower(colnames(kegg_df))) {
	kegg_df$Neg_logP <- -log10(kegg_df$pvalue)
}
if (!'Neg_logP' %in% colnames(kegg_df) && 'p.adjust' %in% colnames(kegg_df)) {
	kegg_df$Neg_logP <- -log10(kegg_df$p.adjust)
}

kegg_top <- kegg_df %>% group_by(exposure) %>% slice_min(order_by = if ('p.adjust' %in% colnames(kegg_df)) p.adjust else -Neg_logP, n = 10)
kegg_top$Description <- stringr::str_wrap(kegg_top$Description, 42)

p_kegg <- ggplot(data = kegg_top, aes(x = exposure, y = Description)) +
	geom_tile(aes(fill = Neg_logP), color = 'white', size = 0.1) +
	scale_fill_gradientn(name = '-log10(P.adjust)', colours = c("white",  "#5875A4")) +
	theme_bw() +
	scale_y_discrete(position = "left") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16.8), plot.title = element_text(size = 16, face = "bold"), axis.text.y = element_text(size = 16.8), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(0.25, 0.25, 0.25, 1, "cm"), legend.key.size = unit(0.4, "cm"), legend.title = element_text(size = 16), legend.text = element_text(size = 16), legend.position = "bottom", panel.grid = element_blank()) +
	labs(title = 'KEGG') +
	scale_y_discrete(position = "right")

# ------------------------------


react_df <- fread(reactome_enrichment_csv) %>% as.data.frame()
react_df$exposure <- gsub('Nitrogen_oxides_2010', 'NOx', react_df$exposure)
react_df$exposure <- gsub('NO2_mean', 'NO2', react_df$exposure)
react_df$exposure <- gsub('PM25_2010', 'PM2.5', react_df$exposure)
react_df$exposure <- gsub('PM10_mean', 'PM10', react_df$exposure)
react_df$exposure <- gsub('ap_score_healthy_stat', 'AP index', react_df$exposure)

if (!'logP' %in% colnames(react_df) && 'p.adjust' %in% colnames(react_df)) {
	react_df$logP <- -log10(react_df$p.adjust)
}


react_df$Description <- gsub('Regulation of Insulin-like Growth Factor \\(?IGF\\)? transport and uptake by Insulin-like Growth Factor Binding Proteins \\(?IGFBPs\\)?',
							  'Regulation of IGF transport and uptake by IGFBPs', react_df$Description)
react_df$Description <- gsub('TNF receptor superfamily \\(?TNFSF\\)? members mediating non\\-canonical NF\\-kB pathway',
							  'TNFSF members mediating non-canonical NF-kB pathway', react_df$Description)
react_df$Description <- gsub('Immunoregulatory interactions between a Lymphoid and a non\\-Lymphoid cell',
							  'Immunoregulatory interactions between Lymphoid and Lymphoid cell', react_df$Description)

react_top <- react_df %>% group_by(exposure) %>% slice_min(order_by = if ('p.adjust' %in% colnames(react_df)) p.adjust else -logP, n = 10)
react_top$logP <- ifelse('logP' %in% colnames(react_top), react_top$logP, -log10(react_top$p.adjust))
react_top$Description <- stringr::str_wrap(react_top$Description, 42)

p_react <- ggplot(data = react_top, aes(x = exposure, y = Description)) +
	geom_tile(aes(fill = logP), color = 'white', size = 0.1) +
	scale_fill_gradientn(name = '-log10(P.adjust)', colours = c("white",  "#5F9E6E")) +
	theme_bw() +
	scale_y_discrete(position = "left") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16.8), plot.title = element_text(size = 16, face = "bold"), axis.text.y = element_text(size = 16.8), axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = margin(0.25, 0.25, 0.25, 1, "cm"), legend.key.size = unit(0.4, "cm"), legend.title = element_text(size = 16), legend.text = element_text(size = 16), legend.position = "bottom", panel.grid = element_blank()) +
	labs(title = 'Reactome') +
	scale_y_discrete(position = "right")

# ------------------------------


df2c <- fread(fig2c_counts_csv) %>% as.data.frame()   # expects columns: pollutants, Common, Specific, All
df2c$Common <- ifelse('Common' %in% colnames(df2c), df2c$Common, df2c$All - df2c$Specific)
data_long <- reshape2::melt(df2c, id.vars = "pollutants", measure.vars = c("Common", "Specific"), variable.name = "Type", value.name = "Value")
data_long$Type <- factor(data_long$Type, levels = c("Common", "Specific"))
data_long$pollutants <- factor(data_long$pollutants, levels = rev(c('PM2.5','PM10','NO2','NOx','AP index','All')))

p_fig2c <- ggplot(data_long, aes(x = pollutants, y = Value, fill = pollutants)) +
	geom_bar(stat = "identity", aes(alpha = Type), width = 0.8) +
	scale_alpha_manual(values = c("Common" = 0.6, "Specific" = 1)) +
	scale_fill_manual(values = c("NOx" = "#5C50A2", "NO2" = "#5875A4", "PM10" = "#5F9E6E", "PM2.5" = "#DAC060", "All" = "#BEBEBE", "AP index" = "#2F7FC1")) +
	labs(y = "Number of Proteins", x = "Pollutants", fill = 'Group') +
	coord_flip() +
	theme_classic() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	scale_y_continuous(expand = c(0, 0)) +
	theme(axis.title.x = element_text(size = 16, margin = margin(t = 10)), axis.title.y = element_text(size = 16, margin = margin(r = 10)), axis.text = element_text(size = 14), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.margin = margin(20, 20, 20, 20), legend.title = element_text(size = 16), legend.text = element_text(size = 14), legend.position = "top", legend.justification = c(0.9, 0.5)) +
	guides(fill = 'none', alpha = 'none')
# ------------------------------


p_top <- p_mht
p_mid <- ggarrange(p_vol, p_fig2c, ncol = 2, widths = c(1, 1.2))
p_enrich <- ggarrange(p_go, p_kegg, p_react, ncol = 3, widths = c(1, 1.05, 1.1))
p_all <- ggarrange(p_top, p_mid, p_enrich, ncol = 1, heights = c(1.2, 1, 1))

ggsave(out_png, plot = p_all, width = 14, height = 11, dpi = 350)
ggsave(out_pdf, plot = p_all, device = "pdf", width = 14, height = 11)


