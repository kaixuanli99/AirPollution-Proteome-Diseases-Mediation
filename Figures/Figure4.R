
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
