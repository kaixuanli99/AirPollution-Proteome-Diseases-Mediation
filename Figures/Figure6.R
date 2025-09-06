
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
enrichment_df = read.xlsx('D:\\DATA\\UKB_AP_Proteome\\PWAS-lwr\\enrichment_df.xlsx',sheet=1)
enrichment_df2 = enrichment_df %>% group_by(exposure) %>% arrange(p.adjust,.by_group = TRUE) %>% dplyr::slice_head(n = 10)
enrichment_df2 = enrichment_df2[enrichment_df2$p.adjust < 0.01,]
"#053264", "#820823"
"#5875A4", "#d2dae2","#B55D60"
exposure_name = c('NOx','NO2','PM2.5','PM10','Common')
plot_list = list()
k=1
for(i in c(3,4,2,1,5)){

  tmp = enrichment_df2[enrichment_df2$exposure == exposure[i],]
  tmp$logP = -log(tmp$p.adjust)
  
  plot_tmp = ggplot(data = tmp, aes(x = Count, y = reorder(Description, logP), fill = logP)) + 
            geom_bar(stat = "identity", width = 0.7) +
            # 修改色标范围和方向，使logP负值越低（绝对值越大）颜色越深
            scale_fill_gradientn(name = '-log10(P.adjust)',
                                 colours = c("#5875A4", "#B55D60")) +
            theme_bw() +
            theme(axis.text.x = element_text(size = 16),
                  plot.title = element_text(size = 16, face = "bold"),
                  axis.text.y = element_text(size = 16),
                  axis.title.x = element_text(size = 16),
                  axis.title.y = element_blank(),
                  plot.margin = margin(0.25, 0.25, 0.25, 1, "cm"),
                  legend.key.size = unit(0.4, "cm"),
                  legend.title = element_text(size = 14),
                  legend.text = element_text(size = 14),
                  legend.position = "right",
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor = element_blank()) +
            labs(title = paste0(exposure_name[i],'-Specific'), x = 'Gene Count') 
  plot_list[[k]] <- plot_tmp
  k=k+1
  
}

library(ggpubr)
combined_plot <- ggarrange(plotlist = plot_list, 
                           ncol = 2,               # 两列排列
                           nrow = ceiling(length(exposure_name)/2), # 自动计算所需行数
                           common.legend = F,   # 共享图例
                           legend = "right",       # 图例位置
                           align = "hv",
                           labels = c("A",'B','C','D','E'), 
                           font.label = list(size = 24, face = "bold"))           # 水平和垂直对齐

ggsave(filename = "./Figure6.pdf", plot = combined_plot, 
       device = "pdf", width = 21, height =12.98701)

