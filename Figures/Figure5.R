
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
data = fread("./mediation_factor_df_FDR.csv") %>% as.data.frame()
protein_counts <- data %>%
  group_by(Disease, exposure) %>%
  summarise(protein_count = n_distinct(mediation_factor)) %>%
  ungroup()

data = data[data$ACME_FDR<0.05 & data$Total_FDR < 0.05 & data$Prop_mediated_ave > 0.04,]

data$mediation_factor

result <- data %>%
  # 按disease和pollutants分组
  group_by(Disease, exposure) %>%
  # 计算每组中不同protein_name的数量
  summarise(protein_count = n_distinct(mediation_factor), .groups = "drop") %>%
  # 转换为宽数据格式，行是disease，列是pollutants
  pivot_wider(
    names_from = exposure,
    values_from = protein_count,
    values_fill = 0  # 如果没有对应的蛋白质，填充0
  )

result = result[,c(1,5,4,3,2,6)]
colnames(result) = c('Disease','PM2.5','PM10','NO2','NOx','AP index')


library(ggradar)
library(RColorBrewer)
c("#9E0142", "#D53E4F", "#F46D43" ,"#FDAE61", "#FEE08B" ,"#FFFFBF" ,"#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
spectral_colors <- brewer.pal(11, "Spectral")
spectral_colors = c(spectral_colors,'darkorchid2','darkblue','darkgreen','bisque3','#053264','#820823','#54A966','#4C4C4C')


plot_radar = function(tmp2,disease,color,seq){
  tmp = ggradar(tmp2,
                grid.min = 0,
                grid.mid = seq[2], grid.max = seq[3],
                values.radar = seq,
                gridline.min.colour = "grey",
                gridline.mid.colour = "#4C4C4C", gridline.max.colour = "#4C4C4C",
                axis.label.size = 5.3, axis.line.colour = "grey",
                plot.title = disease, 
                background.circle.colour = "white",
                background.circle.transparency = 0.1,
                group.point.size = 5,
                group.line.width = 1.8,group.colours = color,
                legend.title = "Disease",
                legend.position = "bottom")+
    theme(legend.title = element_text(size=3),
          plot.title = element_text(size = 24),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  return(tmp)
}
