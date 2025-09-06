
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


plot_PS_ap_cox <- function(tmp,ap,disease,k){

    tmp = tmp[1:5,]
    tmp[6,] = c(1,1,1,'1.00(1.00~1.00)',1,disease,ap)
    tmp = tmp[c(6,1,2,3,4,5),]
    #row.names(tmp)[1] = 'PS_NO2_AVELow PS & Low NO2'
    tmp$group  = rep(c('Low ProteinScore','Intermediate ProteinScore','High ProteinScore'),each = 2)
    tmp$ap_group  = rep(c('Low Exposure','High Exposure'), 3)
    tmp$group  = factor(tmp$group,levels= rev(c('Low ProteinScore','Intermediate ProteinScore','High ProteinScore')))
    tmp$ap_group  = factor(tmp$ap_group,levels= rev(c('Low Exposure','High Exposure')))

    tmp$HR = as.numeric(tmp$HR)
    tmp$Lowerlimit = as.numeric(tmp$Lowerlimit)
    tmp$Upperlimit = as.numeric(tmp$Upperlimit)

    ap = ap

    if(k == 1){
        p <- ggplot(tmp, aes(x = group, y = HR, group = as.factor(ap_group), shape = as.factor(ap_group),colour = as.factor(ap_group))) +
                geom_hline(yintercept = 1, linetype = 2, color = "grey4",size = 1.15) +
                geom_point(position = position_dodge(width = 0.5), size = 5) +
                geom_errorbar(aes(ymax = Upperlimit, ymin = Lowerlimit),
                width = 0.2, alpha = 0.7,
                position = position_dodge(width = 0.5),size=1.5) +
                ylim(0,tmp[6,]$Upperlimit+1)+
                theme_minimal()+
                theme_bw() +
                coord_flip() +labs(title = paste(disease," ",ap,sep=''))+
                guides(shape = guide_legend(title = 'Exposure',order =1 ))+
                theme(axis.title.y = element_blank(),
                    axis.text.y = element_text(size = 18,color="black"),
                    axis.title.x = element_text(size = 18,color="black"),
                    legend.position = "none",
                    plot.title = element_text(size = 18))+
                scale_shape_manual(values = c(17,19))  +
                scale_colour_manual(values = c("Low Exposure" = "lightblue", "High Exposure" = "#C62703FF"))
    }else if(k %in% c(2,3,4)){
        p <- ggplot(tmp, aes(x = group, y = HR, group = as.factor(ap_group), shape = as.factor(ap_group),colour = as.factor(ap_group))) +
                geom_hline(yintercept = 1, linetype = 2, color = "grey4",size = 1.15) +
                geom_point(position = position_dodge(width = 0.5), size = 5) +
                geom_errorbar(aes(ymax = Upperlimit, ymin = Lowerlimit),
                width = 0.2, alpha = 0.7,
                position = position_dodge(width = 0.5),size=1.5) +
                ylim(0,tmp[6,]$Upperlimit+1)+
                theme_minimal()+
                theme_bw() +
                coord_flip() +labs(title = paste(disease," ",ap,sep=''))+
                guides(shape = guide_legend(title = 'Exposure',order =1 ))+
                theme(axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.title.x = element_text(size = 18,color="black"),
                    legend.position = "none",
                    plot.title = element_text(size = 18))+
                scale_shape_manual(values = c(17,19))  +
                scale_colour_manual(values = c("Low Exposure" = "lightblue", "High Exposure" = "#C62703FF"))
    }else{
        p <- ggplot(tmp, aes(x = group, y = HR, group = as.factor(ap_group), shape = as.factor(ap_group),colour = as.factor(ap_group))) +
                geom_hline(yintercept = 1, linetype = 2, color = "grey4",size = 1.15) +
                geom_point(position = position_dodge(width = 0.5), size = 5) +
                geom_errorbar(aes(ymax = Upperlimit, ymin = Lowerlimit),
                width = 0.2, alpha = 0.7,
                position = position_dodge(width = 0.5),size=1.5) +
                ylim(0,tmp[6,]$Upperlimit+1)+
                theme_minimal()+
                theme_bw() +
                coord_flip() +labs(title = paste(disease," ",ap,sep=''))+
                #guides(shape = guide_legend(title = 'Exposure',order =1 ))+
                theme(axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    legend.text = element_text(size=16),
                    axis.title.x = element_text(size = 18,color="black"),
                    plot.title = element_text(size = 18))+
                scale_shape_manual(name = "Exposure",values = c(17,19))  +
                scale_colour_manual(name = "Exposure",values = c("Low Exposure" = "lightblue", "High Exposure" = "#C62703FF"))
        }
    return(p)
}

##### 按照疾病来画
for(j in 1:disease_count){
    outcome = grep(paste(disease_list_name[j],'_2010_2020',sep=''),colnames(UKB_RA_df_clean)[1:184])
    endpoint = colnames(UKB_RA_df_clean)[outcome][1]
    duration = colnames(UKB_RA_df_clean)[outcome][2]
    k=1
    for(i in quantile2_exposure){
        
        j_idx = disease_quantile3[j]

        tmp = ifelse(UKB_RA_df_clean[,j_idx]==1 & UKB_RA_df_clean[,i]==1,1,
                ifelse(UKB_RA_df_clean[,j_idx]==1 & UKB_RA_df_clean[,i]==2,2,
                    ifelse(UKB_RA_df_clean[,j_idx]==2 & UKB_RA_df_clean[,i]==1,3,
                        ifelse(UKB_RA_df_clean[,j_idx]==2 & UKB_RA_df_clean[,i]==2,4,
                            ifelse(UKB_RA_df_clean[,j_idx]==3 & UKB_RA_df_clean[,i]==1,5,
                                ifelse(UKB_RA_df_clean[,j_idx]==3 & UKB_RA_df_clean[,i]==2,6,NA))))))
        UKB_RA_df_clean$tmp = tmp
        UKB_RA_df_clean$tmp = as.factor(UKB_RA_df_clean$tmp)

        UKB_RA_df_clean$tmp = ifelse(UKB_RA_df_clean$tmp==1 ,'Low PS & Low Exposure',
                                    ifelse(UKB_RA_df_clean$tmp==2 ,'Low PS & High Exposure',
                                        ifelse(UKB_RA_df_clean$tmp==3 ,'Intermediate PS & Low Exposure',
                                            ifelse(UKB_RA_df_clean$tmp==4 ,'Intermediate PS & High Exposure',
                                            ifelse(UKB_RA_df_clean$tmp==5 ,'High PS & Low Exposure',
                                                ifelse(UKB_RA_df_clean$tmp==6 ,'High PS & High Exposure',NA))))))
        UKB_RA_df_clean$tmp = factor(UKB_RA_df_clean$tmp,levels = c('Low PS & Low Exposure','Low PS & High Exposure','Intermediate PS & Low Exposure','Intermediate PS & High Exposure','High PS & Low Exposure','High PS & High Exposure'))
        
        a <- c('tmp',covariates)
        mySurv <- Surv(time = UKB_RA_df_clean[,c(duration)], event = UKB_RA_df_clean[,c(endpoint)])
        FML<- as.formula(paste0('mySurv~',paste(a,collapse = "+")))
        fit_cox <- coxph(FML,data = UKB_RA_df_clean) #### cox regression


        res_cox<-exp(cbind(coef(fit_cox),confint.default(fit_cox))) %>% as.data.frame()
        res_cox[,4]<-paste0(round(res_cox[,1],2),'(',round(res_cox[,2],2),'~',round(res_cox[,3],2),')') 
        a<-summary(fit_cox)
        b<-as.data.frame(a$coefficients)
        res_cox<-cbind(res_cox,b[,5])
        colnames(res_cox)<-c("HR","Lowerlimit","Upperlimit","HR_95%CI","P_value")

        exposure = exposure_plot[k]
        res_cox$endpoint = disease_list[j]
        res_cox$exposure = exposure
        endpoint_plot = disease_list[j]
        
        assign(paste('p',k,sep=''),plot_PS_ap_cox(res_cox,exposure,endpoint_plot,k))

        write.csv(res_cox, 
                        file = paste("/home/likaixuan/likaixuan/data/UKB_AP_Proteome/data/result/AP-related_ProteinRiskScore/",endpoint,'_',exposure,'.csv',sep='') ,
                        quote = F
                        )
        k=k+1

    }
    p <- ggarrange(p1, p2, p3, p4, p5, ncol = 5,widths = c(1.45, 0.8, 0.8, 0.8, 1.2))
    ggsave(paste("ProteinScore_",endpoint,'_2.png',sep=''), path = "/home/likaixuan/likaixuan/data/UKB_AP_Proteome/data/result/AP-related_ProteinRiskScore", plot = p, width = 22, height = 4.074, dpi = 350)


}
