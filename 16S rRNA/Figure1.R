
library(ggplot2)
library(phyloseq)
library(ggpubr)
library(dplyr)
library(ggsignif)
library(vegan)
library(readxl)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # working on Rstudio only
getwd()
### Loading data

physeq_rare_mean_rounded = readRDS("ressources/RDS/Rarefied_physeq_mean_rounded.rds")


### Alpha diversity
physeq_rare_mean_rounded@sam_data$Week = ifelse(physeq_rare_mean_rounded@sam_data$Week== 2, "week 2",
                                                ifelse(physeq_rare_mean_rounded@sam_data$Week == 4, "week 4", "week 5"))

alpha_div = estimate_richness(physeq_rare_mean_rounded, measures = c("Shannon", "Chao1")) # compute alpha diversity
alpha_div = cbind(alpha_div, physeq_rare_mean_rounded@sam_data) # add the factor to the new table

alpha_div_2week = alpha_div %>% filter(Week == "week 2")
alpha_div_4week = alpha_div %>% filter(Week == "week 4")
alpha_div_5week = alpha_div %>% filter(Week == "week 5")

alpha_all<-rbind(alpha_div_2week,alpha_div_4week,alpha_div_5week)
alpha_all$Irrigation[alpha_all$Irrigation %in% "Watered"]<-"Control"

symnum.args2 <- list(cutpoints = c(0,  0.001, 0.01, 0.05, Inf), symbols = c( "***", "**", "*", "ns"))
star_size = 8 # change for the size of the stars after test


shannon_plot<-ggplot(alpha_all, aes(x = Week, y = Shannon)) + 
  geom_boxplot(aes( fill=Irrigation), outlier.shape = NA, alpha = 0.5, position = position_dodge2(preserve = "single")) + 
  scale_fill_manual(values=c("#018571","#a6611a"))+
  geom_point(aes( fill=Irrigation),size = 1,alpha = 0.8, shape = 21, position = position_jitterdodge(jitter.width = 0.1))+
  theme_bw()+
  labs(title = "", x = "", y = "Shannon Index")+
  theme(axis.text.x = element_text(size = 10, colour = "black", hjust=0.5,face = "bold"), 
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(face = "bold", size = 12, vjust = 3), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12,face = "bold"),
        panel.grid = element_blank()
  ) + stat_compare_means(aes(group = Irrigation), label = "p.signif",
                         method = "t.test", symnum.args = symnum.args2 , label.y =5.8, size = star_size)+
  scale_y_continuous(limits = c(4.15,6.15),expand = c(0,0))


# ======================= Beta diversity ====================================
plot_pcoa = function(data, color, cluster, Shape, filename){
  list_output = list() # to return the plot and the beta diversity results
  print("Computing PCoA")
  physeq.ord_pcoa = ordinate(data, "PCoA", "bray")
  print("Plotting PCoA")
  
  p3 = plot_ordination(data, physeq.ord_pcoa, shape = Shape) +
    scale_color_manual(values=c("#018571","#a6611a")) +
    geom_point(size = 4, aes(color = !!sym(color))) +  # Use aes with sym
    theme_bw() +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 25), 
          legend.text = element_text(size = 18), legend.title = element_text(size = 20))
  

  
  # ggsave(filename)
  list_output[[1]] = p3
  list_output[[2]] = physeq.ord_pcoa
  return(list_output)
}

sample_data(physeq_rare_mean_rounded)$Irrigation[sample_data(physeq_rare_mean_rounded)$Irrigation %in% "Watered"]<-"Control"
p4 = plot_pcoa(physeq_rare_mean_rounded, "Irrigation", "Group", "Week", "output/PCoA_16S_rRNA_cluster.png")

pcoa_16S<-p4[[1]]


### Plant data
week2d<-read_excel("ressources/2 weeks.xlsx", sheet = "Sheet1")
week4d<-read_excel("ressources/4 weeks.xlsx", sheet = "Sheet1")
week5d<-read_excel("ressources/5 weeks.xlsx", sheet = "Sheet1") 

week2<-week2d %>% select(-Shootdryweight) %>% mutate(Week="Week2")%>% filter(Treatment %in% "Water")
week4<-week4d %>% select(-Treatment1) %>% mutate(Week="Week4")%>% filter(Treatment %in% "Water")
week5<-week5d %>% select(-Treatment1,-Sootdryweight)%>% mutate(Week="Week5")%>% filter(Treatment %in% "Water")

All_data<-rbind(week2,week4,week5)  %>% select(-Treatment)

All_data$Group[All_data$Group %in% "irrigated"]<-"Control"
colnames(All_data)[2]<-"Irrigation"

All_data$Week<-stringr::str_replace(All_data$Week,"Week", "Week ")

#Plot shoot length
shoot_le<-ggplot(All_data, aes(x = Week, y = Shootlength)) + 
  geom_boxplot(aes( fill=Irrigation), outlier.shape = NA, alpha = 0.5, position = position_dodge2(preserve = "single")) + 
  scale_fill_manual(values=c("#018571","#a6611a"))+
  geom_point(aes( fill=Irrigation),size = 1,alpha = 0.8, shape = 21, position = position_jitterdodge(jitter.width = 0.1))+
  theme_bw()+
  labs(title = "", x = "", y = "Shoot Length (cm)")+
  theme(axis.text.x = element_text(size = 10, colour = "black", hjust=0.5,face = "bold"), 
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(face = "bold", size = 12, vjust = 3), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12,face = "bold"),
        panel.grid = element_blank()
  ) + stat_compare_means(aes(group = Irrigation), label = "p.signif",method = "t.test", symnum.args = symnum.args2 , label.y = 39, size = star_size)+
  scale_y_continuous(limits = c(24,40))


#Plot the Shoot fresh weight
shoot_fresh<-ggplot(All_data, aes(x = Week, y = Shootfreshweight)) + 
  geom_boxplot(aes(x = Week, y = Shootfreshweight, fill=Irrigation), outlier.shape = NA, alpha = 0.5, position = position_dodge2(preserve = "single")) + 
  scale_fill_manual(values=c("#018571","#a6611a"))+
  geom_point(aes(x = Week, y = Shootfreshweight, fill=Irrigation),size = 1,alpha = 0.8, shape = 21, position = position_jitterdodge(jitter.width = 0.1))+
  theme_bw()+
  labs(title = "", x = "", y = "Shoot Fresh Weight (g)")+
  theme(axis.text.x = element_text(size = 10, colour = "black", hjust=0.5, face = "bold"), 
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(face = "bold", size = 12, vjust = 3), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12,face = "bold"),
        panel.grid = element_blank()
  )+ stat_compare_means(aes(group = Irrigation), label = "p.signif",method = "t.test", symnum.args = symnum.args2 , label.y = 9.2, size = star_size)+
  scale_y_continuous(limits = c(0,10.22), expand =c(0,0))

#Plot the root dry weight
root_dry<-ggplot(All_data, aes(x = Week, y = Rootdryweight)) + 
  geom_boxplot(aes(x = Week, y = Rootdryweight, fill=Irrigation), outlier.shape = NA, alpha = 0.5, position = position_dodge2(preserve = "single")) + 
  scale_fill_manual(values=c("#018571","#a6611a"))+
  geom_point(aes(x = Week, y = Rootdryweight, fill=Irrigation),size = 1,alpha = 0.8, shape = 21, position = position_jitterdodge(jitter.width = 0.1))+
  theme_bw()+
  labs(title = "", x = "", y = "Root Dry Weight (g)")+
  theme(axis.text.x = element_text(size = 10, colour = "black", hjust=0.5,face = "bold"), 
        axis.text.y = element_text(size = 10,face = "bold"),
        axis.title.y = element_text(face = "bold", size = 12, vjust = 3), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12,face = "bold"),
        panel.grid = element_blank()
  )+ stat_compare_means(aes(group = Irrigation), label = "p.signif",method = "t.test", symnum.args = symnum.args2 ,label.y = 1.8, size = star_size)+
  scale_y_continuous(limits = c(0,2), expand =c(0,0))

#Extract the legend
sha_legend<-get_legend(pcoa_16S)
leg_for_plot<-ggarrange(sha_legend)

#Plot and save
x11(height = 6, width = 7)
ggarrange(shoot_le, shoot_fresh,root_dry,pcoa_16S,shannon_plot,  leg_for_plot, ncol = 3,nrow = 2,
          legend = "none", align = "h", labels = c("A","B","C","D","E"))
ggsave("Fig1_Plant_data.png", dpi = 400, width = 10, height = 7)
dev.off()




### PERMANOVA (Table S1)
permanova = function(physeq, var1, var2 = NULL, var3 = NULL, var4 = NULL, var5 = NULL, matrix = "bray"){
  cat("Preparing data for PERMANOVA\n")
  data = as.data.frame((t(physeq@otu_table))) # convert the otu_table to a dataframe to perform the permANOVA
  
  metdata = physeq@sam_data
  data = cbind(data, physeq@sam_data)
  data$Irrigation = as.factor(data$Irrigation)
  data$Week = as.factor(data$Week)
  data$Compartments = as.factor(data$Compartments)
  data$Treatment = as.factor(data$Treatment)
  data$Group= as.factor(data$Group)
  # remove the 6 last columns (metadata variable)
  data_corr = data[, -((ncol(data)-6):ncol(data))]
  cat("Computing Bray-Curtis dissimilarity matrix\n")
  bray_dist = vegdist(data_corr, method = matrix)
  cat("Computing PERMANOVA\n")
  formula_str = paste("bray_dist ~", var1, if (!is.null(var2)) paste("*", var2) else "", 
                      if (!is.null(var3)) paste("*", var3) else "", 
                      if (!is.null(var4)) paste("*", var4) else "", 
                      if (!is.null(var5)) paste("*", var5) else "")
  formula = as.formula(formula_str)
  
  model = adonis2(formula, data = data, permutations = 999, by = "terms")
  return(model)
}
permanova(physeq = physeq_rare_mean_rounded, var1 = "Irrigation", var2 = "Week")



