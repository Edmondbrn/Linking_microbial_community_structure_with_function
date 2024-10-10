# NRPS Rhizoplane drought - control comparison - alpha and beta diversity

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(phyloseq)
#Descriptive statistics of the data (prior to rarefying)


rhizoplane_NRPS<-readRDS("ressources/RDS/raw_physeq_rhizoplane.rds")
# Subset control
rhizoplane_NRPS_co<-subset_samples(rhizoplane_NRPS, Treatment %in% "Control")
rhizoplane_NRPS_co<-prune_taxa(taxa_sums(rhizoplane_NRPS_co)>0,rhizoplane_NRPS_co)
rhizoplane_NRPS_co # number of ACs
min(sample_sums(rhizoplane_NRPS_co)) 
max(sample_sums(rhizoplane_NRPS_co))
median(sample_sums(rhizoplane_NRPS_co))
sum(sample_sums(rhizoplane_NRPS_co))

#Import data
Rarefied_physeq_mean_rounded_rhizoplane<-readRDS("ressources/RDS/Rarefied_physeq_mean_rounded_rhizoplane.rds")
alpha_div_plane = estimate_richness(Rarefied_physeq_mean_rounded_rhizoplane, measures = c("Shannon", "Chao1")) # compute alpha diversity
alpha_div_plane = cbind(alpha_div_plane, Rarefied_physeq_mean_rounded_rhizoplane@sam_data) # add the factor to the new table
alpha_div_plane$WEEK = gsub("(\\d)(weeks)", "\\1 \\2", alpha_div_plane$WEEK) # add a space between the week number and the word WEEK

alpha_div_plane

alpha_div_plane$Irrigation<-gsub("Watered","Control",alpha_div_plane$Irrigation)
alpha_div_plane$Compartments<-gsub("r","R",alpha_div_plane$Compartments)


symnum.args2 <- list(cutpoints = c(0,  0.001, 0.01, 0.05, Inf), symbols = c( "***", "**", "*", "ns"))


Shan_plot_rp<-alpha_div_plane %>% ggplot(aes(x = Irrigation, y = Shannon, color = Irrigation))+
  geom_point() +
  facet_grid(~ WEEK, scales = "free_x") + #organise the plot by week
  theme_bw() +theme(axis.title.x =element_blank(), legend.position = "none", panel.grid = element_blank(),
                    axis.text = element_text(size = 10), axis.title = element_text (size = 12, face = "bold"),
                    strip.background = element_rect(fill="white", color = "black"),
                    strip.text = element_text(size =12))+
  scale_y_continuous(breaks = c(6,6.5,7,7.5), limits = c(6,7.5), expand = c(0,0))+
  stat_compare_means(aes(group = Irrigation), label = "p.signif",
                     method = "t.test", symnum.args = symnum.args2 , label.y =7.4, label.x = 1.4)+
  scale_color_manual(values=c("#018571","#a6611a"))


# ======================= Beta diversity ====================================
plot_pcoa = function(data, color,  Shape, filename){
  list_output = list() # to return the plot and the beta diversity results
  print("Computing PCoA")
  physeq.ord_pcoa = ordinate(data, "PCoA", "bray")
  print("Plotting PCoA")
  
  p3 = plot_ordination(data, physeq.ord_pcoa, shape = Shape) +
    scale_color_manual(values=c("#a6611a","#018571"),name = "Irrigation", labels = c("Drought","Control")) +
    geom_point(size = 4, aes(color = !!sym(color))) +  # Use aes with sym
    theme_bw() +
    theme(axis.text = element_text(size = 10, face = "bold"), axis.title = element_text(size = 12, face = "bold"), 
          legend.text = element_text(size = 11), legend.title = element_text(size = 13),
          panel.grid = element_blank())
  
   ggsave(filename)
  list_output[[1]] = p3
  list_output[[2]] = physeq.ord_pcoa
  return(list_output)
}

#Import data
p4<- plot_pcoa(Rarefied_physeq_mean_rounded_rhizoplane, "Irrigation",  "WEEK", "output/PCoA_NRPS_rp.png")
pcoa_NRPS_rp<-p4[[1]]

pcoa_NRPS_rp
#Extract the legend
pcoa_legend<-get_legend(pcoa_NRPS_rp)
leg_for_plot<-ggarrange(pcoa_legend)

p3<-pcoa_NRPS_rp+theme(legend.position = "none")
library(patchwork)  

all_plots_t<-Shan_plot_rp+p3+leg_for_plot+patchwork::plot_layout(widths = c(3.6,2,1))
all_plots_t[[3]]<-all_plots_t[[3]]+plot_layout(tag_level ='new')
all_plots_t+patchwork::plot_annotation(tag_levels = ("A"))

x11(width = 8, height = 4)
all_plots_t+patchwork::plot_annotation(tag_levels = ("A"))
ggsave("output/figure3_drought_rp.png", dpi = 300)
dev.off()



#=========================
# Plot richness for suppl


rich_plot_plane<-alpha_div_plane %>% ggplot(aes(x = Irrigation, y = Chao1))+geom_point() +
  facet_grid(~ WEEK,scales = "free_x") + #organise the plot by week
  theme_bw() +theme(axis.title.x =element_blank(), legend.position = "none", panel.grid = element_blank(),
                    axis.text = element_text(size = 10), axis.title = element_text (size = 14),
                    strip.background = element_rect(fill="white", color = "black"),
                    strip.text = element_text(size =12))+
  scale_y_continuous(breaks = c(2000,4000,6000,8000), limits = c(1700,8000), expand = c(0,0))+
  stat_compare_means(aes(group = Irrigation), label = "p.signif",
                     method = "t.test", symnum.args = symnum.args2 , label.y =7700, label.x = 1.4)

x11(height = 4, width =6)
rich_plot_plane

ggsave("output/Suppl_fig_Richness_rp_drought.png",dpi = 300)
dev.off()


# ============= Mantel test
# data imported in line 8
# import 16S rRNA phyloseq object
setwd("..");physeq16S_rare_mean_rounded = readRDS("16S_rRNA/output/Rarefied_physeq_mean_rounded.rds");setwd("NRPS")
Rarefied_physeq_mean_rounded_rhizoplane # 1 sample missing

sample_data(Rarefied_physeq_mean_rounded_rhizoplane)
sample_data(physeq16S_rare_mean_rounded)$Original_label<-gsub("_Wa","",sample_data(physeq16S_rare_mean_rounded)$Original_label)
sample_data(physeq16S_rare_mean_rounded)$Original_label<-gsub("p","P",sample_data(physeq16S_rare_mean_rounded)$Original_label)

sample_data(Rarefied_physeq_mean_rounded_rhizoplane)$Original_label<-gsub("_C","",sample_data(Rarefied_physeq_mean_rounded_rhizoplane)$Original_label)
sample_data(Rarefied_physeq_mean_rounded_rhizoplane)$Original_label<-gsub("_6","",sample_data(Rarefied_physeq_mean_rounded_rhizoplane)$Original_label)

setdiff(sample_data(physeq16S_rare_mean_rounded)$Original_label,sample_data(Rarefied_physeq_mean_rounded_rhizoplane)$Original_label)
#Remove sample 64_RP

new16S<-subset_samples(physeq16S_rare_mean_rounded, !Original_label =="64_RP")

sample_names(new16S)<-sample_data(new16S)$Original_label
sample_names(Rarefied_physeq_mean_rounded_rhizoplane)<-sample_data(Rarefied_physeq_mean_rounded_rhizoplane)$Original_label

library(vegan)

new16S_bc<-vegdist(t(otu_table(new16S)), method = "bray")
NRPS_bc<-vegdist(t(otu_table(Rarefied_physeq_mean_rounded_rhizoplane)), method = "bray")
mantel(new16S_bc,NRPS_bc)
# Mantel statistic r = 0.4531, sign = 0.001