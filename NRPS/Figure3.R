# Script for Figure 3
library(tidyverse)
library(ggplot2)
library(phyloseq)
plot_annotation = function(data, type){
  colors = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3")
  
  ggplot(data = data, aes(x = 2, y = Freq, fill = Var1))+
    geom_bar(stat = "identity", width = 1, color = "black") +
    coord_polar(theta = "y") + # line to remove to have a stacked barplot
    xlim(0.5, 2.5) +
    geom_text(aes(label = paste0(round(Freq, 3), "%")), 
              position = position_stack(vjust = 0.5)) +
    theme_void() +
    theme(legend.position = "right", legend.text = element_text(size = 15)) +
    labs(fill = "ASV")+
    scale_fill_manual(values = colors) # Utilisation de la palette de couleurs
  
  
}

compute_freq = function (data){ # fonction to select and compute the proportion of annotated clusters in the data
  asv_table = otu_table(data)
  asv_table = asv_table[rowSums(asv_table)>0,]
  asv_present = rownames(asv_table)
  taxonomy_table = tax_table(data)
  taxonomy_asv = as.data.frame(taxonomy_table[asv_present, ])
  annot = as.data.frame(prop.table(table(taxonomy_asv$superkingdom)) * 100)
  return(annot)
}
#Import data
physeq_rhizoplane <- readRDS("ressources/RDS/Rarefied_physeq_mean_rounded_rhizoplane.rds")
physeq_rhizoplane_drought<- subset_samples(physeq_rhizoplane, Irrigation %in% "Drought")
physeq_rhizoplane_watered<-subset_samples(physeq_rhizoplane, Irrigation %in% "Watered")

annot_drought<-compute_freq(physeq_rhizoplane_drought)
annot_watered<-compute_freq(physeq_rhizoplane_watered)

# Calculate product
compute_freq_product = function (data){ # function to select and compute the proportion of annotated clusters in the data
  asv_table = otu_table(data)
  asv_table = asv_table[rowSums(asv_table)>0,]
  asv_present = rownames(asv_table)
  taxonomy_table = tax_table(data)
  taxonomy_asv = as.data.frame(taxonomy_table[asv_present, ])
  annot = as.data.frame(prop.table(table(taxonomy_asv$product)) * 100)
  return(annot)
}


compute_freq_taxa_of_assigned = function (data,name){ # function to select and compute the proportion of annotated clusters in the data
  asv_table = otu_table(data)
  asv_table = asv_table[rowSums(asv_table)>0,]
  asv_present = rownames(asv_table)
  taxonomy_table = tax_table(data)
  taxonomy_asv = as.data.frame(taxonomy_table[asv_present, ])
  annot = as.data.frame(prop.table(table(taxonomy_asv$phylum)) * 100)
  colnames(annot) = c("Phylum",name)
  return(annot)
}

physeq_rp_dr_assi<-subset_taxa(physeq_rhizoplane_drought, !product == "unassigned" ) 
physeq_rp_dr_assi<-prune_taxa(taxa_sums(physeq_rp_dr_assi) > 0, physeq_rp_dr_assi)
sample_data(physeq_rp_dr_assi)
physeq_rp_dr_assi_w4<-physeq_rp_dr_assi %>% subset_samples(WEEK %in% "4weeks") 
physeq_rp_dr_assi_w5<-physeq_rp_dr_assi %>% subset_samples(!WEEK %in% "4weeks")
annot_phyl_rp_dr_w4<-compute_freq_taxa_of_assigned(physeq_rp_dr_assi_w4,"W4 D") 
annot_phyl_rp_dr_w5<-compute_freq_taxa_of_assigned(physeq_rp_dr_assi_w5,"W5 D")

physeq_rp_ww_assi<-subset_taxa(physeq_rhizoplane_watered, !product == "unassigned" ) 

physeq_rp_ww_assi <-physeq_rp_ww_assi %>% subset_samples(Treatment %in% "Control") 
physeq_rp_ww_assi <-prune_taxa(taxa_sums(physeq_rp_ww_assi ) > 0, physeq_rp_ww_assi )

physeq_rp_ww_assi_w4<-physeq_rp_ww_assi %>% subset_samples(WEEK %in% "4weeks") 
physeq_rp_ww_assi_w5<-physeq_rp_ww_assi %>% subset_samples(!WEEK %in% "4weeks")
annot_phyl_rp_ww_w4<-compute_freq_taxa_of_assigned(physeq_rp_ww_assi_w4,"W4 C") 
annot_phyl_rp_ww_w5<-compute_freq_taxa_of_assigned(physeq_rp_ww_assi_w5,"W5 C")

#Include week 2
physeq_rp_ww_assi_w2<-physeq_rp_ww_assi %>% subset_samples(WEEK %in% "2weeks") 
annot_phyl_rp_ww_w2<-compute_freq_taxa_of_assigned(physeq_rp_ww_assi_w2,"W2 C") 

#Combined verything for easier plotting
list_dfs<-list(annot_phyl_rp_dr_w4,annot_phyl_rp_dr_w5,annot_phyl_rp_ww_w4,annot_phyl_rp_ww_w5, annot_phyl_rp_ww_w2)
combined_dfs<-list_dfs %>% reduce(full_join, by='Phylum') %>% pivot_longer(!Phylum,names_to = "Treatment",values_to = "Proportion")

combined_dfs # Proportion of NRPS for each phylum

combined_dfs$Phylum<-gsub("Pseudomonadota","Proteobacteria",combined_dfs$Phylum)

annotated_proportions<-ggplot(combined_dfs, aes(x = Treatment, y = Proportion, fill = Phylum)) + geom_bar(position="fill", stat="identity")+
  theme_bw()+
  theme(legend.position = "right", legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8)) +
  scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values =   c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"),
                                   labels = c("Actinomycetota (64-73%)",
                                              "Bacillota (0.32-0.84%)",
                                              "Gemmatimonadota (0.11-0.28%)",
                                              "Myxococcota (8.1-13%)",
                                              "Plantomycetota (0.045%)",
                                              "Proteobacteria (18-24%)",
                                              "Verrucomicrobiota (0.14-0.22%)"))
  


#======================================
  
  # Heatmap
  
physeq_rhizoplane

library(ampvis2)
Rarefied_amp_mean_rp<-amp_load(physeq_rhizoplane)

#Change Treatment variable to control

Rarefied_amp_mean_rp$metadata$Irrigation[Rarefied_amp_mean_rp$metadata$Irrigation %in% "Watered"]<-"C"
Rarefied_amp_mean_rp$metadata$Irrigation[Rarefied_amp_mean_rp$metadata$Irrigation %in% "Drought"]<-"D"
Rarefied_amp_mean_rp$metadata$WEEK[Rarefied_amp_mean_rp$metadata$WEEK %in% "2weeks"]<-"Week 2"
Rarefied_amp_mean_rp$metadata$WEEK[Rarefied_amp_mean_rp$metadata$WEEK %in% "4weeks"]<-"Week 4"
Rarefied_amp_mean_rp$metadata$WEEK[Rarefied_amp_mean_rp$metadata$WEEK %in% "5weeks"]<-"Week 5"

Rarefied_amp_mean_rp$tax$OTU<-gsub("amplicon_cluster","AC",Rarefied_amp_mean_rp$tax$OTU)
heatmap_rp<-Rarefied_amp_mean_rp%>% amp_heatmap(group_by = "Irrigation", facet_by = "WEEK",
                                                tax_show = 25, tax_aggregate = "OTU",
                                                tax_add = "Genus", plot_values_size = 3,
                                                color_vector = c("White", "Brown")) + 
  theme(axis.text.y = element_text(size=10),strip.background = element_rect(fill="white", color = "black"),
        strip.text = element_text(size =12), axis.text.x = element_text(angle = 0, hjust = 0.5))


library(ggpubr)
x11(width = 7, height = 8)
ggarrange(annotated_proportions,heatmap_rp, ncol = 1 , labels = c("A","B"), heights = c(0.4,0.7))

ggsave("output/Fig3.png",dpi = 300)
dev.off()
