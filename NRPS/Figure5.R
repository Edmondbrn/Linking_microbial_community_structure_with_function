# Script for plotting figure 5
library(ggplot2)
library(tidyverse)
library(phyloseq)
library(ampvis2)
# import Mibig blast results

all_mibigs<-read.table("ressources/RDS/full_blast_mibig_cluster2.txt", header = TRUE, sep ="\t", dec = ".")
colnames(all_mibigs)

# import the phyloseq object

#
physeq_raw<- readRDS("ressources/RDS/raw_physeq_rhizoplane.rds")
#Pull out relevant samples

physeq_control<-physeq_raw %>%  subset_samples(Treatment %in% "Control") 

amp_control<-amp_load(physeq_control)
amp_control$metadata

#Mining the data for CLPs based on Zhou et al. (2024) Frontiers in Bioengineering and Biotechnology
#Putisolvin (present), massetolide A (np), surfactin (p), amphisin (np), viscosinamide (p),
# factin (p) - subgroups exist, bananamide (p incl subgroups), mycin (incl subgroups), poaeaamide (np),
#orfamide (p), cocoyamide (np),asplenin (np), entolysin (p), xantholysin (np), tolaasin (p), peptin(p incl subgoups), sclerosin (np)

CLPs<-paste(c("putisolvin","viscosin","massetolide","surfactin","amphisin","factin","banana","orfamide","asplenin","entolysin",
        "tolaasin","peptin","sclerosin","xantholysin","cocoyamide","tensin","lokisin","milkisin", "arthrofactin",
        "anikasin","hodersin","nepenthensin","oakridgin","prosekin","pseudodesmin","pseudophomin",
        "syringotoxin","gacamide","sessilin","cichorinotoxin"),collapse ="|")

CLP_mibig<-all_mibigs[grepl(CLPs,all_mibigs$Compound),]

#Now it needs to be curated, remove compounds not on the list above
all_CLPs_unfil<-CLP_mibig %>% group_by(Compound) %>%
  summarise(n = n()) %>% arrange(desc(n))

#Remove
unique(all_CLPs_unfil$Compound)
CLP_remove<-paste(c("glidopeptin","anabaenopeptin","cinnapeptin","omnipeptin", "empedopeptin", "polyoxypeptin","micropeptin",
              "hypeptin","leupeptin", "incarnatapeptin","sarpentin","stechlisin","pelgipeptin", "livipeptin", "octapeptin"),
              collapse = "|")


CLP_mibig_red<-CLP_mibig[!grepl(CLP_remove,CLP_mibig$Compound),]
red_CLPs<-CLP_mibig_red %>% group_by(Compound) %>%
  summarise(n = n()) %>% arrange(desc(n)) #29 CLPs left in the data

CLP_ACs<-CLP_mibig_red$Query
#Subset the data as above using all of these


heatmap_CLPs_all<-amp_control %>% amp_subset_taxa(tax_vector = CLP_ACs,normalise = TRUE) %>%
  amp_heatmap(group_by = "Irrigation", facet_by ="WEEK",
              tax_aggregate = "Genus", normalise = FALSE, round = 3)
x11(width = 6, height = 5)
heatmap_CLPs_all
ggsave("output/FigS11_heatmap_CLPs_all.png")
dev.off()
#This suggests that CLPs are relatively more abundant in the watered plants 

#Plot Pseudomonas only
amp_control %>% amp_subset_taxa(tax_vector = CLP_ACs,normalise = TRUE) %>%
  amp_subset_samples(Compartments %in% "rhizoplane") %>% amp_subset_taxa(tax_vector = "Pseudomonas", normalise = FALSE) %>%
  amp_heatmap(group_by = "Irrigation", facet_by ="WEEK",
              tax_aggregate = "OTU", normalise = FALSE, round = 3, tax_show = 25) 

# Then I want to pull out the Pseudomonas
CLP_pse_table<-amp_control %>% amp_subset_taxa(tax_vector = CLP_ACs,normalise = TRUE) %>%
  amp_subset_samples(Compartments %in% "rhizoplane") %>% amp_subset_taxa(tax_vector = "Pseudomonas", normalise = FALSE) %>%
  amp_heatmap(group_by = "Irrigation", facet_by ="WEEK",
              tax_aggregate = "OTU", normalise = FALSE, round = 3, tax_show = 105, textmap = TRUE) 
CLP_pse_table$Query<-rownames(CLP_pse_table)

CLPs_abu_mibig<-left_join(CLP_pse_table,CLP_mibig_red,by = "Query")

CLPs_abu_mibig$Compound[CLPs_abu_mibig$Compound %in% "thanafactin A,"]<-"thanafactin A"

CLPs_abu_mibig$Compound[CLPs_abu_mibig$Compound %in% c("putisolvin III, putisolvin VI putisolvin V,")]<-"putisolvin"
CLPs_abu_mibig$Compound[CLPs_abu_mibig$Compound %in% c("bananamide F, bananamide D, bananamide G, bananamide E,")]<-"bananamide"
CLPs_abu_mibig$Compound[CLPs_abu_mibig$Compound %in% c("tolaasin A,","tolaasin I, tolaasin F,")]<-"Tolaasin"
CLPs_abu_mibig$Compound[CLPs_abu_mibig$Compound %in% c("viscosin,")]<-"Viscosin"

#Then plot them
CLP_mean_vals<-CLPs_abu_mibig %>% pivot_longer(cols = c(starts_with("Dr"), starts_with("Wa")),
                         names_to = "Week_treat",
                         values_to ="Rel_abu") %>% select(-E.value,-Info,-Identity...,-Positives...,-Gaps...) %>%
  group_by(Compound, Week_treat) %>%
  summarise(sum = sum(Rel_abu),
            n = n())
#Check that the sum is correct
CLP_mean_vals %>% group_by(Week_treat) %>%
  summarise(total = sum(sum))
#This is fine

CLP_mean_vals$Compound

#Before plotting I round the numbers
CLP_mean_vals$sum<-round(CLP_mean_vals$sum,4)

CLP_mean_vals$Week_treat<-as.character(CLP_mean_vals$Week_treat)
CLP_mean_vals$Week_treat<-factor(CLP_mean_vals$Week_treat, levels = c("Watered week 2","Watered week 4","Drought week 4", "Watered week 5","Drought week 5"))

# add the week numbers as a variable
CLP_means_full<-cbind(CLP_mean_vals,str_split_fixed(CLP_mean_vals$Week_treat," ", 2))
colnames(CLP_means_full)[5:6]<-c("Treatment","Week")
CLP_means_full$Week<-gsub("week","",CLP_means_full$Week)
CLP_means_full$Treatment[CLP_means_full$Treatment %in% "Watered"]<-"C"
CLP_means_full$Treatment[CLP_means_full$Treatment %in% "Drought"]<-"D"  
Pseu_clp_heat<-
ggplot(CLP_means_full, aes(x =Treatment, y = Compound, fill = sum)) +geom_tile(color = "white",
                                                                               lwd = 1.5,
                                                                               linetype = 1)+
  facet_grid(cols = vars(Week), scales = "free_x", space="free")+
  geom_text(aes(label = round(sum,2)), color = "black", size = 3) + 
  scale_fill_gradient(low = "White", high = "Brown")+
  theme(axis.text = element_text(size = 13),                       legend.position = "none",
                       axis.title = element_blank())+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background =  element_blank(),
                       strip.background = element_rect(fill="white", color = "black"))+ggtitle("Week")+
  theme(plot.title = element_text(hjust = 0.5))

x11(height = 4, width = 5.5)
Pseu_clp_heat
ggsave("output/Figure5_heatmap_CLPs_PSeu.png",dpi = 300)
dev.off()
