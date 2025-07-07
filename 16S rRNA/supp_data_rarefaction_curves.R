# ===========================================================


# Library loading


# ===========================================================

library(phyloseq)
library(ampvis2)
library(ggplot2)
library(vegan)
library(dplyr)
library(ggpubr)
library(vegan)


# ===========================================================


# Prepare data


# ===========================================================


# Import data
df_feature = read.csv("ressources/data/otu_table_merged_transpose.csv")
df_meta = read.csv("ressources/data/metadata1.csv")
df_tax = read.csv("ressources/data/tax_table_merged.csv", sep = ";")
   

df_feature = df_feature %>%
  tibble::column_to_rownames("ASV") # format the table for phyloseq

df_tax = df_tax %>% 
    filter(is.na(Family) | Family != "Mitochondria") %>% # remove mitochondria samples (no chloroplasta)
    tibble::column_to_rownames("ASV")

df_meta = df_meta %>% 
  mutate(Week = case_when( # replace number by week string
    Week == 2 ~ "week 2",
    Week == 4 ~ "week 4",
    Week == 5 ~ "week 6", # correct the typo in the metadata file
    TRUE ~ as.character(Week)
  )) %>%
  tibble::column_to_rownames("Sample_ID") 

# convert to phyloseq format
otu_mat = as.matrix(df_feature) # format again
tax_mat = as.matrix(df_tax)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(df_meta)
physeq = phyloseq(OTU, TAX, samples)
physeq = subset_taxa(physeq, Kingdom == "Bacteria") # remove non-bacterial sequences

# the line to filter the data to keep only the control samples
physeq = subset_samples(physeq, Treatment == "Water")
physeq@sam_data$Week = as.factor(physeq@sam_data$Week)
ampvis_data = amp_load(physeq)

# ===========================================================


# rarefaction curves


# ===========================================================

Rarecurve = function(data, col, type) { 
  cat("Creating rarefaction curve for ", type, " samples\n")
  p3 = amp_rarecurve(data, color_by = col, stepsize = 100, facet_by = "Treatment")
  p3 = p3 +
    theme(axis.text = element_text(size = 14), 
             axis.title = element_text(size = 20), 
             legend.text = element_text(size = 15), 
             legend.title = element_text(size = 15),
             strip.text.x = element_blank()) +
    theme_minimal(base_family = "serif") +
    theme(axis.text = element_text(angle = 0, hjust = 0.5, size = 14),
          axis.title = element_text(size = 20), 
          axis.text.y = element_text(angle = 0),
          
          legend.text = element_text(size = 18),  # Agrandir le texte du titre de chaque facette + # Mettre les graduations à l'horizontale et les agrandir
           legend.title = element_text(size = 20))+
    labs(y = "Number of ASVs")+
    geom_line(linewidth = 1.5)
  cat("Saving the plot\n")
 
  return(p3)
}

p1 = Rarecurve(data = ampvis_data, col = "Week", type = "Rhizoplane")
p5<-p1+theme(strip.text.x = element_blank())

ggsave("output/rarefaction_curve_control.png",p5 , width = 10, height = 10)
