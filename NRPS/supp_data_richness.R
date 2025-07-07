
# ===========================================================


# Library loading


# ===========================================================


library(phyloseq)
library(ampvis2)
library(ggplot2)
library(vegan)
library(dplyr)
library(tidyr)
library(ggpubr)
library(readxl)


# ===========================================================


# Prepare data


# ===========================================================

rarefied_physeq = readRDS("/home/be203133/Linking_microbial_community_structure_with_function/NRPS/ressources/RDS/Rarefied_physeq_mean_rounded_rhizoplane.rds")

# ===========================================================


# Richness


# ===========================================================
alpha_div_plot = function(data, x_var,y_var, col, facet, comparizon, nb_test, name){
  y_lab = ifelse(y_var == "Shannon", "Shannon diversity", "Chao1")
  symnum = list(cutpoints = c(0, 0.0001/nb_test, 0.001/nb_test, 0.01/nb_test, 0.05/nb_test, Inf), symbols = c("****", "***", "**", "*", "ns"))
  formule = as.formula(paste0("~", facet))

  p = ggplot(data = data, aes_string(x = x_var, y = y_var, color = col)) +
    geom_point(size = 4) +
    facet_grid(formule, scales = "free") +
    theme_minimal(base_family = "serif") +
    scale_color_manual(values = c("#018571", "#a6611a")) +
    theme(
      axis.text.x = element_text(hjust = 0.5, size = 15),           
      axis.title.x = element_text(hjust = 0.5, size = 18),        
      axis.text.y = element_text(size = 16),
      legend.title = element_text(size = 20),
      legend.text = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      strip.text = element_text(size = 15),
      strip.background = element_rect(fill = "white", color = "black")
    ) +
    ylim(c(min(data[[y_var]]) - 0.1 * min(data[[y_var]]), max(data[[y_var]]) + 0.1 * max(data[[y_var]]))) +
    labs(x = "", y = y_lab) +
    stat_compare_means(comparisons = comparizon, symnum.args = symnum,
                      label = "p.signif", size = 6, bracket.size = 0,
                      tip.length = 0, vjust = 0.5)
    ggsave(paste0("output/alpha_diversity_",name ,"_", y_var, ".png"), plot = p, width = 25, height = 15, units = "cm")
   return(p)
}

comparison = list(c("Control", "Drought")) # add another paired vector in th elist for further comparizon
nb_test = 2
alpha_div_plane = estimate_richness(rarefied_physeq, measures = c("Shannon", "Chao1")) # compute alpha diversity
alpha_div_plane = cbind(alpha_div_plane, rarefied_physeq@sam_data) # add the factor to the new table
alpha_div_plane$Irrigation = ifelse(alpha_div_plane$Irrigation == "Watered", "Control", "Drought")

alpha_div_plane = alpha_div_plane %>%
  mutate(WEEK = case_when(
    WEEK == "2weeks" ~ "week 2",
    WEEK == "4weeks" ~ "week 4",
    WEEK == "5weeks" ~ "week 6",
    TRUE ~ as.character(WEEK)
  ))

alpha_div_plot(data = alpha_div_plane, 
               x_var = "Irrigation",  
               y_var = "Chao1", 
               col = "Irrigation", 
               facet = "WEEK", 
               comparizon = comparison, 
               nb_test = nb_test, 
               name = "rhizoplane")
