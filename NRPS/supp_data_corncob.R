# ==========================================================================================================

# Load the necessary libraries

# ==========================================================================================================

library(corncob)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ampvis2)

# ===========================================================================================================

# Functions

# ==========================================================================================================

binom = function(data, variable){
    # create the formula
    formula_str = paste0("~ ", variable)
    formula = as.formula(formula_str)
    phi_formula = formula
    # perform the analysis
    analysis = differentialTest(formula = formula, phi.formula = phi_formula, data = data,
                                formula_null = ~ 1, phi.formula_null = formula, test = "Wald", boot = FALSE,
                                fdr = "fdr", fdr_cutoff = 0.05, verbose = T)    
    return(analysis)
}


# function to simplify the product name from mibig
simplify_product <- function(product) {
  # divide string by coma
  products = unlist(strsplit(product, ","))
  # extract the ifrst part before a coma
  base_product = trimws(products[1])
  words = unlist(strsplit(base_product, " "))
  # if strin has more than 2 words, return only the two first ones
  if (length(words) > 2 & words[1] != words[2])
    base_product = paste(words[1:2], collapse = " ")
  else if (length(words) > 2 & words[1] == words[2])
    base_product = words[1]
  return(base_product)
}
# put an upper case only to the first character of a string
capitalize_first <- function(word) {
    # function to put a capital lettre at the beginning of a word
    word_list = strsplit(word, "_")[[1]] # Get the first element of the result list
    word_upper = c()
    for (i in 1:length(word_list)){ # browse all the subwords
        if (nchar(word_list[i]) > 0)
            word_upper = c(word_upper, paste0(toupper(substring(word_list[i], 1, 1)), 
                                             tolower(substring(word_list[i], 2)))) 
        else
            word_upper = c(word_upper, word_list[i])
    }
    word = paste(word_upper, collapse = "_") # reconstruct the word
    word = gsub("asv", "ASV", word) # put ASV in capital
    return(word)
}


generate_heatmap = function(model, file_name) {
    good_tax = model$significant_taxa # extract the significant taxa (amplicon_cluster)
    model$data@tax_table = model$data@tax_table[rownames(model$data@tax_table) %in% good_tax,] # filtering
    model$data@otu_table = model$data@otu_table[rownames(model$data@tax_table),]
    mibig_annotation = mibig_annotation[rownames(model$data@tax_table),]
    # add fields to the results
    tax_table_data = as.data.frame(model$data@tax_table)
    tax_table_data$Product = as.vector(mibig_annotation$Compound)
    tax_table_data$Product_clean = sapply(tax_table_data$Product, simplify_product)
    # Create a label column for the y-lab plot
    tax_table_data$Label <- paste(
    ifelse(is.na(tax_table_data$genus), "Unknown", tax_table_data$genus),
    ifelse(is.na(tax_table_data$Product_clean), "Unknown", tax_table_data$Product_clean),
    sep = " | "
    )
    tax_table_data$phylum <- tax_table_data$Label # because ampvis remove other column than classical one
    # updtate the phyloseq object
    model$data@tax_table = tax_table(as.matrix(tax_table_data))
    physeq_mibig = model$data


    # normalisation from the original phyloseq object (with all ASVs)
    otu_physeq = physeq@otu_table
    otu_physeq = transform_sample_counts(otu_physeq, function(x) x / sum(x) * 100) # normalisation to 100%
    physeq_propor = otu_table(otu_physeq)
    physeq_propor = merge_phyloseq(physeq_propor, tax_table(physeq_mibig), sample_data(physeq_mibig))
    ampvis_model = amp_load(physeq_propor)
    # change week string format
    ampvis_model$metadata = ampvis_model$metadata %>%
    mutate(WEEK = case_when(
        WEEK == "2weeks" ~ "week 2",
        WEEK == "4weeks" ~ "week 4",
        WEEK == "5weeks" ~ "week 6",
        TRUE ~ as.character(WEEK)
    ))


    heat_mibig <- ampvis_model %>%
    amp_heatmap(
        group_by = "WEEK",
        tax_aggregate = "Phylum", # dummy column to store product names and organism names
        facet_by = "WEEK",
        round = 2,
        plot_values = TRUE,
        normalise = FALSE, # hand made normalisation before
        tax_show = 30,
        color_vector = c("#018571", "whitesmoke", "#a6611a")
    ) +
    theme_minimal(base_family = "serif") +
    theme(
        axis.text.x = element_text(hjust = 1, size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black", face = "italic"),
        axis.title.x = element_text(size = 18, hjust = 0.5, face = "bold"),
        axis.title.y = element_text(size = 18, hjust = 0.5, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        strip.text = element_text(size = 16, face = "bold", color = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
    ggsave(file_name, heat_mibig, width = 12, height = 10)
    return (heat_mibig)

}

# ==========================================================================================================

# Load data

# ==========================================================================================================

physeq = readRDS("ressources/RDS/raw_physeq_rhizoplane.rds")
mibig_annotation = read.csv("ressources/data/full_cluster_mibig_blast.txt", sep = "\t", row.names = 1)
model_2_4 = readRDS("ressources/RDS/plane_WC_week2_4.rds")
model_4_5 = readRDS("ressources/RDS/plane_WC_week4_5.rds")
model_2_5 = readRDS("ressources/RDS/plane_WC_week2_5.rds")
model_4_5_drought = readRDS("ressources/RDS/plane_DC_week4_5.rds")

p1 = generate_heatmap(model_2_4, "output/S9_figure.png")
p2 = generate_heatmap(model_4_5, "output/S10_figure.png")
# extra plots
p3 = generate_heatmap(model_2_5, "output/week_2_vs_5_abundance.png")
p4 = generate_heatmap(model_4_5_drought, "output/week_4_vs_5_abundance_drought.png")
