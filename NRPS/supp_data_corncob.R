# ==========================================================================================================

# Load the necessary libraries

# ==========================================================================================================

library(corncob)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ampvis2)
library(patchwork)



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

plot_res_mibig = function(model, select_level, title, size){
    if (length(model$significant_taxa) ==  0)
        return("Error, no significant taxa found")
    df = plot(model, total = T, B = 1000, level = select_level, dataonly = T)
    DF = df$data
    DF$title = title
    DF$taxa = tolower(DF$taxa)
    DF = DF %>% arrange(desc(x)) # sort the table before the slicing

    
    if (length(grep("_", DF$taxa)) == 1){ # test if there is no ASV number in the data, it happens some time
        for (i in 1:length(DF$taxa)){
            DF$taxa[[i]] =  gsub(".*", paste0(DF$taxa[[i]] ," (", model$significant_taxa[[i]],")"), DF$taxa[[i]])
        }
    } 
    if (dim(DF)[1]>80){
        DF1 = DF[1:dim(DF)[1] / 2 +1,]
        DF2 = DF[(dim(DF)[1]  / 2 ):dim(DF)[1]+1,]
        # Divide the datframe into two
        # plot according to the treatment
        plot1 = plot_point_mibig(DF1, size)
        plot2 = plot_point_mibig(DF2, size)
        plot_merge = ggpubr::ggarrange(plot2, plot1, ncol = 2, common.legend = TRUE, legend = "right", widths = c(0.4, 0.4))
        print(plot_merge)
        ggsave(paste0("output/", title, ".png"), plot = plot_merge,width = 15, height = 10 )
    }
    else {
        plot = plot_point_mibig(DF, size)
        print(plot)
        ggsave(paste0("output/", title, ".png"), plot = plot, width = 15, height = 10)
    }
    return(DF)
}

plot_point_mibig = function(data, size, title){ # function to plot the graph
    data$taxa = gsub("amplicon_cluster", "AC", data$taxa)
    data$taxa = sapply(data$taxa, function(x) capitalize_first(x))
    data$taxa = gsub("ac", "AC", data$taxa)
    data$taxa = factor(data$taxa, levels = unique(data$taxa)) # sort the taxa for the plot    
    negative_x_index = which(data$x < 0)
    if (length(negative_x_index) > 0)
        y_position = data$taxa[min(negative_x_index)]
    else
        y_position = NA
    plot = ggplot(data, aes(x = x, y = taxa)) +
                geom_errorbar(aes(xmin = xmin, xmax = xmax), width = 0.2, color = "red") + 
                geom_point(color = "black") + 
                geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
                theme_linedraw() +
                # facet_grid(. ~ title)+
                labs(x = "", y = "") + # remove the axis labels
                theme(axis.text.y = element_text(angle = 0, hjust = 1, size = size, family = "serif"))+
                theme(axis.text.x = element_text(angle = 0, hjust = 1, size = size), strip.text = element_text(size = size+2, family = "serif"),
                      plot.margin = unit(c(0.5,0.5,0.5,1.5), "cm"))
    if (!is.na(y_position)) {
        plot = plot + geom_hline(yintercept = as.numeric(y_position) -0.5, color = "black")
    }
    return(plot)
}

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
    good_tax = model_4_5$significant_taxa # extract the significant taxa (amplicon_cluster)
    model_4_5$data@tax_table = model_4_5$data@tax_table[rownames(model_4_5$data@tax_table) %in% good_tax,] # filtering
    model_4_5$data@otu_table = model_4_5$data@otu_table[rownames(model_4_5$data@tax_table),]
    mibig_annotation = mibig_annotation[rownames(model_4_5$data@tax_table),]
    # add fields to the results
    tax_table_data = as.data.frame(model_4_5$data@tax_table)
    tax_table_data$Product = as.vector(mibig_annotation$Compound)
    tax_table_data$Organism = as.vector(mibig_annotation$Organism)
    tax_table_data$Product_clean = sapply(tax_table_data$Product, simplify_product)
    tax_table_data$Label <- paste(
    ifelse(is.na(tax_table_data$Organism), "Unknown", tax_table_data$Organism),
    ifelse(is.na(tax_table_data$Product_clean), "Unknown", tax_table_data$Product_clean),
    sep = " | "
    )
    tax_table_data$phylum <- tax_table_data$Label # because ampvis remove other column than classical one
    model_4_5$data@tax_table = tax_table(as.matrix(tax_table_data))
    physeq_mibig = model_4_5$data


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
        plot_values = TRUE,
        normalise = FALSE,
        tax_show = 25,
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

}

# ==========================================================================================================

# Load data

# ==========================================================================================================

physeq = readRDS("ressources/RDS/raw_physeq_rhizoplane.rds")
mibig_annotation = read.csv("ressources/data/full_cluster_mibig_blast.txt", sep = "\t", row.names = 1)
model_2_4 = readRDS("ressources/RDS/rhizoplanes_week2_4_WEEK.rds")
model_4_5 = readRDS("ressources/RDS/rhizoplanes_week4_5_WEEK.rds")

generate_heatmap(model_2_4, "output/week_2_vs_4_abundance.jpeg")
generate_heatmap(model_4_5, "output/week_4_vs_5_abundance.jpeg")
