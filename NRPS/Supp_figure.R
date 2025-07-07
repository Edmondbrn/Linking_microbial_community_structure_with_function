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
# ==========================================================================================================

# Load the data

# ==========================================================================================================
# setwd("/home/be203133/Linking_microbial_community_structure_with_function/NRPS")
try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path))) # only works on RStudio
model = readRDS("ressources/RDS/rhizoplanes_week5_irrigation.rds")
mibig_annotation = read.csv("ressources/data/full_cluster_mibig_blast.txt", sep = "\t", row.names = 1)
physeq = readRDS("ressources/RDS/raw_physeq_rhizoplane.rds")


# ==========================================================================================================

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


# ============================================================================================================================
# script
# ============================================================================================================================

plot_res_mibig(model, "genus", "Week 5 irrigation_taxa_name", 14)


good_tax = model$significant_taxa # extract the significant taxa (amplicon_cluster)
model$data@tax_table = model$data@tax_table[rownames(model$data@tax_table) %in% good_tax,] # filtering
model$data@otu_table = model$data@otu_table[rownames(model$data@tax_table),]
mibig_annotation = mibig_annotation[rownames(model$data@tax_table),]

tax_table_data = as.data.frame(model$data@tax_table)
tax_table_data$Product = as.vector(mibig_annotation$Compound)
tax_table_data$Organism = as.vector(mibig_annotation$Organism)
tax_table_data$Product_clean = sapply(tax_table_data$Product, simplify_product)


tax_table_data$phylum = tax_table_data$Product_clean # because ampvis remove other column than classical one
model$data@tax_table = tax_table(as.matrix(tax_table_data))

DF_mibig = plot_res_mibig(model, "phylum", "Week 5 irrigation_mibig_product", 14)
physeq_mibig = model$data
ampvis_mibig = amp_load(physeq_mibig)
saveRDS(ampvis_mibig, "ressources/RDS/ampvis2_rhizoplanes_week5_irrigation_mibig.rds")

DF_mibig3 = DF_mibig  # change the ordrer for taxa to x
DF_mibig3$taxa_clean = gsub("amplicon_cluster", "AC",DF_mibig3$taxa) # replace amplicon_cluster by AC
DF_mibig3 = DF_mibig3 %>% mutate(taxa_clean = fct_reorder(taxa_clean, x)) # reorder the taxa by x
DF_mibig3$OTU

#Subset the ampvis object
DF_mibig3$OTU<-str_split_fixed(DF_mibig3$taxa," \\(",2) [,2]
DF_mibig3$OTU<-gsub("\\)","",DF_mibig3$OTU)

DF_mibig3$OTU<-factor(reorder(DF_mibig3$OTU,desc(DF_mibig3$x)))

# DF_mibig3$genus = as.data.frame(physeq4week@tax_table)$genus[rownames(physeq4week@tax_table) %in% DF_mibig3$OTU] # dom2BGC
DF_mibig3$genus = strsplit(mibig_annotation$Organism[rownames(mibig_annotation) %in% DF_mibig3$OTU], " ") # mibig
DF_mibig3$genus = sapply(DF_mibig3$genus, function(x) x[1])


distinct_colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "yellow", "black")


# normalisation from the original phyloseq object (with all ASVs)
otu_physeq = physeq@otu_table
otu_physeq = transform_sample_counts(otu_physeq, function(x) x / sum(x) * 100) # normalisation to 100%
physeq_propor = otu_table(otu_physeq)
physeq_propor = merge_phyloseq(physeq_propor, tax_table(physeq), sample_data(physeq))
ampvis_total = amp_load(physeq_propor)

amp_mibig_un = ampvis_total %>% 
    amp_subset_taxa(tax_vector = DF_mibig3$OTU, normalise = FALSE) # already normalised

amp_mibig_un$abund
heat_mibig <- amp_mibig_un %>% 
    amp_heatmap( 
        group_by = "Irrigation",
        tax_aggregate = "OTU",
        
        facet_by = "Irrigation",
        
        plot_values = T,
        normalise = F,
        tax_show = 25,
        color_vector = c("royalblue4", "whitesmoke", "darkred")) +
        theme(strip.background = element_rect(fill = "darkgray"),
                strip.text = element_text(color = "black", size = 17),
                axis.text = element_text(family = "serif"),
                axis.text.x = element_blank(),
                axis.text.y = element_text(size = 12),
                legend.title = element_text(size = 15),
                legend.text = element_text(size = 12),
                panel.spacing = unit(0.1, "lines"),
                axis.ticks.length = unit(0.1, "cm"),
                plot.margin = margin(1, 10, 1, 10, "cm")) + # adjust the margin to reduce the heatmap
        scale_x_discrete(expand = c(0.01, 0.01)
    )
ggsave("heat_mibig.jpeg", heat_mibig)


taxa_o<-heat_mibig$data$Display
taxa_o<-levels(droplevels(taxa_o))

taxa_re <- factor(levels(DF_mibig3$OTU))
amp_mibig_un$metadata$Irrigation = ifelse(amp_mibig_un$metadata$Irrigation == "Drought", "Drought", "Control")
amp_mibig_un$metadata$Irrigation = factor(amp_mibig_un$metadata$Irrigation, levels = c("Drought", "Control")) # reorder to have control on the right

DF_mibig3$taxa_clean = lapply(as.character(DF_mibig3$taxa_clean), function(x) capitalize_first(x))
DF_mibig3$taxa_clean = gsub("ac", "AC", DF_mibig3$taxa_clean)
DF_mibig3  = DF_mibig3 %>% arrange(x)
DF_mibig3$taxa_clean = factor(DF_mibig3$taxa_clean, levels = unique(DF_mibig3$taxa_clean)) # sort the taxa for the plot

diff_plot<-DF_mibig3  %>% ggplot(aes(x = x, y = taxa_clean)) + geom_point(color = "black")+
  geom_errorbar(aes(xmin = xmin, xmax = xmax, color = genus), width = 0.2) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  geom_hline(yintercept = 11.5, linetype = "solid", color = "black") + 
  theme_linedraw() +
  theme_bw()+
  labs(x = "", y = "") + # remove the axis labels
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 15, family = "serif"),
        panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 15), strip.text = element_text(size = 20, family = "serif"),
        plot.margin = unit(c(0.5,0.5,0.5,1.5), "cm"),
        legend.position = "bottom",
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18))+
    scale_color_manual(values = distinct_colors) +
    guides(color = guide_legend(nrow = 5))



heat_forcom<-amp_mibig_un%>% amp_heatmap( 
  group_by = "Irrigation",
  tax_aggregate = "OTU",
  #tax_add ="Phylum",
  facet_by = "Irrigation",
  order_y_by = rev(taxa_re),
  plot_values = T,
  normalise = F,
  tax_show = 25,
  round = 3,
  color_vector = c("lightblue", "white", "red")) +
  theme(strip.background = element_rect(fill = "lightgray"),
        strip.text = element_text(color = "black", size = 17),
        axis.text = element_text(family = "serif"),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), #element_text(size = 12),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        panel.spacing = unit(0.1, "lines"),
        axis.ticks.length = unit(0.1, "cm"),
        plot.margin = margin(1, 1, 1, 1, "cm")) + # adjust the margin to reduce the heatmap
  scale_x_discrete(expand = c(0.01, 0.01))+
    scale_fill_gradient(low = "White", high = "Brown")


    

diff_plot <- diff_plot + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
heat_forcom <- heat_forcom | theme(plot.margin = margin(0, 0, 0, 0, "cm"))

# Combine plot with patchwork
combined_plot <- diff_plot + heat_forcom + plot_layout(ncol = 2, widths = c(1, 1))

# display combined plot
x11(width = 12, height = 8)
combined_plot
ggsave("output/diff_plot_combin_heat_week5.png", plot = combined_plot, dpi = 300, width = 12, height = 10)
dev.off()
