# ==========================================================================================================

# Load the necessary libraries

# ==========================================================================================================

library(corncob)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(ampvis2)
library(ggpubr)
library(forcats)
library(patchwork)
library(stringr)
# ==========================================================================================================

# Load the data

# ==========================================================================================================
try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path))) # only works on RStudio
physeq = readRDS("ressources/RDS/raw_physeq_rhizoplane.rds")
View(physeq@sam_data)

physeq = subset_samples(physeq, Treatment == "Control") # remove control samples
physeq4week = subset_samples(physeq, WEEK == "week 4") # subset the data to only include 4-week samples


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
    data$taxa = factor(data$taxa, levels = unique(data$taxa)) # sort the taxa for the plot    
    negative_x_index = which(data$x < 0)
    if (length(negative_x_index) > 0)
        y_position = data$taxa[max(negative_x_index)]
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
        plot = plot + geom_hline(yintercept = as.numeric(y_position) +0.5, color = "black")
    }
    return(plot)
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


# Function to extract the cluster sequences from a model for significant taxa (run blastP on the sequences)
extract_cluster = function(model){
    cluster_list = model$significant_taxa # extract the cluster names
    seq_cluster = readAAStringSet("/home/edmond/NRPS_analysis/data/NRPS/amplicon_cluster.txt") # get the association between the cluster and the ASV sequences
    nb_cluster = length(cluster_list)
    top_cluster = AAStringSet() # Initialize an empty DNAStringSet object
    cpt = 1
    for (i in cluster_list){
        top_cluster = c(top_cluster, seq_cluster[names(seq_cluster) == i, ])
    }
    print(top_cluster)
    # ask to tu user the file name
    file_name = readline(prompt = "Enter the file name (without / or special characters): ")
    # # Write the combined ASV sequences to a single fasta file
    writeXStringSet(top_cluster, file = paste0("control_samples_analysis/output_sequence_corncob/", file_name, ".fasta"))
}

# ============================================================================================================================
# script
# ============================================================================================================================


# model = binom(physeq4week, "Irrigation")
# saveRDS(model, "output/rhizoplane_week4_irrigation.rds")

model = readRDS("ressources/RDS/corncob_rhizoplane_week4_irrigation.rds")
df_blast_mibig = read.csv("ressources/data/full_cluster_mibig_blast.txt", sep = "\t", row.names = 1)

# plot the raw results from Dom2BGC
plot_res_mibig(model, "genus" ,"Rhizoplane week 4 irrigation", 11) # plot the results

good_tax = model$significant_taxa # extract the significant taxa (amplicon_cluster)
model$data@tax_table = model$data@tax_table[rownames(model$data@tax_table) %in% good_tax,] # filtering
model$data@otu_table = model$data@otu_table[rownames(model$data@tax_table),]
df_blast_mibig = df_blast_mibig[rownames(model$data@tax_table),]

tax_table_data = as.data.frame(model$data@tax_table)
tax_table_data$Product = as.vector(df_blast_mibig$Compound)
tax_table_data$Organism = as.vector(df_blast_mibig$Organism)
tax_table_data$Product_clean = sapply(tax_table_data$Product, simplify_product)


tax_table_data$phylum = tax_table_data$Product_clean # because ampvis remove other column than classical one
model$data@tax_table = tax_table(as.matrix(tax_table_data))

physeq_mibig = model$data
ampvis_mibig = amp_load(physeq_mibig)
#saveRDS(ampvis_mibig, "ressources/RDS/ampvis2_rhizoplanes_week4_irrigation_mibig.rds")


DF_mibig = plot_res_mibig(model, "Product_clean","Rhizoplane week 4 irrigation", 9) 
DF_mibig$x
model$significant_taxa
tax_table_df = as.data.frame(model$data@tax_table)
products_double = tax_table_df %>% count(Product_clean) %>% filter(n > 1, Product_clean != "NaN") # get the compound with at least two occurences
product_twice = tolower(products_double$Product_clean) # for the selection



loop = dim(DF_mibig)[1]
keep_index = c()
compound = c()
for(i in 1:loop){
    trim_product <- trimws(strsplit(DF_mibig$taxa[i], "\\(")[[1]][1])
    print(trim_product)
    compound = c(compound, trim_product)
    if (trim_product %in% product_twice){
        keep_index = c(keep_index, i) # index with a double compounds
    }
}

DF_mibig$simple_product = compound

DF_mibig2 = DF_mibig[keep_index,] # remove unique compounds


DF_mibig3<-DF_mibig2  # change the ordrer frol taxa to x
DF_mibig3$taxa_clean = gsub("amplicon_cluster", "AC",DF_mibig3$taxa) # replace amplicon_cluster by AC
DF_mibig3 = DF_mibig3 %>% mutate(taxa_clean = fct_reorder(taxa_clean, x)) # reorder the taxa by x
DF_mibig3$OTU

mibig_annotation = read.csv("ressources/data/full_cluster_mibig_blast.txt", sep = "\t")

#Subset the ampvis object
DF_mibig3$OTU<-str_split_fixed(DF_mibig3$taxa," \\(",2) [,2]
DF_mibig3$OTU<-gsub("\\)","",DF_mibig3$OTU)

DF_mibig3$OTU<-factor(reorder(DF_mibig3$OTU,desc(DF_mibig3$x)))

# DF_mibig3$genus = as.data.frame(physeq4week@tax_table)$genus[rownames(physeq4week@tax_table) %in% DF_mibig3$OTU] # dom2BGC
DF_mibig3$genus = strsplit(mibig_annotation$Organism[mibig_annotation$Query %in% DF_mibig3$OTU], " ") # mibig
DF_mibig3$genus = sapply(DF_mibig3$genus, function(x) x[1])


distinct_colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink", "cyan", "yellow", "black")

diff_plot<-DF_mibig3  %>% ggplot(aes(x = x, y = taxa_clean)) + geom_point(color = "black")+
  geom_errorbar(aes(xmin = xmin, xmax = xmax, color = genus), width = 0.2) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
  geom_hline(yintercept = 17.5, linetype = "solid", color = "black") + 
  theme_linedraw() +
  theme_bw()+
  labs(x = "", y = "") + # remove the axis labels
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 15, family = "serif"),
        panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 15), strip.text = element_text(size = 20, family = "serif"),
        plot.margin = unit(c(0.5,0.5,0.5,1.5), "cm"),
        legend.position = "bottom",
        legend.justification = c(1,1),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 18))+
    scale_color_manual(name = "Genus",values = distinct_colors, labels = c("Amycolatopsis","Bacillus","Corallococcus",
                                                                           "Cystobacter","Halomonas","Pseudomonas","Streptomyces","Unknown")) +
    guides(color = guide_legend(nrow = 4)) # 4 row for the legend



# normalisation from the original phyloseq object (with all ASVs)
otu_physeq = physeq@otu_table
otu_physeq = transform_sample_counts(otu_physeq, function(x) x / sum(x) * 100) # normalisation to 100%
physeq_propor = otu_table(otu_physeq)
physeq_propor = merge_phyloseq(physeq_propor, tax_table(physeq), sample_data(physeq))
ampvis_total = amp_load(physeq_propor)

amp_mibig_un = ampvis_total %>% 
    amp_subset_taxa(tax_vector = DF_mibig3$OTU, normalise = FALSE) # already normalised
# amp_mibig_un<-ampvis_mibig %>% 
#   amp_subset_taxa(tax_vector = DF_mibig3$OTU, normalise = TRUE) 

amp_mibig_un$abund
heat_mibig<-amp_mibig_un%>% amp_heatmap( 
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
  scale_x_discrete(expand = c(0.01, 0.01))

taxa_o<-heat_mibig$data$Display
unique(taxa_o)
taxa_o<-levels(droplevels(taxa_o))


levels(DF_mibig3$OTU)
taxa_re[order(ordered(taxa_re, levels =taxa_re))]
taxa_re<-factor(levels(DF_mibig3$OTU))
amp_mibig_un$metadata$Irrigation = ifelse(amp_mibig_un$metadata$Irrigation == "Drought", "Drought", "Control")
amp_mibig_un$metadata$Irrigation = factor(amp_mibig_un$metadata$Irrigation, levels = c("Drought", "Control")) # reorder to have control on the right

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
ggsave("output/diff_plot_combin_heat.png", plot = combined_plot, dpi = 300, width = 12, height = 10)
dev.off()



