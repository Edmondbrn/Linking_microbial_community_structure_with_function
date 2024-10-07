
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
physeq = readRDS("ressources/RDS/physeq_ctrl.rds")

physeq4week = subset_samples(physeq, Week == "4") # subset the data to only include 4-week samples
physeq4week_genus_glom = tax_glom(physeq4week, "Family") # agglommerate the data by genus


physeq_4_5week = subset_samples(physeq, Week %in% c("4", "5")) # subset the data to only include 4-week and 5-week samples
physeq_4_5week_watered = subset_samples(physeq_4_5week, Irrigation == "Watered") # subset the data to only include watered samples
physeq_4_5week_watered_glom = tax_glom(physeq_4_5week_watered, "Family") # glom the data by genus


model = binom(physeq4week, "Irrigation") # perform the analysis
model_glom = binom(physeq4week_genus_glom, "Irrigation") # perform the analysis on the glommed data

model4_5_water = binom(physeq_4_5week_watered, "Week") # perform the analysis on the watered samples
model4_5_water_glom = binom(physeq_4_5week_watered_glom, "Week") # perform the analysis on the glommed watered samples
model4_5 = binom(physeq_4_5week, "Week") # perform the analysis on the 4-week and 5-week samples
saveRDS(model, "ressources/RDS/model_4week_irrigation.rds") # save the model
saveRDS(model_glom, "ressources/RDS/model_4week_irrigation_glom.rds") # save the model
saveRDS(model4_5_water, "ressources/RDS/model_4_5week_watered.rds") # save the model
saveRDS(model4_5_water_glom, "ressources/RDS/model_4_5week_watered_glom.rds") # save the model
saveRDS(model4_5, "ressources/RDS/model_4_5week.rds") # save the model
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

plot_point_16S = function(data, size, title, notaglom){ # function to plot the graph
    if (!notaglom)
        data$taxa = sub("\\(sequence.*", "", data$taxa)
    else
        data$taxa = gsub("sequence", "ASV", data$taxa)
    data = data %>% arrange(x) # sor the value regading x point position
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
                theme(axis.text.y = element_text(angle = 0, hjust = 1, size = size, family = "serif"),
                      axis.title.y = element_blank())+
                theme(axis.text.x = element_text(angle = 0, hjust = 1, size = size), strip.text = element_text(size = size+2, family = "serif"),
                      plot.margin = unit(c(0.5,0.5,0.5,1.5), "cm"))
    if (!is.na(y_position)) {
        plot = plot + geom_hline(yintercept = as.numeric(y_position) +0.5, color = "black")
    }
    return(plot)
}



plot_res_16S = function(model, select_level , title, size, notaglom = TRUE){
    if (length(model$significant_taxa) ==  0)
        return("Error, no significant taxa found")
    df = plot(model, total = T, B = 1000, level = select_level, dataonly = T)
    DF = df$data
    DF$title = title
    DF$taxa = tolower(DF$taxa)
    split_taxa = strsplit(DF$taxa, "-")
    for (i in 1:length(split_taxa)){ # remove unwanted characters after a -
        if (length(split_taxa[[i]]) > 2)
            DF$taxa[[i]] = split_taxa[[i]][length(split_taxa[[i]])]
    }
    
    if (length(grep("_", DF$taxa)) == 1){ # test if there is no ASV number in the data, it happens some time
        for (i in 1:length(DF$taxa)){
            DF$taxa[[i]] =  gsub(".*", paste0(DF$taxa[[i]] ," (", model$significant_taxa[[i]],")"), DF$taxa[[i]])
        }
    } 
    if (dim(DF)[1]>80){
        DF = DF %>% arrange(taxa) # sort the table before the slicing
        DF1 = DF[1:dim(DF)[1] / 2 +1,]
        DF2 = DF[(dim(DF)[1]  / 2 ):dim(DF)[1],]
        # Divide the datframe into two
        # plot according to the treatment
        plot1 = plot_point_16S(DF1, size,title, notaglom)
        plot2 = plot_point_16S(DF2, size,title, notaglom)
        plot_merge = ggpubr::ggarrange(plot1, plot2, ncol = 2, common.legend = TRUE, legend = "right", widths = c(0.4, 0.4))
        print(plot_merge)
        ggsave(paste0("output/", title, ".png"), plot = plot_merge)
    }
    else{
        plot = plot_point_16S(DF, size, title, notaglom)
        print(plot)
        ggsave(paste0("output/", title, ".png"), plot = plot, width = 10, height = 10)
    }
    return(plot)
}


# ============================================================================================================================
# script
# ============================================================================================================================


# model = binom(physeq4week, "Irrigation")
# saveRDS(model, "output/rhizoplane_week4_irrigation.rds")

model4week = readRDS('ressources/RDS/model_4week_irrigation.rds')
model4week_glom = readRDS('ressources/RDS/model_4week_irrigation_glom.rds')
model4_5_water = readRDS('ressources/RDS/model_4_5week_watered.rds')
model4_5_water_glom = readRDS('ressources/RDS/model_4_5week_watered_glom.rds')

plot1 = plot_res_16S(model4week, c("Order","Family","Genus"), "16S_4week", 12)
# plot_res_16S(model4week, c("Family", "Genus"), "16S_4week", 10) for the two levels Family and Genus
plot2 = plot_res_16S(model4week_glom, c("Family","Genus"), "16S_4week_glom", 12)

plot3 = plot_res_16S(model4_5_water, "Family", "16S_4_5week_watered", 12) # Week 5 effect
plot4 = plot_res_16S(model4_5_water_glom, "Family", "16S_4_5week_watered_glom", 12, FALSE) # week 5 effect


# preview for the combined plot
ggarrange(plot1, plot3, ncol = 2, labels = c("A", "B"))
ggarrange(plot2, plot4, ncol = 2, labels = c("A", "B"))


