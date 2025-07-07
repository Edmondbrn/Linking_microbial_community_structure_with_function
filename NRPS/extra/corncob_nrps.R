# ==========================================================================================================

# Load the necessary libraries

# ==========================================================================================================
.libPaths("/home/edmond/R/x86_64-pc-linux-gnu-library/4.3")
library(corncob)
library(phyloseq)
library(ggplot2)
library(dplyr)
library(ggpubr)

# ==========================================================================================================

# Load the data

# ==========================================================================================================
setwd("/home/edmond/NRPS_analysis")
physeq = readRDS("output/carbom.rds")

physeq_rhizoplane = subset_samples(physeq, Compartments == "rhizoplane")
physeq_rhizosphere = subset_samples(physeq, Compartments == "rhizosphere")

physeq_plane_water = subset_samples(physeq_rhizoplane, Irrigation == "Watered")
physeq_plane_water_control = subset_samples(physeq_plane_water, Treatment == "Control")
physeq_plane_WC_2_4 = subset_samples(physeq_plane_water_control, WEEK %in% c("2weeks", "4weeks"))
physeq_plane_WC_4_5 = subset_samples(physeq_plane_water_control, WEEK %in% c("4weeks", "5weeks"))
physeq_plane_WC_2_5 = subset_samples(physeq_plane_water_control, WEEK %in% c("2weeks", "5weeks"))

physeq_sphere_water = subset_samples(physeq_rhizosphere, Irrigation == "Watered")
physeq_sphere_water_control = subset_samples(physeq_sphere_water, Treatment == "Control")
physeq_sphere_WC_2_4 = subset_samples(physeq_sphere_water_control, WEEK %in% c("2weeks", "4weeks"))
physeq_sphere_WC_4_5 = subset_samples(physeq_sphere_water_control, WEEK %in% c("4weeks", "5weeks"))
physeq_sphere_WC_2_5 = subset_samples(physeq_sphere_water_control, WEEK %in% c("2weeks", "5weeks"))

physeq_plane_drought = subset_samples(physeq_rhizoplane, Irrigation == "Drought")
physeq_plane_drought_control = subset_samples(physeq_plane_drought, Treatment == "Control")

physeq_sphere_drought = subset_samples(physeq_rhizosphere, Irrigation == "Drought")
physeq_sphere_drought_control = subset_samples(physeq_sphere_drought, Treatment == "Control")

physeq_sphere_week4 = subset_samples(physeq_rhizosphere, WEEK == "4weeks")
physeq_sphere_week4_ctrl = subset_samples(physeq_sphere_week4, Treatment == "Control")
physeq_sphere_week4_water = subset_samples(physeq_sphere_week4, Irrigation == "Watered")
physeq_sphere_wee4_drought = subset_samples(physeq_sphere_week4, Irrigation == "Drought")

physeq_sphere_W4_SB_Ctrl_water = subset_samples(physeq_sphere_week4_water, Treatment %in% c("Control", "SBW25"))
physeq_sphere_W4_SB_Ctrl_drought = subset_samples(physeq_sphere_wee4_drought, Treatment %in% c("Control", "SBW25"))
physeq_sphere_W4_SB_Mut_water = subset_samples(physeq_sphere_week4_water, Treatment %in% c("Mutant", "SBW25"))
physeq_sphere_W4_SB_Mut_drought = subset_samples(physeq_sphere_wee4_drought, Treatment %in% c("Mutant", "SBW25"))
physeq_sphere_W4_Mut_Ctrl_water = subset_samples(physeq_sphere_week4_water, Treatment %in% c("Control", "Mutant"))
physeq_sphere_W4_Mut_Ctrl_drought = subset_samples(physeq_sphere_wee4_drought, Treatment %in% c("Control", "Mutant"))


physeq_sphere_week5 = subset_samples(physeq_rhizosphere, WEEK == "5weeks")
physeq_sphere_week5_ctrl = subset_samples(physeq_sphere_week5, Treatment == "Control")
physeq_sphere_week5_water = subset_samples(physeq_sphere_week5, Irrigation == "Watered")
physeq_sphere_week5_drought = subset_samples(physeq_sphere_week5, Irrigation == "Drought")

physeq_sphere_W5_SB_Ctrl_water = subset_samples(physeq_sphere_week5_water, Treatment %in% c("Control", "SBW25"))
physeq_sphere_W5_SB_Ctrl_drought = subset_samples(physeq_sphere_week5_drought, Treatment %in% c("Control", "SBW25"))
physeq_sphere_W5_SB_Mut_water = subset_samples(physeq_sphere_week5_water, Treatment %in% c("Mutant", "SBW25"))
physeq_sphere_W5_SB_Mut_drought = subset_samples(physeq_sphere_week5_drought, Treatment %in% c("Mutant", "SBW25"))
physeq_sphere_W5_Mut_Ctrl_water = subset_samples(physeq_sphere_week5_water, Treatment %in% c("Control", "Mutant"))
physeq_sphere_W5_Mut_Ctrl_drought = subset_samples(physeq_sphere_week5_drought, Treatment %in% c("Control", "Mutant"))


physeq_plane_week4 = subset_samples(physeq_rhizoplane, WEEK == "4weeks")
physeq_plane_week4_control = subset_samples(physeq_plane_week4, Treatment == "Control")
physeq_plane_week4_water = subset_samples(physeq_plane_week4, Irrigation == "Watered")
physeq_plane_week4_drought = subset_samples(physeq_plane_week4, Irrigation == "Drought")

physeq_plane_W4_SB_Ctrl_water = subset_samples(physeq_plane_week4_water, Treatment %in% c("Control", "SBW25"))
physeq_plane_W4_SB_Ctrl_drought = subset_samples(physeq_plane_week4_drought, Treatment %in% c("Control", "SBW25"))
physeq_plane_W4_SB_Mut_water = subset_samples(physeq_plane_week4_water, Treatment %in% c("Mutant", "SBW25"))
physeq_plane_W4_SB_Mut_drought = subset_samples(physeq_plane_week4_drought, Treatment %in% c("Mutant", "SBW25"))
physeq_plane_W4_Mut_Ctrl_water = subset_samples(physeq_plane_week4_water, Treatment %in% c("Control", "Mutant"))
physeq_plane_W4_Mut_Ctrl_drought = subset_samples(physeq_plane_week4_drought, Treatment %in% c("Control", "Mutant"))



physeq_plane_week5 = subset_samples(physeq_rhizoplane, WEEK == "5weeks")
physeq_plane_week5_control = subset_samples(physeq_plane_week5, Treatment == "Control")
physeq_plane_week5_water = subset_samples(physeq_plane_week5, Irrigation == "Watered")
physeq_plane_week5_drought = subset_samples(physeq_plane_week5, Irrigation == "Drought")

physeq_plane_W5_SB_Ctrl_water = subset_samples(physeq_plane_week5_water, Treatment %in% c("Control", "SBW25"))
physeq_plane_W5_SB_Ctrl_drought = subset_samples(physeq_plane_week5_drought, Treatment %in% c("Control", "SBW25"))
physeq_plane_W5_SB_Mut_water = subset_samples(physeq_plane_week5_water, Treatment %in% c("Mutant", "SBW25"))
physeq_plane_W5_SB_Mut_drought = subset_samples(physeq_plane_week5_drought, Treatment %in% c("Mutant", "SBW25"))
physeq_plane_W5_Mut_Ctrl_water = subset_samples(physeq_plane_week5_water, Treatment %in% c("Control", "Mutant"))
physeq_plane_W5_Mut_Ctrl_drought = subset_samples(physeq_plane_week5_drought, Treatment %in% c("Control", "Mutant"))


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

binom_variation = function(data, variable){
    # create the formula
    formula_str = paste0("~ ", variable)
    formula = as.formula(formula_str)
    phi_formula = formula
    # perform the analysis
    analysis = differentialTest(formula = formula, phi.formula = phi_formula, data = data,
                                formula_null = ~ formula, phi.formula_null = 1, test = "LRT", boot = FALSE,
                                fdr = "fdr", fdr_cutoff = 0.05, verbose = T)    
    return(analysis)
}

plot_point = function(data, size){ # function to plot the graph
    plot = ggplot(data, aes(x = x, y = taxa, color = color)) +
                geom_point(color = "black") + 
                geom_errorbar(aes(xmin = xmin, xmax = xmax), width = 0.2) + 
                geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
                theme_linedraw() +
                facet_grid(. ~ title)+
                labs(x = "", y = "Genus") +
                theme(axis.text.y = element_text(angle = 0, hjust = 1, size = size))
                theme(axis.text.x = element_text(angle = 0, hjust = 1, size = size))+
    return(plot)
}

plot_res_unique = function(model, title, size){
    if (length(model$significant_taxa) ==  0){
        return("Error, no significant taxa found")
    }
    # extract data

    df = plot(model, total = T, B = 1000, level = c("genus", "family"), dataonly = T)
    DF = df$data
    DF$color = sign(DF$xmax) == sign(DF$xmin) 
    DF$title = title
    
    if (length(grep("_", DF$taxa)) == 1){ # test if there is no ASV number in the data, it happens some time
        for (i in 1:length(DF$taxa)){
            DF$taxa[[i]] =  gsub(".*", paste0(DF$taxa[[i]] ," (", model$significant_taxa[[i]],")"), DF$taxa[[i]])
        }
    } 
    if (dim(DF)[1]>80){
        DF2 = DF[1:dim(DF)[1] / 2,]
        DF1 = DF[(dim(DF)[1] / 2 + 1):dim(DF)[1],]
        # Divide the datframe into two
        # plot according to the treatment
        plot1 = plot_point(DF1, size)
        plot2 = plot_point(DF2, size)
        print(ggpubr::ggarrange(plot1, plot2, ncol = 2, common.legend = TRUE, legend = "right"))
    }
    else{
        plot = plot_point(DF, size)
        print(plot)
    }
    return(DF)

}

plot_proportion = function(data){
    propor = as.data.frame(prop.table(table(data))*100)
    ggplot(propor, aes(x = data, y = Freq, fill = data)) +
        geom_bar(stat = "identity") +
        theme_linedraw() +
        labs(x = "Taxa", y = "Proportion (%)", fill = "Taxa colour", size = 15)+
        # scale_fill_manual("rainbow") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12))+
        theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12))
}

# Main function to plot the data by comparing another ASVs
plot_res_unique_common = function(model, title, size, common_table){
    if (length(model$significant_taxa) ==  0){
        return("Error, no significant taxa found")
    }
    # extract data
    df = plot(model, total = T, B = 1000, level = c("genus", "family"), dataonly = T)
    DF = df$data
    DF$color = sign(DF$xmax) == sign(DF$xmin) 
    DF$title = title
    
    DF$common = ifelse(DF$taxa %in% common_table$taxa, "common", "new")
    

    if ( length(grep("_", DF$taxa)) == 0){ # test if there is no ASV number in the data, happens some time
        for (i in 1:length(DF$taxa)){
            DF$taxa[[i]] =  gsub(".*", paste0(DF$taxa[[i]] ," (", model$significant_taxa[[i]],")"), DF$taxa[[i]])
        }
    } 
    # Divide the datframe into two
    DF2 = DF[1:dim(DF)[1] / 2,]
    DF1 = DF[(dim(DF)[1] / 2 + 1):dim(DF)[1],]
    # plot according to the treatment
    plot1 = plot_point(DF1, size)
    plot2 = plot_point(DF2, size)
    print(ggpubr::ggarrange(plot1, plot2, ncol = 2, common.legend = TRUE, legend = "right"))
    return(DF)

}

# plot to fill the data according to another one to see common ASVs
plot_point_common = function(data, size){ # function to plot the graph
    plot = ggplot(data, aes(x = x, y = taxa, color = common)) +
                geom_point(color = "black") + 
                geom_errorbar(aes(xmin = xmin, xmax = xmax), width = 0.2) + 
                geom_vline(xintercept = 0, linetype = "dashed", color = "black") + 
                theme_linedraw() +
                facet_grid(. ~ title)+
                labs(x = "", y = "Taxa") +
                theme(axis.text.y = element_text(angle = 0, hjust = 1, size = size))
                theme(axis.text.x = element_text(angle = 0, hjust = 1, size = size))+
    return(plot)
}

# ==========================================================================================================

# Analysis
# Rhizoplane
# model1 = binom(physeq_plane_WC_2_4, "WEEK")
# saveRDS(model1, "RDS/corncob/plane_WC_week2_4.rds")
# model2 = binom(physeq_plane_WC_2_5, "WEEK")
# saveRDS(model2, "RDS/corncob/plane_WC_week2_5.rds")
# model3 = binom(physeq_plane_WC_4_5, "WEEK")
# saveRDS(model3, "RDS/corncob/plane_WC_week4_5.rds")
# # Rhizosphere
# model4 = binom(physeq_sphere_WC_2_4, "WEEK")
# saveRDS(model4, "RDS/corncob/sphere_WC_week2_4.rds")
# model5 = binom(physeq_sphere_WC_2_5, "WEEK")
# saveRDS(model5, "RDS/corncob/sphere_WC_week2_5.rds")
# model6 = binom(physeq_sphere_WC_4_5, "WEEK")
# saveRDS(model6, "RDS/corncob/sphere_WC_week4_5.rds")

# # Rhizoplane
# model7 = binom(physeq_plane_drought_control, "WEEK")
# saveRDS(model7, "RDS/corncob/plane_DC_week4_5.rds")

# # Rhizosphere
# model8 = binom(physeq_sphere_drought_control, "WEEK")
# saveRDS(model8, "RDS/corncob/sphere_DC_week2_4.rds")

# model9 = binom(physeq_sphere_week4, "Irrigation")
# saveRDS(model9, "RDS/corncob/rhizosphere_week4.rds")
# model10 = binom(physeq_sphere_week5, "Irrigation")
# saveRDS(model10, "RDS/corncob/rhizosphere_week5.rds")
# model11 = binom(physeq_plane_week4, "Irrigation")
# saveRDS(model11, "RDS/corncob/rhizoplane_week4.rds")
# model12 = binom(physeq_plane_week5, "Irrigation")
# saveRDS(model12, "RDS/corncob/rhizoplane_week5.rds")

# model13 = binom(physeq_sphere_week4_ctrl, "Irrigation")
# saveRDS(model13, "RDS/corncob/SBW/rhizosphere_week4_irrigation.rds")
# model14 = binom(physeq_sphere_week4_ctrl, "Irrigation")
# saveRDS(model14, "RDS/corncob/SBW/rhizosphere_week5_irrigation.rds")
# model15 = binom(physeq_plane_week4_control, "Irrigation")
# saveRDS(model15, "RDS/corncob/SBW/rhizoplane_week4_irrigation.rds")
# model16 = binom(physeq_plane_week5_control, "Irrigation")
# saveRDS(model16, "RDS/corncob/SBW/rhizoplane_week5_irrigation.rds")

# model17 = binom(physeq_sphere_W4_SB_Ctrl_water, "Treatment")
# saveRDS(model17, "RDS/corncob/SBW/rhizosphere_week4_SB_Ctrl_water.rds")
# model18 = binom(physeq_sphere_W4_SB_Ctrl_drought, "Treatment")
# saveRDS(model18, "RDS/corncob/SBW/rhizosphere_week4_SB_Ctrl_drought.rds")
# model19 = binom(physeq_sphere_W4_SB_Mut_water, "Treatment")
# saveRDS(model19, "RDS/corncob/SBW/rhizosphere_week4_SB_Mut_water.rds")
# model20 = binom(physeq_sphere_W4_SB_Mut_drought, "Treatment")
# saveRDS(model20, "RDS/corncob/SBW/rhizosphere_week4_SB_Mut_drought.rds")
# model21 = binom(physeq_sphere_W4_Mut_Ctrl_water, "Treatment")
# saveRDS(model21, "RDS/corncob/SBW/rhizosphere_week4_Mut_Ctrl_water.rds")
# model22 = binom(physeq_sphere_W4_Mut_Ctrl_drought, "Treatment")
# saveRDS(model22, "RDS/corncob/SBW/rhizosphere_week4_Mut_Ctrl_drought.rds")

# model23 = binom(physeq_sphere_W5_SB_Ctrl_water, "Treatment")
# saveRDS(model23, "RDS/corncob/SBW/rhizosphere_week5_SB_Ctrl_water.rds")
# model24 = binom(physeq_sphere_W5_SB_Ctrl_drought, "Treatment")
# saveRDS(model24, "RDS/corncob/SBW/rhizosphere_week5_SB_Ctrl_drought.rds")
# model25 = binom(physeq_sphere_W5_SB_Mut_water, "Treatment")
# saveRDS(model25, "RDS/corncob/SBW/rhizosphere_week5_SB_Mut_water.rds")
# model26 = binom(physeq_sphere_W5_SB_Mut_drought, "Treatment")
# saveRDS(model26, "RDS/corncob/SBW/rhizosphere_week5_SB_Mut_drought.rds")
# model27 = binom(physeq_sphere_W5_Mut_Ctrl_water, "Treatment")
# saveRDS(model27, "RDS/corncob/SBW/rhizosphere_week5_Mut_Ctrl_water.rds")
# model28 = binom(physeq_sphere_W5_Mut_Ctrl_drought, "Treatment")
# saveRDS(model28, "RDS/corncob/SBW/rhizosphere_week5_Mut_Ctrl_drought.rds")

# model29 = binom(physeq_plane_W4_SB_Ctrl_water, "Treatment")
# saveRDS(model29, "RDS/corncob/SBW/rhizoplane_week4_SB_Ctrl_water.rds")
# model30 = binom(physeq_plane_W4_SB_Ctrl_drought, "Treatment")
# saveRDS(model30, "RDS/corncob/SBW/rhizoplane_week4_SB_Ctrl_drought.rds")
# model31 = binom(physeq_plane_W4_SB_Mut_water, "Treatment")
# saveRDS(model31, "RDS/corncob/SBW/rhizoplane_week4_SB_Mut_water.rds")
# model32 = binom(physeq_plane_W4_SB_Mut_drought, "Treatment")
# saveRDS(model32, "RDS/corncob/SBW/rhizoplane_week4_SB_Mut_drought.rds")
# model33 = binom(physeq_plane_W4_Mut_Ctrl_water, "Treatment")
# saveRDS(model33, "RDS/corncob/SBW/rhizoplane_week4_Mut_Ctrl_water.rds")
# model34 = binom(physeq_plane_W4_Mut_Ctrl_drought, "Treatment")  
# saveRDS(model34, "RDS/corncob/SBW/rhizoplane_week4_Mut_Ctrl_drought.rds")

# model35 = binom(physeq_plane_W5_SB_Ctrl_water, "Treatment")
# saveRDS(model35, "RDS/corncob/SBW/rhizoplane_week5_SB_Ctrl_water.rds")
# model36 = binom(physeq_plane_W5_SB_Ctrl_drought, "Treatment")
# saveRDS(model36, "RDS/corncob/SBW/rhizoplane_week5_SB_Ctrl_drought.rds")
# model37 = binom(physeq_plane_W5_SB_Mut_water, "Treatment")
# saveRDS(model37, "RDS/corncob/SBW/rhizoplane_week5_SB_Mut_water.rds")
# model38 = binom(physeq_plane_W5_SB_Mut_drought, "Treatment")
# saveRDS(model38, "RDS/corncob/SBW/rhizoplane_week5_SB_Mut_drought.rds")
# model39 = binom(physeq_plane_W5_Mut_Ctrl_water, "Treatment")
# saveRDS(model39, "RDS/corncob/SBW/rhizoplane_week5_Mut_Ctrl_water.rds")
# model40 = binom(physeq_plane_W5_Mut_Ctrl_drought, "Treatment")
# saveRDS(model40, "RDS/corncob/SBW/rhizoplane_week5_Mut_Ctrl_drought.rds")

physeq_plane_2week = subset_samples(physeq_rhizoplane, WEEK == "2weeks")
physeq_plane_2week_SBW_ctrl = subset_samples(physeq_plane_2week, Treatment %in% c("Control", "SBW25"))
physeq_plane_2week_SBW_mut = subset_samples(physeq_plane_2week, Treatment %in% c("Mutant", "SBW25"))
physeq_plane_2week_mut_ctrl = subset_samples(physeq_plane_2week, Treatment %in% c("Control", "Mutant"))


physeq_sphere_2_weeks = subset_samples(physeq_rhizosphere, WEEK == "2weeks")
physeq_sphere_2weeks_SBW_ctrl = subset_samples(physeq_sphere_2_weeks, Treatment %in% c("Control", "SBW25"))
physeq_sphere_2weeks_SBW_mut = subset_samples(physeq_sphere_2_weeks, Treatment %in% c("Mutant", "SBW25"))
physeq_sphere_2weeks_mut_ctrl = subset_samples(physeq_sphere_2_weeks, Treatment %in% c("Control", "Mutant"))

model41 = binom(physeq_plane_2week_SBW_ctrl, "Treatment")
saveRDS(model41, "RDS/corncob/rhizoplane_week2_SBW_ctrl.rds")
model42 = binom(physeq_plane_2week_SBW_mut, "Treatment")
saveRDS(model42, "RDS/corncob/rhizoplane_week2_SBW_mut.rds")
model43 = binom(physeq_plane_2week_mut_ctrl, "Treatment")
saveRDS(model43, "RDS/corncob/rhizoplane_week2_mut_ctrl.rds")

model44 = binom(physeq_sphere_2weeks_SBW_ctrl, "Treatment")
saveRDS(model44, "RDS/corncob/rhizosphere_week2_SBW_ctrl.rds")
model45 = binom(physeq_sphere_2weeks_SBW_mut, "Treatment")
saveRDS(model45, "RDS/corncob/rhizosphere_week2_SBW_mut.rds")
model46 = binom(physeq_sphere_2weeks_mut_ctrl, "Treatment")
saveRDS(model46, "RDS/corncob/rhizosphere_week2_mut_ctrl.rds")

# ==========================================================================================================

print("Touts'est bien passé")
q()

# physeq_plane_week5

# model1 = readRDS("RDS/corncob/plane_WC_week2_4.rds")
# model2 = readRDS("RDS/corncob/plane_WC_week2_5.rds")
# model3 = readRDS("RDS/corncob/plane_WC_week4_5.rds")
# model4 = readRDS("RDS/corncob/sphere_WC_week2_4.rds")
# model5 = readRDS("RDS/corncob/sphere_WC_week2_5.rds")
# model6 = readRDS("RDS/corncob/sphere_WC_week4_5.rds")
# model7 = readRDS("RDS/corncob/plane_DC_week4_5.rds")
# model8 = readRDS("RDS/corncob/sphere_DC_week2_4.rds")
model9 = readRDS("RDS/corncob/rhizosphere_week4.rds")
model10 = readRDS("RDS/corncob/rhizosphere_week5.rds")
model11 = readRDS("RDS/corncob/rhizoplane_week4.rds")
model12 = readRDS("RDS/corncob/rhizoplane_week5.rds")

model13 = readRDS("NRPS_analysis/RDS/corncob/SBW/rhizoplane_week4_Mut_Ctrl_water.rds")
model14 = readRDS("NRPS_analysis/RDS/corncob/SBW/rhizoplane_week4_SB_Ctrl_water.rds")
model15 = readRDS("NRPS_analysis/RDS/corncob/SBW/rhizoplane_week4_SB_Mut_water.rds")

model16 = readRDS("NRPS_analysis/RDS/corncob/SBW/rhizoplane_week5_Mut_Ctrl_water.rds")
model17 = readRDS("NRPS_analysis/RDS/corncob/SBW/rhizoplane_week5_SB_Ctrl_water.rds")
model18 = readRDS("NRPS_analysis/RDS/corncob/SBW/rhizoplane_week5_SB_Mut_water.rds")


model19 = readRDS("RDS/corncob/SBW/rhizoplane_week4_Mut_Ctrl_drought.rds")
model20 = readRDS("RDS/corncob/SBW/rhizoplane_week4_SB_Ctrl_drought.rds")
model21 = readRDS("RDS/corncob/SBW/rhizoplane_week4_SB_Mut_drought.rds")

model22 = readRDS("RDS/corncob/SBW/rhizoplane_week5_Mut_Ctrl_drought.rds")
model23 = readRDS("RDS/corncob/SBW/rhizoplane_week5_SB_Ctrl_drought.rds")
model24 = readRDS("RDS/corncob/SBW/rhizoplane_week5_SB_Mut_drought.rds")

model25 = readRDS("RDS/corncob/SBW/rhizosphere_week4_Mut_Ctrl_drought.rds")
model26 = readRDS("RDS/corncob/SBW/rhizosphere_week4_SB_Ctrl_drought.rds")
model27 = readRDS("RDS/corncob/SBW/rhizosphere_week4_SB_Mut_drought.rds")

model28 = readRDS("RDS/corncob/SBW/rhizosphere_week5_Mut_Ctrl_drought.rds")
model29 = readRDS("RDS/corncob/SBW/rhizosphere_week5_SB_Ctrl_drought.rds")
model30 = readRDS("RDS/corncob/SBW/rhizosphere_week5_SB_Mut_drought.rds")

model31 = readRDS("RDS/corncob/SBW/rhizosphere_week4_Mut_Ctrl_water.rds")
model32 = readRDS("RDS/corncob/SBW/rhizosphere_week4_SB_Ctrl_water.rds")
model33 = readRDS("RDS/corncob/SBW/rhizosphere_week4_SB_Mut_water.rds")

model34 = readRDS("RDS/corncob/SBW/rhizosphere_week5_Mut_Ctrl_water.rds")
model35 = readRDS("RDS/corncob/SBW/rhizosphere_week5_SB_Ctrl_water.rds")
model36 = readRDS("RDS/corncob/SBW/rhizosphere_week5_SB_Mut_water.rds")

# plot_res_unique(model1, "Rhizoplane watered week 2 vs 4", 10)
# plot_res_unique(model2, "Rhizoplane watered week 2 vs 5", 10)
# plot_res_unique(model3, "Rhizoplane watered week 4 vs 5", 10)
# plot_res_unique(model4, "Rhizosphere watered week 2 vs 4", 10)
# plot_res_unique(model5, "Rhizosphere watered week 2 vs 5", 10)
# plot_res_unique(model6, "Rhizosphere watered week 4 vs 5", 10)
# plot_res_unique(model7, "Rhizoplane drought week 4 vs 5", 10)
# plot_res_unique(model8, "Rhizosphere drought week 4 vs 5", 10)
plot_res_unique(model9, "Rhizosphere week 4 irrigation effect", 10)
plot_res_unique(model10, "Rhizosphere week 5 irrigation effect", 10)
plot_res_unique(model11, "Rhizoplane week 4 irrigation effect", 10)
plot_res_unique(model12, "Rhizoplane week 5 irrigation effect", 10)

plot_res_unique(model13, "Rhizoplane week 4 Mutant vs Ctrl", 10)
plot_res_unique(model14, "Rhizoplane week 4 SBW25 vs Ctrl", 10)
plot_res_unique(model15, "Rhizoplane week 4 SBW25 vs Mutant", 10)


plot_res_unique(model16, "Rhizoplane week 5 Mutant vs Ctrl", 10)
plot_res_unique(model17, "Rhizoplane week 5 SBW25 vs Ctrl", 10)
plot_res_unique(model18, "Rhizoplane week 5 SBW25 vs Mutant", 10)

plot_res_unique(model19, "Rhizoplane week 4 Mutant vs Ctrl", 10)
plot_res_unique(model20, "Rhizoplane week 4 SBW25 vs Ctrl", 10)
plot_res_unique(model21, "Rhizoplane week 4 SBW25 vs Mutant", 10)

plot_res_unique(model22, "Rhizoplane week 5 Mutant vs Ctrl", 10)
plot_res_unique(model23, "Rhizoplane week 5 SBW25 vs Ctrl", 10)
plot_res_unique(model24, "Rhizoplane week 5 SBW25 vs Mutant", 10)

plot_res_unique(model19, "Rhizoplane week 4 Mutant vs Ctrl", 10)
plot_res_unique(model20, "Rhizoplane week 4 SBW25 vs Ctrl", 10)
plot_res_unique(model21, "Rhizoplane week 4 SBW25 vs Mutant", 10)

plot_res_unique(model22, "Rhizoplane week 5 Mutant vs Ctrl", 10)
plot_res_unique(model23, "Rhizoplane week 5 SBW25 vs Ctrl", 10)
plot_res_unique(model24, "Rhizoplane week 5 SBW25 vs Mutant", 10)


plot_res_unique(model25, "Rhizosphere week 4 Mutant vs Ctrl", 10)
plot_res_unique(model26, "Rhizosphere week 4 SBW25 vs Ctrl", 10)
plot_res_unique(model27, "Rhizosphere week 4 SBW25 vs Mutant", 10)

plot_res_unique(model28, "Rhizosphere week 5 Mutant vs Ctrl", 10)
plot_res_unique(model29, "Rhizosphere week 5 SBW25 vs Ctrl", 10)
plot_res_unique(model30, "Rhizosphere week 5 SBW25 vs Mutant", 10)


plot_res_unique(model31, "Rhizosphere week 4 Mutant vs Ctrl", 10)
plot_res_unique(model32, "Rhizosphere week 4 SBW25 vs Ctrl", 10)
plot_res_unique(model33, "Rhizosphere week 4 SBW25 vs Mutant", 10)

plot_res_unique(model34, "Rhizosphere week 5 Mutant vs Ctrl", 10)
plot_res_unique(model35, "Rhizosphere week 5 SBW25 vs Ctrl", 10)
plot_res_unique(model36, "Rhizosphere week 5 SBW25 vs Mutant", 10)







