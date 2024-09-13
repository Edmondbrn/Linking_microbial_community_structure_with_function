
library(ggplot2)
library(phyloseq)
library(ggpubr)
library(dplyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # working on Rstudio only
getwd()
### Loading data

physeq_rare_median = readRDS("ressources/RDS/Rarefied_physeq_median.rds")
physeq_rare_mean_rounded = readRDS("ressources/RDS/Rarefied_physeq_mean_rounded.rds")



### Alpha diversity


alpha_box_plot = function(data, x_var,y_var, col, comparizon, nb_test = 1, name){
  method = ifelse(shapiro.test(data[[y_var]])$p > 0.05, "t.test", "wilcox.test") # set the test for the comparizon
  symnum = list(cutpoints = c(0, 0.0001/nb_test, 0.001/nb_test, 0.01/nb_test, 0.05/nb_test, Inf), symbols = c("****", "***", "**", "*", "")) # format the output according to the correction
  p = ggplot(data = data, aes(x = !!sym(x_var), y = !!sym(y_var)))+
    scale_fill_manual(values=c("#018571","#a6611a"))+
    geom_boxplot(aes(fill = !!sym(x_var)))+
    geom_point(aes(x = !!sym(x_var), y = !!sym(y_var), fill= !!sym(x_var)),size = 2,alpha = 0.8, shape = 21, position = position_jitterdodge())+
    theme_bw()+
    # facet_grid(~ Week, scales = "free") + # organise the plot by week
    theme(axis.text.x = element_text(hjust = 1, size = 20, angle = 0),  
          axis.text.y = element_text(size = 20),  
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.title.x = element_text(size = 24),  
          axis.title.y = element_text(size = 24),
          strip.text = element_text(size = 15)) +  
    ylim(c(min(data[[y_var]]) - 0.1 * min(data[[y_var]]) , max(data[[y_var]]) + 0.1 *max(data[[y_var]]))) +
    labs(x = "", y = "")+
    stat_compare_means(method = method ,symnum.args = symnum, label = "p.signif", size = 12, label.x = 1.5, label.y = max(data[[y_var]]) )
  ggsave(paste0("output/alpha_diversity_",name ,"_", y_var, ".png"), plot = p, width = 25, height = 15, units = "cm")
  return(p)
}

physeq_rare_mean_rounded@sam_data$Week = ifelse(physeq_rare_mean_rounded@sam_data$Week== 2, "week 2",
                                          ifelse(physeq_rare_mean_rounded@sam_data$Week == 4, "week 4", "week 5"))

alpha_div = estimate_richness(physeq_rare_mean_rounded, measures = c("Shannon", "Chao1")) # compute alpha diversity
alpha_div = cbind(alpha_div, physeq_rare_mean_rounded@sam_data) # add the factor to the new table

comparison = list(c("Watered", "Drought")) # change this according to the variable you want to compare
nb_test = 2  # for Bonferonni correction


alpha_div_2week = alpha_div %>% filter(Week == "week 2")
alpha_div_4week = alpha_div %>% filter(Week == "week 4")
alpha_div_5week = alpha_div %>% filter(Week == "week 5")


p1 = alpha_box_plot(data = alpha_div_2week, x_var = "Irrigation", 
                    y_var = "Shannon", 
                    col = "Irrigation", 
                    comparizon = comparison , 
                    nb_test = nb_test, 
                    name = "control")

p2 = alpha_box_plot(data = alpha_div_4week, x_var = "Irrigation", 
                    y_var = "Shannon", 
                    col = "Irrigation", 
                    comparizon = comparison , 
                    nb_test = nb_test, 
                    name = "control")

p3 = alpha_box_plot(data = alpha_div_5week, x_var = "Irrigation", 
                    y_var = "Shannon", 
                    col = "Irrigation", 
                    comparizon = comparison , 
                    nb_test = nb_test, 
                    name = "control")

x11(height = 10, width = 10)
(pfinal = ggarrange(p2, p3, ncol = 1, nrow = 2, 
                   labels = c("Week 4", "Week 5"), label.x = -0.025,
                   common.legend = TRUE, legend = "none"))


# Normality test
ggqqplot(alpha_div_5week$Shannon)
ggqqplot(alpha_div_4week$Shannon[alpha_div_4week$Irrigation == "Watered"])
ggqqplot((alpha_div_4week$Shannon[alpha_div_4week$Irrigation == "Drought"]))
# looks okay, confirmed by shapiro test (not normal with week 2 and 5 I think but not sure)

shapiro.test(alpha_div_5week$Shannon)
shapiro.test(alpha_div_4week$Shannon[alpha_div_4week$Irrigation == "Watered"])
shapiro.test(alpha_div_4week$Shannon[alpha_div_4week$Irrigation == "Drought"])


t.test(alpha_div_4week$Shannon[alpha_div_4week$Irrigation == "Watered"],
       alpha_div_4week$Shannon[alpha_div_4week$Irrigation == "Drought"])

# ======================= Beta diversity ====================================
plot_pcoa = function(data, color, cluster, Shape, filename){
  list_output = list() # to return the plot and the beta diversity results
  print("Computing PCoA")
  physeq.ord_pcoa = ordinate(data, "PCoA", "bray")
  print("Plotting PCoA")
  
  p3 = plot_ordination(data, physeq.ord_pcoa, shape = Shape) +
    scale_color_manual(values=c("#018571","#a6611a")) +
    geom_point(size = 4, aes(color = !!sym(color))) +  # Use aes with sym
    theme_bw() +
    theme(axis.text = element_text(size = 15), axis.title = element_text(size = 25), 
          legend.text = element_text(size = 18), legend.title = element_text(size = 20))
  
  # add the elipse on the plot
  res = readline("Would you like to plot the elipse ? (Y/N)")
  if (res %in% "Yy"){
    p3 = p3 + stat_ellipse(aes(color = !!sym(cluster), fill = !!sym(cluster)), geom = "polygon", show.legend = T, linetype = "dashed", linewidth = 1, alpha = 0.1) + 
      labs(color = "Groups", fill = "Groups") +
      scale_color_manual(values = c("Week 2 control" = "#F8766D", # change the color of the groups points
                                    "Week 4 drought" = "#B79F00", 
                                    "Week 4 control" = "#00BA38", 
                                    "Week 5 control" = "#F564E3", 
                                    "Week 5 drought" = "#619CFF")) +
      scale_fill_manual(values = c("Week 2 control" = "#F8766D", # change the color of the groups ellipse
                                   "Week 4 drought" = "#B79F00", 
                                   "Week 4 control" = "#00BA38", 
                                   "Week 5 control" = "#F564E3", 
                                   "Week 5 drought" = "#619CFF")) 
  }
  
  ggsave(filename)
  list_output[[1]] = p3
  list_output[[2]] = physeq.ord_pcoa
  return(list_output)
}

p4 = plot_pcoa(physeq_rare_mean_rounded, "Irrigation", "Group", "Week", "output/PCoA_16S_rRNA_cluster.png")
p4[[1]]

physeq_rare_mean_rounded@sam_data$Group = ifelse(physeq_rare_mean_rounded@sam_data$Group == "2Watered" , "Week 2 control",
                                                 ifelse(physeq_rare_mean_rounded@sam_data$Group == "4Drought", "Week 4 drought", 
                                                        ifelse(physeq_rare_mean_rounded@sam_data$Group == "4Watered", "Week 4 control",
                                                               ifelse(physeq_rare_mean_rounded@sam_data$Group == "5Watered", "Week 5 control","Week 5 drought"))))

physeq_rare_mean_rounded@sam_data$Week = as.factor(physeq_rare_mean_rounded@sam_data$Week)
saveRDS(p4[[2]], file = "output/physeq_ord_pcoa.rds") # save beta diversity results


# for week 4 and 5 samples
x11(height = 10, width = 20)
plot_comb = ggarrange(p4[[1]], p2,  ncol = 2, 
          # labels , label.x = 0.5, label.y = 0,
          widths = c(2,1), common.legend = TRUE, legend = "bottom")

# Put a title for each plots
annotate_figure(plot_comb,
                top = text_grob("\t\t\t\t\t\t\t\t Beta diversity (PCoA) \t\t\t\t\t\t\t\t\t\t\t\t\t   Alpha diversity (Shannon)", 
                size = 20))

### PERMANOVA
permanova = function(physeq, var1, var2 = NULL, var3 = NULL, var4 = NULL, var5 = NULL, matrix = "bray"){
  cat("Preparing data for PERMANOVA\n")
  data = as.data.frame((t(physeq@otu_table))) # convert the otu_table to a dataframe to perform the permANOVA
  
  metdata = physeq@sam_data
  data = cbind(data, physeq@sam_data)
  data$Irrigation = as.factor(data$Irrigation)
  data$Week = as.factor(data$Week)
  data$Compartments = as.factor(data$Compartments)
  data$Treatment = as.factor(data$Treatment)
  data$Group= as.factor(data$Group)
  # remove the 6 last columns (metadata variable)
  data_corr = data[, -((ncol(data)-6):ncol(data))]
  cat("Computing Bray-Curtis dissimilarity matrix\n")
  bray_dist = vegdist(data_corr, method = matrix)
  cat("Computing PERMANOVA\n")
  formula_str = paste("bray_dist ~", var1, if (!is.null(var2)) paste("*", var2) else "", 
                      if (!is.null(var3)) paste("*", var3) else "", 
                      if (!is.null(var4)) paste("*", var4) else "", 
                      if (!is.null(var5)) paste("*", var5) else "")
  formula = as.formula(formula_str)
  
  model = adonis2(formula, data = data, permutations = 999, by = "terms")
  return(model)
}
permanova(physeq = physeq_rare_mean_rounded, var1 = "Irrigation", var2 = "Week")



