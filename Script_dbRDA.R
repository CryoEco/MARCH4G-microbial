setwd("/Users/...")

# load betadiv data
set.seed(123)

library(phyloseq)
library(vegan)
library(tidyverse)
library(microViz)
library(dplyr)


#import metadata
#Extract metadata
metadata_16S_2 <- readr::read_tsv("/.../metadata.tsv")
metadata_16S_2
str(metadata_16S_2)

otu_table(ps_hell)

#Log-transform concentrations
#EC, NO3, NH4, PO4, CH4, SPM, DOC, DN, SO4, Fe, POC, 16SCps_L, Avg_cps_gSPM, Avg_cps_gPOC

met_env <- decostand(metadata_16S_2[,c(17, 19:22, 25:29, 31:34)], method = "log", logbase = 10)
met_env

met_env_all <- cbind(met_env, metadata_16S_2[,c(1, 2, 6:8, 11:16, 18, 23,24)])
met_env_all

#Standardize all variables
met_env_std <- decostand(met_env_all[, c(1:14,21,22,24:26)], method = "standardize")
met_env_std

met_env_std_all <- cbind(met_env_std, met_env_all[, c(15:20,23, 27,28)])
met_env_std_all


#### model including NO3, PO4, CH4, DOC, SO4, pH, Fe, CatchSize, latitude + Cps_L, Avg_cps_gPOC + POC + SPM + Avg_cps_gSPM + Condition (DOY)
##with this model, samples to be removed (NA) are: KA-02, NA-05, NA-03, KA-04-b, KA-05-a, IL-02-b, QA-2-3, KA-01-d, IL-02-D

ps_model1 <- subset_samples(ps_hell, OrSiteName != "NA05" & ConsSiteName != "IL-02-b" & ConsSiteName != "NA-03" & ConsSiteName != "KA-02" & ConsSiteName != "KA-04-b" & ConsSiteName != "KA-05-a" & ConsSiteName != "KA-01-d" & OrSiteName != "QA-2-3" & OrSiteName != "IL-02-D")

#model with 47 samples
ps_model1
ps_model1 <- phyloseq_validate(ps_model1, remove_undetected = TRUE, verbose = TRUE)
ps_model1

#Create distance matrix from ps_model1
dist <- phyloseq::distance(ps_model1, method = "bray")
dist

env_dbrda1 <- met_env_std_all %>% select(-NH4, -EC, -DN, -Q)
env_dbrda1 <- env_dbrda1 %>% na.omit()


partialdbRDA1 <- dbrda(dist ~ NO3+CH4+DOC+SO4+pH+Fe+latitudeN+PO4+CatchmentSize+Cps_L+Avg_cps_gSPM+Avg_cps_gPOC + POC + SPM + Condition(DOY), env_dbrda1, distance = "bray") #Some constraints or conditions were aliased because they were redundant. This can happen if terms are linearly dependent (collinear):‘POC’, ‘SPM’

summary(partialdbRDA1) #the model with all variables explains 61,41% of the variation
vif.cca(partialdbRDA1) # NO3, Cps_L, pH and Avg_cps_gSPM with VIF > 5, POC and SPM aliased

#run model again, removing NO3, pH and Cps_L

env_dbrda2 <- env_dbrda1 %>% select(-NO3, -Cps_L, -pH)


partialdbRDA2 <- dbrda(dist ~ CH4+DOC+SO4+Fe+latitudeN+PO4+CatchmentSize+Avg_cps_gSPM+Avg_cps_gPOC + POC + SPM + Condition(DOY), env_dbrda2, distance = "bray") #Some constraints or conditions were aliased because they were redundant. This can happen if terms are linearly dependent (collinear):‘SPM’

summary(partialdbRDA2) #the model with all variables explains 58,46% of the variation
vif.cca(partialdbRDA2) #Avg_cps_gSPM, Avg_cps_gPOC and POC with VIF > 5

##run model again, removing Avg_cps_gSPM and POC

env_dbrda3 <- env_dbrda2 %>% select(-Avg_cps_gSPM, -POC)

partialdbRDA3 <- dbrda(dist ~ CH4+DOC+SO4+Fe+latitudeN+PO4+CatchmentSize+Avg_cps_gPOC + SPM + Condition(DOY), env_dbrda3, distance = "bray") 

summary(partialdbRDA3) #the model with all variables explains 55,79% of the variation
vif.cca(partialdbRDA3) #all variables with VIF < 5

#fwd selection
fwd_partialdbRDA <- ordiR2step(dbrda(dist ~ 1, data = env_dbrda3, distance = "bray"), 
                                 scope = formula(partialdbRDA3), 
                                 direction = "forward",
                                 R2scope = TRUE,
                                 pstep = 1000, 
                                 trace = TRUE)

fwd_partialdbRDA$call #PO4 + CH4 + SPM + latitudeN + CatchmentSize + DOC + Fe

vif.cca(fwd_partialdbRDA)

final_model <- dbrda(dist ~ PO4 + CH4 + SPM + latitudeN + DOC + CatchmentSize + Fe + Condition(DOY), data = env_dbrda3, distance = "bray")
summary(final_model) ##50,55%
capture.output(summary(final_model), file = "final_model_dbRDA.txt")
vif.cca(final_model)

R2_final_model <- RsquareAdj(final_model)$r.squared
R2_final_model

#adjusted R-squared (corrected for the number of explanatory variables)
R2adj_final_model <- RsquareAdj(final_model)$adj.r.squared
R2adj_final_model ## correcting for the number of variables, the % variation explained is now 43.67

#Significance testing
#whole model
anova(final_model, step = 1000, set.seed(123))

#significance of each variable

anova_final_model_margin <- anova(final_model, step = 1000, by = "margin", set.seed(123))
anova_final_model_margin
anova_final_model_margin$`Pr(>F)`<- p.adjust(anova_final_model_margin$`Pr(>F)`, method = 'fdr')
anova_final_model_margin 

#significance of each canonical axis
anova(final_model, step = 1000, by = "axis", set.seed(123)) 


#Plotting the final dbRDA

#extract the constrained eigenvalues
eig_constrained <- final_model$CCA$eig

#compute % explained by each axis
var1 <- eig_constrained[1] / sum(eig_constrained) * 100
var2 <- eig_constrained[2] / sum(eig_constrained) * 100

dbRDA_plot <- plot_ordination(
  physeq = ps_model25, 
  ordination = final_model,
  color = "ConsSiteName",
  shape = "Location",
  axes = c(1,2)
) +
  geom_point(aes(colour = ConsSiteName), alpha = 1, size = 5, stroke = 1.5) +
  scale_colour_paletteer_d("ggsci::default_ucscgb")

#Add the environmental variables as arrows
arrowmat <- vegan::scores(final_model, display = "bp")
arrowmat
## rename

#Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map <- aes(xend = dbRDA1, 
                 yend = dbRDA2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * dbRDA1, 
                 y = 1.3 * dbRDA2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

#Make a new graphic
dbRDA_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = 0.5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4.3, 
    #colour = "black", 
    fontface = "bold",
    data = arrowdf, 
    show.legend = FALSE
  ) +
  labs(
    x = paste0("dbRDA1 (", round(var1, 1), "%)"),
    y = paste0("dbRDA2 (", round(var2, 1), "%)")
  ) + labs(col = "Site Name")

#new plot
dbRDA_plot_2 <- plot_ordination(
  physeq = ps_model25, 
  ordination = final_model,
  color = "Location",
  shape = "Methane_conc",
  axes = c(1,2)
) +
  geom_point(aes(colour = Location), alpha = 1, size = 5, stroke = 1.5) +
  scale_colour_paletteer_d("ggsci::default_ucscgb")

#Add the environmental variables as arrows
arrowmat_2 <- vegan::scores(final_model, display = "bp")
arrowmat_2
## rename

#Add labels, make a data.frame
arrowdf_2 <- data.frame(labels = rownames(arrowmat_2), arrowmat_2)
# Define the arrow aesthetic mapping
arrow_map_2 <- aes(xend = dbRDA1, 
                   yend = dbRDA2, 
                   x = 0, 
                   y = 0, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)

label_map_2 <- aes(x = 1.3 * dbRDA1, 
                   y = 1.3 * dbRDA2, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)

arrowhead_2 = arrow(length = unit(0.02, "npc"))

arrowdf_2 <- data.frame(labels = rownames(arrowmat_2), arrowmat_2)

arrowdf_2$labels <- dplyr::recode(
  arrowdf_2$labels,
  CatchmentSize = "Catchment size",
  latitudeN     = "Latitude"
)

#Make a new graphic

dbRDA_plot_2 + 
  geom_segment(
    mapping = arrow_map_2, 
    size = 0.5, 
    data = arrowdf_2, 
    color = "gray", 
    arrow = arrowhead_2
  ) + 
  geom_text(
    mapping = label_map_2, 
    size = 4.3, 
    #colour = "black", 
    fontface = "bold",
    data = arrowdf_2, 
    show.legend = FALSE
  ) +
  labs(col = "Location", shape = "Methane concentration")+ theme_minimal()


ps_model1@sam_data[["Location"]] <- factor(ps_model1@sam_data[["Location"]], levels = c("Qaanaaq", "Upernavik", "Ilulissat", "Kangerlussuaq", "Nuuk", "Narsarsuaq"))

colors <- c("Qaanaaq" = "royalblue", "Upernavik" = "gold", "Ilulissat" = "red4", "Kangerlussuaq" = "#1F9698", "Nuuk" = "magenta", "Narsarsuaq" = "#783FC1")

dbRDA_plot_3 <- plot_ordination(
  physeq = ps_model25, 
  ordination = final_model,
  color = "Location",
  axes = c(1,2)
) 

#Add the environmental variables as arrows
arrowmat_2 <- vegan::scores(final_model, display = "bp")
arrowmat_2
## rename

#Add labels, make a data.frame
arrowdf_2 <- data.frame(labels = rownames(arrowmat_2), arrowmat_2)
# Define the arrow aesthetic mapping
arrow_map_2 <- aes(xend = dbRDA1, 
                   yend = dbRDA2, 
                   x = 0, 
                   y = 0, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)

label_map_2 <- aes(x = 1.3 * dbRDA1, 
                   y = 1.3 * dbRDA2, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)

arrowhead_2 = arrow(length = unit(0.02, "npc"))

arrowdf_2 <- data.frame(labels = rownames(arrowmat_2), arrowmat_2)


arrowdf_2$labels <- dplyr::recode(
  arrowdf_2$labels,
  CH4           = "CH[4]",
  PO4           = "PO[4]^\"3-\"",
  CatchmentSize = "Catchment~size",
  latitudeN     = "Latitude"
)

final_plot_dbRDA <- dbRDA_plot_3 + 
  geom_segment(
    mapping = arrow_map_2, 
    size = 0.5, 
    data = arrowdf_2, 
    color = "gray", 
    arrow = arrowhead_2
  ) + 
  geom_text(
    mapping = label_map_2, 
    size = 7, 
    #colour = "black", 
    fontface = "bold",
    data = arrowdf_2, 
    show.legend = FALSE,
    parse = TRUE
  ) +
  labs(col = "Location")+
  geom_point(aes(fill = Location),size = 6,shape = 21,stroke = 1.5,color = "black")+ 
  scale_color_manual(values = colors, name = "Region")+ 
  scale_fill_manual(values = colors, name = "Region")+
  theme_bw()+
  theme(
    legend.title = element_text(size = 16),  
    legend.text = element_text(size = 14), 
    legend.key.size = unit(0.75, "cm"),
    legend.position = "right",
    axis.title = element_text(size = 14),    
    axis.text = element_text(size = 12), 
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )

final_plot_dbRDA
