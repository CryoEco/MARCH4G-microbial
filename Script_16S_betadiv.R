setwd("/Users/...")

# load phyloseq object

##Beta-diversity - dissimilarities among samples (use rarefied dataset, and create distance matrix from hellinger transformed count data)

ps_rare

ps_hell <- microbiome::transform(ps_rare, 'hell')
ps_hell
otu_table(ps_rare)
otu_table(ps_hell)

#PCoA with bray-curtis dissimilarity matrix using MicroViz
library(RColorBrewer)

distinct_palette(n=NA, pal = "brewerPlus", add = "lightgrey")
brewerPlus <- distinct_palette()
RColorBrewer::brewer.pal(41, brewerPlus)
library(paletteer)

ps_rare %>%
  tax_transform("hell", rank = "Genus") %>% 
  dist_calc("bray") %>%
  ord_calc("PCoA") %>%
  ord_plot(color = "Site", shape = "Location", size = 5) + scale_colour_paletteer_d("ggsci::default_ucscgb")

bc <- phyloseq::distance(ps_hell, method = "bray")
bc.pcoa <- ordinate(ps_hell, method = "MDS", distance = bc)
bc.pcoa2 <- ordinate(ps_hell, method = "PCoA", distance = bc)
bc.egvals <- bc.pcoa$values$Eigenvalues
bc.egvals2 <- bc.pcoa2$values$Eigenvalues

plot_ordination(ps_hell, bc.pcoa, color = "ConsSiteName", shape = "Location")+
  labs(col = "Site Name")+
  coord_fixed(sqrt(bc.egvals[2]/bc.egvals[1]))+
  geom_point(size = 5)

bray_pcoa <- plot_ordination(ps_hell, bc.pcoa2, color = "ConsSiteName", shape = "Location")+
  labs(col = "Site Name")+
  coord_fixed(sqrt(bc.egvals2[2]/bc.egvals2[1]))+
  geom_point(size = 5)+ scale_colour_paletteer_d("ggsci::default_ucscgb")

bray_pcoa
ggsave("pcoa-bray-rare.jpeg", dpi = 300)

plot_ordination(ps_hell, bc.pcoa2, color = "ConsSiteName", shape = "Location")+
  labs(col = "Site Name")+
  coord_fixed(sqrt(bc.egvals2[2]/bc.egvals2[1]))+
  geom_point(size = 6)+ scale_colour_paletteer_d("ggsci::default_ucscgb")+ theme_bw()


plot_ordination(ps_hell, bc.pcoa2, color = "Location", shape = "Methane_conc")+
  labs(col = "Location", shape = "Methane concentration")+
  coord_fixed(sqrt(bc.egvals2[2]/bc.egvals2[1]))+
  geom_point(size = 6)+ scale_colour_brewer(palette = "Dark2")+ theme_bw()


##Plot PCoA as regions

#Order locations by decreasing latitude

ps_hell@sam_data[["Location"]] <- factor(ps_hell@sam_data[["Location"]], levels = c("Qaanaaq", "Upernavik", "Ilulissat", "Kangerlussuaq", "Nuuk", "Narsarsuaq"))

colors <- c("Qaanaaq" = "royalblue", "Upernavik" = "gold", "Ilulissat" = "red4", "Kangerlussuaq" = "#1F9698", "Nuuk" = "magenta", "Narsarsuaq" = "#783FC1")

final_pcoa <- plot_ordination(ps_hell, bc.pcoa2, color = "Location")+
  labs(col = "Location")+
  coord_fixed(sqrt(bc.egvals2[2]/bc.egvals2[1]))+
  geom_point(aes(fill = Location),size = 6,shape = 21,stroke = 1.5,color = "black")+ 
  scale_color_manual(values = colors, name = "Region")+ 
  scale_fill_manual(values = colors, name = "Region")+
  theme_bw()+
  theme(
    legend.title = element_text(size = 16),  
    legend.text = element_text(size = 14),   
    axis.title = element_text(size = 14),    
    axis.text = element_text(size = 12),     
  ) 

final_pcoa

#Extract metadata
metadata_16S <- data.frame(sample_data(ps_hell))
metadata_16S

#ADONIS test (PERMANOVA)
region <- vegan::adonis2(bc ~ phyloseq::sample_data(ps_hell)$Location, set.seed(123))
region
region$R2

samplingSite <- vegan::adonis2(bc ~ phyloseq::sample_data(ps_hell)$ConsSiteName, set.seed(123))
samplingSite


#Dispersion test and plot
dispr <- vegan::betadisper(bc, phyloseq::sample_data(ps_hell)$Location)
dispr
plot(dispr, main = "Ordination Centroids and Dispersion Labeled: Bray Curtis Distance", sub = "")
boxplot(dispr, main = "", xlab = "")
set.seed(123)
permutest(dispr)

dispr2 <- vegan::betadisper(bc, phyloseq::sample_data(ps_hell)$ConsSiteName)
dispr2
set.seed(123)
permutest(dispr2)
