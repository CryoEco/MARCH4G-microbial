setwd("/Users/.../")

#load phyloseq object

library(ANCOMBC)
library(phyloseq)
library(dplyr)
library(ggplot2)

ps_rare
sample_data(ps_rare)

##### 1. Phylum #####
ancom.16s.phylum <- ancombc(data = ps_rare, tax_level = "Phylum", meta_data = NULL,
                            p_adj_method = "holm", prv_cut = 0.10,
                            lib_cut = 0, formula = "CH4", 
                            struc_zero = FALSE, neg_lb = FALSE, 
                            alpha = 0.05, n_cl = 2, verbose = TRUE); ancom.16s.phylum

ancom.16s.res.phylum <- ancom.16s.phylum$res; ancom.16s.res.phylum

df_lfc.phylum = data.frame(ancom.16s.res.phylum$lfc[, -1] * ancom.16s.res.phylum$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = ancom.16s.res.phylum$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se.phylum = data.frame(ancom.16s.res.phylum$se[, -1] * ancom.16s.res.phylum$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = ancom.16s.res.phylum$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se.phylum)[-1] = paste0(colnames(df_se.phylum)[-1], "SE")

df_fig_.phylum = df_lfc.phylum %>% 
  dplyr::left_join(df_se.phylum, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, CH4, CH4SE) %>%
  dplyr::filter(CH4 != 0) %>% 
  dplyr::arrange(desc(CH4)) %>%
  dplyr::mutate(direct = ifelse(CH4 > 0, "Positive LFC", "Negative LFC"))
df_fig_.phylum$taxon_id = factor(df_fig_.phylum$taxon_id, levels = df_fig_.phylum$taxon_id)
df_fig_.phylum$direct = factor(df_fig_.phylum$direct, 
                               levels = c("Positive LFC", "Negative LFC"))

p_.phylum = ggplot(data = df_fig_.phylum, 
                   aes(x = taxon_id, y = CH4, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = CH4 - CH4SE, ymax = CH4 + CH4SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of CH4 (Phylum)") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_.phylum

##### 2. Class #####
ancom.16s.class <- ancombc(data = ps_rare, tax_level = "Class", meta_data = NULL,
                           p_adj_method = "holm", prv_cut = 0.10,
                           lib_cut = 0, formula = "CH4", 
                           struc_zero = FALSE, neg_lb = FALSE, 
                           alpha = 0.05, n_cl = 2, verbose = TRUE); ancom.16s.class

ancom.16s.res.class <- ancom.16s.class$res; ancom.16s.res.class

df_lfc.class = data.frame(ancom.16s.res.class$lfc[, -1] * ancom.16s.res.class$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = ancom.16s.res.class$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se.class = data.frame(ancom.16s.res.class$se[, -1] * ancom.16s.res.class$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = ancom.16s.res.class$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se.class)[-1] = paste0(colnames(df_se.class)[-1], "SE")

df_fig_.class = df_lfc.class %>% 
  dplyr::left_join(df_se.class, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, CH4, CH4SE) %>%
  dplyr::filter(CH4 != 0) %>% 
  dplyr::arrange(desc(CH4)) %>%
  dplyr::mutate(direct = ifelse(CH4 > 0, "Positive LFC", "Negative LFC"))
df_fig_.class$taxon_id = factor(df_fig_.class$taxon_id, levels = df_fig_.class$taxon_id)
df_fig_.class$direct = factor(df_fig_.class$direct, 
                              levels = c("Positive LFC", "Negative LFC"))

p_.class = ggplot(data = df_fig_.class, 
                  aes(x = taxon_id, y = CH4, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = CH4 - CH4SE, ymax = CH4 + CH4SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of CH4 (Class)") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_.class


##### 3. Order #####
ancom.16s.order <- ancombc(data = ps_rare, tax_level = "Order", meta_data = NULL,
                           p_adj_method = "holm", prv_cut = 0.10,
                           lib_cut = 0, formula = "CH4", 
                           struc_zero = FALSE, neg_lb = FALSE, 
                           alpha = 0.05, n_cl = 2, verbose = TRUE); ancom.16s.order

ancom.16s.res.order <- ancom.16s.order$res; ancom.16s.res.order

df_lfc.order = data.frame(ancom.16s.res.order$lfc[, -1] * ancom.16s.res.order$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = ancom.16s.res.order$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se.order = data.frame(ancom.16s.res.order$se[, -1] * ancom.16s.res.order$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = ancom.16s.res.order$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se.order)[-1] = paste0(colnames(df_se.order)[-1], "SE")

df_fig_.order = df_lfc.order %>% 
  dplyr::left_join(df_se.order, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, CH4, CH4SE) %>%
  dplyr::filter(CH4 != 0) %>% 
  dplyr::arrange(desc(CH4)) %>%
  dplyr::mutate(direct = ifelse(CH4 > 0, "Positive LFC", "Negative LFC"))
df_fig_.order$taxon_id = factor(df_fig_.order$taxon_id, levels = df_fig_.order$taxon_id)
df_fig_.order$direct = factor(df_fig_.order$direct, 
                              levels = c("Positive LFC", "Negative LFC"))

p_.order = ggplot(data = df_fig_.order, 
                  aes(x = taxon_id, y = CH4, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = CH4 - CH4SE, ymax = CH4 + CH4SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of CH4 (Order)") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_.order

##### 4. Family #####
ancom.16s.family <- ancombc(data = ps_rare, tax_level = "Family", meta_data = NULL,
                            p_adj_method = "holm", prv_cut = 0.10,
                            lib_cut = 0, formula = "CH4", 
                            struc_zero = FALSE, neg_lb = FALSE, 
                            alpha = 0.05, n_cl = 2, verbose = TRUE); ancom.16s.family

ancom.16s.res.family <- ancom.16s.family$res; ancom.16s.res.family

df_lfc.family = data.frame(ancom.16s.res.family$lfc[, -1] * ancom.16s.res.family$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = ancom.16s.res.family$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se.family = data.frame(ancom.16s.res.family$se[, -1] * ancom.16s.res.family$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = ancom.16s.res.family$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se.family)[-1] = paste0(colnames(df_se.family)[-1], "SE")

df_fig_.family = df_lfc.family %>% 
  dplyr::left_join(df_se.family, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, CH4, CH4SE) %>%
  dplyr::filter(CH4 != 0) %>% 
  dplyr::arrange(desc(CH4)) %>%
  dplyr::mutate(direct = ifelse(CH4 > 0, "Positive LFC", "Negative LFC"))
df_fig_.family$taxon_id = factor(df_fig_.family$taxon_id, levels = df_fig_.family$taxon_id)
df_fig_.family$direct = factor(df_fig_.family$direct, 
                               levels = c("Positive LFC", "Negative LFC"))

p_.family = ggplot(data = df_fig_.family, 
                   aes(x = taxon_id, y = CH4, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = CH4 - CH4SE, ymax = CH4 + CH4SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of CH4 (Family)") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_.family

##### 5. Genus #####
ancom.16s.genus <- ancombc(data = ps_rare, tax_level = "Genus", meta_data = NULL,
                           p_adj_method = "holm", prv_cut = 0.10,
                           lib_cut = 0, formula = "CH4", 
                           struc_zero = FALSE, neg_lb = FALSE, 
                           alpha = 0.05, n_cl = 2, verbose = TRUE); ancom.16s.genus

ancom.16s.res.genus <- ancom.16s.genus$res; ancom.16s.res.genus

df_lfc.genus = data.frame(ancom.16s.res.genus$lfc[, -1] * ancom.16s.res.genus$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = ancom.16s.res.genus$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se.genus = data.frame(ancom.16s.res.genus$se[, -1] * ancom.16s.res.genus$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = ancom.16s.res.genus$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se.genus)[-1] = paste0(colnames(df_se.genus)[-1], "SE")

df_fig_.genus = df_lfc.genus %>% 
  dplyr::left_join(df_se.genus, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, CH4, CH4SE) %>%
  dplyr::filter(CH4 != 0) %>% 
  dplyr::arrange(desc(CH4)) %>%
  dplyr::mutate(direct = ifelse(CH4 > 0, "Positive LFC", "Negative LFC"))
df_fig_.genus$taxon_id = factor(df_fig_.genus$taxon_id, levels = df_fig_.genus$taxon_id)
df_fig_.genus$direct = factor(df_fig_.genus$direct, 
                              levels = c("Positive LFC", "Negative LFC"))

p_.genus = ggplot(data = df_fig_.genus, 
                  aes(x = taxon_id, y = CH4, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = CH4 - CH4SE, ymax = CH4 + CH4SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of CH4 (Genus)") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_.genus

##### For better interpretation, log-transform CH4 concentrations, so you can relate it better to the log fold changes in abundance

sample_data(ps_rare)$logCH4 <- log10(sample_data(ps_rare)$CH4)
sample_data(ps_rare)

##### 4. Family #####
ancom.16s.family2 <- ancombc(data = ps_rare, tax_level = "Family", meta_data = NULL,
                             p_adj_method = "holm", prv_cut = 0.10,
                             lib_cut = 0, formula = "logCH4", 
                             struc_zero = FALSE, neg_lb = FALSE, 
                             alpha = 0.05, n_cl = 2, verbose = TRUE); ancom.16s.family

ancom.16s.res.family2 <- ancom.16s.family2$res; ancom.16s.res.family2

write.table(ancom.16s.res.family2, "/Users/.../ancombc-CH4groups-family.tsv", sep = "\t", quote = F, col.names = NA)

df_lfc.family2 = data.frame(ancom.16s.res.family2$lfc[, -1] * ancom.16s.res.family2$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = ancom.16s.res.family2$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se.family2 = data.frame(ancom.16s.res.family2$se[, -1] * ancom.16s.res.family2$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = ancom.16s.res.family2$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se.family2)[-1] = paste0(colnames(df_se.family2)[-1], "SE")

df_fig_.family2 = df_lfc.family2 %>% 
  dplyr::left_join(df_se.family2, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, logCH4, logCH4SE) %>%
  dplyr::filter(logCH4 != 0) %>%
  dplyr::arrange(desc(logCH4)) %>%
  dplyr::mutate(direct = ifelse(logCH4 > 0, "Positive LFC", "Negative LFC"))
df_fig_.family2$taxon_id = factor(df_fig_.family2$taxon_id, levels = df_fig_.family2$taxon_id)
df_fig_.family2$direct = factor(df_fig_.family2$direct, 
                                levels = c("Positive LFC", "Negative LFC"))

##Plot only the top 10 positive and negative LFC
df_fig_filtered <- df_fig_.family2 %>%
  top_n(10, wt = logCH4) %>%
  bind_rows(df_fig_.family2 %>% top_n(-10, wt = logCH4)) %>%
  mutate(taxon_id = factor(taxon_id, levels = taxon_id[order(logCH4, decreasing = TRUE)]))

p_.family2 = ggplot(data = df_fig_filtered, 
                    aes(x = taxon_id, y = logCH4, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = logCH4 - logCH4SE, ymax = logCH4 + logCH4SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of CH4 (Family)") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_.family2

##### 5. Genus #####
ancom.16s.genus2 <- ancombc(data = ps_rare, tax_level = "Genus", meta_data = NULL,
                            p_adj_method = "holm", prv_cut = 0.10,
                            lib_cut = 0, formula = "logCH4", 
                            struc_zero = FALSE, neg_lb = FALSE, 
                            alpha = 0.05, n_cl = 2, verbose = TRUE); ancom.16s.genus

ancom.16s.res.genus2 <- ancom.16s.genus2$res; ancom.16s.res.genus2

write.table(ancom.16s.res.genus2, "/Users/.../ancombc-CH4groups-genus.tsv", sep = "\t", quote = F, col.names = NA)

df_lfc.genus2 = data.frame(ancom.16s.res.genus2$lfc[, -1] * ancom.16s.res.genus2$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = ancom.16s.res.genus2$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se.genus2 = data.frame(ancom.16s.res.genus2$se[, -1] * ancom.16s.res.genus2$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = ancom.16s.res.genus2$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se.genus2)[-1] = paste0(colnames(df_se.genus2)[-1], "SE")
df_q.genus2 = data.frame(ancom.16s.res.genus2$q_val[, -1] * ancom.16s.res.genus2$diff_abn[, -1], check.names = FALSE) %>% ##include fdr corrected p-values
  mutate(taxon_id = ancom.16s.res.genus2$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_q.genus2)[-1] = paste0(colnames(df_q.genus2)[-1], "qval")

df_fig_.genus2 = df_lfc.genus2 %>% 
  dplyr::left_join(df_se.genus2, by = "taxon_id") %>%
  dplyr::left_join(df_q.genus2, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, logCH4, logCH4SE, logCH4qval) %>%
  dplyr::filter(logCH4 != 0) %>% 
  dplyr::arrange(desc(logCH4)) %>%
  dplyr::mutate(direct = ifelse(logCH4 > 0, "Positive LFC", "Negative LFC"))
df_fig_.genus2$taxon_id = factor(df_fig_.genus2$taxon_id, levels = df_fig_.genus2$taxon_id)
df_fig_.genus2$direct = factor(df_fig_.genus2$direct, 
                               levels = c("Positive LFC", "Negative LFC"))

p_.genus2 = ggplot(data = df_fig_.genus2, 
                   aes(x = taxon_id, y = logCH4, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = logCH4 - logCH4SE, ymax = logCH4 + logCH4SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of CH4 (Genus)") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_.genus2

##only top 10 pos and neg
df_fig_filtered_genus <- df_fig_.genus2 %>%
  top_n(10, wt = logCH4) %>%
  bind_rows(df_fig_.genus2 %>% top_n(-10, wt = logCH4)) %>%
  mutate(taxon_id = factor(taxon_id, levels = taxon_id[order(logCH4, decreasing = TRUE)]))

p_.genus3 = ggplot(data = df_fig_filtered_genus, 
                   aes(x = taxon_id, y = logCH4, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = logCH4 - logCH4SE, ymax = logCH4 + logCH4SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of CH4 (Genus)") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_.genus3

##Filter based on q-value (p-value fdr corrected)
df_fig_filtered_genus_q <- df_fig_.genus2 %>%
  filter(logCH4qval < 0.05) %>%
  top_n(10, wt = logCH4) %>%
  bind_rows(df_fig_.genus2 %>% filter(logCH4qval < 0.05) %>% top_n(-10, wt = logCH4)) %>%
  mutate(taxon_id = factor(taxon_id, levels = taxon_id[order(logCH4, decreasing = TRUE)]))

p_.genus4 = ggplot(data = df_fig_filtered_genus_q, 
                   aes(x = taxon_id, y = logCH4, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = logCH4 - logCH4SE, ymax = logCH4 + logCH4SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Top 10 genera positively and negatively associated with CH4") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_.genus4



df_fig_filtered_genus_q2 <- df_fig_.genus2 %>%
  filter(logCH4qval < 0.05) %>%
  bind_rows(df_fig_.genus2 %>% filter(logCH4qval < 0.05) %>% 
              mutate(taxon_id = factor(taxon_id, levels = taxon_id[order(logCH4, decreasing = TRUE)])))

p_.genus5 = ggplot(data = df_fig_filtered_genus_q2, 
                   aes(x = taxon_id, y = logCH4, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = logCH4 - logCH4SE, ymax = logCH4 + logCH4SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Taxon Log Fold Change", 
       title = expression("Log-fold changes per 10-fold increase in CH"[4(aq)]*" concentration")) + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1, size = 14, vjust = 1),
        axis.title.y = element_text(size = 14),
        plot.margin = margin(10, 10, 30, 25),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.3, "lines"))  #the order is top, right, bottom, left
p_.genus5
ggsave("ancom-bc-allgenera.jpg", p_.genus5, dpi = 300, width = 12)



signif_taxa_lfc <- ancom.16s.res.genus2$q_val %>%
  filter(logCH4 < 0.05)

signif_taxa_lfc

diff_abund <- ancom.16s.res.genus2$diff_abn$logCH4 
diff_abund
str(diff_abund)
diff_abund <- as.data.frame(diff_abund)
str(diff_abund)
diff_abund_true <- diff_abund %>%
  filter(diff_abund == "TRUE")
diff_abund_true
