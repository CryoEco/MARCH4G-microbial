setwd("/.../")

#load table with site names and their respective CH4 concentrations

library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

###Separate CH4 concentrations and sites
df_CH4 <- df_pd %>%
  select(ConsSiteName, Methane_conc, CH4) 
str(df_CH4)

df_CH4$logCH4 <- log(df_CH4$CH4)

CH4log <- decostand(df_CH4[,3], method = "log", logbase = 10)
CH4log

df_CH4_log <- cbind(df_CH4, CH4log)
df_CH4_log

df_CH4_log$CH4_sqrt <- sqrt(df_CH4_log$CH4)


#Plot distribution of CH4 with red lines separating each range

# Define your boundaries (log-transformed)
boundaries <- c(log(4.9), log(20), log(40), log(80), log(490))
boundary_labels <- c("4.9", "20", "40", "80", "490")

p <- ggplot(df_CH4, 
            aes(x = logCH4,
                y = reorder(ConsSiteName, logCH4))) +
  geom_point(size = 5) +
  
  # Add vertical lines
  geom_vline(xintercept = boundaries, 
             color = "red", 
             linetype = "dashed", 
             linewidth = 1) +
  
  # Add labels for the boundaries
  geom_text(data = data.frame(x = boundaries, 
                              y = Inf, 
                              label = boundary_labels),
            aes(x = x, y = y, label = label),
            vjust = 2,        # Position above the line
            color = "red",
            size = 4) +
  
  labs(y = "Site") +
  theme_bw() +
  xlab(bquote("Log-transformed"~CH[4]~"concentrations"))+
  theme(
    # Increase legend text size
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    
    # Increase legend key (color box) size
    legend.key.size = unit(0.75, "cm"),
    
    # Increase axis label size
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    
    # Increase axis tick label size
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    
    # Adjust legend position needed
    legend.position = "right"
  )

p
