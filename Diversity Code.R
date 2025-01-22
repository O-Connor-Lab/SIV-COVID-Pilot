rm(list = ls())
dir.create(file.path('Figures', 'Diversity'), showWarnings = FALSE)

#Install Packages from CRAN:
# install.packages(c("tidyverse", "reshape2", "ggnewscale", "RColorBrewer", "ggtext", "ggrepel", "gghighlight", "cowplot", "ggpubr", "rstatix"))

#Install qiime2R:
#if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R")

#load libraries
library(tidyverse)
library(reshape2)
library(ggnewscale)
library(RColorBrewer)
library(ggtext)
library(ggrepel)
library(gghighlight)
library(cowplot)
library(ggpubr)
library(rstatix)

library(qiime2R)

#Reading in core file types
FeatureTable <- read_qza("Data/featuretable.qza")$data #Feature Table

metadata<-read.table("Data/Sample-Metadata-Numeric.txt", header = TRUE) #Read in metadata
colnames(metadata) <-  c('SampleID', 'patient', 'location', 'time') #make sure the column names are readable
metadata$location[metadata$location == "NasalSwab"] <- "Nasal Swab" #Adjust the location labels
metadata$location[metadata$location == "RectalSwab"] <- "Rectal Swab"
metadata$location[metadata$location == "ThroatSwab"] <- "Tracheal Swab"
metadata$time <- as.factor(metadata$time) #Turns numeric values into categorical values
metadata$location <- fct_relevel(metadata$location, "Nasal Swab", "Tracheal Swab", "Rectal Swab", "Stool")
metadata$time <- fct_relevel(metadata$time, "-7", "3", "5", "7", "10", "14")

#Alpha Diversity

#Reads in shannon diversity from QIIME2
shannon <- read_qza("Outputs/Diversity/rarefied-shannon.qza")$data %>%
  rownames_to_column(var='SampleID') %>%
  inner_join(metadata) %>%
  filter(location != 'control')
write.csv(shannon,"Figures/Supplemental/Shannon-Diversity.csv")

#Kruskal Wallis test for different compartments
kruskal_stool <- kruskal.test(shannon_entropy ~ time, data = shannon[which(shannon$location == "Stool"),])

kruskal_Rectal <- kruskal.test(shannon_entropy ~ time, data = shannon[which(shannon$location == "Rectal Swab"),])

kruskal_Nasal <- kruskal.test(shannon_entropy ~ time, data = shannon[which(shannon$location == "Nasal Swab"),])

kruskal_Tracheal <- kruskal.test(shannon_entropy ~ time, data = shannon[which(shannon$location == "Tracheal Swab"),])
pairwise.wilcox.test(shannon[which(shannon$location == "Tracheal Swab"),]$shannon_entropy, shannon[which(shannon$location == "Tracheal Swab"),]$time, p.adjust.method = "BH")

#Friedman test for different compartments
friedman_stool <-
  shannon[which(shannon$location == "Stool"),] %>%
  friedman_test(shannon_entropy ~ time | patient)

friedman_Rectal <-
  shannon[which(shannon$location == "Rectal Swab"),] %>%
  friedman_test(shannon_entropy ~ time | patient)

friedman_Nasal <-
  shannon[which(shannon$location == "Nasal Swab"),] %>%
  friedman_test(shannon_entropy ~ time | patient)

friedman_Tracheal <-
  shannon[which(shannon$location == "Tracheal Swab"),] %>%
  friedman_test(shannon_entropy ~ time | patient)

# Post Hoc test for Tracheal swabs
dunn_Tracheal <- 
  dunn_test(shannon[which(shannon$location == "Tracheal Swab"),], 
            shannon_entropy ~ time, 
            p.adjust.method = "none")
dunn_Tracheal

dunn_stool <- 
  dunn_test(shannon[which(shannon$location == "Stool"),], 
            shannon_entropy ~ time, 
            p.adjust.method = "none")
dunn_stool

#Plotting shannon diversity for different compartments
shannon_rectal <- ggplot(data = shannon[which(shannon$location == "Rectal Swab"),], aes(x = time, y = shannon_entropy, group = time)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  facet_wrap(~location, nrow = 1, scales = "free",  axes = 'all') +
  xlab("Days post-SARS-CoV-2") +
  ylab("Shannon Index") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 18, face = "bold"),
        text = element_text(family = "Arial"),
        axis.title = element_text(size = 16),
        axis.line = element_line(linewidth = 0.353),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text = element_text(size = 13)) 
shannon_rectal
ggsave("Figures/Diversity/rectal_shannon.png", plot = shannon_rectal, scale = 1, width = 5, height = 5) 

shannon_stool <- ggplot(data = shannon[which(shannon$location == "Stool"),], aes(x = time, y = shannon_entropy, group = time)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  facet_wrap(~location, nrow = 1, scales = "free",  axes = 'all') +
  xlab("Days post-SARS-CoV-2") +
  ylab("Shannon Index") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 18, face = "bold"),
        text = element_text(family = "Arial"),
        axis.title = element_text(size = 16),
        axis.line = element_line(linewidth = 0.353),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        aspect.ratio = 1,
        axis.text = element_text(size = 13),
        axis.title.y = element_blank()) 
shannon_stool

ggsave("Figures/Diversity/stool_shannon.png", plot = shannon_stool, scale = 1, width = 5, height = 5) 

shannon_nasal <- ggplot(data = shannon[which(shannon$location == "Nasal Swab"),], aes(x = time, y = shannon_entropy, group = time)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  facet_wrap(~location, nrow = 1, scales = "free",  axes = 'all') +
  xlab("Days post-SARS-CoV-2") +
  ylab("Shannon Index") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 18, face = "bold"),
        text = element_text(family = "Arial"),
        axis.title = element_text(size = 16),
        axis.line = element_line(linewidth = 0.353),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        axis.text = element_text(size = 13),
        aspect.ratio = 1) 
shannon_nasal
ggsave("Figures/Diversity/nasal_shannon.png", plot = shannon_nasal, scale = 1, width = 5, height = 5) 

shannon_Tracheal <- ggplot(data = shannon[which(shannon$location == "Tracheal Swab"),], aes(x = time, y = shannon_entropy, group = time)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  geom_bracket(
    xmin = "-7", xmax = "5", y.position = 7.75,
    label = "*", label.size = 8, size = 0.5) +
  facet_wrap(~location, nrow = 1, scales = "free",  axes = 'all') +
  xlab("Days post-SARS-CoV-2") +
  ylab("Shannon Index") +
  ylim(2, 8) +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 18, face = "bold"),
        text = element_text(family = "Arial"),
        axis.title = element_text(size = 16),
        axis.line = element_line(linewidth = 0.353),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        aspect.ratio = 1,
        axis.title.y = element_blank(),
        axis.text = element_text(size = 13)) 
shannon_Tracheal
ggsave("Figures/Diversity/Tracheal_shannon.png", plot = shannon_Tracheal, scale = 1, width = 5, height = 5) 


p_shannon <- plot_grid(shannon_stool, shannon_rectal, shannon_nasal, shannon_Tracheal, nrow = 1, align = 'v')
ggsave("Figures/Diversity/shannon_all.png", plot = p_shannon, scale = 1, width = 15, height = 5) 


# Observed Features between compartments and controls
observed_features <- read_qza("Outputs/Diversity/unrarefied-observed-features.qza")$data %>%
  rownames_to_column(var='SampleID') %>%
  inner_join(metadata) %>%
  mutate(location = if_else(location == 'control', if_else(grepl('S_Aureus', patient), "Positive Control", "Negative Control"), location)) %>%
  filter(location != 'Positive Control')
observed_features$location <- fct_relevel(observed_features$location, "Stool", "Rectal Swab", "Tracheal Swab", "Nasal Swab", "Negative Control") #Controls the order of the plots

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))

observed_features <- read_qza("Outputs/Diversity/unrarefied-observed-features.qza")$data %>%
  rownames_to_column(var='SampleID') %>%
  inner_join(metadata) %>%
  mutate(location = if_else(location == 'control', if_else(grepl('S_Aureus', patient), "Positive Control", "Negative Control"), location)) %>%
  filter(location != 'Positive Control')
observed_features$location <- fct_relevel(observed_features$location, "Stool", "Rectal Swab", "Tracheal Swab", "Nasal Swab", "Negative Control") #Controls the order of the plots

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))

p_observed_features_location <- ggplot(data = observed_features, aes(x = location, y = observed_features)) +
  geom_boxplot() +
  geom_point() +
  geom_pwc(stat = "pwc", method = "dunn_test",label = "{p.adj.signif}", vjust = 0.6, hide.ns = TRUE, p.adjust.method = "BH", 
           label.size = 6, size = 0.5, bracket.nudge.y = .10) +
  labs(y = "Observed ASVs", x = "Sampling Location") +
  theme_custom +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1), axis.title.x = element_blank()) +
  guides(fill= "none")

p_observed_features_location

ggsave("Figures/Diversity/observed_features_location.png", plot = p_observed_features_location, scale = 1, width = 5, height = 5) 

observed_features %>%
  group_by(location) %>%
  get_summary_stats()

kruskal.test(observed_features ~ location, data = observed_features)
pairwise.wilcox.test(observed_features$observed_features, observed_features$location, p.adjust.method = 'BH')

#Observed Features over time
observed_features_no_control <- 
  observed_features %>%
  filter(location != "control")
  
p_observed_features_time <- ggplot(data = observed_features_no_control, aes(x = time, y = observed_features)) +
  geom_boxplot() +
  geom_point() +
  labs(y = "Observed Features", x = "Days post-SARS-CoV-2") +
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 14, face = "bold"),
        text = element_text(family = "Arial"),
        axis.title = element_text(size = 12),
        axis.line = element_line(linewidth = 0.353),
        legend.text = element_text(size = 12),
        legend.title = element_blank(),
        aspect.ratio = 1,
        panel.spacing = unit(0.5, "cm", data = NULL)) 

kruskal.test(observed_features ~ time, data = observed_features_no_control)
pairwise.wilcox.test(observed_features_no_control$observed_features, observed_features_no_control$time, p.adjust.method = 'BH')


p_observed_features_time

# Beta Diversity - Unweighted unifrac
unweighted_unifrac<-read_qza("Outputs/Diversity/unweighted_unifrac_pcoa_results.qza") #Reads inputed PCoA plot using the qiime2r package
  
unweighted_unifrac_vector <- unweighted_unifrac$data$Vectors %>%
  select("SampleID", PC1, PC2) %>% #Selects the data from the pcoa plot
  left_join(metadata) %>%
  filter(SampleID != 'TDP_Broncoscope_BAL') %>%
  mutate(location = if_else(location == 'control', if_else(grepl('S_Aureus', patient), "Positive Control", "Negative Control"), location))
unweighted_unifrac_vector$location <- fct_relevel(unweighted_unifrac_vector$location, "Stool", "Rectal Swab", "Tracheal Swab", "Nasal Swab")

p_beta_diversity <- ggplot(data = unweighted_unifrac_vector, aes(x=PC1, y=PC2, color=location)) + #Plots the PC1 and PC2. Colors the dots based on the category of intrest
  geom_point(alpha = 0.75) + #Makes the dots semi-transparent to better see the plot
  scale_color_manual(values = c('#6699CC', '#004488', '#EE99AA', '#994455', "orange", "grey")) +
  xlab(paste("PC1: ", round(100*unweighted_unifrac$data$ProportionExplained[1]), "%")) + #Calculates the percentage covered by the x axis
  ylab(paste("PC2: ", round(100*unweighted_unifrac$data$ProportionExplained[2]), "%")) + #Calculates the percentage covered by the y axis
  stat_ellipse() + #Creates 95% confidence ellipses around the different categories
  theme_classic() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size = 14, face = "bold"),
        text = element_text(family = "Arial"),
        axis.title = element_text(size = 12),
        axis.line = element_line(linewidth = 0.353),
        legend.text = element_text(size = 12, family = "Arial"),
        legend.title = element_blank(),
        aspect.ratio = 1,
        panel.spacing = unit(0.5, "cm", data = NULL)) 
p_beta_diversity
ggsave("Figures/Diversity/unweighted_unifrac.png", plot = p_beta_diversity, scale = 1, width = 5, height = 5) 


