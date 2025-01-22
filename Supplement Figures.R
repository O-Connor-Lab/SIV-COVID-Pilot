library(tidyverse)
library(reshape2)
library(ggnewscale)
library(RColorBrewer)
library(ggtext)
library(ggrepel)
library(gghighlight)
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


#Reads in shannon diversity from QIIME2
shannon <- read_qza("Outputs/Diversity/rarefied-shannon.qza")$data %>%
  rownames_to_column(var='SampleID') %>%
  inner_join(metadata) %>%
  filter(location != 'control')
write.csv(shannon,"Figures/Supplemental/Shannon-Diversity.csv")

#Read in taxonomy
taxonomy<-read_qza("Data/taxonomy_2.qza") #read in taxonomy
write.csv(taxonomy$data,"Figures/Supplemental/Taxonomy.csv")
taxonomy<-parse_taxonomy(taxonomy$data) #Removes the confidence from the taxonomy and breaks up the taxonomy

#Relative Abundance (Genus Level)
taxasums <- summarize_taxa(FeatureTable, taxonomy)$Genus
write.csv(taxasums, "Figures/Supplemental/Feature_Table.csv", row.names = TRUE)

col <- c("f_Lachnospiraceae" = '#DDD8EF',
         "f_Rhodocyclaceae" = '#D1C1E1', 
         "f_Ruminococcaceae" = '#C3A8D1',
         "g_[Eubacterium]_coprostanoligenes_group" = '#B58FC2',
         "g_Acinetobacter" = '#A778B4',
         "g_Bacteroidales_RF16_group" = '#9B62A7', 
         "g_Blautia" = '#8C4E99', 
         "g_Bradymonadales" = '#6F4C9B', 
         "g_Campylobacter" = '#6059A9', 
         "g_Christensenellaceae_R-7_group" = '#5568B8', 
         "g_Clostridia_UCG-014" = '#4E79C5', 
         "g_Clostridia_vadinBB60_group" = '#4D8AC6', 
         "g_Escherichia-Shigella" = '#4E96BC', 
         "g_Fibrobacter" = '#549EB3', 
         "g_Helicobacter" = '#59A5A9', 
         "g_Intestinibacter" = '#60AB9E',
         "g_Lactobacillus" = '#69B190', 
         "g_Muribaculaceae" = '#77B77D', 
         "g_Myroides" = '#8CBC68', 
         "g_NK4A214_group" = '#A6BE54', 
         "g_Prevotella" = '#BEBC48', 
         "g_Rikenellaceae_RC9_gut_group" = '#D1B541',
         "g_Ruminococcus" = '#DDAA3C', 
         "g_Sarcina" = '#E49C39',
         "g_Streptococcus" = '#E78C35',
         "g_Succinivibrio" = '#E67932',
         "g_UCG-002" = '#E4632D',
         "g_UCG-005" = '#DA2222', 
         "f_Moraxellaceae" = '#88CCEE', 
         "f_Pasteurellaceae" = '#99DDFF',
         "g_Alloprevotella" = '#77AADD', 
         "g_Corynebacterium" = '#44AA99', 
         "g_Dolosigranulum" = '#CCDDAA',
         "g_Filobacterium" = '#BBCC33',
         "g_Fusobacterium" = '#AAAA00',
         "g_Gemella" = '#DDCC77',
         "g_Haemophilus" = '#FFCCCC',
         "g_Porphyromonas" = '#FFAABB',
         "g_Psychrobacter" = '#EE8866',
         "g_Staphylococcus" = '#CC6677',
         "g_Veillonella" = '#95211B',
         "f_Enterobacteriaceae" = "#C79070",
         "g_Adhaeribacter" = "#9075BF",
         "g_Arthrobacter" = '#90C9BC',
         "g_Cutibacterium" = '#C32D6F',
         "g_Sphingomonas" = "#729C33",
         "Other" = "#777777")

theme_custom <- list(theme_classic() + theme(text = element_text(family = "Arial", colour = 'black'),
                                             strip.background = element_blank(),
                                             strip.text = element_text(size = 14, face = "bold"),
                                             plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                                             axis.text = element_text(size = 14, family = "Arial"),
                                             axis.title = element_text(size = 14, family = "Arial"),
                                             aspect.ratio = 1))

#11A - Observed Features
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

ggsave("Figures/Diversity/observed_features_location.png", plot = p_observed_features_location, scale =1.5, height = 3, width = 3, unit = "in") 

observed_features %>%
  group_by(location) %>%
  get_summary_stats()

kruskal.test(observed_features ~ location, data = observed_features)
pairwise.wilcox.test(observed_features$observed_features, observed_features$location, p.adjust.method = 'BH')


#11B - RA controls
taxa_long_genus_control <-
  taxasums %>%
  as.data.frame() %>% #Makes it a dataframe
  rownames_to_column("Taxon") %>% #Extracts the row names into a new column
  melt(id.vars = "Taxon", variable.name = "SampleID", value.name = "Abundance", na.rm = TRUE) %>% #Makes the dataframe long
  merge(metadata, by.x = 'SampleID', all = TRUE) %>%#Adds in the metadata for each sample
  filter(Abundance != 0) %>% # Filters out 0s to speed things up
  filter(location == 'control') %>% #Deletes the controls
  mutate(location = if_else(grepl('S_Aureus', patient), "Positive Control", "Negative Control")) %>%
  filter(location == 'Negative Control') %>% #Deletes the Postive controls
  group_by(SampleID) %>%
  mutate(Total_Abundance = sum(Abundance)) %>%
  ungroup()  %>%
  mutate(Relative_Abundance = Abundance/Total_Abundance) %>%
  mutate(LTC = if_else(Relative_Abundance > 0.05, Taxon, "Other")) %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";") %>%
  mutate(LTC = if_else(LTC == "Other", "Other", paste0("g_", str_trim(Genus)))) %>%
  mutate(LTC = if_else(LTC =="g_NA", paste0("f_", str_trim(Family)), LTC))

relative_abundance_control <- ggplot(taxa_long_genus_control, aes(x = patient,
                                                                  y = Abundance, fill = LTC)) +
  geom_bar(stat = 'identity', position = "fill") +
  labs(y = " Relative Abundance", x = "Negative Controls") +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00), 
                     labels = c("0%", "25%", "50%", "75%", "100%")) +
  scale_fill_manual(values = col) +
  theme_custom +
  theme(axis.text.x = element_blank(), legend.key.size = unit(12, "pt"), legend.text = element_markdown(size = 10), legend.position = "bottom") +
  #guides(fill = "none")
  guides(fill=guide_legend(title=""))
relative_abundance_control 
ggsave("Figures/Relative Abundance/relative_abundnace_control.png", plot = relative_abundance_control, scale =1.5, height = 3, width = 7.5, unit = "in")


most_abundant_taxa_genus_control <- 
  taxa_long_genus_control %>%
  group_by(SampleID) %>%
  filter(Relative_Abundance == max(Relative_Abundance))

#11C - Beta diversity
unweighted_unifrac<-read_qza("Outputs/Diversity/unweighted_unifrac_pcoa_results.qza") #Reads inputed PCoA plot using the qiime2r package

unweighted_unifrac_vector <- unweighted_unifrac$data$Vectors %>%
  select("SampleID", PC1, PC2) %>% #Selects the data from the pcoa plot
  left_join(metadata) %>%
  filter(SampleID != 'TDP_Broncoscope_BAL') %>%
  mutate(location = if_else(location == 'control', if_else(grepl('S_Aureus', patient), "Positive Control", "Negative Control"), location))
unweighted_unifrac_vector$location <- fct_relevel(unweighted_unifrac_vector$location, "Rectal Swab", "Stool", "Nasal Swab", "Tracheal Swab")

p_beta_diversity <- ggplot(data = unweighted_unifrac_vector, aes(x=PC1, y=PC2, color=location)) + #Plots the PC1 and PC2. Colors the dots based on the category of intrest
  geom_point(alpha = 0.75) + #Makes the dots semi-transparent to better see the plot
  scale_color_manual(values = c('#6699CC', '#004488', '#EE99AA', '#994455', "orange", "grey")) +
  xlab(paste("PC1: ", round(100*unweighted_unifrac$data$ProportionExplained[1]), "%")) + #Calculates the percentage covered by the x axis
  ylab(paste("PC2: ", round(100*unweighted_unifrac$data$ProportionExplained[2]), "%")) + #Calculates the percentage covered by the y axis
  stat_ellipse() + #Creates 95% confidence ellipses around the different categories
  theme_custom +
  theme(legend.key.size = unit(12, "pt"), legend.text = element_markdown(size = 10)) +
  #guides(color= guide_legend(title=""))
  guides(color = "none")
p_beta_diversity
ggsave("Figures/Diversity/unweighted_unifrac.png", plot = p_beta_diversity, scale =1.5, height = 3, width = 3, unit = "in") 

#12 Relative Abundance (Phylum Level)
taxasums_phylum <- summarize_taxa(FeatureTable, taxonomy)$Phylum

taxa_long_phylum <-
  taxasums_phylum %>%
  as.data.frame() %>% #Makes it a dataframe
  rownames_to_column("Taxon") %>% #Extracts the row names into a new column
  melt(id.vars = "Taxon", variable.name = "SampleID", value.name = "Abundance", na.rm = TRUE) %>% #Makes the dataframe long
  merge(metadata, by.x = 'SampleID', all = TRUE) %>%#Adds in the metadata for each sample
  filter(Abundance != 0) %>% # Filters out 0s to speed things up
  filter(location != 'control') %>% #Deletes the controls
  mutate(Taxon = str_remove(Taxon, "^[^;]*;")) %>%
  group_by(SampleID) %>%
  mutate(Total_Abundance = sum(Abundance)) %>%
  ungroup()  %>%
  mutate(Relative_Abundance = Abundance/Total_Abundance) %>%
  mutate(Taxon = if_else(Relative_Abundance > 0.01, Taxon, "Other"))

col_p <- c(" Actinobacteriota" = '#CAACCB',
           " Bacteroidota" = '#AE76A3', 
           " Campilobacterota" = '#882E72', 
           " Cyanobacteria" = '#1965B0',
           " Desulfobacterota" = '#5289C7',	
           " Elusimicrobiota" =  '#7BAFDE', 
           " Fibrobacterota" = '#4EB265', 
           " Firmicutes" = '#90C987', 
           " Fusobacteriota" = '#CAE0AB', 
           " Patescibacteria" = '#F4A736', 
           " Proteobacteria" = '#E8601C', 
           " Spirochaetota" = '#E65518',
           " WPS-2" = '#DC050C', 
           "Other" = "#777777")

relative_abundance_phylum <- ggplot(taxa_long_phylum, aes(x = time,
                                                                y = Abundance, fill = Taxon)) +
  geom_bar(stat = 'identity', position = "fill") +
  facet_grid(location ~ patient) +
  labs(x = "Days post-SARS-CoV-2", y = " Relative Abundance (%)") +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00), 
                     labels = c("0", "25", "50", "75", "100")) +
  scale_fill_manual(values = col_p) +
  theme_custom +
  guides(fill=guide_legend(title=""))
relative_abundance_phylum

ggsave("Figures/Relative Abundance/relative_abundnace_phylum.png", plot = relative_abundance_phylum, scale =1.5, height = 5, width = 9, unit = "in")

# Differential Abundance Load from ANCOMBC script and transform to log2
DA_stool<-read.csv("Outputs/Ancombc2/ancombc2-primary-results-stool.csv") %>%
  mutate(across(lfc_timeD3:lfc_timeD14, function(x) exp(x))) %>%
  mutate(across(lfc_timeD3:lfc_timeD14, function(x) log2(x)))

DA_Tracheal<-read.csv("Outputs/Ancombc2/ancombc2-primary-results-throat.csv") %>%
  mutate(across(lfc_timeD3:lfc_timeD14, function(x) exp(x))) %>%
  mutate(across(lfc_timeD3:lfc_timeD14, function(x) log2(x))) 

DA_nasal<-read.csv("Outputs/Ancombc2/ancombc2-primary-results-nasal.csv") %>%
  mutate(across(lfc_timeD3:lfc_timeD14, function(x) exp(x))) %>%
  mutate(across(lfc_timeD3:lfc_timeD14, function(x) log2(x))) 

#Volcano plots
DA_stool_to_plot <- DA_stool %>%  
  select(taxon, lfc_timeD3:lfc_timeD14, q_timeD3:q_timeD14) %>% 
  melt("taxon") %>%
  separate(variable, into = c("variable", "time"), sep = "_") %>%
  dcast(taxon + time ~ variable) %>%
  separate(taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = "(?<=[a-zA-Z])__") %>%
  mutate(across(Kingdom:Family,function(x) str_remove(x, pattern = "_.$"))) %>%
  rowwise() %>%
  mutate(LTC = if_else(!is.na(Genus),paste0("g_", Genus), if_else(!is.na(Family), paste0("f_",Family), if_else(!is.na(Order), Order, if_else(!is.na(Class), Class, if_else(!is.na(Phylum), Phylum, Kingdom))))))

DA_stool_to_plot$time = str_replace(DA_stool_to_plot$time, "timeD", "Day ")
DA_stool_to_plot$time = fct_relevel(DA_stool_to_plot$time, "Day 3", "Day 5", "Day 7", "Day 10", "Day 14")

DA_Tracheal_to_plot<-DA_Tracheal %>%
  select(taxon, lfc_timeD3:lfc_timeD14, q_timeD3:q_timeD14) %>% 
  melt("taxon") %>%
  separate(variable, into = c("variable", "time"), sep = "_") %>%
  dcast(taxon + time ~ variable) %>%
  separate(taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = "(?<=[a-zA-Z])__") %>%
  mutate(across(Kingdom:Family,function(x) str_remove(x, pattern = "_.$"))) %>%
  rowwise() %>%
  mutate(LTC = if_else(!is.na(Genus),paste0("g_", Genus), if_else(!is.na(Family), paste0("f_",Family), if_else(!is.na(Order), Order, if_else(!is.na(Class), Class, if_else(!is.na(Phylum), Phylum, Kingdom))))))


DA_Tracheal_to_plot$time = str_replace(DA_Tracheal_to_plot$time, "timeD", "Day ")
DA_Tracheal_to_plot$time = fct_relevel(DA_Tracheal_to_plot$time, "Day 3", "Day 5", "Day 7", "Day 10", "Day 14")
DA_nasal_to_plot<-DA_nasal %>%
  select(taxon, lfc_timeD3:lfc_timeD14, q_timeD3:q_timeD14) %>% 
  melt("taxon") %>%
  separate(variable, into = c("variable", "time"), sep = "_") %>%
  dcast(taxon + time ~ variable) %>%
  separate(taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = "(?<=[a-zA-Z])__") %>%
  mutate(across(Kingdom:Family,function(x) str_remove(x, pattern = "_.$"))) %>%
  rowwise() %>%
  mutate(LTC = if_else(!is.na(Genus),paste0("g_", Genus), if_else(!is.na(Family), paste0("f_",Family), if_else(!is.na(Order), Order, if_else(!is.na(Class), Class, if_else(!is.na(Phylum), Phylum, Kingdom))))))

DA_nasal_to_plot$time = str_replace(DA_nasal_to_plot$time, "timeD", "Day ")
DA_nasal_to_plot$time = fct_relevel(DA_nasal_to_plot$time, "Day 3", "Day 5", "Day 7", "Day 10", "Day 14")

p_stool_volcano <- ggplot(data = DA_stool_to_plot) +
  geom_point(aes(x = lfc, y = -log(q), color = time), size = 2) +
  geom_hline(yintercept=-log(0.05), linetype='dotted') +
  xlim(-10,10) +
  ylim(0, 15) +
  gghighlight(q < 0.05, use_direct_label = FALSE) +
  geom_label_repel(aes(x = lfc, y = -log(q), label = if_else(Genus == 'Succinivibrio' | Genus == 'Streptococcus', Genus, "")), box.padding = 0.75, size = 5)+
  #geom_label_repel(aes(x = lfc, y = -log(q), label = LTC, segment.size = 0.5, fontface = if_else(LTC == "g_Streptococcus" | LTC == "g_Succinivibrio", "bold", "plain")), force_pull = 0, size = 3, max.overlaps = Inf, box.padding = 1, force = 3)+
  labs(x = expression("log"[2]*" Fold Change"), y = expression("-log"[10]*"FDR"), title = "Stool") +
  scale_color_manual("Timepoint", values = c("Day 3" = '#CC3311',"Day 5" = '#0077BB',"Day 7" = '#009988','Day 10' = '#EE7733', 'Day 14' = '#DDCC77')) +
  theme_custom +
  guides(color = "none") 
p_stool_volcano
ggsave(file="Figures/Differential Abundance/DA_stool_volcano.png", plot=p_stool_volcano, scale =1.5, height = 3, width = 3, unit = "in")

p_nasal_volcano <- ggplot(data = DA_nasal_to_plot) +
  geom_point(aes(x = lfc, y = -log(q), color = time),size = 2) +
  geom_hline(yintercept=-log(0.05), linetype='dotted') +
  xlim(-10,10) +
  ylim(0, 15) +
  gghighlight(q < 0.05, label_key = Genus, label_params = list(color = "black", size = 5, box.padding = 2)) +
  labs(x = expression("log"[2]*" Fold Change"), y = expression("-log"[10]*"FDR"), title = "Nasal Swab") +
  scale_color_manual("Timepoint", values = c("Day 3" = '#CC3311',"Day 5" = '#0077BB',"Day 7" = '#009988','Day 10' = '#EE7733', 'Day 14' = '#DDCC77')) +
  theme_custom +
  guides(color = "none") 
p_nasal_volcano
ggsave(file="Figures/Differential Abundance/DA_nasal_volcano.png", plot=p_nasal_volcano, scale =1.5, height = 3, width = 3, unit = "in")

p_Tracheal_volcano <- ggplot(data = DA_Tracheal_to_plot) +
  geom_point(aes(x = lfc, y = -log(q), color = time), size = 2) +
  geom_hline(yintercept=-log(0.05), linetype='dotted') +
  ylim(0, 15) +
  xlim(-10,10) +
  gghighlight(q < 0.05, label_key = Genus, label_params = list(color = "black", size = 5, box.padding = 2.5, force = 3)) +
  labs(x = expression("log"[2]*" Fold Change"), y = expression("-log"[10]*"FDR"), title = "Tracheal Swab") +
  scale_color_manual("Timepoint", values = c("Day 3" = '#CC3311',"Day 5" = '#0077BB',"Day 7" = '#009988','Day 10' = '#EE7733', 'Day 14' = '#DDCC77')) +
  theme_custom +
  guides(color = "none") 
p_Tracheal_volcano
ggsave(file="Figures/Differential Abundance/DA_Tracheal_volcano.png", plot=p_Tracheal_volcano, scale =1.5, height = 3, width = 3, unit = "in")

