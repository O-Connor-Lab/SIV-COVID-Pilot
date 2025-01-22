#Load in packages
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
metadata$location[metadata$location == "RectalSwab"] <- "Rectal Swab" #Adjust the location labels
metadata$location[metadata$location == "ThroatSwab"] <- "Tracheal Swab" #Adjust the location labels
metadata$time <- as.factor(metadata$time) #Turns numeric values into categorical values
metadata$location <- fct_relevel(metadata$location, "Nasal Swab", "Tracheal Swab", "Rectal Swab", "Stool") #Relevels location column
metadata$time <- fct_relevel(metadata$time, "-7", "3", "5", "7", "10", "14") #Relevels time to be chronological


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

#Saves relative abundance with taxonomical information
taxasums <- summarize_taxa(FeatureTable, taxonomy)$Genus
write.csv(taxasums, "Figures/Supplemental/Feature_Table.csv", row.names = TRUE)

#Keeps colors consistant across figures
col <- c("f_Lachnospiraceae" = '#DDD8EF',
         "f_Rhodocyclaceae" = '#D1C1E1', 
         "f_Ruminococcaceae" = '#C3A8D1',
         "g_[Eubacterium]_coprostanoligenes_group" = '#B58FC2',
         "g_[Eubacterium]_ruminantium_group" = '#B58B97',
         "g_Acinetobacter" = '#A778B4',
         "g_Agathobacter" = '#b946aa', 
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
         "g_Intestinibacter" = '#679698',
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

# Load differential abundance analysis from ANCOM Differential Abundance script and transform into to log2
DA_rectal<-read.csv("Outputs/Ancombc2/ancombc2-primary-results-rectal.csv") %>%
  mutate(across(lfc_timeD3:lfc_timeD14, function(x) exp(x))) %>%
  mutate(across(lfc_timeD3:lfc_timeD14, function(x) log2(x)))

#Keep theme consistant across graphs
theme_custom <- list(theme_classic() + theme(text = element_text(family = "Arial", colour = 'black'),
                                             strip.background = element_blank(),
                                             strip.text = element_text(size = 14, face = "bold"),
                                             plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                                             axis.text = element_text(size = 14),
                                             axis.title = element_text(size = 14, family = "Arial"),
                                             aspect.ratio = 1))
#Shannon Diversity (5A)
shannon_rectal <- ggplot(data = shannon[which(shannon$location == "Rectal Swab"),], aes(x = time, y = shannon_entropy, group = time)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  facet_wrap(~location, nrow = 1, scales = "free",  axes = 'all') +
  xlab("Days post-SARS-CoV-2") +
  ylab("Shannon Index") +
  theme_custom
shannon_rectal
ggsave("Figures/Diversity/rectal_shannon.png", plot = shannon_rectal, scale = 1.5, width = 2.5, height = 2.5, unit = "in") 

shannon_stool <- ggplot(data = shannon[which(shannon$location == "Stool"),], aes(x = time, y = shannon_entropy, group = time)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  facet_wrap(~location, nrow = 1, scales = "free",  axes = 'all') +
  xlab("Days post-SARS-CoV-2") +
  ylab("") +
  theme_custom

shannon_stool
ggsave("Figures/Diversity/stool_shannon.png", plot = shannon_stool, scale = 1.5, width = 2.5, height = 2.5, unit = "in")

#Statistical tests for shannon diversity
friedman_stool <-
  shannon[which(shannon$location == "Stool"),] %>%
  friedman_test(shannon_entropy ~ time | patient)

dunn_stool <- 
  dunn_test(shannon[which(shannon$location == "Stool"),], 
            shannon_entropy ~ time, 
            p.adjust.method = "none")
dunn_stool

friedman_Rectal <-
  shannon[which(shannon$location == "Rectal Swab"),] %>%
  friedman_test(shannon_entropy ~ time | patient)

#Relative Abundance (5B)
taxa_long_genus_GI <-
  taxasums %>%
  as.data.frame() %>% #Makes it a dataframe
  rownames_to_column("Taxon") %>% #Extracts the row names into a new column
  melt(id.vars = "Taxon", variable.name = "SampleID", value.name = "Abundance", na.rm = TRUE) %>% #Makes the dataframe long
  merge(metadata, by.x = 'SampleID', all = TRUE) %>%#Adds in the metadata for each sample
  filter(Abundance != 0) %>% # Filters out 0s
  filter(location != 'control') %>% #Deletes the controls
  filter(location == 'Rectal Swab' | location == 'Stool') %>% #Selects GI samples
  group_by(SampleID) %>%
  mutate(Total_Abundance = sum(Abundance)) %>% #Calculates total abundance for each sample
  ungroup()  %>%
  mutate(Relative_Abundance = Abundance/Total_Abundance) %>% #Calculates relative abundance
  mutate(LTC = if_else(Relative_Abundance > 0.05, "", "Other")) %>% #Designates any microbe that has a relative abundance of less than 5% as "Other"
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";") %>%
  mutate(LTC = if_else(LTC == "Other", "Other", paste0("g_", str_trim(Genus)))) %>%
  mutate(LTC = if_else(LTC =="g_NA", paste0("f_", str_trim(Family)), LTC)) %>%
  mutate(LTC = if_else(LTC =="f_NA", "Other", LTC)) %>%
  mutate(LTC = if_else(LTC =="g_NA", "g_uncultured", LTC))

relative_abundance_genus_GI <- ggplot(taxa_long_genus_GI, aes(x = time,
                                                              y = Abundance, fill = LTC)) +
  geom_bar(stat = 'identity', position = "fill", width = 0.9) +
  facet_grid(location ~ patient) +
  labs(x = "Days post-SARS-CoV-2", y = " Relative Abundance (%)") +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00), 
                     labels = c("0", "25", "50", "75", "100")) +
  scale_fill_manual(values = col) +
  theme_custom +
  theme(legend.key.size = unit(12, "pt"), legend.position = "bottom", legend.text = element_markdown(size = 10)) +
  #guides(fill=guide_legend(title=""))
  guides(fill= "none")

relative_abundance_genus_GI
ggsave("Figures/Relative Abundance/relative_abundnace_GI_genus.png", plot = relative_abundance_genus_GI, scale =1.5, height = 2.5, width = 8.33, unit = "in")

#Volcano plot for differential abundance (5C)

#Clean up data names
DA_rectal_to_plot <- DA_rectal %>% 
  select(taxon, lfc_timeD3:lfc_timeD14, q_timeD3:q_timeD14) %>% 
  melt("taxon") %>%
  separate(variable, into = c("variable", "time"), sep = "_") %>%
  dcast(taxon + time ~ variable) %>%
  separate(taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = "(?<=[a-zA-Z])__") %>%
  mutate(across(Kingdom:Family,function(x) str_remove(x, pattern = "_.$"))) %>%
  rowwise() %>%
  mutate(LTC = if_else(!is.na(Genus),paste0("g_", Genus), if_else(!is.na(Family), paste0("f_",Family), if_else(!is.na(Order), Order, if_else(!is.na(Class), Class, if_else(!is.na(Phylum), Phylum, Kingdom))))))

DA_rectal_to_plot$time = str_replace(DA_rectal_to_plot$time, "timeD", "Day ")
DA_rectal_to_plot$time = fct_relevel(DA_rectal_to_plot$time, "Day 3", "Day 5", "Day 7", "Day 10", "Day 14")

p_rectal_volcano <- ggplot(data = DA_rectal_to_plot) +
  geom_point(aes(x = lfc, y = -log(q), color = time), size = 2.5) +
  geom_hline(yintercept=-log(0.05), linetype='dotted') +
  xlim(-10,10) +
  ylim(0, 15) +
  gghighlight(q < 0.05, use_direct_label = FALSE) +
  geom_label_repel(aes(x = lfc, y = -log(q), label = if_else(Genus == 'Succinivibrio' | Genus == 'Streptococcus', Genus, ""), family = 'Arial'), box.padding = 1, max.overlaps = Inf, size = 5)+
  #geom_label_repel(aes(x = lfc, y = -log(q), label = LTC, segment.size = 0.5, fontface = if_else(LTC == "g_Streptococcus" | LTC == "g_Succinivibrio", "bold", "plain")), force_pull = 0, size = 3, max.overlaps = Inf, box.padding = 1, force = 3)+
  labs(x = expression("log"[2]*" Fold Change"), y = expression("-log"[10]*"FDR"), title = "Rectal Swab") +
  scale_color_manual("Timepoint", values = c("Day 3" = '#CC3311',"Day 5" = '#0077BB',"Day 7" = '#009988','Day 10' = '#EE7733', 'Day 14' = '#DDCC77')) +
  #guides(color = guide_legend()) +
  guides(color= "none") +
  theme_custom

p_rectal_volcano
ggsave(file="Figures/Differential Abundance/DA_rectal_volcano.png", plot=p_rectal_volcano, scale = 1.5, width=5.82, height=4.82, unit = "in")
ggsave(file="Figures/Differential Abundance/DA_rectal_volcano.png", plot=p_rectal_volcano, scale = 2, width=4, height=4, unit = "in")

# Relative abundance of significant taxa (5D-E)
plot_abundance_rectal <-
  taxasums %>%
  as.data.frame() %>% #Makes it a dataframe
  rownames_to_column("Taxon") %>% #Extracts the row names into a new column
  melt(id.vars = "Taxon", variable.name = "SampleID", value.name = "Abundance", na.rm = TRUE) %>% #Makes the dataframe long
  merge(metadata, by.x = 'SampleID', all = TRUE) %>%#Adds in the metadata for each sample
  filter(location != 'control') %>% #Deletes the controls
  filter(location == 'Rectal Swab') %>%
  group_by(SampleID) %>%
  mutate(Total_Abundance = sum(Abundance)) %>%
  ungroup()  %>%
  mutate(Relative_Abundance = Abundance/Total_Abundance) %>%
  filter(Taxon %in% c('d__Bacteria; Firmicutes; Bacilli; Lactobacillales; Streptococcaceae; Streptococcus',
                      'd__Bacteria; Proteobacteria; Gammaproteobacteria; Aeromonadales; Succinivibrionaceae; Succinivibrio')) #Selects taxa of intrest

strep_rectal <-
  plot_abundance_rectal %>%
  filter(Taxon == 'd__Bacteria; Firmicutes; Bacilli; Lactobacillales; Streptococcaceae; Streptococcus') 

succini_rectal <-
  plot_abundance_rectal %>%
  filter(Taxon == 'd__Bacteria; Proteobacteria; Gammaproteobacteria; Aeromonadales; Succinivibrionaceae; Succinivibrio')

strep_plot_rectal <- ggplot(data = strep_rectal, aes(x = time, y = Relative_Abundance, group = patient)) +
  geom_point()+
  geom_line() +
  geom_bracket(
    xmin = c("-7", "-7", "-7"), xmax = c("5", "7", "10"), #Significance comes from ANCOMBC2 calculations
    y.position = 0.2, label = c("*"), label.size = 6, size = 0.5, vjust = 0.6,
    step.increase = 0.1
  ) +
  labs(x = "Days post-SARS-CoV-2", y =" Relative Abundance (%)", title = "Streptococcus") +
  scale_y_continuous(breaks = c(0.00, 0.1, 0.2, 0.3),
                     labels = c("0", "10", "20", "30"), 
                     limits = c(0.00, 0.3)) +
  theme_custom

strep_plot_rectal
ggsave("Figures/Differential Abundance/strep_plot_rectal.png", plot = strep_plot_rectal, scale = 1.5, width = 2.5, height = 2.5, unit = "in")


succini_plot_rectal <- ggplot(data = succini_rectal, aes(x = time, y = Relative_Abundance, group = patient)) +
  geom_point() +
  geom_line() +
  geom_bracket(
    xmin = c("-7", "-7", "-7", "-7", "-7"), xmax = c("3", "5", "7", "10", "14"),  #Significance comes from ANCOMBC2 calculations
    y.position = 0.2, label = c("****", "**", "***", "****", "****"), label.size = 6, size = 0.5, vjust = 0.6,
    step.increase = 0.1
  ) +
  labs(x = "Days post-SARS-CoV-2", y =" Relative Abundance (%)", title = "Succinivibrio") +
  scale_y_continuous(breaks = c(0.00, 0.1, 0.2, 0.3),
                     labels = c("0", "10", "20", "30"), 
                     limits = c(0.00, 0.3)) +
  theme_custom

succini_plot_rectal
ggsave("Figures/Differential Abundance/succini_plot_rectal.png", plot = succini_plot_rectal, scale = 1.5, width = 2.5, height = 2.5, unit = "in")
