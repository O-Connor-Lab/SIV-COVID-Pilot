#Load in Packages
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
metadata$location <- fct_relevel(metadata$location, "Nasal Swab", "Tracheal Swab", "Rectal Swab", "Stool") #Relevels location so they are in a specific order
metadata$time <- fct_relevel(metadata$time, "-7", "3", "5", "7", "10", "14") #Relevels time so it is in chronological order


#Reads in shannon diversity from QIIME2
shannon <- read_qza("Outputs/Diversity/rarefied-shannon.qza")$data %>%
  rownames_to_column(var='SampleID') %>%
  inner_join(metadata) %>%
  filter(location != 'control')
write.csv(shannon,"Figures/Supplemental/Shannon-Diversity.csv")

#Read in taxonomy assigned in QIIME2
taxonomy<-read_qza("Data/taxonomy_2.qza") #read in taxonomy
write.csv(taxonomy$data,"Figures/Supplemental/Taxonomy.csv")
taxonomy<-parse_taxonomy(taxonomy$data) #Removes the confidence from the taxonomy and breaks up the taxonomy

#Save feature table with taxonomical information (Genus Level)
taxasums <- summarize_taxa(FeatureTable, taxonomy)$Genus
write.csv(taxasums, "Figures/Supplemental/Feature_Table.csv", row.names = TRUE)

#Assign consistant colors with bacteria of intrest
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

#Assigns theme for all graphs
theme_custom <- list(theme_classic() + theme(text = element_text(family = "Arial", colour = 'black'),
                                             strip.background = element_blank(),
                                             strip.text = element_text(size = 14, face = "bold"),
                                             plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                                             axis.text = element_text(size = 14),
                                             axis.title = element_text(size = 14, family = "Arial"),
                                             aspect.ratio = 1))

#Shannon diversity (4A)
shannon_nasal <- ggplot(data = shannon[which(shannon$location == "Nasal Swab"),], aes(x = time, y = shannon_entropy, group = time)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  facet_wrap(~location, nrow = 1, scales = "free",  axes = 'all') +
  xlab("Days post-SARS-CoV-2") +
  ylab("Shannon Index") +
  theme_custom

shannon_nasal
ggsave("Figures/Diversity/nasal_shannon.png", plot = shannon_nasal, scale = 1.5, width = 2.5, height = 2.5, unit = "in") 

shannon_Tracheal <- ggplot(data = shannon[which(shannon$location == "Tracheal Swab"),], aes(x = time, y = shannon_entropy, group = time)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  geom_bracket(
    xmin = "-7", xmax = "5", y.position = 7.65,   #Label taken from Dunn's test results
    label = "*", label.size = 6, size = 0.5, vjust = 0.6) + 
  facet_wrap(~location, nrow = 1, scales = "free",  axes = 'all') +
  xlab("Days post-SARS-CoV-2") +
  ylab("Shannon Index") +
  ylim(2, 8) +
  theme_custom +
  theme(axis.title.y = element_blank())

#Statistical tests for Shannon Diversity
friedman_Nasal <-
  shannon[which(shannon$location == "Nasal Swab"),] %>%
  friedman_test(shannon_entropy ~ time | patient)

friedman_Tracheal <-
  shannon[which(shannon$location == "Tracheal Swab"),] %>%
  friedman_test(shannon_entropy ~ time | patient)

dunn_Tracheal <- 
  dunn_test(shannon[which(shannon$location == "Tracheal Swab"),], 
            shannon_entropy ~ time, 
            p.adjust.method = "none")
dunn_Tracheal

shannon_Tracheal
ggsave("Figures/Diversity/Tracheal_shannon.png", plot = shannon_Tracheal, scale = 1.5, width = 2.5, height = 2.5, unit = "in") 

#Relative Abundance (4B)
taxa_long_genus_RS <-
  taxasums %>%
  as.data.frame() %>% #Makes it a dataframe
  rownames_to_column("Taxon") %>% #Extracts the row names into a new column
  melt(id.vars = "Taxon", variable.name = "SampleID", value.name = "Abundance", na.rm = TRUE) %>% #Makes the dataframe long
  merge(metadata, by.x = 'SampleID', all = TRUE) %>%#Adds in the metadata for each sample
  filter(Abundance != 0) %>% # Filters out 0 values
  filter(location != 'control') %>% #Filters the controls
  filter(location == 'Nasal Swab' | location == 'Tracheal Swab') %>% #Filters for respiratory samples
  group_by(SampleID) %>%
  mutate(Total_Abundance = sum(Abundance)) %>% #Find total abundance of each sample
  ungroup()  %>%
  mutate(Relative_Abundance = Abundance/Total_Abundance) %>% #Calculate relative abundance
  mutate(LTC = if_else(Relative_Abundance > 0.1, "", "Other")) %>% #Designate any Genus that represents less than 10% of a sample as other
  #Clean names
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = ";") %>%
  mutate(LTC = if_else(LTC == "Other", "Other", paste0("g_", str_trim(Genus)))) %>%
  mutate(LTC = if_else(LTC =="g_NA", paste0("f_", str_trim(Family)), LTC)) %>%
  mutate(LTC = if_else(LTC =="g_NA", "g_uncultured", LTC))

relative_abundance_genus_RS <- ggplot(taxa_long_genus_RS, aes(x = time,
                                                              y = Abundance, fill = LTC)) +
  geom_bar(stat = 'identity', position = "fill") +
  facet_grid(location ~ patient) +
  labs(x = "Days post-SARS-CoV-2", y = " Relative Abundance (%)") +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 1.00), 
                     labels = c("0", "25", "50", "75", "100")) +
  scale_fill_manual(values = col) +
  theme_custom +
  theme(legend.key.size = unit(10, "pt"), legend.position = "bottom", legend.text = element_markdown(size = 8), axis.title = element_text(20)) +
  #guides(fill=guide_legend(title=""))
  guides(fill = "none")


relative_abundance_genus_RS 
ggsave("Figures/Relative Abundance/relative_abundnace_RS_genus.png", plot = relative_abundance_genus_RS, scale =1.5, height = 2.5, width = 8.33, unit = "in")
