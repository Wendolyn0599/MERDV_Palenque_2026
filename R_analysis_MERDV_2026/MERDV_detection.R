# MERDV Palenque 2026
# Script 1: MERDV_detection.R
# Description: Processing of MERDV detection data from BLAST and Bowtie2 
#              results. Generation of positivity tables and preparation for 
#              statistical analyses.
# Author: Wendolyn Estrada
# Last update: August 2026

# Set working directory (use here::here() for reproducibility)
# If using RStudio, ensure the project is open.
# Otherwise, set the path manually:
# setwd("/path/to/your/project")

################################################################################

# 0. Load required libraries====================================================
library(dplyr)      # Data manipulation
library(purrr)      # Functional programming (map_chr)
library(stringr)    # String manipulation
library(rentrez)    # NCBI Entrez API for taxonomic annotation
library(ggplot2)    # Data visualization
library(patchwork)  # Combining multiple plots
library(scales)     # Plot formatting and scales
library(writexl)    # Export files

#1. LOAD DATA===================================================================
getwd()#Veryfy project path
# Read Bowtie2 positivity table against Rhabdoviridae.
#ID: pool identifier; Total_reads, Reads_rhabdo: reads that apped to contings, Contigs:count.
#It was generated from  MapBBEnsamSpadeDupMagBwt_list.sh outputs in the Bash terminal
bowtie_results <- read.table("../DATA/positiveSEQ_BowtieRhabdo.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#Read BLASTn results from BLASTN_contings.sh 
blast_results <- read.table("../DATA/rhabdo_blastn_unificada.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE) 

#Read contings results 
contings_results <- read.table("../DATA/contings_magnitude.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) 

#Read sample metadata
METADATA <- read.table("../DATA/METADATA.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Clean genus: extract only the genus (remove subgenus)
METADATA$Genus <- sub(" .*", "", METADATA$Genus)

#2. CONTIG VALIDATION WITH BLASTN ==============================================

# Clean sample IDs (remove "DNA_" and "_S*" patterns)
blast_results$sample <- sub("^DNA_?", "", blast_results$sample)
dim(blast_results) #Size of original table

# Taxonomic anotation
get_name <- function(acc) {
  tryCatch({
    res <- entrez_summary(db = "nuccore", id = acc)
    return(res$organism)
  }, error = function(e) {
    return(NA)
  })
}

blast_results$Hit_name <- map_chr(blast_results$sseqid, get_name)

# Filter blast_results to keep only validated contigs
blast_results_filt <- blast_results %>%filter(Hit_name == "Merida virus")
dim(blast_results_filt) #Size after filtering

# Join with contig length and abundance (magnitude)
summary_table <- blast_results_filt %>%
  left_join(contings_results %>% select(qseqid,Conting_length, Magnitude), by = "qseqid")

# Add metadata
summary_table$sample <- as.numeric(summary_table$sample)
summary_table <- summary_table %>%
  left_join(METADATA , by = c("sample" = "ID"))

# Select and order columns
summary_table <- summary_table %>%
  select(sample, Place, Genus, Species, Num_mosquitoes, Set, qseqid, Conting_length, Magnitude, sseqid, pident, 
         length, qcovs, evalue, bitscore, Hit_name)

#4. IDENTIFICATION OF MERDV-POSITIVE POOLS =====================================
# Get IDs of validated positive samples
validIDs <- unique(summary_table$sample)
# Filter Bowtie results to keep only validated positive samples
positiveSEQ <- bowtie_results[bowtie_results$ID %in% validIDs, ]
# Calculate reads per million (RPM)
positiveSEQ <- positiveSEQ %>% mutate(
    RPM = (Reads_rhabdo / Total_reads) * 1e6,
    log_RPM = log10(RPM + 1))

#5. CREATE POSITIVITY TABLE WITH METADATA======================================= 
# Create complete table with all pools, number of reads and their positivity status
set1_positiveRhab <- METADATA %>% filter(Set == 1 ) %>%
  left_join(bowtie_results %>% select(ID, Total_reads,Contings), by ="ID") %>%
  left_join(positiveSEQ %>% select(ID, Reads_rhabdo, RPM, log_RPM), by = "ID") %>%
  mutate(MERDV_status = ifelse(ID %in% positiveSEQ$ID, 1, 0))

#Replace NA by 0 innegative pools
set1_positiveRhab <- set1_positiveRhab %>% mutate(
  Reads_rhabdo = ifelse(is.na(Reads_rhabdo), 0, Reads_rhabdo),
  RPM = ifelse(is.na(RPM), 0, RPM),
  log_RPM = ifelse(is.na(log_RPM), 0, log_RPM),
  log_TotalReads = log10(Total_reads),
  # Convert categorical variables to factors
  Place = as.factor(Place),
  Genus = as.factor(Genus),
  MERDV_status = as.factor(MERDV_status)) %>%
  select( # Select and order columns
    ID, Place, Genus, Num_mosquitoes,
    Total_reads, log_TotalReads, MERDV_status,
    Reads_rhabdo,Contings, RPM, log_RPM)
#Save as excel
write_xlsx(set1_positiveRhab,"../OUTPUTS/supplementary_table_S1.xlsx" )

#6. DESCRIPTIVE STATISTICS======================================================

desc_stat <- set1_positiveRhab %>%
  group_by(MERDV_status) %>%
  summarise(
    Mediana = median(Num_mosquitoes, na.rm = TRUE),
    IQR = IQR(Num_mosquitoes, na.rm = TRUE),
    Q1 = quantile(Num_mosquitoes, 0.25, na.rm = TRUE),   
    Q3 = quantile(Num_mosquitoes, 0.75, na.rm = TRUE),   
    .groups = "drop")

# Positivity frequency by habitat
freq_habitat <- set1_positiveRhab %>% group_by(Place) %>%
  summarise(
    Negative_pools = sum(MERDV_status == 0),
    Positive_pools = sum(MERDV_status == 1),
    Total = n(),
    Percentage = round((Positive_pools / Total) * 100, 2))

# Positivity frequency by genus
freq_host <- set1_positiveRhab %>% group_by(Genus) %>%
  summarise(
    Negative_pools = sum(MERDV_status == 0),
    Positive_pools = sum(MERDV_status == 1),
    Total = n(),
    Percentage = round((Positive_pools / Total) * 100, 2))

# 7.  CONFIDENCE INTERVALS======================================================
# For sampling site
freq_habitat <- freq_habitat %>% rowwise() %>% mutate(
  CI = paste0(round(binom.test(Positive_pools, Total)$conf.int[1] * 100, 2), "-", 
              round(binom.test(Positive_pools, Total)$conf.int[2] * 100, 2)))
# For host genus
freq_host <- freq_host %>% rowwise() %>% mutate(
  CI = paste0(round(binom.test(Positive_pools, Total)$conf.int[1] * 100, 2), "-", 
              round(binom.test(Positive_pools, Total)$conf.int[2] * 100, 2)))

#8. FISHER TESTS================================================================
#For sampling site
place <- table(set1_positiveRhab$Place, set1_positiveRhab$MERDV_status)
fisher_hab <- fisher.test(place)
print(fisher_hab)

#For host genus
genus <- table(set1_positiveRhab$Genus, set1_positiveRhab$MERDV_status)
fisher_genus <- fisher.test(genus)
print(fisher_genus)

#9. STATISTICAL TESTS ==========================================================

# Kruskal-Wallis to compare log_TotalReads 
kruskal_totalReads <- kruskal.test(log_TotalReads ~ Place, data = set1_positiveRhab)
kruskal_totalReads

# Kruskal-Wallis to compare pool size between positive pools
kruskal_numMosco <- kruskal.test(Num_mosquitoes~ Place, data = set1_positiveRhab)
kruskal_numMosco

# Wilcoxon log_TotalReads reads by positivity status
wilcoxon_TotaReads_PosPools <- wilcox.test(log_TotalReads ~ MERDV_status, data = set1_positiveRhab)
wilcoxon_TotaReads_PosPools

# Convert MERDV_status to factor for Wilcoxon test
set1_positiveRhab <- set1_positiveRhab %>% mutate(MERDV_status = as.factor(MERDV_status))
# Wilcoxon relacion Totalres vs pools positivos
wilcoxon_moscos_PosPools <- wilcox.test(Num_mosquitoes ~ MERDV_status, data = set1_positiveRhab)
wilcoxon_moscos_PosPools

#Kruskal-Wallis abundance (RPM) by sampling site
positiveSEQ <- positiveSEQ %>% left_join(METADATA, by = "ID")
kruskal_hab <- kruskal.test(RPM ~ Place, data = positiveSEQ)
kruskal_hab
#Spearman correlation to evaluate RPM vs total reads
spe_cor_rpm <- cor.test(~ RPM + Total_reads, data = positiveSEQ, method = "spearman")
spe_cor_rpm
#10 VISUALIZATION SETTINGS======================================================
# Custom theme
my_theme <- theme_minimal() +
  theme(
    axis.title = element_text(size = 12, family = "sans"),
    axis.text = element_text(size = 12, family = "sans"),
    legend.text = element_text(size = 12, family = "sans"),
    legend.title = element_text(size = 12, family = "sans"),
    strip.text = element_text(size = 12, family = "sans"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.tag = element_text(size = 12, face = "bold")
  )

# Color palette
my_colors <- c( "#E63946", "#FF9F1C","#0057B8","#2E294E","#00A878","#9D4EDD", "#C5A000", "#FF80FF"  )
# Jitter settings
my_jitter <- geom_jitter(aes(color = Genus), width = 0.2, height = 0.1, size = 2.5 )

# Function to italicize genus labels
italic_labels <- function(x) {
  sapply(x, function(g) parse(text = paste0("italic('", g, "')")))
}

#14. PLOTING====================================================================
# Figure 1A: Distribution of total reads by site
  fig_1a <- ggplot(set1_positiveRhab, aes(x=Place, y=log_TotalReads)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  my_jitter +
  facet_wrap(~Place,scales = "free_x", space = "free_x", labeller = labeller( Place = c( "Conserved-area" = "Conserved\narea",
                                                                                         "Cultural-center" = "Cultural\ncenter","Eco-Hotel" = "Eco-\nHotel")))+
  scale_color_manual(values = my_colors, labels = italic_labels) + 
  scale_x_discrete(labels = c( "Cultural-center" = "Cultural\ncenter","Eco-Hotel" = "Eco-\nHotel","Conserved-area" = "Conserved\narea" )) +
  labs( y = "Log10 (total reads)",x="", color = "Host genus") +
  coord_cartesian(ylim = c(6,8))+
  my_theme

# Figure 1B: Distribution of pool size by site
fig_1b <- ggplot(set1_positiveRhab, aes(x=Place, y=Num_mosquitoes)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA)+
  my_jitter +
  facet_wrap(~Place,scales = "free_x", space = "free_x", labeller = labeller( Place = c( "Conserved-area" = "Conserved\narea",
                                                                                         "Cultural-center" = "Cultural\ncenter","Eco-Hotel" = "Eco-\nHotel")))+
  scale_color_manual(values = my_colors, labels = italic_labels) + 
  scale_x_discrete(labels = c( "Cultural-center" = "Cultural\ncenter","Eco-Hotel" = "Eco-\nHotel","Conserved-area" = "Conserved\narea" )) +
  labs(y = "Pool size",x="",color = "Host genus") +
  coord_cartesian(ylim = c(0, 60))+
  my_theme 

# Figure 1C: Total reads vs detection by site
fig_1c <- ggplot(set1_positiveRhab, aes( x=MERDV_status, y=log_TotalReads)) +
  geom_boxplot(aes(fill = MERDV_status), alpha = 0.3,outlier.shape = NA)+
  my_jitter + 
  facet_wrap(~Place, labeller = labeller( Place = c( "Conserved-area" = "Conserved\narea",
                                                     "Cultural-center" = "Cultural\ncenter","Eco-Hotel" = "Eco-\nHotel")))+
  scale_x_discrete( labels = c("0" = "Negative", "1" = "Positive")) +
  scale_color_manual(values = my_colors, labels = italic_labels) + 
  scale_fill_manual(values = c("0" = "white", "1" = "blue"),labels = c( "0" = "Negative", "1" = "Positive"),name = "MERDV detection status") +
  labs(x="", y= "Log10(total reads)")+
  coord_cartesian(ylim = c(6, 8))+
  my_theme 

# Figure 1D: Pool size vs detection by site
fig_1d<- ggplot(set1_positiveRhab, aes( x=MERDV_status, y=Num_mosquitoes)) +
  geom_boxplot(aes(fill = MERDV_status) , alpha = 0.3,outlier.shape = NA)+
  my_jitter + 
  facet_wrap(~Place, labeller = labeller( Place = c( "Conserved-area" = "Conserved\narea",
                                                     "Cultural-center" = "Cultural\ncenter","Eco-Hotel" = "Eco-\nHotel")))+
  scale_x_discrete( labels = c("0" = "Negative", "1" = "Positive")) +
  scale_color_manual(values = my_colors, labels = italic_labels) + 
  scale_fill_manual(values = c("0" = "white", "1" = "blue"),labels = c( "0" = "Negative", "1" = "Positive"),name = "MERDV detection status") +
  labs(x="", y= "Pool size")+
  coord_cartesian(ylim = c(0, 60))+
  my_theme

#Save figure 1
FIG1 <- (fig_1a| fig_1b)/(fig_1c | fig_1d)+
  plot_annotation(tag_levels = "a")  
ggsave(
  "../OUTPUTS/FIG1.png",
  FIG1,
  width = 3300,
  height = 2500,
  units = "px",
  dpi = 300,
)

# Figure 2A: Stacked bar plot of detection by genus and site
fig_2a <- ggplot(set1_positiveRhab, aes(x=MERDV_status,  fill = Genus)) +
  geom_bar(position = "stack") +
  facet_wrap(~Place, labeller = labeller( Place = c( "Conserved-area" = "Conserved\narea",
                                                     "Cultural-center" = "Cultural\ncenter","Eco-Hotel" = "Eco-\nHotel")))+
  scale_x_discrete( labels = c("0" = "Negative", "1" = "Positive")) +
  scale_fill_manual(values = my_colors, name = "Host genus",labels = italic_labels) +
  labs(x = "MERDV detection status", y = "Number of pools") +
  geom_text(stat = "count",  aes(label = after_stat(count)), position = position_stack(vjust = 0.5),size = 3.2,color = "white",) +
  my_theme 

# Figure 2B: Pool size by genus (positive samples only)
fig_2b <- ggplot(positiveSEQ, aes(x = Genus, y = Num_mosquitoes, color = Genus)) +
  geom_boxplot(width = 1,outlier.shape = NA)+
  my_jitter +
  facet_wrap(~Place,scales = "free_x", space = "free_x", labeller = labeller( Place = c( "Conserved-area" = "Conserved\narea",
                                                                                         "Cultural-center" = "Cultural\ncenter","Eco-Hotel" = "Eco-\nHotel")))+
  scale_color_manual(values = c("#E63946","#0057B8"), labels = italic_labels) + 
  scale_x_discrete(labels = italic_labels) +
  labs( y = "Pool size",x= "", color = "Host genus") + 
  coord_cartesian(ylim = c(0, 60))+
  my_theme

# Figure 2C: Abundance by site
fig_2c <- ggplot(positiveSEQ, aes(x = Place, y = log_RPM, color = Genus)) +
  geom_violin(trim = TRUE)+
  geom_boxplot(width = 0.07,outlier.shape = NA)+
  my_jitter +
  facet_wrap(~Place,scales = "free_x", space = "free_x", labeller = labeller( Place = c( "Conserved-area" = "Conserved\narea",
                                                                                         "Cultural-center" = "Cultural\ncenter","Eco-Hotel" = "Eco-\nHotel")))+
  scale_x_discrete(labels = c( "Cultural-center" = "Cultural\ncenter","Eco-Hotel" = "Eco-\nHotel","Conserved-area" = "Conserved\narea" )) +
  scale_color_manual(values = c("#E63946","#0057B8"), labels = italic_labels) +
  coord_cartesian(ylim = c(0, 6))+
  labs(x = "", y = "Log10 (RPM)", color = "Host genus")+
  my_theme


# Figure 2D: Abundance vs sequencing depth
fig_2d <-ggplot(positiveSEQ, aes(x = log10(Total_reads) , y = log_RPM, color = Genus)) +
  geom_point() +
  facet_wrap(~Place, labeller = labeller( Place = c( "Conserved-area" = "Conserved\narea", 
                                                     "Cultural-center" = "Cultural\ncenter","Eco-Hotel" = "Eco-\nHotel")))+
  scale_color_manual(values = c("#E63946","#0057B8"), labels = italic_labels) +
  labs(x = "Log10 (Total reads)", y = "Log10 (RPM)", color = "Host genus") + 
  coord_cartesian(xlim = c(6.5, 7.5),ylim = c(0, 6))+
  my_theme

#Save FIG 2 
#Save figure 1
FIG2 <- (fig_2a| fig_2b)/(fig_2c | fig_2d)+
  plot_annotation(tag_levels = "a")  
ggsave(
  "../OUTPUTS/FIG2.png",
  FIG1,
  width = 3300,
  height = 2500,
  units = "px",
  dpi = 300,
)
