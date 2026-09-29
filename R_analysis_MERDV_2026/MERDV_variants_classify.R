# MERDV Palenque 2026
# Script 3: MERDV_variants_classify.R
# Description: Descriptive classification of genomic variants in 
#              MERDV consensus genomes predicted by iVAR variants
# Author: Wendolyn Estrada
# Last update: August 2026

# Set working directory (use here::here() for reproducibility)
# If using RStudio, ensure the project is open.
# Otherwise, set the path manually:
# setwd("/path/to/your/project")

################################################################################

#0. Load required libraries====================================================
library(readr)          # For reading TSV files
library(purrr)          # Iterates over a list or vector 
library(stringr)        # Detect and extract string patterns
library(dplyr)          # Work with DataFrames 
library(tidyr)          # Work with DataFrames 
library(writexl)        # Export excel files 

#1. DEFINITION OF FUNCTIONS====================================================

# Function to read iVar variants reports________________________________________
#Receives: Path to the directory containing individual iVar reports
#Returns: An unified table with all iVar calls
read_individual_table <- function(file_path){
  
  file_name <- basename(file_path) # Obtener nombre archivo
  sample_id <- str_extract( file_name, "(?<=DNA_?)\\d+(?=_S)" )# Get ID
  sample_id <- as.numeric(sample_id)
  df <- read_tsv( file_path, show_col_types = FALSE )  # Read TSV
  print(paste0("Sample:", sample_id))
  print(paste0("Total variants (unfiltered): ", nrow(df)))
  df <- df %>% mutate(ID = sample_id, .before = 1) # Add ID as first column
  
  return(df)
}

# Function to classify variants_________________________________________________
#Receives: DataFrame of variants (unified iVar reports)
#          Threshold used for consensus building
#Return: An updated DataFrame filtered by allelic frecuency with classification of variants 

classify_variants <- function (df, t){
  filt_df <- df %>% filter(ALT_FREQ >= t) # Filter by allelic frecuency
  print(paste0("Variants over threshold: ", nrow(filt_df)))
  filt_df <- filt_df %>% mutate(Mutation_Type = case_when(
    str_detect(ALT, "^-") ~ "Deletion",
    str_detect(ALT, "^\\+") ~ "Insertion",
    is.na(GFF_FEATURE) ~ "Intergenic",
    REF_AA == ALT_AA ~ "Synonymous",
    REF_AA != ALT_AA ~ "Nonsynonymous",
    TRUE ~ "Undetermined")) 
  return(filt_df)
}

#Function to calculate depth statistis__________________________________________
#Receives: DataFrame with filtered and classified variants
#Returns: Statistics by genome
calculate_depth_stats <- function(filt_df) {
  stats <- filt_df %>%
    group_by(ID) %>%
    summarise(
      Mean_depth = mean(TOTAL_DP, na.rm = TRUE),
      Median_depth = median(TOTAL_DP, na.rm = TRUE),
      Pct_ge_10x = mean(TOTAL_DP >= 10, na.rm = TRUE) * 100,
      Pct_ge_20x = mean(TOTAL_DP >= 20, na.rm = TRUE) * 100,
      N_variant_positions = n())
  return(stats)
}

#2. MAIN========================================================================

#1. READ AND CONCATENATE IVAR REPORTS___________________________________________
path_files <- "../DATA/MERDV_iVArVariants"
files <- list.files(path = path_files, pattern = "_TRUE.tsv$",full.names = TRUE)
files
DFivar <- map_dfr(files, read_individual_table)

#2. ADD METADATA________________________________________________________________
DFivar <- DFivar %>% left_join(METADATA , by = "ID") %>%
  select(-Pool,-Species,-Set) %>% relocate(Place, Genus, Num_mosquitoes, .after = ID)

#3. ANNOTATION WITH REFERENCE GFF FILE__________________________________________
# Read gff
gff <- read_tsv( "../DATA/MERDV_NC_040599_1.gff3", comment = "#", col_names = FALSE, show_col_types = FALSE)
# Rename columns
colnames(gff) <- c("seqid","source","type","start","end", "score","strand","phase","attributes")
# Filter CDS
gff_cds <- gff %>% filter(type == "CDS")
# Get GeneID y product
gff_cds <- gff_cds %>% mutate(
  GFF_FEATURE = str_extract(attributes, "(?<=ID=)[^;]+"),
  GeneID = str_extract(attributes,"(?<=GeneID:)\\d+"),
  REF_CDS = paste(start, end, sep = "-"),
  Protein = str_extract(attributes, "(?<=product=)[^;]+" ),
  ProteinID = str_extract(attributes, "(?<=protein_id=)[^;]+")) %>%
  select(GFF_FEATURE,GeneID,REF_CDS, Protein, ProteinID,)

#Add annotation 
DFivar <- DFivar %>% left_join(gff_cds, by = "GFF_FEATURE")
print(dim(DFivar))

#3. CLASSIFY VARIANTS___________________________________________________________
filt_DFivar <- classify_variants(DFivar, 0.5)
print(dim(filt_DFivar))

#4 CALCULATE DEPHT STATISTICS_________________________________________________
depth_stats <- calculate_depth_stats(filt_DFivar)

#7.SUMMARY OF VARIANTS PER GENOME_______________________________________________
#Group by distic position in reference genome and alternative base
variants_per_pool <- filt_DFivar %>% 
  distinct(ID, POS, REF, ALT, .keep_all = TRUE) %>%
  group_by(ID, Protein) %>% 
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Protein, values_from = Count, values_fill = 0)
#Add depth statistics
variants_per_pool <- variants_per_pool %>% left_join(depth_stats, by = "ID")
#Save excel
write_xlsx(variants_per_pool, "../OUTPUTS/supplementary_table_S4.xlsx")  

#8. DIFFERENT VARIANTS ACROSS MERDV GENOMES (ALT_FREQ > 0.5)____________________
diff_variants <- filt_DFivar %>% group_by(POS, REF, ALT,) %>%
  summarise(
    n = n_distinct(ID),
    REF_CODON = dplyr::first(REF_CODON),
    ALT_CODON = dplyr::first(ALT_CODON),
    REF_AA = dplyr::first(REF_AA),
    ALT_AA = dplyr::first(ALT_AA),
    TYPE = dplyr::first(Mutation_Type),
    PROTEIN = dplyr::first(Protein),
    ALT_FREQ = mean(ALT_FREQ, na.rm = TRUE),  # promedio
    .groups = "drop"
  ) %>% arrange(desc(n))

print(table(diff_variants$TYPE))
print(table(diff_variants$PROTEIN))

#Split the 'Product' column into separate columns and indicate
#whether or not the mutation falls within that gene region.
diff_variants_wide <- diff_variants %>%  mutate(
  PROTEIN = ifelse(is.na(PROTEIN)| PROTEIN == "NA","Non-coding region", PROTEIN)) %>%    
  mutate(Count = 1) %>% #Conteo, 1 si la mutacion pertenece a esa trgion
  pivot_wider(names_from = PROTEIN, values_from = Count, #pasar a formato ancho
              values_fill = 0) #Si la mutacion no corresponde a ese producto
#Save excel
write_xlsx(diff_variants_wide,  "../OUTPUTS/supplementary_table_S5.xlsx")  

#9. SHARED VARIANTS SUMMARY_____________________________________________________
summary_sharedV <- diff_variants %>% mutate(
  Proportion = n / 18,
  Category = case_when( #Calcular numero de mutaciones segun la proporcion de pools
    n == 18 ~ "All genomes (n = 18)",
    Proportion >= 0.75 ~ "75% of genomes (n ≥ 14)",
    Proportion >= 0.50 ~ "50% of genomes (n ≥ 9)",
    Proportion >= 0.25 ~ "25% of genomes (n ≥ 5)",
    n == 1 ~ "One genome",
    TRUE ~ "2-4 genomes",
  )) %>% group_by(PROTEIN, Category) %>% summarise(Count = n(), .groups = 'drop')  %>%
  pivot_wider(names_from = PROTEIN, values_from = Count, values_fill = 0) #Desglosar por gen

#CALCULATE VARIANTS PER KB BY GENE AND TYPE_____________________________________
#Add reference gene length and normalization factor
diff_variants <- diff_variants %>% mutate(
  nuc_length = case_when(
    PROTEIN == "nucleoprotein"~ 1437,
    PROTEIN == "phosphoprotein"~ 1203,
    PROTEIN == "matrix protein"~ 564,
    PROTEIN == "glycoprotein"~ 1530,
    PROTEIN == "RNA-dependent RNA polymerase"~ 6411,
    TRUE ~ NA_real_
  ),
  norm_factor = 1000/nuc_length
)

#Symmary by gene
summary_ByGene <- diff_variants %>% group_by(PROTEIN,TYPE) %>%
  summarise(Count = n(),                           
            nuc_length = max(nuc_length),   # Como es constante por gen, max da el mismo valor
            norm_factor = max(norm_factor),
            .groups = 'drop') %>% 
  pivot_wider(names_from = TYPE,values_from = Count,values_fill =  0) 

#Normalization
summary_ByGene<-summary_ByGene %>% mutate(
  NoSyns_perkb = ifelse(!is.na(nuc_length), Nonsynonymous * norm_factor, NA_real_),
  Syns_perkb = ifelse(!is.na(nuc_length), Synonymous * norm_factor, NA_real_))

#NOTE: Intergenic, Deletion, and Intergenic were not normalized; they are located in non-coding regions.
#Guardar tabla como excel
write_xlsx(summary_ByGene, "../DATA/MERDV_summaryByGene_normalized.xlsx")
