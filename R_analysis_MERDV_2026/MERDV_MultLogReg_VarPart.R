# MERDV Palenque 2026
# Script 4: MERDV_MultLogReg_VarPart.R
# Description: Multivariable logistic regression and variance partitioning 
# Author: Wendolyn Estrada
# Last update: August 2026

################################################################################

#0. Load required libraries=====================================================
library(dplyr)

# MULTIVARIABLE LOGISTIC REGRESSION============================================

#Use the compleate table with metadata and positiviy status genarated 
#in MERDV_detection as input

#1. DEFINE MODELS__________________________________________________________________
#1.1 Full set model: Genus (Others vs Culex) + Site
set1_positiveRhab$Genus_binary <- ifelse(set1_positiveRhab$Genus == "Culex", "Culex", "Other")
set1_positiveRhab$Genus_binary <- factor(set1_positiveRhab$Genus_binary,levels = c("Other", "Culex"))
fullSet_model <- glm(MERDV_status ~ Genus_binary + Place, data = set1_positiveRhab, family = binomial)

# 1.2 Only Culex model
culex <- set1_positiveRhab %>% filter(Genus == "Culex")
culex_model <- glm(MERDV_status ~ Place, data = culex, family = binomial)

#1.3 Null model (intercept-only) - reference for R² (full dataset, n=80)
null_model <- glm(MERDV_status ~ 1,data = set1_positiveRhab,family = binomial)

#1.4 Null model  (only Culex dataset, n=52)
null_model_culex <- glm(MERDV_status ~ 1,  data = culex,  family = binomial)

# 2. GENERATE DATAFRAME FOR FULL SET MODEL______________________________________

full_coef <- coef(summary(fullSet_model))
full_OR <- exp(full_coef[, "Estimate"])
full_CI <- exp(confint.default(fullSet_model))
full_p <- full_coef[, "Pr(>|z|)"]

full_df <- data.frame(
  Predictor = c("Intercept", 
                "Genus (Culex vs Other)", 
                "Place (Cultural-center vs Conserved)", 
                "Place (Eco-Hotel vs Conserved)"),
  OR = round(full_OR, 2),
  CI_lower = round(full_CI[,1], 2),
  CI_upper = round(full_CI[,2], 2),
  P_value = round(full_p, 4)
)

#3. GENERATE DATAFRAME FOR CULEX MODEL_________________________________________

culex_coef <- coef(summary(culex_model))
culex_OR <- exp(culex_coef[, "Estimate"])
culex_CI <- exp(confint.default(culex_model))
culex_p <- culex_coef[, "Pr(>|z|)"]

culex_df <- data.frame(
  Predictor = c("Intercept", 
                "Place (Cultural-center vs Conserved)", 
                "Place (Eco-Hotel vs Conserved)"),
  OR = round(culex_OR, 2),
  CI_lower = round(culex_CI[,1], 2),
  CI_upper = round(culex_CI[,2], 2),
  P_value = round(culex_p, 4)
)

#4. MCFADDEN R SQUERE___________________________________________________________

#  Rsquere = 1 - (logLik(model) / logLik(nul_model))
r2_full   <- 1 - as.numeric(logLik(fullSet_model) / logLik(null_model))
r2_culex   <- 1 - as.numeric(logLik(culex_model) / logLik(null_model_culex))

print("McFaddeN pseudoR (Multivariable logistic regression):")
print(paste0("Full model: ", r2_full))
print(paste0("Culex model: ", r2_culex))

# 4. VARIANCE PARTITIONING =====================================================

# The partitioning is performed using the models that utilize ALL the data (80 pools)
# It does NOT include the *Culex* model because it uses a subset of the data
# The Rsquere of ful_model is decomposed into:
#   - Unique to the genus (what the genus alone contributes)
#   - Unique to the site (what the site alone contributes)
#   - Shared (what cannot be attributed to just one)

#1. ADITIONAL MODELS____________________________________________________________
#1.1 Genus-only model (without site)
only_genus_model <- glm(MERDV_status ~ Genus_binary, data = set1_positiveRhab, family = binomial)

#1.2 Site-only model (without genus)
only_site_model <- glm(MERDV_status ~ Place,data = set1_positiveRhab, family = binomial)

#2. MCFADDEN R SQUERE___________________________________________________________
r2_full_part <- 1 - as.numeric(logLik(fullSet_model) / logLik(null_model))
r2_only_genus <- 1 - as.numeric(logLik(only_genus_model) / logLik(null_model))
r2_only_site  <- 1 - as.numeric(logLik(only_site_model) / logLik(null_model))

#3 Variant partition____________________________________________________________

# Unique fraction explained only by genus (genus after removing site effect)
genus_part <- r2_full_part - r2_only_site

# Unique fraction explained only by sampling site (site after removing genus effect)
site_part <- r2_full_part - r2_only_genus

# Shared fraction, attributed to both
shared <- r2_only_genus + r2_only_site - r2_full_part

#4. VERIFICATION (must yield exactly r2_full_part)______________________________
sum(genus_part, site_part, shared)  
df_part <- data.frame(
  Component = c("Genus-specific",
                "Site-specific",
                "Shared (genus + site)",
                "Total explained"),
  R2_McFadden = round(c(genus_part, site_part, shared, r2_full_part), 4),
  Percentage = round(c(genus_part/r2_full_part*100,
                       site_part/r2_full_part*100,
                       shared/r2_full_part*100,
                       100), 1))
write_xlsx(df_part, "../OUTPUTS/MERDV_variancePartition.xlsx")