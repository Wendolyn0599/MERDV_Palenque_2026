# MERDV Palenque 2026
# Script 5: MERDV_MLE_rates.R
# Description: Maximum Likelihood Estimation (MLE) of infection rates 
#              by site and genus group (Culex, Aedes, Other)
# Author: Wendolyn Estrada
# Last update: September 2026

################################################################################

# 0. Load required libraries====================================================
library(dplyr)
library(writexl)
#===============================================================================

# This script uses the positivity table from MERDV_detection.R:
#  set1_positiveRhab, it must be charged in the environment

# 1. DEFINE MLE FUNCTION========================================================
#The function calculates Maximum Likelihood Estimation per 1000 mosquitoes.
#-Receives: pool_data; A column from the "summary_data" table, and correspond to a list of tables
#          that indicates the size (Num_mosquitoes) and the MERDV detection status for each 
#          individual pool within the evaluated group (MERDV_status: 1 if positive, 0 if negative).
#-Return: MLE/1000, the infection rate per 1000 mosquitoes; CI_lower and CI_upper 

#The log-likelihood for a set of pools with varying sizes is given by the formula:
#   l(p) = Σ_i[ x_i * log(1 - (1-p)^m_i)] + log(1-p) * Σ_i[m_i(n_i - x_i)]
# Where:
# - l(p) = log-likelihood, given p
# - p = probability of infection
# - m_i = pool size
# - x_i = number of positive pools of size m_i
# - n_i = total number of pools of size m_i

# The MLE (p_hat) is the value of p that maximizes l(p).
# 95% confidence intervals are obtained by likelihood profiling, solving:
#    -2*(l(p) - l(p_hat)) = qchisq(0.95, 1) for p.

calc_mle <- function(pool_data) {
  #Get individual pool information and convert to dataFrame
  df <- as.data.frame(pool_data)
  # Group by pool size to obtain m_i, x_i, n_i
  df_sum <- df %>% group_by(Num_mosquitoes) %>%
    summarise( x = sum(MERDV_status),           # positive pools of a certain size
               n = n(),.groups = "drop" ) %>%   # The total pools of that size
    arrange(Num_mosquitoes)
  
  m <- df_sum$Num_mosquitoes   # vector with pool sizes of each pool (m_i)
  x <- df_sum$x                # vector of positive pools by size (x_i)
  n <- df_sum$n                # vector of total pools by size (n_i)
  total_pools <- sum(n)        #Total pools
  positive_pools <- sum(x)     #Positive pools
  
  # Case 1: Non positive pools
  if (positive_pools == 0) {
    return(data.frame(MLE = 0, CI_lower = 0, CI_upper = 0))
  } else if (positive_pools == total_pools) { 
    # Case 2: All pools are positive
    return(data.frame(MLE = 1000, CI_lower = 0, CI_upper = 1000))
  } else {
    #Case 3: General:
    
    #Internal helper function log_lik(p): calculates the log-likelihood for a given value of p. 
    # -Since these functions are defined within calc_mle, it has access to the values of m, x, and n 
    # without needing to pass them as extra arguments. 
    # -This functions take p as input and when is called by optimize(),
    #  it test different values of p within the interval and return the final value.
    # As p is a variable parameter that change while is tested until optimized, it not defined before. 
    log_lik <- function(p) {
      if (p <= 0 || p >= 1) return(-Inf)
      log1_p <- log1p(-p)  # log(1-p)
      sum(x * log1p(-exp(m * log1_p)) + (n - x) * m * log1_p)}
    
    #1. Find the MLE by minimizing -log_lik (equivalent to maximizing log_lik)
    opt <- optimize(function(p) -log_lik(p), interval = c(1e-6, 1 - 1e-6))
    p_hat <- opt$minimum
    max_logL <- log_lik(p_hat)
    
    #2. Calculate confidence intervals by likelihood profiling 
    q <- qchisq(0.95, 1)
    # Target function: -2*(l(p) - l(p_hat)) - q = 0   
    target <- function(p) {
      if (p <= 0 || p >= 1) return(1e10)
      -2 * (log_lik(p) - max_logL) - q}
    
    # Lower limit: search in (0, p_hat)
    lower <- tryCatch({
      uniroot(target, interval = c(1e-6, p_hat), extendInt = "downX")$root
    }, error = function(e) {
      # If uniroot fails, fallback to grid search
      grid <- seq(1e-6, p_hat, length.out = 1000)
      idx <- which.min(abs(target(grid)))
      grid[idx]})
    
    # Upper limit: search in (p_hat, 1)
    upper <- tryCatch({
      uniroot(target, interval = c(p_hat, 1 - 1e-6), extendInt = "upX")$root
    }, error = function(e) {
      grid <- seq(p_hat, 1 - 1e-6, length.out = 1000)
      idx <- which.min(abs(target(grid)))
      grid[idx]})
    
    # Convert to a rate per 1,000 mosquitoes
    data.frame(
      MLE = round(p_hat * 1000, 2),
      CI_lower = round(max(0, lower) * 1000, 2),
      CI_upper = round(min(1, upper) * 1000, 2))
    
  }
}

# 2. CREATE GENUS GROUPS (Culex vs Others)=====================================
data <- set1_positiveRhab %>%
  mutate(
    Genus_group = ifelse(Genus == "Culex", "Culex", "Other"),
    MERDV_status = as.numeric(as.character(MERDV_status)))

# 3. SUMMARIZE BY SITE AND GENUS GROUP==========================================
summary_data <- data %>%
  group_by(Place, Genus_group) %>%
  summarise(
    Pools = n(),
    Mosquitoes = sum(Num_mosquitoes),
    Positive = sum(MERDV_status),
    pool_data = list(data.frame(
      Num_mosquitoes = Num_mosquitoes,
      MERDV_status = MERDV_status )),
    .groups = 'drop')

# 4. CALCULATE MIR AND MLE======================================================
result <- summary_data %>%
  rowwise() %>%
  mutate(
    MIR = round((Positive / Mosquitoes) * 1000, 2),
    MLE_result = list(calc_mle(pool_data)),
    MLE = MLE_result$MLE,
    CI_lower = MLE_result$CI_lower,
    CI_upper = MLE_result$CI_upper
  ) %>%
  select(-MLE_result, -pool_data) %>%
  ungroup()

# 5. ADD TOTALS=================================================================
totals_by_genus <- data %>%
  group_by(Genus_group) %>%
  summarise(
    Place = "Total",
    Pools = n(),
    Mosquitoes = sum(Num_mosquitoes),
    Positive = sum(MERDV_status),
    pool_data = list(data.frame(
      Num_mosquitoes = Num_mosquitoes,
      MERDV_status = MERDV_status)),
    .groups = 'drop')

# Total by site (all genus)
totals_by_site <- data %>%
  group_by(Place) %>%
  summarise(
    Genus_group = "Total",
    Pools = n(),
    Mosquitoes = sum(Num_mosquitoes),
    Positive = sum(MERDV_status),
    pool_data = list(data.frame(
      Num_mosquitoes = Num_mosquitoes,
      MERDV_status = MERDV_status
    )), .groups = 'drop')

# All sites and genera 
final_total <- data %>%
  summarise(
    Place = "Total",
    Genus_group = "All",
    Pools = n(),
    Mosquitoes = sum(Num_mosquitoes),
    Positive = sum(MERDV_status),
    pool_data = list(data.frame(
      Num_mosquitoes = Num_mosquitoes,
      MERDV_status = MERDV_status))
  )

# Calculate MLE for totals
totals_combined <- bind_rows(final_total, totals_by_genus, totals_by_site) %>%
  rowwise() %>%
  mutate(
    MIR = round((Positive / Mosquitoes) * 1000, 2),
    MLE_result = list(calc_mle(pool_data)),
    MLE = MLE_result$MLE,
    CI_lower = MLE_result$CI_lower,
    CI_upper = MLE_result$CI_upper
  ) %>%
  select(-MLE_result, -pool_data) %>%
  ungroup()

# 6. FINAL TABLE================================================================
final_table <- bind_rows(totals_combined, result) %>%
  arrange(
    factor(Place, levels = c("Total", "Conserved-area", "Cultural-center", "Eco-Hotel")),
    factor(Genus_group, levels = c("All", "Total", "Culex", "Other"))) %>%
  mutate(
    `CI (95%)` = paste0(CI_lower, "–", CI_upper)) %>%
  select(
    Site = Place,
    Genus = Genus_group,
    Pools,
    n = Mosquitoes,
    Positive,
    `MIR/1000` = MIR,
    `MLE/1000` = MLE,
    `CI (95%)`
  )
#Export
write_xlsx(final_table, "../OUTPUTS/MERDV_MLE_infection_rates.xlsx")