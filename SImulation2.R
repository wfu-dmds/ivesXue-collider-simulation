# Assuming the original setup is already defined
# Modifying the data generation function in Box 3 to include early measurements

generateData <- function(n, seed){
  set.seed(seed)
  # Early measurements
  Sodium_gr_0 <- rnorm(n, mean = 140, sd = 15) # Sodium intake at time 0
  SBP_1 <- 1.05 * Sodium_gr_0 + rnorm(n) # SBP at time 1
  Proteinuria_2 <- 2.00 * SBP_1 + 2.80 * Sodium_gr_0 + rnorm(n) # Proteinuria at time 2
  
  # Original measurements
  Age_years <- rnorm(n, 65, 5)
  Sodium_gr_3 <- Age_years / 18 + rnorm(n) # Sodium intake at time 3
  sbp_in_mmHg_4 <- 1.05 * Sodium_gr_3 + 2.00 * Age_years + rnorm(n) # SBP at time 4
  Proteinuria_5 <- 2.00 * sbp_in_mmHg_4 + 2.80 * Sodium_gr_3 + rnorm(n) # Proteinuria at time 5
  
  data.frame(Sodium_gr_0, SBP_1, Proteinuria_2, Age_years, Sodium_gr_3, sbp_in_mmHg_4, Proteinuria_5)
}

ObsData <- generateData(n = 1000, seed = 777)

# Adjusting regression models in Box 4 to use the newly simulated variables
fit0 <- lm(sbp_in_mmHg_4 ~ Sodium_gr_3, data = ObsData)
fit1 <- lm(sbp_in_mmHg_4 ~ Sodium_gr_3 + Age_years, data = ObsData)
# Using Proteinuria_2 (early measured) instead of Proteinuria_5
fit2 <- lm(sbp_in_mmHg_4 ~ Sodium_gr_3 + Age_years + Proteinuria_2, data = ObsData)

# Visualizations and further analysis would follow as per the original structure,
# but focusing on the adjusted models and data.

