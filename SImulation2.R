# Step 1: Install and Load Necessary Libraries
# You may need to install packages using install.packages() if you haven't already.
library(visreg)     # To visualize regression output
library(ggplot2)    # For general plotting
library(gridExtra)  # To combine plots
library(forestplot) # To make forest plots

library(ggdag)


# Step 2: Original Data Generation and Analysis
# Box 1 and Box 2 from original script remain unchanged

# Box 1: data generation, crude and adjusted models
N <- 1000
set.seed(777)
W <- rnorm(N)
A <- 0.5 * W + rnorm(N)
Y <- 0.3 * A + 0.4 * W + rnorm(N)
fit1 <- lm(Y ~ A)
fit2 <- lm(Y ~ A + W)

## Y: outcome
## X: exposure
## W: confounder
dag_1 <- dagify(
  A ~ W,
  Y ~ A + W,
  exposure = "A",
  outcome = "Y",
  coords = time_ordered_coords()
)

ggdag(dag_1)


g1 <- visreg(fit1, "A", gg = TRUE, line = list(col = "blue"),
             points = list(size = 2, pch = 1, col = "black")) + 
  theme_classic() +
  coord_cartesian(ylim = c(-4, 4)) 
plot(g1) # Not plotted in paper
g2 <- visreg(fit2, "A", gg = TRUE, line = list(col = "blue"),
             points = list(size = 2, pch = 1, col = "black")) + 
  theme_classic() +
  coord_cartesian(ylim = c(-4, 4)) +
  ggtitle("Figure 2A")
plot(g2) 

# Box 2: data generation, with collider
A <- rnorm(N)
Y <- 0.3 * A + rnorm(N)
C <- 1.2 * A + 0.9 * Y + rnorm(N)
fit3 <- lm(Y ~ A)
fit4 <- lm(Y ~ A + C)

## A: exposure
## Y: outcome
## C: collider
dag_2 <- dagify(
  Y ~ A,
  C ~ A + Y,
  exposure = "A",
  outcome = "Y",
  coords = time_ordered_coords()
)
ggdag(dag_2)

g4 <- visreg(fit4, "A", gg = TRUE, line = list(col = "red"),
             points = list(size = 2, pch = 1, col = "black")) + 
  theme_classic() +
  coord_cartesian(ylim = c(-4, 4)) +
  ggtitle("Figure 2B") 
plot(g4)

# Figure  2
grid.arrange(g2, g4, ncol = 2)

# Box 3: Data generation for proteinuria example (Original)

dag_3 <- dagify(
  Sodium_gr ~ Age_years,
  sbp_in_mmHg ~ Sodium_gr + Age_years,
  Proteinureia_in_mg ~ sbp_in_mmHg + Sodium_gr,
  exposure = "Sodium_gr",
  outcome = "sbp_in_mmHg",
  coords = time_ordered_coords()
)
ggdag(dag_3)

generateData <- function(n, seed) {
  set.seed(seed)
  Age_years <- rnorm(n, 65, 5)
  Sodium_gr <- Age_years / 18 + rnorm(n)
  sbp_in_mmHg <- 1.05 * Sodium_gr + 2.00 * Age_years + rnorm(n)
  hypertension <- ifelse(sbp_in_mmHg >= 140, 1, 0)
  Proteinuria_in_mg <- 2.00 * sbp_in_mmHg + 2.80 * Sodium_gr + rnorm(n)
  data.frame(sbp_in_mmHg, hypertension, Sodium_gr, Age_years, Proteinuria_in_mg)
}
ObsData <- generateData(n = 1000, seed = 777)
head(ObsData)

# Original Regression Models and Visualization
# Box 4: regression models: crude, adjusted for confounder, inclusion of  collider
#        Note that we use the object names "fit1, fit2" as above, for consistency with the paper
fit0 <- lm(sbp_in_mmHg ~ Sodium_gr, data = ObsData)
fit1 <- lm(sbp_in_mmHg ~ Sodium_gr + Age_years , data = ObsData)
fit2 <- lm(sbp_in_mmHg ~ Sodium_gr + Age_years + Proteinuria_in_mg, data = ObsData)

m1 <- visreg(fit0, "Sodium_gr", gg = TRUE, xlab = "Sodium (gr)", ylab = "SBP (mmHg)",
             line = list(col = "blue"),
             points = list(size = 2, pch = 1, col = "black"), bty = "n") + 
  theme_classic() +
  coord_cartesian(ylim = c(110, 165)) +
  ggtitle("Figure 4A") 
m2 <- visreg(fit1, "Sodium_gr", gg = TRUE, xlab = "Sodium (gr)", ylab = "SBP (mmHg)",
             line = list(col = "blue"),
             points = list(size = 2, pch = 1, col = "black"), bty = "n") + 
  theme_classic() +
  coord_cartesian(ylim = c(129, 140)) +
  ggtitle("Figure 4B") 
m3 <- visreg(fit2, "Sodium_gr", gg = TRUE, xlab = "Sodium (gr)", ylab = "SBP (mmHg)",
             line = list(col = "red"),
             points = list(size = 2, pch = 1, col = "black"), bty = "n") + 
  theme_classic() +
  coord_cartesian(ylim = c(129, 140)) +
  ggtitle("Figure 4C") 

# Figure 4
# Note: figures combined slightly differently as in paper
grid.arrange(m1, m2, m3, ncol = 3)

# Box 5: proteinuria example as in Box 4; but with hypertension as binary outcome
#        again: crude, adjusted for confounder, inclusion of collider 
fit3 <- glm(hypertension ~ Sodium_gr, family = binomial(link = "logit"), data = ObsData)
or <- round(exp(fit3$coef)[2], 3)
ci95 <- exp(confint(fit3))[-1, ]
lci <- round(ci95[1], 3)
uci <- round(ci95[2], 3)
model <- c("Crude")
result1 <- data.frame(model, or, lci, uci, stringsAsFactors = FALSE)

fit4 <- glm(hypertension ~ Sodium_gr + Age_years, family = binomial(link = "logit"), data = ObsData) 
or <- round(exp(fit4$coef)[2], 3)
ci95 <- exp(confint(fit4))[2, ]
lci <- round(ci95[1], 3)
uci <- round(ci95[2], 3)
model <- c("Adjusted")
result2 <- data.frame(model, or, lci, uci, stringsAsFactors = FALSE)

fit5 <- glm(hypertension ~ Sodium_gr + Age_years + Proteinuria_in_mg,
            family = binomial(link = "logit"), data = ObsData) 
or <- round(exp(fit5$coef)[2], 3)
ci95 <- exp(confint(fit5))[2, ]
lci <- round(ci95[1], 3)
uci <- round(ci95[2], 3)
model <- c("Collider")
result3 <- data.frame(model, or, lci, uci, stringsAsFactors = FALSE)

# Models fit visualization (Forest plot function and plot) to depict the collider effect
or_graph <- function(fp){
  
  tabla <- cbind(c("Model", paste(fp$model)),
                 c("Odds ratio", fp$or), 
                 c("95%CI", paste0("(", fp$lci, " - ", fp$uci, ")")))
  
  forestplot(labeltext = tabla,
             graph.pos = 3,
             mean = c(NA, fp$or),
             is.summary = c(TRUE, rep(FALSE, nrow(fp))),
             lower = c(NA, fp$lci),
             upper = c(NA, fp$uci),
             xlab = "Odds ratio",
             txt_gp = fpTxtGp(label = gpar(cex = 1.25),
                              ticks = gpar(cex = 1.1),
                              xlab  = gpar(cex = 1.2),
                              title = gpar(cex = 1.2)),
             col = fpColors(box = "blue", lines = "blue", zero = "black"),
             cex = 0.9,
             clip = c(0, 10),
             zero = 1,
             boxsize = 0.05,
             lwd.ci = 2,
             ci.vertices = TRUE,
             lineheight = "auto",
             xticks = seq(0, 10, 1),
             ci.vertices.height = 0.1,
             grid = TRUE
  )
}

# Figure 5
fp <- rbind(result1,result2,result3)
or_graph(fp)

# Box 6:  Monte Carlo Simulation
set.seed(050472) # Seed for reproducibility
R <- 1000        # Number of simulation runs
# Vectors to store results
true <- rep(NA, R)    
collider <- rep(NA, R)
se <- rep(NA, R)

for(r in 1:R) {
  if (r%%10 == 0) cat(paste("This is simulation run number", r, "\n"))
  # Function to generate data 
  generateData <- function(n){
    Age_years <- rnorm(n, 65, 5)
    Sodium_gr <- Age_years / 18 + rnorm(n)
    sbp_in_mmHg <- 1.05 * Sodium_gr + 2.00 * Age_years + rnorm(n)
    Proteinuria_in_mg <- 2.00 * sbp_in_mmHg + 2.80 * Sodium_gr + rnorm(n)
    data.frame(sbp_in_mmHg, Sodium_gr, Age_years, Proteinuria_in_mg)
  }
  ObsData <- generateData(n = 10000) 
  # Estimated effect using the correct model (including Sodium_gr and Age_years)
  true[r] <- summary(lm(sbp_in_mmHg ~ Sodium_gr + Age_years, data = ObsData))$coef[2, 1]
  # Estimated effect using the incorrect model (including the collider)
  collider[r] <- summary(lm(sbp_in_mmHg ~ Sodium_gr + Age_years + Proteinuria_in_mg, data = ObsData))$coef[2, 1]
  # Standard error
  se[r]       <- summary(lm(sbp_in_mmHg ~ Sodium_gr + Age_years + Proteinuria_in_mg, data = ObsData))$coef[2, 2]
}


mean(true)       # Average estimated effect using the correct model (including Sodium_gr and Age_years)
# Correct!
mean(collider)   # Average estimated effect using the incorrect model (including the collider)
# Biased!
# Bias 
Bias <- true - collider
mean(Bias)
relBias <- (true - collider) / true
mean(relBias) * 100

# Pragmatic confidence interval based on average standard error
lci <- (mean(collider) - 1.96 * mean(se))
uci <- (mean(collider) + 1.96 * mean(se))

# Step 3: Adapted Data Generation with New Variables

## time 0: age
## time 1: sodium
## time 2: sbp
## time 3: proteinuria

generateData <- function(n, seed) {
  set.seed(seed)
  Age_years <- rnorm(n, 65, 5)
  Sodium_gr <- Age_years / 18 + rnorm(n)
  sbp_in_mmHg <- 1.05 * Sodium_gr + 2.00 * Age_years + rnorm(n)
  hypertension <- ifelse(sbp_in_mmHg >= 140, 1, 0)
  Proteinuria_in_mg <- 2.00 * sbp_in_mmHg + 2.80 * Sodium_gr + rnorm(n)
  data.frame(sbp_in_mmHg, hypertension, Sodium_gr, Age_years, Proteinuria_in_mg)
}

new_dag <- dagify(
  sodium_0 ~ age,
  sbp_1 ~ age + sodium_0,
  sodium_1 ~ age + sodium_0,
  proteinuria_2 ~ sodium_1 + sbp_1,
  sodium_2 ~ age + sodium_1,
  sbp_2 ~ age + sbp_1 + sodium_1,
  proteinuria_3 ~ sodium_2 + sbp_2,
  sbp_3 ~ age + sbp_2 + sodium_2,
  exposure = "sodium_1",
  outcome = "sbp_2",
  coords = time_ordered_coords()
)

ggdag(new_dag)
ggdag_adjustment_set(adjust_for(new_dag, "age"))

new_dag <- dagify(
  sodium_neg1 ~ age,
  sbp_0 ~ age + sodium_neg1,
  proteinuria_1 ~ sodium_neg1 + sbp_0,
  sodium_0 ~ sodium_neg1,
  sbp_1 ~ age + sodium_0,
  proteinuria_2 ~ sodium_0 + sbp_1,
  exposure = "sodium_0",
  outcome = "sbp_1",
  coords = time_ordered_coords()
)

ggdag(new_dag)
ggdag_adjustment_set(new_dag)
generateDataExtended <- function(n, seed) {
  set.seed(seed)
  Age_years <- rnorm(n, 64)
  Sodium_gr_neg1 <- Age_years / 18 + rnorm(n, sd = sqrt(0.5))
  sbp_in_mmHg_0 <-  1.05 * Sodium_gr_neg1 + 2.00 * Age_years + rnorm(n)
  Proteinuria_in_mg_1 <- 2.00 * sbp_in_mmHg_0 + 2.80 * Sodium_gr_neg1 + rnorm(n)
  Sodium_gr_0 <- Sodium_gr_neg1 + rnorm(n, sd = sqrt(0.5))
  sbp_in_mmHg_1 <- 1.05 * Sodium_gr_0 + 2.00 * Age_years + rnorm(n)
  Proteinuria_in_mg_2 <- 2.00 * sbp_in_mmHg_1 + 2.80 * Sodium_gr_0 + rnorm(n)
  hypertension <- ifelse(sbp_in_mmHg_1 >= 140, 1, 0)
  data.frame(sbp_in_mmHg_1, hypertension, Sodium_gr_0, Age_years, Sodium_gr_neg1,
             sbp_in_mmHg_0, Proteinuria_in_mg_1, Proteinuria_in_mg_2)
}

# Step 4: Adapted Regression Analyses
ObsDataExtended <- generateDataExtended(n = 1000, seed = 777)

fit6 <- lm(sbp_in_mmHg_1 ~ Sodium_gr_0 + Age_years + Proteinuria_in_mg_2, data = ObsDataExtended)
fit6
fit7 <- lm(sbp_in_mmHg_1 ~ Sodium_gr_0 + Age_years + Proteinuria_in_mg_1, data = ObsDataExtended)
fit7
# Step 5: Visualizations and Comparisons
# Adapted Model Visualization
m6 <- visreg(fit6, "Sodium_gr_0", gg = TRUE, xlab = "Sodium at Time 0 (gr)", ylab = "SBP at Time 1 (mmHg)",
             line = list(col = "red"),
             points = list(size = 2, pch = 1, col = "black"), bty = "n") + 
  theme_classic() +
  ggtitle("Adapted Model with Post-Exposure/Outcome Proteinuria")

# Visualization of the adapted model
print(m6)

m7 <- visreg(fit7, "Sodium_gr_0", gg = TRUE, xlab = "Sodium at Time 0 (gr)", ylab = "SBP at Time 1 (mmHg)",
             line = list(col = "red"),
             points = list(size = 2, pch = 1, col = "black"), bty = "n") + 
  theme_classic() +
  ggtitle("Adapted Model with Pre-Exposure/Outcome Proteinuria")

# Visualization of the adapted model
print(m7)
