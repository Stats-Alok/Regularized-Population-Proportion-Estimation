library(glmnet)
library(ggplot2)
library(reshape2)
library(ggplot2)
library(ggridges)
library(psych)
library(dplyr)
library(openxlsx)


set.seed(123)

# -----------------------------
# Parameters
# -----------------------------
N <- 400        # sample size
p <- 15          # number of predictors (try 3, 6, 9)
rho <- 0.95      # correlation (0.80 to 0.99)
sigma <- 1      # error variance

# -----------------------------
# Step 1: Generate Z
# -----------------------------
Z <- matrix(rnorm(N * (p + 1)), N, p + 1)

# -----------------------------
# Step 2: Generate correlated X
# -----------------------------
X <- matrix(0, N, p)

for(i in 1:N){
  for(j in 1:p){
    X[i, j] <- sqrt(1 - rho^2) * Z[i, j] + rho * Z[i, p + 1]
  }
}

# -----------------------------
# Step 3: Standardize X
# -----------------------------
X <- scale(X)   # ensures X'X is correlation-like

# -----------------------------
# Step 4: Generate beta (normalized)
# -----------------------------
beta_raw <- runif(p, 0.5, 1.5)   # random positive values
true_beta <- beta_raw / sqrt(sum(beta_raw^2))  # enforce beta'β = 1

# -----------------------------
# Step 5: Generate error term
# -----------------------------
epsilon <- rnorm(N, mean = 0, sd = sigma)

# -----------------------------
# Step 6: Generate response y
# -----------------------------
prob <-( 1/(1+exp(X %*% true_beta + epsilon)))
Y <- rbinom(N,1,prob)

colnames(X) <- c("X1", "X2", "X3","X4", "X5", "X6","X7", "X8", "X9","X10", "X11", "X12","X13", "X14", "X15")


pairs.panels(
  X,
  gap = 0,
  pch = 21,
  lm = TRUE,
  labels = c("X1", "X2", "X3","X4", "X5", "X6","X7", "X8", "X9","X10", "X11", "X12","X13", "X14", "X15"),
  main = "Distributions and Correlations",
  cex.labels = 1,
  cex.main = 1,
  cex.cor = 0.7   # 👈 reduce correlation text size
)



my_data <- data.frame(Y, X)

# -----------------------------
# Parameters
# -----------------------------
m <- 10         # set size
r <- 10         # number of cycles
n <- m * r     # RSS sample size
rep <- 1000

pop_prop <- mean(Y)

# -----------------------------
# RSS function with cycles
# -----------------------------
rss_sample <- function(data, m, r){
  
  X_all <- matrix(0, m*r, p)
  Y_all <- numeric(m*r)
  
  counter <- 1
  
  for(cycle in 1:r){
    
    # draw m sets each of size m
    idx <- sample(1:nrow(data), m*m, replace = FALSE)
    temp <- data[idx, ]
    
    sets <- split(temp, rep(1:m, each = m))
    
    for(i in 1:m){
      s <- sets[[i]]
      
      # ranking using X1 (perfect ranking proxy)
      s <- s[order(s$X1), ]
      
      # select ith order statistic
      X_all[counter, ] <- as.matrix(s[i, 2:16])
      Y_all[counter] <- s$Y[i]
      
      counter <- counter + 1
    }
  }
  
  list(X = X_all, Y = Y_all)
}

# -----------------------------
# Storage
# -----------------------------
SRS_prop <- numeric(rep)
RSS_prop <- numeric(rep)
CRSS_prop <- numeric(rep)

# -----------------------------
# Simulation
# -----------------------------
for(iter in 1:rep){
  
  # -----------------
  # SRS
  # -----------------
  srs <- my_data[sample(1:N, n), ]
  SRS_prop[iter] <- mean(srs$Y)
  
  # -----------------
  # RSS
  # -----------------
  rss <- rss_sample(my_data, m, r)
  
  X_rss <- rss$X
  y_rss <- rss$Y
  
  RSS_prop[iter] <- mean(y_rss)
  
  # -----------------
  # Ridge Logistic Regression
  # -----------------
  cv_fit <- cv.glmnet(X_rss, y_rss,
                      alpha = 0,
                      family = "binomial",
                      standardize = FALSE)
  
  fit <- glmnet(X_rss, y_rss,
                alpha = 0,
                family = "binomial",
                lambda = cv_fit$lambda.min,
                standardize = FALSE)
  
  beta_hat <- as.vector(fit$beta)
  intercept <- fit$a0
  
  # -----------------
  # Construct CRSS estimator
  # -----------------
  CRSS_values <- numeric(m*m*r)
  k <- 1
  
  for(cycle in 1:r){
    
    idx <- sample(1:nrow(my_data), m*m, replace = FALSE)
    temp <- my_data[idx, ]
    
    sets <- split(temp, rep(1:m, each = m))
    
    for(i in 1:m){
      s <- sets[[i]]
      s <- s[order(s$X1), ]
      
      for(j in 1:m){
        
        if(i == j){
          CRSS_values[k] <- s$Y[j]
        } else {
          x_vec <- as.matrix(s[j, 2:16])
          linpred <- intercept + x_vec %*% beta_hat
          prob_hat <- 1 / (1 + exp(-linpred))
          
          CRSS_values[k] <- as.numeric(prob_hat)
        }
        
        k <- k + 1
      }
    }
  }
  
  CRSS_prop[iter] <- mean(CRSS_values)
}

# -----------------------------
# Performance Metrics
# -----------------------------
ARB_SRS <- mean(abs(SRS_prop - pop_prop))
ARB_RSS <- mean(abs(RSS_prop - pop_prop))
ARB_Con.RSS <- mean(abs(CRSS_prop - pop_prop))

MSE_SRS <- mean((SRS_prop - pop_prop)^2)
MSE_RSS <- mean((RSS_prop - pop_prop)^2)
MSE_Con.RSS <- mean((CRSS_prop - pop_prop)^2)

SRS_RE <- MSE_SRS / MSE_Con.RSS
RSS_RE <- MSE_RSS / MSE_Con.RSS

# -----------------------------
# Output
# -----------------------------
cat("ARB_SRS =", ARB_SRS, "\n")
cat("ARB_RSS =", ARB_RSS, "\n")
cat("ARB_CRSS =", ARB_Con.RSS, "\n\n")

cat("MSE_SRS =", MSE_SRS, "\n")
cat("MSE_RSS =", MSE_RSS, "\n")
cat("MSE_CRSS =", MSE_Con.RSS, "\n\n")

cat("Relative Efficiency (CRSS vs SRS) =", SRS_RE, "\n")
cat("Relative Efficiency (CRSS vs RSS) =", RSS_RE, "\n")



# -----------------------------
# Combine results into one data frame
# -----------------------------
plot_data <- data.frame(
  value = c(SRS_prop, RSS_prop, CRSS_prop),
  Method = factor(rep(c("SRS", "RSS", "Con.RSS"), each = length(SRS_prop)))
)

# -----------------------------
# Calculate 95% CI for each Method
# -----------------------------
ci_df <- plot_data %>%
  group_by(Method) %>%
  summarise(
    ci_lower = quantile(value, 0.05),
    ci_upper = quantile(value, 0.95),
    mean_value = mean(value),
    .groups = "drop"
  )

# -----------------------------
# Density Plot
# -----------------------------
ggplot(plot_data, aes(x = value, fill = Method))   +
  geom_density(alpha = 0.7, adjust = 2) +
  
  # True value line
  geom_vline(xintercept = pop_prop, color = "black") +
  
  # Mean lines
  geom_vline(data = ci_df,
             aes(xintercept = mean_value, color = Method),
             linetype = "dashed",
             show.legend = FALSE) +
  
  # CI bands
  geom_rect(data = ci_df,
            aes(xmin = ci_lower, xmax = ci_upper,
                ymin = 0, ymax = Inf, fill = Method),
            alpha = 0.2,
            inherit.aes = FALSE) +
  
  labs(
    #title = "Distribution of Estimated Proportions",
    x = "Estimated Value",
    y = "Density"
  ) +
  
  theme(plot.title = element_text(hjust = 0.5))



results <- data.frame(
  Estimator = c("SRS","RSS","CRSS"),
  
  MSE = c(MSE_SRS, MSE_RSS, MSE_Con.RSS),
  
  RE = c(1, SRS_RE, RSS_RE),
  
  ARB = c(ARB_SRS, ARB_RSS, ARB_Con.RSS)
)
results
