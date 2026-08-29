library(glmnet)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(psych)
library(corrplot)

set.seed(123)

# -----------------------------
# Load Data
# -----------------------------
my_data <- read.csv("Dataset Path")
my_data <- my_data[ ,-1]          # Auxiliary Informations 

# -----------------------------
#  Binary Study Variable
# -----------------------------
my_data$diagnosis <- ifelse(
  my_data$diagnosis == "M",
  0, 1
)

# -----------------------------
# Create X and Y
# -----------------------------

X <- as.matrix(my_data[,-c(1,32)])
Y <- my_data[,1]

e = eigen(t(X) %*% X)$val
CI=sqrt(max(e)/min(e)); CI

# -----------------------------
# Correlation Plot
# -----------------------------

df_corr <- cor(X)

corrplot(
  df_corr,
  method = "circle",
  type = "upper",
  order = "hclust",
  tl.col = "black",
  tl.srt = 45,
  addrect = 5,
  col = colorRampPalette(c("#2c7bb6", "#ffffbf", "#d7191c"))(200)
)

# -----------------------------
# Create Dataset
# -----------------------------

my_data <- data.frame(Y, X)

# -----------------------------
# Parameters
# -----------------------------
N <- length(Y)
p <- 30
m <- 5         # set size
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
      s <- s[order(s$radius_mean), ]
      
      # select ith order statistic
      X_all[counter, ] <- as.matrix(s[i, 2:31])
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
      s <- s[order(s$radius_mean), ]
      
      for(j in 1:m){
        
        if(i == j){
          CRSS_values[k] <- s$Y[j]
        } else {
          x_vec <- as.matrix(s[j, 2:31])
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
  
  PRE = c(100, SRS_RE, RSS_RE),
  
  ARB = c(ARB_SRS, ARB_RSS, ARB_Con.RSS)
)
results

