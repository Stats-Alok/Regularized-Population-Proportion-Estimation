library(glmnet)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(psych)

set.seed(123)

# -----------------------------
# Load Data
# -----------------------------
my_data <- read.csv("/Volumes/ALOK (M&C)/Dataset/BodyFat.csv")
my_data <- my_data[,-1]          # Auxilairy Information 
# -----------------------------
# Create Binary Response
# -----------------------------
my_data$Y <- ifelse(
  my_data$Percent.body.fat.from.Siri.s..1956..equation > 20,
  1, 0
)

# -----------------------------
# Define Predictors
# -----------------------------
X_names <- setdiff(names(my_data),
                   c("Y", "Percent.body.fat.from.Siri.s..1956..equation"))

X <- as.matrix(my_data[,c(-1,-15)])

e = eigen(t(X) %*% X)$val
CI=sqrt(max(e)/min(e)); CI

colnames(X) <- c(
  "Age", "Weight", "Height", "Neck",
  "Chest", "Abdomen", "Hip", "Thigh",
  "Knee", "Ankle", "Biceps", "Forearm", "Wrist"
)


pairs.panels(
  X,
  gap = 0,
  pch = 21,
  lm = TRUE,
  main = "Distributions and Correlations",
  cex.labels = 1,
  cex.main = 1,
  cex.cor = 0.7   # 👈 reduce correlation text size
)

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


# Standardize predictors
my_data[X_names] <- scale(my_data[X_names])

# -----------------------------
# Parameters
# -----------------------------
m <- 10
r <- 10
n <- m * r
rep <- 10000

N <- nrow(my_data)
pop_prop <- mean(my_data$Y)

# -----------------------------
# RSS Function
# -----------------------------
rss_sample <- function(data, m, r, X_names){
  
  X_all <- matrix(0, m*r, length(X_names))
  Y_all <- numeric(m*r)
  
  counter <- 1
  
  for(cycle in 1:r){
    
    idx <- sample(1:nrow(data), m*m, replace = FALSE)
    temp <- data[idx, ]
    
    sets <- split(temp, rep(1:m, each = m))
    
    for(i in 1:m){
      s <- sets[[i]]
      
      # Ranking variable (important!)
      s <- s[order(s$Abdomen.2.circumference..cm), ]
      
      selected <- s[i, ]
      
      X_all[counter, ] <- as.matrix(selected[, X_names])
      Y_all[counter] <- selected$Y
      
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
# Simulation Loop
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
  rss <- rss_sample(my_data, m, r, X_names)
  
  X_rss <- rss$X
  y_rss <- rss$Y
  
  RSS_prop[iter] <- mean(y_rss)
  
  # -----------------
  # Check for class issue
  # -----------------
  if(length(unique(y_rss)) < 2){
    CRSS_prop[iter] <- NA
    next
  }
  
  # -----------------
  # Ridge Logistic Regression
  # -----------------
  cv_fit <- cv.glmnet(X_rss, y_rss,
                      alpha = 0, nfolds = 10,
                      family = "binomial",
                      standardize = TRUE)
  
  fit <- glmnet(X_rss, y_rss,
                alpha = 0,
                family = "binomial",
                lambda = cv_fit$lambda.min,
                standardize = TRUE)
  
  beta_hat <- as.vector(fit$beta)
  intercept <- fit$a0
  
  # -----------------
  # CRSS Estimator
  # -----------------
  CRSS_values <- numeric(m*m*r)
  k <- 1
  
  for(cycle in 1:r){
    
    idx <- sample(1:nrow(my_data), m*m, replace = FALSE)
    temp <- my_data[idx, ]
    
    sets <- split(temp, rep(1:m, each = m))
    
    for(i in 1:m){
      s <- sets[[i]]
      
      # SAME ranking variable
      s <- s[order(s$Abdomen.2.circumference..cm), ]
      
      for(j in 1:m){
        
        if(i == j){
          CRSS_values[k] <- s$Y[j]
        } else {
          x_vec <- as.matrix(s[j, X_names])
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
# Remove NA values
# -----------------------------
valid_idx <- !is.na(CRSS_prop)

SRS_prop <- SRS_prop[valid_idx]
RSS_prop <- RSS_prop[valid_idx]
CRSS_prop <- CRSS_prop[valid_idx]

# -----------------------------
# Performance Metrics
# -----------------------------
ARB_SRS <- mean(abs(SRS_prop - pop_prop))
ARB_RSS <- mean(abs(RSS_prop - pop_prop))
ARB_CRSS <- mean(abs(CRSS_prop - pop_prop))

MSE_SRS <- mean((SRS_prop - pop_prop)^2)
MSE_RSS <- mean((RSS_prop - pop_prop)^2)
MSE_CRSS <- mean((CRSS_prop - pop_prop)^2)

SRS_RE <- MSE_SRS / MSE_CRSS
RSS_RE <- MSE_RSS / MSE_CRSS

# -----------------------------
# Output
# -----------------------------
cat("ARB_SRS =", ARB_SRS, "\n")
cat("ARB_RSS =", ARB_RSS, "\n")
cat("ARB_CRSS =", ARB_CRSS, "\n\n")

cat("MSE_SRS =", MSE_SRS, "\n")
cat("MSE_RSS =", MSE_RSS, "\n")
cat("MSE_CRSS =", MSE_CRSS, "\n\n")

cat("Relative Efficiency (CRSS vs SRS) =", SRS_RE, "\n")
cat("Relative Efficiency (CRSS vs RSS) =", RSS_RE, "\n")

# -----------------------------
# Plot
# -----------------------------
plot_data <- data.frame(
  value = c(SRS_prop, RSS_prop, CRSS_prop),
  Method = factor(rep(c("SRS", "RSS", "CRSS"),
                      each = length(SRS_prop)))
)

ci_df <- plot_data %>%
  group_by(Method) %>%
  summarise(
    ci_lower = quantile(value, 0.05),
    ci_upper = quantile(value, 0.95),
    mean_value = mean(value),
    .groups = "drop"
  )

ggplot(plot_data, aes(x = value, fill = Method)) +
  geom_density(alpha = 0.7, adjust = 2) +
  geom_vline(xintercept = pop_prop, color = "black") +
  geom_vline(data = ci_df,
             aes(xintercept = mean_value, color = Method),
             linetype = "dashed",
             show.legend = FALSE) +
  geom_rect(data = ci_df,
            aes(xmin = ci_lower, xmax = ci_upper,
                ymin = 0, ymax = Inf, fill = Method),
            alpha = 0.2,
            inherit.aes = FALSE) +
  labs(x = "Estimated Value", y = "Density")

# -----------------------------
# Save Results
# -----------------------------
results <- data.frame(
  Estimator = c("SRS","RSS","CRSS"),
  MSE = c(MSE_SRS, MSE_RSS, MSE_CRSS),
  RE = c(1, SRS_RE, RSS_RE),
  ARB = c(ARB_SRS, ARB_RSS, ARB_CRSS)
)
results
