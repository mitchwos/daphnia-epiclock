# CODE EDITED FROM EPICLOCK CREATED BY EAMONN MALLON

# Remove unused objects and free up memory
rm(list = ls())
gc()  # Force garbage collection
options(max.print=1000)  # Limit console print output
options(stringsAsFactors = FALSE)  # Set default for character data to avoid memory issues with factors

# Load necessary libraries -----------------------------------------------------
library(dplyr)
library(purrr)
library(tidyverse)
library(ggfortify)
library(glmnet)
library(ggplot2)
library(ggpubr)
library(matrixStats)
library(vegan)
library(modeest)
library(Metrics)
library(vip)
library(viridis)
library(broom)
library(lme4)
library(emmeans)
library(data.table)
library(glmmTMB)
library(RColorBrewer)
library(scales)
library(doParallel)
library(foreach)
library(caret)
library(report)


set.seed(123)  # Set a seed for reproducibility

setwd('~/Documents/IRP/R/clock/C32')

# Load full dataset using fread for efficiency
methyl_data <- fread("~/Documents/IRP/R/clock/C32/C32_daphnia_percent")

# Efficiently convert percentages to beta values and replace NA with 0
methyl_data[, (1:24) := lapply(.SD, function(x) replace_na(x / 100, 0)), .SDcols = 1:24]

# Free up memory
gc()  # Garbage collection to reclaim memory

# Filtering on diff cpgs
time_related <- read.delim("~/Documents/IRP/R/clock/C32/C32_Significant_Age_CpG.txt")

# Filter rows where fdrs < 0.05
time_related <- time_related[time_related$fdrs < 0.05, ]

time_related <- time_related %>%
  mutate(cpg = paste(chr, pos, sep = "."))

# Ensure both are data.tables
setDT(methyl_data)
setDT(time_related)

# Perform a semi-join using the 'cpg' column
#filtered_data <- diapause_data[(cpg %in% time_related$cpg) | (cpgtwo %in% time_related$cpg)]
filtered_data <- methyl_data[(Row_ID %in% time_related$cpg)]
#filtered_data <- diapause_data

# Transpose the filtered data -------------------------------------------------
# Transpose the data and convert it to a data.table
# Add the column names of f_diapause_data as the first row
filtered_data_with_colnames <- rbindlist(list(as.list(colnames(filtered_data)), filtered_data))

# Transpose the data including the column names
trans_data <- as.data.table(t(filtered_data_with_colnames))
rm(filtered_data, filtered_data_with_colnames)

# Set row 25 as the column names
setnames(trans_data, as.character(trans_data[25, ]))  # Row 25 becomes the column names

# Step 2: Remove rows 81
trans_data <- trans_data[-c(25)]#change to just 81

# Add 'diet' column: Extract text between first and second underscore
trans_data$diet <- sub("^[^_]*_([^_]*)_.*", "\\1", trans_data$Row_ID)

# Add 'day' column: Extract the old/young and converts to 50/10 days
trans_data <- trans_data %>%
  mutate(day = case_when(
    str_detect(Row_ID, "old") ~ 50,
    str_detect(Row_ID, "young") ~ 10,
    TRUE ~ NA_real_
  ))

# Move diet and day columns up in data table
trans_data <- trans_data %>%
  select(Row_ID, diet, day, everything())

fwrite(trans_data, "~/Documents/IRP/R/clock/C32/C32_trans_data.csv")


#==============================Epigenetic clock=========================================

# Free up memory
gc()  # Garbage collection to reclaim memory

# Filter the data 
train_data <- trans_data
numeric_data <- train_data[, -c("Row_ID", "day", "diet")]
numeric_data[] <- lapply(numeric_data, as.numeric)
x <- as.matrix(numeric_data)
y <- train_data$day  # 'day' is the response variable

# Step 1: Split the data into training and testing sets
train_index <- createDataPartition(y, p = 0.8, list = FALSE)  # 80% for training
x_train <- x[train_index, ]
x_test <- x[-train_index, ]
y_train <- y[train_index]
y_test <- y[-train_index]

# Step 2: Define the train control (can be cross-validation or bootstrapping)
# Define refined train control with more rigorous cross-validation
train_control <- trainControl(
  method = "loocv") # leave one out cross validation

tuneGrid <- expand.grid(
  alpha = c(0.5),  # Ridge (0), ElasticNet (0.5), and Lasso (1)
  lambda = 10^seq(-4, 1, length = 100)
)

# Step 3: Set up parallel backend for training (optional)
num_cores <- detectCores() - 1  
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Step 4: Train the model on the training set
cv_model <- train(
  x = x_train,
  y = y_train,
  method = "glmnet",
  preProcess = c("center", "scale"),  # PCA included , "pca"
  thresh = 0.8,  # Retain 80% of variance
  trControl = train_control,
  tuneGrid = tuneGrid,
  maxit = 500  # Specify maxit directly here for glmnet
)

plot(cv_model$finalModel)

# Step 5: Evaluate the model on the test set
y_pred <- predict(cv_model, newdata = x_test)

# Step 6: Calculate performance metrics on the test set
# (Use different metrics based on the type of model)
results <- postResample(pred = y_pred, obs = y_test)

# Print results
print(cv_model)

# Create a data frame for the actual vs predicted values
results_df <- data.frame(
  Actual = y_test,
  Predicted = y_pred
)

# Scatter plot with regression line
ggplot(results_df, aes(x = Actual, y = Predicted)) +
  geom_point(color = "blue") +        # Scatter plot of actual vs predicted
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Regression line
  labs(
    title = "Actual vs Predicted Values",
    x = "Actual Values",
    y = "Predicted Values"
  ) +
  theme_minimal()

# Retrieve the best parameters
best_lambda <- cv_model$bestTune$lambda
best_alpha <- cv_model$bestTune$alpha

# Step 5: Fit the final model on the 'Control' data with the best parameters
final_model <- train(
  x = x,
  y = y,
  method = "glmnet",
  tuneGrid = expand.grid(alpha = best_alpha, lambda = best_lambda),
  trControl = trainControl(method = "none")
)

# Step 6: Prepare the full dataset predictor matrix (x_full) for prediction
numeric_data_full <- trans_data[, -c("Row_ID", "day","diet")]
numeric_data_full[] <- lapply(numeric_data_full, as.numeric)
x_full <- as.matrix(numeric_data_full)

# Step 7: Make predictions on the full dataset using the trained model
y_pred <- predict(final_model, newdata = x_full)

# Clean up parallel backend
stopCluster(cl)
registerDoSEQ()

# Step 8: Add predictions to the full dataset for visualization
graph_data <- trans_data %>%
  mutate(pred_age = y_pred)

# Step 9: Visualize predicted vs actual age
ggplot(graph_data, aes(x = day, y = pred_age)) +
  geom_abline(intercept = 0, slope = 1, color='grey') +
  geom_jitter(width = 0.4, height = 0.4, alpha = 0.7, size = 2.5, color = '#00ADFF', fill = '#94DEFF') +  # Increase jitter width/height, reduce size
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", linewidth = 1.2, color = '#00ADFF') +  # Use linewidth for line thickness
  labs(
    title = "C32 Epigenetic Clock",
    x = "Chronological Age (Days)",
    y = "Epigenetic Age (Days)",
  ) +
  theme_classic(base_size = 14) +  # Use a classic theme with base font size adjustment
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),  # Center and bold title
    axis.title = element_text(size = 14),  # Increase axis title font size
    axis.text = element_text(size = 12),   # Increase axis text font size
  ) 

ggsave("~/Documents/IRP/R/clock/C32/C32_model_figure.pdf")

# Extract Spearman's rho
cor_val <- cor.test(graph_data$day, graph_data$pred_age, method ='spearman')
print(cor_val)

# Extract slope and intercept of model
fit <- lm(pred_age ~ day, data = graph_data)
coef(fit)

# Extract RSME
rmse <- sqrt(mean(fit$residuals^2))
print(rmse)
                               
# Extract coefficients from the final glmnet model within the caret object
final_glmnet <- final_model$finalModel

# Extract non-zero coefficients for the selected lambda value
non_zero_coefs <- as.data.frame(as.matrix(coef(final_glmnet, s = best_lambda)))
non_zero_coefs <- non_zero_coefs[non_zero_coefs != 0, , drop = FALSE]  # Only keep non-zero coefficients

# Exporting to CSV with row names
write.csv(non_zero_coefs, "~/Documents/IRP/R/clock/C32/C3S2_clock_loci.csv", row.names = TRUE)

