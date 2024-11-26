#Time series Analysis of MB Microsporidia 
# Load required libraries
library(zoo)
library(tseries)
library(timeSeries)
library(tsibble)
library(xts)
library(TTR) #For decomposing time series
library(lubridate)
library(dplyr)
library(tidyr)
library(readxl)
library(stats)
library(ggplot2)
library(astsa) 
library(gridExtra)


##loading data
ndata=read_excel("micro_daily.xlsx")
View(ndata)
data=ndata[,c(1,3,4,6,10:21)]
View(data)

#changing y to 0-1
data$MB_Prevalence=data$MB_Prevalence/100

# Grouping data to weekly clusters
# Creating a complete sequence of weeks from min to max timestamp
all_weeks <- data.frame(
  week_start = seq(floor_date(min(data$Date), unit="weeks"),
                   floor_date(max(data$Date), unit="weeks"),
                   by = "1 week")
)

# Generating a cluster ID for each week_start value
week_clusters <- all_weeks %>%
  mutate(cluster_id = paste0("week_", row_number()))

# Adding week_start to the data and then join with week_clusters to assign cluster IDs
data_with_clusters <- data %>%
  mutate(week_start = floor_date(data$Date, unit="weeks")) %>%
  left_join(week_clusters, by = "week_start")

# Rearrange columns: Move week_start and cluster_id to the front
data_reordered <- data_with_clusters %>%
  select(week_start, cluster_id, everything())

# Create a Reference Table:
cluster_date_reference <- data_reordered %>%
  group_by(cluster_id) %>%
  summarise(min_date = min(week_start))
names(data_reordered)
column_names = c("MB_Prevalence","Oviposited_screened","No_MB_positives","temperature","humidity","evaporation","solar_radiation",
                 "windspeed","dewpoint","cloud_cover","solarenergy","uvindex","rainfall","wind_direction")

# Aggregate data based on 'cluster_id'
aggregated_by_cluster <- data_reordered %>%
  group_by(cluster_id) %>%
  summarise(across(all_of(column_names), mean, .names = "mean_{.col}"))

# Join the Reference to the Aggregated Data:
aggregated_by_cluster <- aggregated_by_cluster %>%
  left_join(cluster_date_reference, by = "cluster_id")

# Order by the Minimum Date:
aggregated_by_cluster <- aggregated_by_cluster %>%
  arrange(min_date)
View(aggregated_by_cluster)

# Adding week numbers to all_weeks
all_weeks$week_number <- 1:nrow(all_weeks)


# Join the Reference to the Aggregated Data:
aggregated_with_week_start <- aggregated_by_cluster %>%
  left_join(cluster_date_reference, by = "cluster_id")

# Check the column names
View(aggregated_with_week_start)

# Rename the week_start column:
aggregated_with_week_start <- aggregated_with_week_start %>%
  rename(week_start = min_date.x)

complete_aggregated_data <- left_join(all_weeks, aggregated_with_week_start, by = "week_start")

# Plotting function
plot_trends <- function(df, colname) {
  aggregated_colname <- paste0("mean_", colname)
  
p <- ggplot(df, aes(x = week_number, y = get(aggregated_colname))) +
    geom_line(color = "lightblue", linewidth= 1, na.rm = TRUE) +  # Using na.rm to handle NA values
    geom_point(aes(x = week_number, y = get(aggregated_colname)), na.rm = TRUE) +  # Add points
    geom_smooth(se = FALSE, color = 'maroon', na.rm = TRUE) +  # Using na.rm to handle NA values
    ggtitle(paste("Trends for", colname)) +
    xlab("Week") + ylab(paste(colname)) +
    scale_x_continuous(
      breaks = seq(1, max(df$week_number), by = 2), 
      labels = paste0("Week ", seq(1, max(df$week_number), by = 2))
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis text for better readability
  
  print(p)
}

# Looping through the specified column names
column_names1 = c("MB_Prevalence")
for (colname in column_names1) {
  plot_trends(complete_aggregated_data, colname)
}


#Renaming to weekly data
weekly_data<-aggregated_by_cluster
weekly_data

# Rename columns
names(weekly_data)
weekly_data <- weekly_data %>%
  rename(
    Oviposited_screened=mean_Oviposited_screened,
    MB_Prevalence=mean_MB_Prevalence,
    No_MB_positives=mean_No_MB_positives,
    temperature = mean_temperature,
    humidity=mean_humidity,
    evaporation=mean_evaporation,
    solar_radiation=mean_solar_radiation,
    windspeed=mean_windspeed,
    dewpoint=mean_dewpoint,
    cloud_cover=mean_cloud_cover,
    solarenergy=mean_solarenergy,
    uvindex=mean_uvindex,
    rainfall=mean_rainfall,
    wind_direction=mean_wind_direction
  )


View(weekly_data)

#BINOMIAL GAM
gam_data=weekly_data

names(gam_data)
# Categorize UV index
gam_data$cat_uvindex <- cut(gam_data$uvindex,
                            breaks = c(0, 2, 5, 7, 10, Inf),
                            labels = c("Low", "Moderate", "High", "Very High", "Extreme"),
                            right = TRUE)
 
failures=as.integer(gam_data$Oviposited_screened-gam_data$No_MB_positives)
successes=as.integer(gam_data$No_MB_positives)
gam_data=cbind(gam_data,successes,failures)
View(gam_data)

# FITTING BINOMIAL GAM MODEL WITHOUT SPLINES
model=gam(cbind(successes,failures) ~ temperature + humidity + evaporation + solar_radiation+ 
            windspeed + dewpoint + cloud_cover + solarenergy + 
            uvindex+ rainfall+ wind_direction,family=binomial(link="logit"), data = gam_data)
summary(model)


library(utils)
cov <- colnames(gam_data)[5:15]
min_aic <- Inf
best_model1 <- NULL


# Looping over the non-empty different combinations of covariates
for(i in 1:11){
  combos <- combn(cov, i)
  for(j in 1:ncol(combos)){
    current_covariates <- combos[, j]
    formula_string <- paste("cbind(successes, failures) ~", paste(sapply(current_covariates, function(cov) paste("s(", cov, ", k=3)", sep="")), collapse = "+"))
    gam_model1 <- gam(as.formula(formula_string), family=binomial(link="logit"), data = gam_data)
    model_aic <- AIC(gam_model1)
    if(model_aic < min_aic){
      min_aic <- model_aic
      best_model1 <- gam_model1
    }
  }
}

print(best_model1)
#best_model11=gam(cbind(successes, failures) ~ s(temperature) + s(evaporation) + 
#s(windspeed) + s(dewpoint) + s(cloud_cover) + s(solarenergy) + 
  #s(uvindex) + s(rainfall) + s(wind_direction)+family=binomial(link="logit"), data = gam_data)

summary(best_model1)
AIC(best_model1)


####GAM MODELS WITH LAGGED COVARIATES#####
covariate_names <- c("MB_Prevalence","temperature","humidity","evaporation","solar_radiation",
                     "windspeed","uvindex","dewpoint","cloud_cover","solarenergy","rainfall","wind_direction" )

# Create a copy of the original dataset
lagged_data <- gam_data

# Loop through the covariate names and create lagged versions in lagged_data
for (col_name in covariate_names) {
  lag_name <- paste0(col_name, "_lag1")
  lagged_data[[lag_name]] <- c(NA, head(lagged_data[[col_name]], -1))
}

library(car)

# Fit a linear regression model
lm_model <- lm(MB_Prevalence ~ temperature + humidity + evaporation + 
                 solar_radiation + windspeed + dewpoint + cloud_cover + 
                 solarenergy + rainfall + wind_direction + 
                 MB_Prevalence_lag1 + temperature_lag1 + humidity_lag1 + 
                 evaporation_lag1 + windspeed_lag1 + solar_radiation_lag1+
                 dewpoint_lag1 + cloud_cover_lag1 + 
                 rainfall_lag1 + wind_direction_lag1, 
               data = lagged_data)
# Calculate VIF for each predictor
vif_values <- vif(lm_model)

# Print the VIF values
print(vif_values)

library(car)

# Continue dropping the variable with the highest VIF until all VIFs are less than 10
#while(TRUE) {
#m_model <- lm(MB_Prevalence ~., data = lagged_data)
#vif_values <- vif(lm_model)

# Check if the highest VIF is greater than 10
#if(max(vif_values) > 10) {
# Drop the variable with the highest VIF
#drop_var <- names(vif_values)[which.max(vif_values)]
#lagged_data <- lagged_data[, !(names(lagged_data) %in% drop_var)]
#} else {
# break
#}
#}





# Load the mgcv package
library(mgcv)
# Your existing covariates
existing_covariates2 <- c("temperature", "evaporation", "windspeed","dewpoint",
                          "cloud_cover","rainfall","solarenergy","uvindex","wind_direction")

# List of additional covariates to consider
additional_covariates2 <- c("temperature_lag1","humidity_lag1","evaporation_lag1","solar_radiation_lag1",
                            "windspeed_lag1","uvindex_lag1","dewpoint_lag1","cloud_cover_lag1",
                            "rainfall_lag1","wind_direction_lag1")

# Initialize variables to store the best model and its AIC
best_model2 <- NULL
best_aic2 <- Inf

# Loop through different combinations of additional covariates
for (i in 0:length(additional_covariates2)) {
  # Create a combination of covariates to include
  covariates_to_include2 <- c(existing_covariates2, head(additional_covariates2, i))
  
  # Create the formula for the GAM with smoothed covariates
  formula_string <- paste("cbind(successes, failures) ~", paste(sapply(covariates_to_include2, function(cov) paste("s(", cov, ", k=3)", sep="")), collapse = "+"))
  # Fit the GAM
  gam_model2 <- gam(as.formula(formula_string), family=binomial(link="logit"), data = lagged_data)
  
  # Calculate AIC
  aic2 <- AIC(gam_model2)
  
  # Check if this model has a lower AIC than the current best
  if (aic2 < best_aic2) {
    best_aic2 <- aic2
    best_model2 <- gam_model2
  }
}

# Print the best model and its AIC
print(best_model2)
summary(best_model2)
#best_model2=cbind(successes, failures) ~ s(temperature) + s(evaporation) + 
  #s(windspeed) + s(dewpoint) + s(cloud_cover) + s(rainfall) + 
  #s(solarenergy) + s(uvindex) + s(wind_direction) + s(temperature_lag1)


####GAM MODELS WITH INTERACTIONS#####
# Compute interaction terms
inter_data=gam_data[,5:15]
interaction_data <- as.data.frame(apply(combn(names(inter_data), 2), 2, function(cols) {
  inter_data[, cols[1]] * inter_data[, cols[2]]
}))

# Rename the interaction columns
colnames(interaction_data) <- apply(combn(names(inter_data), 2), 2, paste, collapse="_x_")

interaction_data <- cbind(gam_data, interaction_data)
interaction_data
names(interaction_data)

# Your existing covariates
existing_covariates3 <-c("temperature", "evaporation", "windspeed","dewpoint",
                         "cloud_cover","solarenergy","rainfall","uvindex","wind_direction")

# List of additional covariates to consider
additional_covariates3 <- c(
  "temperature_x_humidity",
  "temperature_x_evaporation",
  "temperature_x_solar_radiation",
  "temperature_x_windspeed",
  "temperature_x_dewpoint",
  "temperature_x_cloud_cover",
  "temperature_x_solarenergy",
  "temperature_x_uvindex",
  "temperature_x_rainfall",
  "temperature_x_wind_direction",
  "humidity_x_evaporation",
  "humidity_x_solar_radiation",
  "humidity_x_windspeed",
  "humidity_x_dewpoint",
  "humidity_x_cloud_cover",
  "humidity_x_solarenergy",
  "humidity_x_uvindex",
  "humidity_x_rainfall",
  "humidity_x_wind_direction",
  "evaporation_x_solar_radiation",
  "evaporation_x_windspeed",
  "evaporation_x_dewpoint",
  "evaporation_x_cloud_cover",
  "evaporation_x_solarenergy",
  "evaporation_x_uvindex",
  "evaporation_x_rainfall",
  "evaporation_x_wind_direction",
  "solar_radiation_x_windspeed",
  "solar_radiation_x_dewpoint",
  "solar_radiation_x_cloud_cover",
  "solar_radiation_x_solarenergy",
  "solar_radiation_x_uvindex",
  "solar_radiation_x_rainfall",
  "solar_radiation_x_wind_direction",
  "windspeed_x_dewpoint",
  "windspeed_x_cloud_cover",
  "windspeed_x_solarenergy",
  "windspeed_x_uvindex",
  "windspeed_x_rainfall",
  "windspeed_x_wind_direction",
  "dewpoint_x_cloud_cover",
  "dewpoint_x_solarenergy",
  "dewpoint_x_uvindex",
  "dewpoint_x_rainfall",
  "dewpoint_x_wind_direction",
  "cloud_cover_x_solarenergy",
  "cloud_cover_x_uvindex",
  "cloud_cover_x_rainfall",
  "cloud_cover_x_wind_direction",
  "solarenergy_x_uvindex",
  "solarenergy_x_rainfall",
  "solarenergy_x_wind_direction",
  "uvindex_x_rainfall",
  "uvindex_x_wind_direction",
  "rainfall_x_wind_direction")


# Initialize variables to store the best model and its AIC
best_model3 <- NULL
best_aic3 <- Inf
names(interaction_data)
# Loop through different combinations of additional covariates
for (i in 0:length(additional_covariates3)) {
  # Create a combination of covariates to include
  covariates_to_include3 <- c(existing_covariates3, head(additional_covariates3, i))
  
  # Create the formula for the GAM with smoothed covariates
  formula_string <- paste("cbind(successes, failures) ~", paste(sapply(covariates_to_include3, function(cov) paste("s(", cov, ")", sep="")), collapse = "+"))
  # Fit the GAM
  gam_model3 <- gam(as.formula(formula_string), family=binomial(link="logit"), data = interaction_data)
  
  # Calculate AIC
  aic3 <- AIC(gam_model3)
  
  # Check if this model has a lower AIC than the current best
  if (aic3 < best_aic3) {
    best_aic3 <- aic3
    best_model3 <- gam_model3
  }
}

# Print the best model and its AIC
print(best_model3)
summary(best_model3)
AIC(best_model3)



#GAM MODELS WITH INTERACTIONS AND LAGGED COVARIATES
# Compute interaction terms

lag_interact_data <- cbind(lagged_data, interaction_data)
names(lag_interact_data)
# Your existing covariates
existing_covariates4 <-c("temperature", "evaporation", "windspeed","dewpoint",
                         "cloud_cover","solarenergy","rainfall","uvindex","wind_direction")

# List of additional covariates to consider
additional_covariates4 <- c("temperature_lag1","humidity_lag1","temperature_x_humidity",
                            "temperature_x_evaporation","temperature_x_solarenergy")

# Initialize variables to store the best model and its AIC
best_model4 <- NULL
best_aic4 <- Inf

# Loop through different combinations of additional covariates
for (i in 0:length(additional_covariates4)) {
  # Create a combination of covariates to include
  covariates_to_include4 <- c(existing_covariates4, head(additional_covariates4, i))
  
  # Create the formula for the GAM with smoothed covariates
  formula4 <- as.formula(paste("cbind(successes,failures)~", paste(sapply(covariates_to_include4, function(cov) paste("s(", cov, ")", sep="")), collapse = "+")))
  
  # Fit the GAM
  gam_model4 <- gam(formula4,family = binomial(link="logit"), data = lag_interact_data)
  
  # Calculate AIC
  aic4 <- AIC(gam_model4)
  
  # Check if this model has a lower AIC than the current best
  if (aic4 < best_aic4) {
    best_aic4 <- aic4
    best_model4 <- gam_model4
  }
}

# Print the best model and its AIC
print(best_model4)
summary(best_model4)
AIC(best_model4)

summary(model)
summary(best_model1)
summary(best_model2)
summary(best_model3)
summary(best_model4)

AIC(model)
AIC(best_model1)
AIC(best_model2)
AIC(best_model3)
AIC(best_model4)

BIC(model)
BIC(best_model1)
BIC(best_model2)
BIC(best_model3)
BIC(best_model4)


library(mgcv)
library(ggplot2)
model_summary=summary(best_model2)
model_summary


sig_terms_indices <- c(1,3,7,8)
# Determine the number of rows needed for the layout
num_rows <- ceiling(length(sig_terms_indices) / 2)
png("significant_smooth_terms.png", width=800, height=400*num_rows) 
par(mfrow=c(num_rows,2))
# Loop through and plot each significant term
for (i in sig_terms_indices) {
  plot(best_model2, select = i, shade = TRUE, shade.col = "lightblue", seWithMean = TRUE)
}
dev.off()





#CROSS VALIDATION

#Cross validation of model without splines
# Split data into k folds
k <- 10
set.seed(123)
folds1 <- sample(1:k, nrow(gam_data), replace = TRUE)

# Initialize variables to store results
mae_values1 <- numeric(k)
rmse_values1 <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:k) {
  # Split data into training and testing sets
  train_data1 <- gam_data[folds1 != i, ]
  test_data1 <- gam_data[folds1 == i, ]
  
  # Fit the model on the training data
  model_train1=gam(cbind(successes,failures) ~ temperature + humidity + evaporation + solar_radiation + 
                     windspeed + dewpoint + cloud_cover + solarenergy + cat_uvindex + 
                     rainfall + wind_direction,family=binomial(link = "logit"),data=train_data1)
  
  # Make predictions on the test data
  predictions1 <- predict(model_train1, newdata = test_data1, type = "response")
  
  actual_values1=test_data1$MB_Prevalence
  
  mae1 <- mean(abs(actual_values1 - predictions1))
  rmse1 <- sqrt(mean((actual_values1 - predictions1)^2))
  
  print(paste("MAE: ", mae1))
  print(paste("RMSE: ", rmse1))
  
}

# Calculate average MAE and RMSE
mean_mae1 <- mean(mae1)
mean_rmse1 <- mean(rmse1)


# Print results
print(paste("Average MAE1: ", round(mean_mae1, 5)))
print(paste("Average RMSE1: ", round(mean_rmse1, 5)))


#cross validation of best model with splines
# Split data into k folds
k <- 10
set.seed(123)
folds2 <- sample(1:k, nrow(gam_data), replace = TRUE)

# Initialize variables to store results
mae_values2 <- numeric(k)
rmse_values2 <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:k) {
  # Split data into training and testing sets
  train_data2 <- gam_data[folds2 != i, ]
  test_data2 <- gam_data[folds2 == i, ]
  
  # Fit the model on the training data
  model_train2=gam(cbind(successes,failures) ~ s(temperature,k=3) + s(windspeed,k=3) + s(evaporation,k=3)+
                     s(dewpoint,k=3) + s(cloud_cover,k=3) + s(solarenergy,k=3) + s(rainfall,k=3) + 
                     s(uvindex,k=3) +s(wind_direction,k=3),family=binomial(link = "logit"),data = train_data2)
  
  # Make predictions on the test data
  predictions2 <- predict(model_train2, newdata = test_data2, type = "response")
  
  actual_values2 =test_data2$MB_Prevalence
  mae2 <- mean(abs(actual_values2 - predictions2))
  rmse2 <- sqrt(mean((actual_values2 - predictions2)^2))
  
  print(paste("MAE: ", mae2))
  print(paste("RMSE: ", rmse2))
  
}

# Calculate average MAE and RMSE
mean_mae2 <- mean(mae2)
mean_rmse2 <- mean(rmse2)

# Print results
print(paste("Average MAE2: ", round(mean_mae2, 5)))
print(paste("Average RMSE2: ", round(mean_rmse2, 5)))



#cross validation model with lagged effects
# Split data into k folds
k <- 10
set.seed(123)
folds3 <- sample(1:k, nrow(lagged_data), replace = TRUE)

# Initialize variables to store results
mae_values3 <- numeric(k)
rmse_values3 <- numeric(k)

# Perform k-fold cross-validation
for (i in 1:k) {
  # Split data into training and testing sets
  train_data3 <- lagged_data[folds3 != i, ]
  test_data3 <- lagged_data[folds3 == i, ]
  
  # Fit the model on the training data
  model_train3=gam(cbind(successes,failures) ~ s(temperature,k=3) + s(windspeed,k=3) + s(evaporation,k=3)+
                     s(dewpoint,k=3) + s(cloud_cover,k=3) + s(solarenergy,k=3) + s(rainfall,k=3) + 
                     s(uvindex,k=3) +s(wind_direction,k=3) +s(temperature_lag1,k=3),family=binomial(link = "logit"),data = train_data3)
  
  # Make predictions on the test data
  predictions3 <- predict(model_train3, newdata = test_data3, type = "response")
  
  actual_values3=test_data3$MB_Prevalence
  mae3 <- mean(abs(actual_values3 - predictions3))
  rmse3 <- sqrt(mean((actual_values3 - predictions3)^2))
  
  print(paste("MAE: ", mae3))
  print(paste("RMSE: ", rmse3))
  
}

# Calculate average MAE and RMSE
mean_mae3 <- mean(mae3)
mean_rmse3 <- mean(rmse3)



# Print results
print(paste("Average MAE1: ", round(mean_mae1, 5)))
print(paste("Average RMSE1: ", round(mean_rmse1, 5)))
print(paste("Average MAE2: ", round(mean_mae2, 5)))
print(paste("Average RMSE2: ", round(mean_rmse2, 5)))
print(paste("Average MAE3: ", round(mean_mae3, 5)))
print(paste("Average RMSE3: ", round(mean_rmse3, 5)))


###RESIDUAL ANALYSIS###
library(mgcv)
library(ggplot2)
library(gridExtra)

residuals <- residuals(best_model2)
# 1. Residuals vs. Fitted
p1 <- ggplot(data.frame(Fitted = fitted(best_model2), Residuals = residuals(best_model2)), aes(x = Fitted, y = Residuals)) +
  geom_point(alpha = 0.6) +
  geom_smooth(se = FALSE, color = "red") +
  labs(title = "Residuals vs. Fitted") +
  theme_minimal()

# 2. Histogram of Residuals
resid_mean <- mean(residuals(best_model2))
resid_sd <- sd(residuals(best_model2))
p2 <- ggplot(data.frame(Residuals = residuals(best_model2)), aes(x = Residuals)) +
  geom_histogram(aes(y = ..density..), fill = "steelblue", color = "white", bins = 20) +
  geom_density(color = "red", size = 1) +  # Overlay the actual density
  stat_function(fun = dnorm, args = list(mean = resid_mean, sd = resid_sd), color = "blue", size = 1) +  # Normal curve
  labs(title = "Histogram of Residuals with Normal Curve") +
  theme_minimal()
# 3. Q-Q Plot
p3 <- ggplot(data.frame(Residuals = residuals(best_model2)), aes(sample = Residuals)) +
  geom_qq() +
  geom_qq_line(color = "red") +
  labs(title = "Q-Q Plot of Residuals") +
  theme_minimal()

# 4. ACF Plot
resid_acf <- acf(residuals(best_model2), plot = FALSE)
df_acf <- data.frame(Lag = resid_acf$lag[-1], ACF = resid_acf$acf[-1])  # Exclude the lag=0 value
n <- nrow(df_acf)


p4 <- ggplot(df_acf, aes(x = Lag, y = ACF)) +
  geom_bar(stat = 'identity', fill = "steelblue") +
  geom_hline(yintercept = c(1.96/sqrt(n), -1.96/sqrt(n)), linetype = "dashed", color = "red") +
  labs(title = "ACF of Residuals") +
  theme_minimal()

# Combine plots
grid.arrange(p1, p2, p3, p4, ncol = 2)

View(lagged_data)

#PLOT FOR ACTUAL VS PREDICTED MB_PREVALANCE
fitted_values=best_model2$fitted.values


if (!is.na(fitted_values[1])) {
  fitted_values <- c(NA, fitted_values)
}

fitted_values
#add fitted values to data
lagged_data$fitted <- fitted_values
Time=weekly_data$min_date
lagged_data=cbind(Time,lagged_data)

library(ggplot2)

# Assuming your data is in lagged_data dataframe
# Sort by Time just to be safe
lagged_data <- lagged_data[order(lagged_data$Time), ]

# Convert the Date to relative week number
lagged_data$Week <- as.numeric(difftime(lagged_data$Time, min(lagged_data$Time), units = "weeks")) + 1

# Colors for actual and fitted values
color_actual <- "blue"
color_fitted <- "red"

# Range for the y-axis
range_y <- range(c(lagged_data$MB_Prevalence, lagged_data$fitted), na.rm = TRUE)

# Plot using ggplot2
plot <- ggplot(lagged_data, aes(x = Week)) + 
  geom_line(aes(y = MB_Prevalence, color = "Actual")) + 
  geom_line(aes(y = fitted, color = "Fitted"), linetype = "dashed") +
  scale_color_manual(values = c("Actual" = color_actual, "Fitted" = color_fitted)) +
  labs(title = "Actual vs. Fitted Time Series",
       x = "Time (in weeks)", y = "MB_Prevalance") +
  theme_minimal() +
  theme(legend.title = element_blank()) +  # Remove the legend title
  scale_x_continuous(breaks = seq(1, max(lagged_data$Week), by = 10), 
                     labels = paste("Week", seq(1, max(lagged_data$Week), by = 10)))  # Custom x-axis labels

print(plot)

