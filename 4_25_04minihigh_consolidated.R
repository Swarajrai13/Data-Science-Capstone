# ===============================
# CAPSTONE ANALYSIS SCRIPT (UPDATED)
# ===============================
# Covers: SUMMARY STATS → OLS → PANEL → CLUSTER → ANOVA → FORECASTING

# ----------
# 1. LIBRARIES
# ----------
library(MASS)
library(dplyr)        # data manipulation
library(skimr)        # summary stats
library(broom)        # tidy regression output
library(corrr)        # correlation matrices
library(plm)          # panel models
library(lmtest)       # hypothesis tests
library(sandwich)     # robust SEs
library(ggplot2)      # plotting
library(cluster)      # clustering
library(factoextra)   # clustering viz
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)          # spatial data
library(fmsb)         # radar chart
library(glmnet)    #LASSO Model
library(lme4)     #Mixed Effects
library(sjPlot)  #Models Comparison
library(viridis)
install.packages("dtwclust")
library(dtwclust)
library(tidyverse)
library(tibble)
library(purrr)
library(forcats)



# ========================================================
# 2. EXPLORATORY + CLEANING WORKFLOW
# ========================================================
setwd("~/Desktop")
final_data <- read.csv("final_data_with_region_4_24.csv") %>%
  janitor::clean_names()

# Check missingness
colSums(is.na(final_data))

# ========================================================
# 3. SUMMARY STATISTICS (DV & IVs)
# ========================================================
summary_data <- final_data %>%
  dplyr::select(new_business_density_rate,
                account_age_15, wefs, estimate,
                gdp_constant_2015usd, gov_edu_spending_pct_gdp,
                working_age_pop, population_growth) %>%
  na.omit()

# Display summaries
t(summary(summary_data))
skim(summary_data)

# ========================================================
# 4. BASELINE OLS REGRESSIONS
# ========================================================
ols_data <- summary_data

# Univariate models
wefs   <- lm(new_business_density_rate ~ wefs, data = ols_data)
wgi_estimate <- lm(new_business_density_rate ~ estimate, data = ols_data)
gdp <-  lm(new_business_density_rate ~ gdp_constant_2015usd, data = ols_data)
edu    <- lm(new_business_density_rate ~ gov_edu_spending_pct_gdp, data = ols_data)
account<- lm(new_business_density_rate ~ account_age_15, data = ols_data)
work  <- lm(new_business_density_rate ~ working_age_pop, data = ols_data)
popg <- lm(new_business_density_rate ~ population_growth, data = ols_data)

tab_model(
  wefs,
  wgi_estimate,
  gdp,
  edu,
  account,
  work,
  popg,
  dv.labels = c(
    "WEFS","WGI Estimate", "GDP","GOV SPEND EDU","FIN ACCOUNT AGE 15+",
    "WORKING AGE POP","POP GROWTH"
  ),
  show.se  = TRUE,
  show.aic = TRUE
)


# Multivariate OLS
ols_multi <- lm(new_business_density_rate ~
                  account_age_15 + wefs + estimate +
                  gdp_constant_2015usd + gov_edu_spending_pct_gdp +
                  working_age_pop + population_growth,
                data = ols_data)
summary(ols_multi)

# ========================================================
#   Cook's Distance (Outlier Plotting)
# ========================================================
# Step 1: Load model (already created as ols_multi)
# Check it exists
summary(ols_multi)

# Step 2: Calculate Cook’s Distance
cooksd <- cooks.distance(ols_multi)

# Step 3: Plot Cook's Distance
plot(cooksd, 
     pch = 20, 
     cex = 1.5,
     main = "Cook's Distance for Multivariate OLS",
     ylab = "Cook's Distance",
     col = ifelse(cooksd > (4/length(cooksd)), "red", "black"))

# Add threshold line
abline(h = 4/length(cooksd), col = "red", lty = 2)

# Step 4: Identify influential points (above threshold)
influential_points <- which(cooksd > (4/length(cooksd)))

# Step 5: Print the rows of concern
ols_data[influential_points, ]

# ========================================================
#   EXCLUDING OUTLIERS & COMPARING (OLS_MULTI)
# ========================================================
# Step 1: Store influential indices (already done)
influential_points <- which(cooksd > (4/length(cooksd)))

# Step 2: Create a cleaned version of ols_data
ols_data_clean <- ols_data[-influential_points, ]

# Step 3: Re-run the multivariate OLS on cleaned data
ols_multi_clean <- lm(new_business_density_rate ~
                        account_age_15 + wefs + estimate +
                        gdp_constant_2015usd + gov_edu_spending_pct_gdp +
                        working_age_pop + population_growth,
                      data = ols_data_clean)

# Step 4: Compare original vs. cleaned model side-by-side
tab_model(
  ols_multi,
  ols_multi_clean,
  dv.labels = c("Full OLS Model", "OLS w/o Influential Points"),
  show.aic = TRUE,
  show.se  = TRUE,
  show.r2  = TRUE
)

#Mapping Outliers

# Step 1: Make sure the influential indices are stored
influential_points <- which(cooks.distance(ols_multi) > (4 / length(ols_multi$fitted.values)))

# Step 2: Retrieve country and year for those rows
# (Assuming `country` and `year` columns are in ols_data or accessible via final_data)
influential_info <- final_data[influential_points, c("country", "year")]

# Step 3: View the result
print(influential_info)

# ========================================================
#   STEPWISE AIC Model
# ========================================================

step_model <- stepAIC(ols_multi, direction = "both",
                      trace = FALSE)

summary(step_model)

# ========================================================
# 5. MULTICOLLINEARITY CHECK
# ========================================================
numeric_predictors <- ols_data %>%
  dplyr::select(new_business_density_rate, account_age_15, wefs,
                estimate, gdp_constant_2015usd,
                gov_edu_spending_pct_gdp, working_age_pop,
                population_growth)
correlate(numeric_predictors) %>% rplot()

# ========================================================
# 6. POOLED PANEL OLS
# ========================================================
pdata <- pdata.frame(final_data, index = c("country","year"))

pooling_model <- plm(new_business_density_rate ~ 
                       account_age_15 + wefs + estimate + 
                       gov_edu_spending_pct_gdp,
                     data = pdata,
                     model = "pooling")
summary(pooling_model)


# ========================================================
# 6. PANEL MODELING (RE/FE/HAUSMAN)
# ========================================================
# Random effects (dropping collinear variables)
re_model <- plm(new_business_density_rate ~
                  wefs + account_age_15 + estimate +
                  gov_edu_spending_pct_gdp,
                data = pdata, model = "random")
summary(re_model)
coeftest(re_model, vcovHC(re_model, type = "HC1", cluster = "group"))

# Fixed effects
fe_model <- plm(new_business_density_rate ~
                  wefs + account_age_15 + estimate +
                  gov_edu_spending_pct_gdp,
                data = pdata, model = "within")
summary(fe_model)

# Hausman test
phtest(fe_model, re_model)

# ========================================================
# 7. CLUSTERING & REGION ANALYSIS
# ========================================================
final_data_clean <- final_data %>%
  filter(complete.cases(select(., where(is.numeric))))
normalized <- scale(select(final_data_clean, where(is.numeric)))
set.seed(42)
kmeans_result <- kmeans(normalized, centers = 4, nstart = 25)
final_data_clean$cluster <- factor(kmeans_result$cluster)

cluster_summary <- final_data_clean %>%
  dplyr::group_by(cluster, region_2, sub_region) %>%
  dplyr::summarise(across(where(is.numeric), ~mean(.x, na.rm=TRUE)), .groups="drop")
print(cluster_summary)

# ========================================================
#  MIXED EFFECTS (HIERARCHICAL) 
# ========================================================

mixed_model <- lmer(new_business_density_rate ~
                      account_age_15 + wefs + estimate +
                      gov_edu_spending_pct_gdp +
                      (1 | country),
                    data = final_data_clean)
summary(mixed_model)

# ========================================================
# 8. OLS with Region‐Level Fixed Effects
# ========================================================

region_fe <- lm(new_business_density_rate ~
                  account_age_15 + wefs + estimate +
                  gov_edu_spending_pct_gdp +
                  factor(sub_region),
                data = final_data_clean)

summary(region_fe)

# ========================================================
# 9. ANOVA + POST-HOC REGION EFFECTS
# ========================================================
final_data_clean$region <- factor(final_data_clean$sub_region)
anova_model <- aov(new_business_density_rate ~ sub_region, data=final_data_clean)
summary(anova_model)
TukeyHSD(anova_model)

# ========================================================
# 10. COUNTRY CLUSTERS (RADAR CHART)
# ========================================================

# Re‐create the feature‐summary data frame:
features_df <- final_data %>%
  filter(!is.na(new_business_density_rate)) %>%
  group_by(country) %>%
  summarise(
    slope = coef(lm(new_business_density_rate ~ year))[2],
    mean = mean(new_business_density_rate),
    variance = var(new_business_density_rate),
    range = max(new_business_density_rate) - min(new_business_density_rate),
    first_diff_var = var(diff(new_business_density_rate)),
    autocorr = cor(head(new_business_density_rate, -1), tail(new_business_density_rate, -1)),
    .groups = "drop"
  )

features_df <- na.omit(features_df)
features_scaled <- features_df %>%
  dplyr::select(-country) %>%
  scale()

dist_matrix <- dist(features_scaled, method = "euclidean")

fviz_nbclust(features_scaled, kmeans, method = "wss") +
  labs(title = "Elbow Method for Optimal k")

set.seed(123)



gap_stat <- clusGap(features_scaled, FUN = kmeans, nstart = 25, K.max = 10, B = 50)



fviz_gap_stat(gap_stat)



# K-means clustering (choose k = 5)

set.seed(123)

kmeans_result <- kmeans(features_scaled, centers = 5)



# Add cluster label
features_df$cluster <- kmeans_result$cluster

table(features_df$cluster) #number of obsv in each cluster
#(note clusters 2 & 4) only have 2.

fviz_cluster(kmeans_result, data = features_scaled,
             geom = "point", ellipse.type = "norm",
             main = "K-Means Clustering of Countries")


world <- ne_countries(scale = "medium", returnclass = "sf")

world_clusters <- world %>%
  
  left_join(features_df, by = c("name" = "country"))



ggplot(data = world_clusters) +
  geom_sf(aes(fill = factor(cluster)), color = "white", size = 0.1) +
  scale_fill_brewer(palette = "Set1", na.value = "gray90") +
  labs(title = "Clusters of Countries by Business Density Features",
       fill = "Cluster") +
  theme_minimal()


# Join cluster labels with features
cluster_profiles <- features_df %>%
  pivot_longer(cols = -c(country, cluster), names_to = "feature", values_to = "value") %>%
  group_by(cluster, feature) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

# Standardize within features for visual comparison
cluster_profiles <- cluster_profiles %>%
  group_by(feature) %>%
  mutate(std_score = scale(mean_value)[,1]) %>%
  ungroup()

# Plot
ggplot(cluster_profiles, aes(x = feature, y = std_score, fill = factor(cluster))) +
  geom_col(position = position_dodge()) +
  facet_wrap(~ cluster, ncol = 1) +
  coord_flip() +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Mean Cluster Profiles (Standardized Feature Scores)",
    x = NULL, y = "Standardized Score"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

# ========================================================
# 11. K-Medoid Clusters (5/3)
# ========================================================

#11a. for n=4
actual_years <- c(2011, 2014, 2017, 2021)


ts_long <- final_data_clean %>%
  select(country, year, new_business_density_rate) %>%
  arrange(country, year)

country_year_counts <- ts_long %>%
  group_by(country) %>%
  summarise(n_years = n()) %>%
  arrange(desc(n_years))

View(country_year_counts) #n = 45 for n_year >= 4 (all = 4 - data limitation)

valid_countries <- country_year_counts %>%
  filter(n_years >= 4)

ts_filtered <- ts_long %>%
  filter(country %in% valid_countries$country)

ts_list <- ts_filtered %>%
  group_by(country) %>%
  summarise(ts = list(new_business_density_rate), .groups = "drop") %>%
  pull(ts)

lengths(ts_list)

set.seed(123)
kmedoid_dtw <- tsclust(
  ts_list,
  type = "partitional",
  k = 4,
  distance = "dtw_basic",
  centroid = "pam",
  trace = TRUE
)
# Extract medoid time series
medoid_series <- kmedoid_dtw@centroids

medoid_df <- map2_dfr(
  medoid_series,
  .y = seq_along(medoid_series),
  .f = function(vec, cluster_id) {
    tibble(
      year_index = 1:length(vec),
      business_density = vec,
      cluster = paste0("Cluster ", cluster_id)
    )
  }
)
# Replace generic year_index with actual years
medoid_df$year <- rep(actual_years, times = length(medoid_series))

# Plot
ggplot(medoid_df, aes(x = year, y = business_density, color = cluster)) +
  geom_line(size = 1.3) +
  scale_color_manual(
    values = c(
      "Cluster 1" = "#e41a1c",
      "Cluster 2" = "#377eb8",
      "Cluster 3" = "#4daf4a",
      "Cluster 4" = "#984ea3"
    )
  ) +
  labs(
    title = "Representative Time Series (Medoids) for Each DTW Cluster",
    x = "Year",
    y = "New Business Density Rate"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")


# Get country names for ts_list (same order as clustering output)
country_names <- ts_filtered %>%
  group_by(country) %>%
  summarise(n = n(), .groups = "drop") %>%
  pull(country)

# Pair with cluster results
dtw_clusters_updated <- tibble(
  country = country_names,
  dtw_cluster = factor(kmedoid_dtw@cluster)
)
world <- ne_countries(scale = "medium", returnclass = "sf")

# Join country-cluster data to map
map_data_dtw_updated <- left_join(world, dtw_clusters_updated, by = c("name" = "country"))

ggplot(map_data_dtw_updated) +
  geom_sf(aes(fill = dtw_cluster), color = "white", size = 0.15) +
  scale_fill_manual(
    name = "DTW Cluster",
    values = c(
      "1" = "#e41a1c",
      "2" = "#377eb8",
      "3" = "#4daf4a",
      "4" = "#984ea3"
    ),
    labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")
  ) +
  labs(
    title = "Country Clusters by Time Series of Business Density",
    subtitle = "Updated DTW Clustering (K-Medoids, 4-Year Minimum)"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

#11b. for n>=3

# STEP 1: Filter countries with ≥ 3 years of data
ts_long <- final_data_clean %>%
  select(country, year, new_business_density_rate) %>%
  arrange(country, year)

country_year_counts <- ts_long %>%
  group_by(country) %>%
  summarise(n_years = n(), .groups = "drop")

valid_countries <- country_year_counts %>%
  filter(n_years >= 3)

ts_filtered <- ts_long %>%
  filter(country %in% valid_countries$country)

# STEP 2: Create time series list for DTW
ts_list <- ts_filtered %>%
  group_by(country) %>%
  summarise(ts = list(new_business_density_rate), .groups = "drop") %>%
  pull(ts)

# STEP 3: Run DTW-based K-Medoids clustering
set.seed(123)
kmedoid_dtw <- tsclust(
  ts_list,
  type = "partitional",
  k = 4,
  distance = "dtw_basic",
  centroid = "pam",
  trace = TRUE
)

# STEP 6: Count number of countries in each DTW cluster
dtw_clusters_updated %>%
  count(dtw_cluster) %>%
  arrange(dtw_cluster)


# STEP 4: Plot DTW medoid time series (relative year index)
medoid_series <- kmedoid_dtw@centroids

medoid_df <- map2_dfr(
  medoid_series,
  .y = seq_along(medoid_series),
  .f = function(vec, cluster_id) {
    tibble(
      year_index = 1:length(vec),
      business_density = vec,
      cluster = paste0("Cluster ", cluster_id)
    )
  }
)

ggplot(medoid_df, aes(x = year_index, y = business_density, color = cluster)) +
  geom_line(size = 1.3) +
  scale_color_manual(
    values = c(
      "Cluster 1" = "#e41a1c",
      "Cluster 2" = "#377eb8",
      "Cluster 3" = "#4daf4a",
      "Cluster 4" = "#984ea3"
    )
  ) +
  labs(
    title = "Representative Time Series (Medoids) for Each DTW Cluster",
    x = "Year Index (Relative)",
    y = "New Business Density Rate"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# STEP 5: Map DTW cluster assignments to countries
country_names <- ts_filtered %>%
  group_by(country) %>%
  summarise(n = n(), .groups = "drop") %>%
  pull(country)

dtw_clusters_updated <- tibble(
  country = country_names,
  dtw_cluster = factor(kmedoid_dtw@cluster)
)

world <- ne_countries(scale = "medium", returnclass = "sf")

map_data_dtw_updated <- left_join(world, dtw_clusters_updated, by = c("name" = "country"))

ggplot(map_data_dtw_updated) +
  geom_sf(aes(fill = dtw_cluster), color = "white", size = 0.15) +
  scale_fill_manual(
    name = "DTW Cluster",
    values = c(
      "1" = "#e41a1c",
      "2" = "#377eb8",
      "3" = "#4daf4a",
      "4" = "#984ea3"
    ),
    labels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")
  ) +
  labs(
    title = "Country Clusters by Time Series of Business Density",
    subtitle = "DTW K-Medoids Clustering (≥3-Year Sequences)"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )



# Step 1: Ensure sub_region is a factor
dtw_geo <- dtw_clusters_updated %>%
  left_join(final_data_clean %>% select(country, sub_region) %>% distinct(), by = "country") %>%
  mutate(sub_region = factor(sub_region))  # convert to factor before dropping levels

# Step 2: Plot with dropped levels directly inside aes()
ggplot(dtw_geo, aes(x = dtw_cluster, fill = fct_drop(sub_region))) +
  geom_bar(position = "stack", color = "black") +
  labs(
    title = "Sub-Regional Representation of DTW Clusters",
    x = "DTW Cluster",
    y = "Number of Countries",
    fill = "Sub-Region"
  ) +
  theme_minimal()




#Mean Cluster Analysis

# STEP 6: Add DTW cluster assignments to full dataset
dtw_merged <- final_data_clean %>%
  left_join(dtw_clusters_updated, by = c("country"))

# STEP 7: Select relevant features and standardize
dtw_profile_data <- dtw_merged %>%
  select(dtw_cluster, account_age_15, wefs, estimate, 
         gdp_constant_2015usd, gov_edu_spending_pct_gdp, 
         working_age_pop, population_growth) %>%
  filter(!is.na(dtw_cluster)) %>%
  mutate(dtw_cluster = factor(dtw_cluster)) %>%
  pivot_longer(-dtw_cluster, names_to = "feature", values_to = "value") %>%
  group_by(feature) %>%
  mutate(scaled_value = scale(value)) %>%
  ungroup()

# STEP 8: Calculate cluster-wise average of standardized values
dtw_cluster_summary <- dtw_profile_data %>%
  group_by(dtw_cluster, feature) %>%
  summarise(mean_scaled = mean(scaled_value, na.rm = TRUE), .groups = "drop")

# STEP 9: Plot mean profile per cluster
ggplot(dtw_cluster_summary, aes(x = reorder(feature, mean_scaled), y = mean_scaled, fill = dtw_cluster)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ dtw_cluster, nrow = 1, strip.position = "top") +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "1" = "#e41a1c",
      "2" = "#377eb8",
      "3" = "#4daf4a",
      "4" = "#984ea3"
    )
  ) +
  labs(
    title = "Mean Cluster Profiles — DTW K-Medoids",
    x = NULL,
    y = "Standardized Mean of Feature"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 11),
    panel.spacing = unit(2, "lines"),
    legend.position = "none"
  )

# ========================================================
# 12. DYNAMIC POOLED OLS FORECASTING (LOG SCALE)
# ========================================================
# Commenting out per-country ARIMA/Granger/GARCH due to short T
# … ARIMA/Granger/GARCH sections removed …

# 12a. Build dynamic panel data with lag and log DV
dynamic_data <- final_data %>%
  arrange(country, year) %>%
  group_by(country) %>%
  mutate(
    lag_density = lag(new_business_density_rate),
    log_density = log(new_business_density_rate + 1e-3)
  ) %>%
  ungroup() %>%
  na.omit()

# 12b. Fit  Dynamic Pooled OLS (Log-Linear with Lag) + region FEs
model_log_region <- lm(
  log_density ~ lag_density + wefs + estimate +
    gdp_constant_2015usd + gov_edu_spending_pct_gdp + factor(sub_region),
  data = dynamic_data
)
summary(model_log_region)
coeftest(
  model_log_region,
  vcov = vcovCL(model_log_region, cluster = dynamic_data$country)
)

# 12c. One-year ahead forecasts (back-transformed)
last_obs <- dynamic_data %>%
  group_by(country) %>%
  filter(year == max(year)) %>%
  ungroup()

pred_log <- predict(
  model_log_region,
  newdata = last_obs,
  interval = "prediction",
  level = 0.95
)

forecast_df <- last_obs %>%
  select(country, year) %>%
  bind_cols(as.data.frame(pred_log)) %>%
  rename(log_fit = fit, log_lwr = lwr, log_upr = upr) %>%
  mutate(
    fit = exp(log_fit) - 1e-3,
    lwr = exp(log_lwr) - 1e-3,
    upr = exp(log_upr) - 1e-3
  )

print(forecast_df)

#One-Year-Ahead Forecasts (+-2 SD Filter)

# 1. Z-score & filter extremes
forecast_df_filtered <- forecast_df %>%
  mutate(fit_z = as.numeric(scale(fit))) %>%   # standardize
  filter(abs(fit_z) <= 2) %>%                  # drop |z|>2
  select(-fit_z)                               # we don't need z‐score after filtering

cat("Kept", nrow(forecast_df_filtered), "of", nrow(forecast_df),
    "observations (|z| ≤ 2 filter)\n")

# 2. Basic forest‐plot (no labels, tiny y‐text)
ggplot(forecast_df_filtered, aes(
  x = fit,
  y = reorder(country, fit),
  xmin = lwr,
  xmax = upr
)) +
  geom_errorbarh(height = 0.3, alpha = 0.7) +
  geom_point(size = 2) +
  coord_cartesian(expand = FALSE) +
  scale_color_brewer("Region", palette = "Set2") +
  labs(
    x       = "Forecast New Business Density",
    y       = NULL,
    title   = "One-Year-Ahead Forecasts (±2 SD Filter)",
    subtitle= "Extreme outliers removed",
    caption = "Capstone dynamic pooled-OLS"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 6),
    legend.position = "bottom"
  )



#Geomapping
library(scales)    # for squish()


# 1. Read in world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# 2. Join your forecast and create a clamped version of the fit
map_forecast <- world %>%
  left_join(forecast_df %>% select(country, fit), by = c("name" = "country")) %>%
  mutate(fit_clamped = pmin(fit, 25))   # anything >25 becomes 25

# 3. Plot with a 0–25 color scale, squishing anything outside that range
ggplot(map_forecast) +
  geom_sf(aes(fill = fit_clamped), color = "white", size = 0.1) +
  scale_fill_viridis(
    option   = "plasma",
    name     = "Forecast\nDensity",
    limits   = c(0, 25),            # force the legend to 0–25
    oob      = squish,              # squish values >25 into the top color
    na.value = "grey90"
  ) +
  labs(
    title    = "Global One-Year-Ahead Forecasts of New Business Density",
    subtitle = paste0("Clamped at 0–25 (", max(dynamic_data$year)+1, ")"),
    caption  = "Capstone dynamic pooled-OLS"
  ) +
  theme_minimal() +
  theme(
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    panel.grid   = element_blank(),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 8)
  )

# ===============================
# 13. Interaction Model
# ===============================

int_mod <- lm(new_business_density_rate ~ 
                wefs * account_age_15 + estimate +
                gov_edu_spending_pct_gdp,
              data = ols_data)

tab_model(int_mod,
          dv.labels = "WEFS × Account Interaction")
#WEFS * Financial Inclusion is insignificant (p= 0.737)



# ===============================
# 14. SJPLOT Models Plot
# ===============================

tab_model(
  ols_multi,
  pooling_model,
  re_model,
  fe_model,
  region_fe,
  mixed_model,
  model_log_region,
  step_model,      # <-- no lasso_model
  dv.labels = c(
    "Multivar OLS","Pooled OLS","Random Effects","Fixed Effects",
    "Region Fixed Effects","Mixed Effects","Dynamic Pooled OLS","Stepwise AIC"
  ),
  show.se  = TRUE,
  show.aic = TRUE
)

