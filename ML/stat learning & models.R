
# ---- Install and load packages ----

pkgs <- c(
  "sf",
  "terra",
  "dplyr",
  "future",
  "lgr",
  "mlr3",
  "mlr3learners",
  "mlr3extralearners",
  "mlr3proba",
  "mlr3spatiotempcv",
  "mlr3tuning",
  "mlr3viz",
  "progressr",
  "pROC",
  "spDataLarge"
)

# Install missing packages
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(to_install) > 0) {
  install.packages(
    to_install,
    repos = "https://packagemanager.posit.co/cran/latest",
    Ncpus = parallel::detectCores()
  )
}
install.packages("spDataLarge", repos = "https://geocompr.r-universe.dev")

# Load packages
invisible(
  lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  })
)


# ---- Introduction ----

data("lsl", "study_mask", package = "spDataLarge")
ta = terra::rast(system.file("raster/ta.tif", package = "spDataLarge"))
library(terra)
library(sf)
library(ggplot2)


# 1. Extract elevation raster and mask to study area

elev_rast <- ta[["elev"]]        # Use the correct layer
elev_masked <- mask(elev_rast, study_mask)


# 2. Create hillshade

slope <- terrain(elev_masked, v = "slope", unit = "radians")
aspect <- terrain(elev_masked, v = "aspect", unit = "radians")
hill <- shade(slope, aspect, angle = 45, direction = 315)


# 3. Convert hillshade raster to data.frame for ggplot

hill_df <- as.data.frame(hill, xy = TRUE)
colnames(hill_df) <- c("x", "y", "hillshade")


# 4. Convert landslide points to sf

lsl_sf <- st_as_sf(lsl, coords = c("x", "y"), crs = st_crs(study_mask))


# 5. Plot

ggplot() +
  # Hillshade background
  geom_raster(data = hill_df, aes(x = x, y = y, fill = hillshade)) +
  scale_fill_gradient(low = "black", high = "white") +
  
  # Study area outline
  geom_sf(data = study_mask, fill = NA, color = "black", size = 0.5) +
  
  # Landslide points
  geom_sf(data = lsl_sf, aes(color = lslpts), size = 2) +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red")) +
  
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Landslide Occurrence on Elevation Map",
    fill = "Hillshade",
    color = "Landslide"
  )

head(lsl)



# ---- Conventional modelling approach in R ----

fit = glm(lslpts ~ slope + cplan + cprof + elev + log10_carea,
          family = binomial(),
          data = lsl)
summary(fit)

library(effects)
plot(allEffects(fit))


pred_glm = predict(object = fit, type = "response")
head(pred_glm)

# making the prediction
pred = terra::predict(ta, model = fit, type = "response")

library(terra)
library(sf)
library(ggplot2)
library(dplyr)


# 1. Mask to study area

pred_masked <- mask(pred, study_mask)


# 2. Classify into 5 categories (0-0.2, 0.2-0.4, ..., 0.8-1)

breaks <- seq(0, 1, by = 0.2)
labels <- c("Very Low", "Low", "Moderate", "High", "Very High")
pred_class <- classify(pred_masked, rcl = cbind(breaks[-length(breaks)], breaks[-1], 1:5))
levels(pred_class) <- data.frame(id = 1:5, class = labels)


# 3. Convert to data.frame for ggplot

pred_df <- as.data.frame(pred_class, xy = TRUE)
colnames(pred_df) <- c("x", "y", "class")
pred_df$class <- factor(pred_df$class, levels = labels)


# 4. Plot with ggplot2

ggplot() +
  geom_raster(data = pred_df, aes(x = x, y = y, fill = class)) +
  scale_fill_manual(
    values = c("Very Low" = "#ffffcc",
               "Low"      = "#a1dab4",
               "Moderate" = "#41b6c4",
               "High"     = "#2c7fb8",
               "Very High"= "#253494"),
    name = "Landslide Susceptibility"
  ) +
  geom_sf(data = study_mask, fill = NA, color = "black", size = 0.5) +
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Landslide Susceptibility Map (5-class scheme)",
    x = "Easting (m)",
    y = "Northing (m)"
  )

pROC::auc(pROC::roc(lsl$lslpts, fitted(fit)))




# --- Simple model evaluation
set.seed(1)

n <- nrow(lsl)
id_train <- sample(seq_len(n), size = 0.7 * n)

train <- lsl[id_train, ]
test  <- lsl[-id_train, ]

fit <- glm(
  lslpts ~ slope + cplan + cprof + elev + log10_carea,
  family = binomial(),
  data = train
)

summary(fit)
p_test <- predict(fit, newdata = test, type = "response")
library(pROC)

auc_test <- pROC::auc(pROC::roc(test$lslpts, p_test, quiet = TRUE))
auc_test
p_train <- fitted(fit)

auc_train <- pROC::auc(pROC::roc(train$lslpts, p_train, quiet = TRUE))

auc_train

#We train the model on 70 % of the data and evaluate its predictive performance on the remaining 30 %, which mimics prediction at new, unseen locations.


# ---- CV ----

#Define evaluation function based on AUC
auc01 <- function(y, p) {
  pROC::auc(pROC::roc(y, p, quiet = TRUE))
}

# --- Your formula once ---
form <- lslpts ~ slope + cplan + cprof + elev + log10_carea

set.seed(1)
k <- 5
fold_id <- sample(rep(1:k, length.out = nrow(lsl)))

auc_classic <- numeric(k)

for (i in 1:k) {
  train <- lsl[fold_id != i, ]
  test  <- lsl[fold_id == i, ]
  
  fit <- glm(form,
             family = binomial(), data = train)
  
  p <- predict(fit, newdata = test, type = "response")
  auc_kfoldcv[i] <- auc01(test$lslpts, p)
}

auc_kfoldcv
mean(auc_kfoldcv)

# LOOCV

p_loocv <- numeric(nrow(lsl))

for (i in 1:nrow(lsl)) {
  fit <- glm(form, family = binomial(), data = lsl[-i, ])
  p_loocv[i] <- predict(fit, newdata = lsl[i, ], type = "response")
}

auc_loocv <- pROC::auc(pROC::roc(lsl$lslpts, p_loocv))



# ---- Tidymodel ----
install.packages("tidymodels")
library(tidymodels)
set.seed(234589)
# split the data into trainng (75%) and testing (25%)
lsl_split <- initial_split(lsl, prop = 3/4)
lsl_split

# extract training and testing sets
lsl_train <- training(lsl_split)
lsl_test <- testing(lsl_split)

# create CV object from training data
# define the recipe
lsl_recipe <- 
  # which consists of the formula (outcome ~ predictors)
  recipe(lslpts ~ slope + cplan + cprof + elev + log10_carea, 
         data = lsl) %>%
  # and some pre-processing steps
  step_normalize(all_numeric()) %>% # normalise all the numeric predictors
  step_nzv(all_predictors())  # remove any near zero variance features

lr_model <- 
  # specify that the model is a logistic regression
  logistic_reg() %>%
  # select the engine/package that underlies the model
  set_engine("glm") %>%
  # choose either the continuous regression or binary classification mode
  set_mode("classification") 

# set the workflow
lr_workflow <- workflow() %>%
  # add the recipe
  add_recipe(lsl_recipe) %>%
  # add the model
  add_model(lr_model)

lr_fit <- lr_workflow %>%
  # fit on the training set and evaluate on test set
  last_fit(lsl_split)
test_performance <- lr_fit %>% collect_metrics()
test_performance

# --- CV folds on training set ---
set.seed(1)
ls_cv <- vfold_cv(lsl_train, v = 5, strata = lslpts)

# --- CV resampling fit ---
lr_cv_fit <- lr_workflow %>%
  fit_resamples(
    resamples = ls_cv,
    metrics   = metric_set(roc_auc),     # AUC like your pROC
    control   = control_resamples(save_pred = TRUE)
  )

# --- mean AUC across folds (+ SE) ---
lr_cv_fit %>% collect_metrics()

final_fit <- lr_workflow %>%
  fit(lsl_train)

glm_fit <- extract_fit_engine(final_fit)
pred <- terra::predict(
  terra::scale(ta), #in not scaled its wrong prediction because we pre-processed the data
  model = glm_fit,
  type  = "response"
)

plot(pred)


# ---- Spatial CV