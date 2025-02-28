library(sf)         # Spatial vector data handling  
library(amt)        # Animal movement analysis  
library(dplyr)      # Data manipulation and transformation  
library(raster)     # Raster data processing  
library(spatstat)   # Spatial point pattern analysis  
library(sp)         # Spatial data classes and methods  
library(terra)      # Raster and vector processing  
library(tibble)     # Improved data frames handling  
library(ggplot2)    # Data visualization and plotting  
library(survival)   # Survival analysis modeling  
library(lubridate)  # Date and time manipulation  


#Define work directory
setwd("Y:/Home/esguerrl/Temp")
# Load northern lapwing positions
position <- st_read("Y:/Home/esguerrl/locations_breeding/locations_points.shp")
# Load raster layers: Land Use /Land cover (LULC), Diversity Shannon Index, Mean Patch Area Index, and Contagion Index
landcover <- raster("Y:/Home/esguerrl/dw_2022_mosaic/LandcoverDW2022.tif")
shannon <- raster("Y:/Home/esguerrl/Indices_Final/Shannon_Total.tif")
mean_patch_area <- raster("Y:/Home/esguerrl/Indices_Final/MPA_Resampled.tif")
contagio <- raster("Y:/Home/esguerrl/Indices_Final/Contagion_Resampled.tif")
elevation <- raster("Y:/Home/esguerrl/Indices_Final/dem30m_merged.tif")

# Convert the required data types of the fields x and y coordinates, timestamp (t) and individu-al id (id)
position$x <- as.numeric(position$x)
position$y <- as.numeric(position$y)
position$t <- as.POSIXct(position$t, format="%Y-%m-%d %H:%M:%S")
position$id <- as.factor(position$id)
position$geometry <- NULL  # Delete the field geometry to avoid errors during track creation

# Filter positions for year. It's necessary to set the year of analysis
position <- position %>% filter(year(t) == 2021)
year_filter <- unique(year(position$t))
cat("Filtered positions for the year... ", year_filter,"\n")
# Validate the number of positions remaining after filtering
cat("Number of positions for the year ", year_filter,": ", nrow(position), "\n")
# Unique individual IDs
individuos <- unique(position$id)
head(individuos)

#Uncomment the rates and tolerances for the year required for the analysis.
############## Define rates, tolerance per individual year 2021###################
config <- data.frame(
  individuo = c("5181", "5264", "5268", "5270", "5277", "5278", "5282", "5284", "5285", "5286", "5301", "5303", "5308", "5439"),
  resample_rate = c(7, 10, 2, 2, 2, 2, 2, 2, 30, 2, 2, 2, 2, 2),  # Resample time in minutes
  resample_tolerance = c(42, 60, 12, 12, 12, 12, 12, 12, 180, 12, 12, 12, 12, 12)  # Resample tolerance in seconds
)

############## Define rates, tolerance per individual year 2022###################
#config <- data.frame(
#  individuo = c("5181", "5264", "5285", "5301", "5303", "5439" ),
#  resample_rate = c(60, 2, 2, 2, 2, 10),  # Resample time in minutes
#  resample_tolerance = c(360, 10, 10, 10, 10, 60)      # Resample tolerance in seconds
#)
# Loop to create subsets, tracks, bursts, and steps for each individual
for (individuo in individuos) {
  cat("\nProcessing individual: ", individuo, "\n")
  
  # Set the individual-specific resampling parameters according to the data frame “config”.
  rate <- config$resample_rate[config$individuo == individuo]
  tolerance <- config$resample_tolerance[config$individuo == individuo]
  cat("Rate (resample): ", rate, "minutes\n")
  cat("Tolerance (resample): ", tolerance, "seconds\n")

  # Filter data for the current individual and store it locally in position_ind_xxxx.rds
  cat("Filtering data for individual: ", individuo, "\n")
  position_ind <- position %>% 
    filter(id == individuo) %>% 
    filter(!is.na(t))  # Remove row with NA in the timestamps column  `t`
  saveRDS(position_ind, file = paste0("position_ind_", individuo, "_", year_filter, ".rds"))
  cat("Number of points for individual ", individuo, ": ", nrow(position_ind), "\n")
  rm(position_ind)
  gc()
  
  # From the individual data stored locally (rds file), create the tracks and bursts
  track_ind <- readRDS(paste0("position_ind_", individuo, "_", year_filter, ".rds")) %>%
    # Creating track (track_xyt format)
    make_track(.x = x, .y = y, .t = t, id = id, crs = 3035) %>%
    # Resample data according to the config per individual
    track_resample(rate = minutes(rate), tolerance = seconds(tolerance))
   head(track_ind)
  # Save the data in local directory 
  saveRDS(track_ind, file = paste0("track_ind_", individuo, "_", year_filter, ".rds"))
  rm(track_ind)
  gc()
  
  
  # Adjust overlapped consecutive positions (step length=0), adding 0.01 m to the (x,y) coordi-nates.
  # Stored as bursts "burst_ind_xxxx.rds".
  burst_ind <- readRDS(paste0("track_ind_", individuo,  "_", year_filter, ".rds")) %>%
    mutate(
      x_ = ifelse(lead(x_) == x_ & lead(y_) == y_, x_ + 0.01, x_),
      y_ = ifelse(lead(x_) == x_ & lead(y_) == y_, y_ + 0.01, y_)
    ) 
  print(head(burst_ind))  # Debugging: Print first rows
  saveRDS(burst_ind, file = paste0("burst_ind_", individuo,  "_", year_filter, ".rds"))
  rm(burst_ind)
  gc()
  
  # Create used and available steps and store them.
  cat("Creating used and available steps for individual: ", individuo, "\n")
  burst_ind <- readRDS(paste0("burst_ind_", individuo,  "_", year_filter, ".rds"))
  # Validar si hay suficientes bursts para procesar
  if (n_distinct(burst_ind$burst_) < 2) {
    cat("Skipping individual", individuo, "as there is only one burst.\n")
    next  # Saltar al siguiente individuo
  }
    # Delete bursts with less than 3 positions, otherwise it will show an error.
  burst_ind <- burst_ind %>%
    group_by(burst_) %>%
    filter(n() >= 3) %>%
    ungroup()
  
    # Debugging: Print the structure and unique values of burst_
  print("Structure of burst_:")
  print(str(burst_ind$burst_))
  print("Unique values of burst_:")
  print(unique(burst_ind$burst_))
  burst_ind$burst_ <- as.integer(burst_ind$burst_)  # Ensure burst_ is treated as an integer
    
  #From tracks and bursts, create the random steps. 10 random step by 1 used.
  steps_ind <- burst_ind %>%
    steps_by_burst(burst_ = burst_) %>%
    random_steps(n = 10) %>%
    group_by(step_id_) %>%
    mutate(
      log_sl_ = log(sl_),
      cos_ta_ = cos(ta_)
    ) %>% 
    ungroup()  %>% 
    filter(!is.na(x2_) & !is.na(y2_) & !is.na(log_sl_) & !is.na(cos_ta_)) 
  saveRDS(steps_ind, file = paste0("steps_ind_", individuo,  "_", year_filter, ".rds"))
  rm(steps_ind)
  gc()
  
  # Convert steps to sf (simple feature) for raster extraction and store it
  cat("Converting steps to sf for individual: ", individuo, "\n")
  steps_sf <- readRDS(paste0("steps_ind_", individuo,  "_", year_filter, ".rds")) %>%
    filter(!is.na(x2_) & !is.na(y2_)) %>%  # Remove rows with "na" missing coordinates before conversion
    st_as_sf(coords = c("x2_", "y2_"), crs = 3035, remove = FALSE) 
# remove = FALSE to keep the x2 and y2 columns in the original dataframe
  
  # Extract raster values for each step
  cat("Extracting raster values for individual: ", individuo, "\n")
  steps_sf$landcover <- raster::extract(landcover, steps_sf)
  steps_sf$shannon <- raster::extract(shannon, steps_sf)
  steps_sf$mean_patch_area <- raster::extract(mean_patch_area, steps_sf)
  steps_sf$contagio <- raster::extract(contagio, steps_sf)
  steps_sf$elevation <- raster::extract(elevation, steps_sf)
  
    cat("Reclassifying landcover and removing NAs for individual: ", individuo, "\n")
  # Specify the landcover classes of interest as numerical codes, where 0: water; 1: trees; 2: grass; 3: flooded vegetation; 4: crops, 5: Shrub & Scrub, 6: Built area. The class 6 is used as a reference level for the other classes
  selected_classes <- c(6, 0, 1, 2, 3, 4, 5)
    # Filter and clean the data for the selected landcover classes
  steps_sf <- steps_sf %>%
    filter(landcover %in% selected_classes) %>%  # Retain only rows where landcover is in the selected classes
    mutate(
      # Reclassify landcover into a factor with meaningful labels: "Class_0", "Class_1", etc.
      landcover = factor(landcover, levels = selected_classes, labels = paste0("Class_", select-ed_classes)) # convert landcover as a categorical variable type factor
    ) %>%
    na.omit()  # Remove rows with missing values in any column
    # Save the processed data for the individual to an rds file steps_sf_xxxx.rds
  saveRDS(steps_sf, file = paste0("steps_sf_", individuo,  "_", year_filter, ".rds"))
  rm(steps_sf)
  gc()
 }

# Loop to fit multiple SSF models for each individual
for (individuo in individuos) {
    # Get the steps for the current individual
  steps_sf <- readRDS(paste0("steps_sf_", individuo,  "_", year_filter, ".rds"))
  print(head(steps_sf))
    
  #  Fit a SSF model for each landscape metric
  models <- list(
    model_1 = steps_sf %>% fit_issf(case_ ~ landcover + shannon + elevation +  log_sl_ + cos_ta_ + strata(step_id_), model = TRUE),
    model_2 = steps_sf %>% fit_issf(case_ ~ landcover + mean_patch_area + elevation + log_sl_ + cos_ta_ + strata(step_id_), model = TRUE),
    model_3 = steps_sf %>% fit_issf(case_ ~ landcover + contagio + elevation + log_sl_ + cos_ta_ + strata(step_id_), model = TRUE)
  )
  
  # Store each fitted model
  for (model_name in names(models)) {
    saveRDS(models[[model_name]], file = paste0(model_name, "_", individuo,  "_", year_filter, ".rds"))
  }
  rm(models)
  gc()
}
# Summarize the models for each individual
for (individuo in individuos) {
  for (model_name in c("model_1", "model_2", "model_3")) {
    cat("\nSummary for Individual: ", individuo, " - Model: ", model_name, "\n")
    model_ind <- readRDS(paste0(model_name, "_", individuo,  "_", year_filter, ".rds"))
    print(summary(model_ind))
    rm(model_ind)
    gc()
  }
}
