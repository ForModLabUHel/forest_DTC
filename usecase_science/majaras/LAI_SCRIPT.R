
library(sf)
library(terra)
library(ncdf4)
library(lubridate)
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)

#1-set directory
base_dir<- "C:/Users/iacom/Desktop/COPERNICUS_LAI"
aoi_path <- file.path(base_dir, "INPUT/AOI/Bounding_Box.shp")
out_dir <- file.path(base_dir, "OUTPUT")  
dir.create(out_dir, showWarnings = FALSE)

#setting with resolution and LAI_scale factors
my_resolution <- 300
LAI_scaling_factor <- 30

#---------------------------------------------------------------------
#1.1-create folders
nc_dir<- file.path(base_dir, "INPUT/nc_files")
setwd(nc_dir)
#LAI
dir.create(file.path(nc_dir, "LAI"), showWarnings = FALSE)
all_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
file_LAI <- all_files[grepl("__LAI\\.nc$", basename(all_files))]
for (file in file_LAI) {
  file.rename(file, file.path(nc_dir, "LAI", basename(file)))
}
lai_dir<- file.path(nc_dir, "INPUT/LAI")

#create LENGTH AFTER folder
dir.create(file.path(nc_dir, "LENGTH_AFTER"), showWarnings = FALSE)
all_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
file_LENGTH_AFTER <- all_files[grepl("__LENGTH_AFTER\\.nc$", basename(all_files))]
for (file in file_LENGTH_AFTER) {
  file.rename(file, file.path(nc_dir, "LENGTH_AFTER", basename(file)))
}
LENGTH_AFTER_dir<- file.path(nc_dir, "INPUT/LENGTH_AFTER")

#create LENGTH BEFORE folder
dir.create(file.path(nc_dir, "LENGTH_BEFORE"), showWarnings = FALSE)
all_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
file_LENGTH_BEFORE <- all_files[grepl("__LENGTH_BEFORE\\.nc$", basename(all_files))]
for (file in file_LENGTH_BEFORE) {
  file.rename(file, file.path(nc_dir, "LENGTH_BEFORE", basename(file)))
}
LENGTH_BEFORE_dir<- file.path(nc_dir, "INPUT/LENGTH_BEFORE")

#create NOBS folder
dir.create(file.path(nc_dir, "NOBS"), showWarnings = FALSE)
all_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
file_NOBS <- all_files[grepl("__NOBS\\.nc$", basename(all_files))]
for (file in file_NOBS) {
  file.rename(file, file.path(nc_dir, "NOBS", basename(file)))
}
NOBS_dir<- file.path(nc_dir, "INPUT/NOBS")

#create QFLAG folder
dir.create(file.path(nc_dir, "QFLAG"), showWarnings = FALSE)
all_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
file_QFLAG <- all_files[grepl("__QFLAG\\.nc$", basename(all_files))]
for (file in file_QFLAG) {
  file.rename(file, file.path(nc_dir, "QFLAG", basename(file)))
}
QFLAG_dir<- file.path(nc_dir, "INPUT/QFLAG")

#create RMSE folder
dir.create(file.path(nc_dir, "RMSE"), showWarnings = FALSE)
all_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
file_RMSE <- all_files[grepl("__RMSE\\.nc$", basename(all_files))]
for (file in file_RMSE) {
  file.rename(file, file.path(nc_dir, "RMSE", basename(file)))
}
RMSE_dir<- file.path(nc_dir, "INPUT/RMSE")

#--------------------------------------------------------------------

#2-Variables folders
var_dirs <- list(
  LAI          = file.path(base_dir, "INPUT/nc_files/LAI"),  # the foldere where you have .nc files
  LENGTH_AFTER = file.path(base_dir, "INPUT/nc_files/LENGTH_AFTER"),
  LENGTH_BEFORE= file.path(base_dir, "INPUT/nc_files/LENGTH_BEFORE"),
  NOBS         = file.path(base_dir, "INPUT/nc_files/NOBS"),
  QFLAG        = file.path(base_dir, "INPUT/nc_files/QFLAG"),
  RMSE         = file.path(base_dir, "INPUT/nc_files/RMSE")
)

#3-move the directory to OUTPUT_LAI where I'll save the results
setwd(out_dir)

#load AOI and convert into vect per terra
#4-Read AOI and take to EPSG:4326
aoi_sf <- st_read(aoi_path, quiet = TRUE)
aoi_sf <- st_transform(aoi_sf, 4326)
aoi_vect <- terra::vect(aoi_sf)

#5-function to extract the date from file name estrarre la data dal nome file (supports 12,8,6 digits)
extract_date_from_filename <- function(fname){
  bn <- basename(as.character(fname))
  m12 <- regmatches(bn, gregexpr("\\d{12}", bn))[[1]]
  m8 <- regmatches(bn, gregexpr("\\d{8}", bn))[[1]]
  m6 <- regmatches(bn, gregexpr("\\d{6}", bn))[[1]]
  if (length(m12) >0 && nzchar(m12[[1]])) return(as.Date(substr(m12[1],1,8), "%Y%m%d"))
  if (length(m8) > 0 && nzchar(m8[1])) return(as.Date(m8[1], "%Y%m%d"))
  if (length(m6) > 0 && nzchar(m6[1])) return(as.Date(paste0(m6[1], "01"), "%Y%m%d"))
  return(NA)
}

#6-Generic function to read one NetCDF variable
process_variable_nc_to_dt <- function(nc_path, varname, aoi_vect = NULL){
  fname <- basename(nc_path)
  date_val <- extract_date_from_filename(fname)
  
  nc <- try(nc_open(nc_path), silent = TRUE)
  if(inherits(nc, "try-error")) return(NULL)
  #if varname isn't exact, try case-intensitive to find variable into file
  vars <- names(nc$var)
  alt <- vars[grepl(varname, vars, ignore.case = TRUE)]
  if(length(alt) == 0) { nc_close(nc); 
    cat(" - WARNING:Variable '", varname, "' not found in ", basename(nc_path), "\n", sep = "") 
    return(NULL)}
  var_actual <- alt[1] #use the variable name found in the file
  
  scale_att <- try(ncatt_get(nc, var_actual, "scale_factor"), silent = TRUE)
  scale_factor <- if(is.list(scale_att) && !is.null(scale_att$value)) as.numeric(scale_att$value) else NA
  fill_att <- try(ncatt_get(nc, var_actual, "_FillValue"), silent = TRUE)
  fillv <- if(is.list(fill_att) && !is.null(fill_att$value)) as.numeric(fill_att$value) else NA
  nc_close(nc)
  #read raster
  r <- try(rast(nc_path, subds = var_actual), silent = TRUE)
  if(inherits(r, "try-error")) {
    r <- try(rast(nc_path), silent = TRUE)
    if(inherits(r, "try-error")) return(NULL)
  }
  #apply scale e fill
  if(!is.na(scale_factor)) r <- r * scale_factor
  if(!is.na(fillv)) r[r == fillv] <- NA
  #if date isn't extract from name, try to read time from the object SpatRaster
  if(is.na(date_val)) {
    tt <- time(r)
    if(length(tt) >= 1 && !is.na(tt[1])) date_val <- as.Date(tt[1])
  }
  #crop/mask to the AOI 
  if(!is.null(aoi_vect)) {
    r <- crop(r, aoi_vect, snap = "out")
    r <- mask(r, aoi_vect)
  }
  #convert into data.frame and remove NA
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  val_col <- setdiff(names(df), c("x", "y"))[1]
  setDT(df)
  dt <- df[, .(x = get("x"), y = get("y"), value = get(val_col))]
  dt <- dt[!is.na(value)]
  if(nrow(dt) == 0) return(NULL)
  dt[, date := date_val]
  dt[, variable := varname] #it retains the canonic variable name
  return(dt)
}
#7-loop on variables and files, save chunk.rds for each file to avoid OOM
chunk_files <- character()
for(var in names(var_dirs)){
  folder <- var_dirs[[var]]
  if(!dir.exists(folder)) next
  files <- list.files(folder, pattern = "\\.nc$", full.names = TRUE)
  #Test A: if you want testing only 1 file per variable, 
  for(f in files) {
    cat("Processing", var, "->", basename(f), "...\n")
    dt <- process_variable_nc_to_dt(f, var, aoi_vect)
    if(is.null(dt)) { cat(" - SKIPPED:", basename(f), "(No data or variable not found)\n"); next}
    #save chunk for file
    rds_name <- file.path(out_dir, paste0(var, "_", tools::file_path_sans_ext(basename(f)), ".rds"))
    saveRDS(dt, rds_name, compress = "xz") 
    chunk_files <- c(chunk_files, rds_name)
    cat(" - SAVED:", basename(rds_name), "\n")
    # free memory
    rm(dt); gc()
  } 
}

#8-merge all chunk in a unique table (data.table)
cat("Merging all chunk into a single data.table...\n")
dt_list <- lapply(chunk_files, readRDS)
big_dt <- data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE)
cat("Merge complete. Total rows:", nrow(big_dt), "\n")

#9-save output
cat("Saving final output...\n")
#check the NOBS, QFLAG, RMSE are including
vars_in_output <- unique(big_dt$variable)
cat("Variables included in final table:", paste(vars_in_output, collapse = ", "), "\n")
if(all(c("NOBS", "QFLAG", "RMSE") %in% vars_in_output)) {
  cat("SUCCESS: NOBS, QFLAG, and RMSE chunk were successfully processed and merged.\n")
} else {
  cat("WARNING: Not all variables (NOBS, QFLAG, RMSE) are present in the output table.\n")
}

saveRDS(big_dt, file = file.path(out_dir, "combined_timeseries.rds"), compress = "xz")
data.table::fwrite(big_dt, file = file.path(out_dir, "combined_timeseries.csv"))
cat("Final files saved in:", out_dir, "\n")

#10-load unique table with all the variables
raw_dt <- readRDS("C:/Users/iacom/Desktop/COPERNICUS_LAI/OUTPUT/combined_timeseries.rds")
#first visualization
cat("Names in big_dt:", names(raw_dt), "\n")
cat("Variables present:", paste(unique(raw_dt$variable), collapse = ", "), "\n")

#11-Divide variables for join
lai_dt <- raw_dt[variable == "LAI", .(x, y, date, LAI_value = value)]
qflag_dt <- raw_dt [variable == "QFLAG", .(x, y, date, QFLAG = value)]
rmse_dt <- raw_dt[variable == "RMSE", .(x, y, date, RMSE = value)]
nobs_dt <- raw_dt[variable == "NOBS", .(x, y, date, NOBS = value)]
length_before_dt <- raw_dt[variable == "LENGTH_BEFORE", .(x, y, date, LENGTH_BEFORE = value)]
length_after_dt <- raw_dt[variable == "LENGTH_AFTER", .(x, y, date, LENGTH_AFTER = value)]

#remove duplicates
lai_dt <- unique(lai_dt, by = c("x", "y", "date"))
qflag_dt <- unique(qflag_dt, by = c("x", "y", "date"))
rmse_dt <- unique(rmse_dt, by = c("x", "y", "date"))
nobs_dt <- unique(nobs_dt, by = c("x", "y", "date"))
length_before_dt <- unique(length_before_dt, by = c("x", "y", "date"))
length_after_dt <- unique(length_after_dt, by = c("x", "y", "date"))
cat(sprintf("   - Righe LAI (unique): %d\n", nrow(lai_dt)))

#merge LAI con the other variables
lai_combined_dt <- merge(lai_dt, qflag_dt, by = c("x", "y", "date"), all.x = TRUE)
lai_combined_dt <- merge(lai_combined_dt, rmse_dt, by = c("x", "y", "date"), all.x = TRUE)
lai_combined_dt <- merge(lai_combined_dt, nobs_dt, by = c("x", "y", "date"), all.x = TRUE)
lai_combined_dt <- merge(lai_combined_dt, length_before_dt, by = c("x", "y", "date"), all.x = TRUE)
lai_combined_dt <- merge(lai_combined_dt, length_after_dt, by = c("x", "y", "date"), all.x = TRUE)

cat("Righe totali dopo il merge:", nrow(lai_combined_dt), "\n")
cat("Valori unici di QFLAG:", paste(unique(lai_combined_dt$QFLAG), collapse = ", "), "\n")
cat("Min/Max NOBS:", min(lai_combined_dt$NOBS, na.rm = TRUE), "/", max(lai_combined_dt$NOBS, na.rm = TRUE), "\n")
cat("Min/Max RMSE:", min(lai_combined_dt$RMSE, na.rm = TRUE), "/", max(lai_combined_dt$RMSE, na.rm = TRUE), "\n")

#apply quality filter
lai_combined_dt_filtered <- lai_combined_dt[
  QFLAG %in% c(0, 1) &
    LAI_value <= 6 & 
    RMSE < 1,
]

#rename value
lai_dt <- lai_combined_dt_filtered
#scaling
lai_dt[, value := 30*LAI_value]
#ID
lai_dt$ID <- match(paste(lai_dt$x, lai_dt$y), unique(paste(lai_dt$x, lai_dt$y)))
lai_dt[, LAI_value := NULL]
lai_dt[, variable := "LAI"]

#12-plot
lai_plot <- ggplot(lai_dt, aes(x = date, y = value, group = ID, color = factor(ID))) +
  geom_line(linewidth = 0.8) +
  labs(title = "Temporal series filtered LAI (per Pixel)",
       x = "Data",
       y = "scaled LAI",
       color = "ID pixel") +
  theme_minimal()

lai_plot
#graphic corrected
lai_dt_filtered <- lai_dt[value<6]
lai_plot_filtered <- ggplot(lai_dt_filtered, aes(x = date, y = value, group = ID, color = factor(ID))) +
  geom_line(size = 0.8) +
  labs(title = "Temporal series filtered LAI (per Pixel)",
       x = "Data",
       y = "scaled LAI",
       color = "ID pixel") +
  theme_minimal()

lai_plot_filtered

#13-Daily interpolation
#convert date
lai_dt_filtered[, date := as.Date(date)]
#order and remove final duplicates per ID/date
setorder(lai_dt_filtered, ID, date)
lai_dt_filtered <- lai_dt_filtered[, .SD[!duplicated(date)], by = ID]
# Create a daily date range for each ID
daily_dates <- lai_dt_filtered[, .(date = seq(min(date), max(date), by = "1 day")), by = ID]

#option 1 Interpolate LAI for each ID with approxfun
daily_lai <- daily_dates[, {
  od <- lai_dt_filtered[ID == .BY$ID]
  if(nrow(od) < 2) {
    .(date = date, value = NA_real_)
  } else {
    f <- approxfun(x = as.numeric(od$date), y = od$value, rule = 1)
    .(date = date, value = f(as.numeric(date)))
  }
}, by = ID]

#add coordinates
coords <- unique(lai_dt_filtered[, .(ID, x, y)])
daily_lai <- merge(daily_lai, coords, by = "ID", all.x = TRUE)

#14-View result
lai_plot_daily <- ggplot(daily_lai, aes(x = date, y = value, group = ID, color = factor(ID))) +
  geom_line(linewidth = 0.8) +
  labs(title = "Time series interpolate LAI (daily per Pixel)",
       x = "Date",
       y = "LAI (interpolate)",
       color = "ID Pixel") +
  theme_minimal()

lai_plot_daily



#TIME GRAPH (trend time of LAI)
lai_time <- daily_lai %>%
  group_by(date) %>%
  summarise(mean_LAI=mean(value, na.rm = TRUE))
ggplot(lai_time, aes (x = date, y = mean_LAI)) +
  geom_line(color="forestgreen") +
  geom_point(color="darkgreen") +
  labs(title = "Time Trend of LAI",
       x = "Year",
       y = "Mean LAI") +
  theme_minimal()
#it shows how LAI changes over time over area

#seasonal or annual Boxplot
lai_box <- daily_lai %>%
  mutate(year = format(date, "%Y"),
         month = format(date, "%m"))

ggplot(lai_box, aes(x = month, y = value, group = month)) +
  geom_boxplot(fill = "lightgreen", color = "darkgreen") +
  labs(title = "Monthly mean LAI values",
       x = "Month",
       y = "LAI") +
  theme_minimal()


# Find the closest pixel(ID):
# Target coordinates
target_lat <- 39.94033
target_lon <- -5.77465

# Compute distance and find closest ID
closest <- lai_dt_filtered[, .(ID, x, y,
                               dist = sqrt((x - target_lon)^2 + (y - target_lat)^2))][
                                 order(dist)][1]

# View result
print(closest)