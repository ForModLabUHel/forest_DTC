###############################################################################
# Script: process_local_nc_lai.R
# Objective: process file .nc locali (CGLS / PROBA-V LAI) - crop on AOI,
#           apply scale/fill, save GeoTIFF and create metadata for QGIS.
# Packages: sf, terra, ncdf4, lubridate, dplyr
###############################################################################

library(sf)
library(terra)
library(ncdf4)
library(lubridate)
library(dplyr)

#set directory
base_dir<- "C:/Users/iacom/Desktop/COPERNICUS_LAI"
setwd(base_dir)

aoi_path <- file.path(base_dir, "INPUT/AOI/Bounding_Box.shp")
out_dir <- file.path(base_dir, "OUTPUT")  # the folder where you save GeoTIFF and metadata
dir.create(out_dir, showWarnings = FALSE)

#setting with resolution and LAI_scale factors
my_resolution <- 300
LAI_scaling_factor <- 30

#---------------------------------------------------------------------
#create a LAI folder
nc_dir<- file.path(base_dir, "INPUT/nc_files")
setwd(nc_dir)

dir.create(file.path(nc_dir, "LAI"), showWarnings = FALSE)
all_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
file_LAI <- all_files[grepl("__LAI\\.nc$", basename(all_files))]
for (file in file_LAI) {
  file.rename(file, file.path(nc_dir, "LAI", basename(file)))
}
lai_dir<- file.path(nc_dir, "INPUT/LAI")

#create LENGTH AFTER folder
setwd(nc_dir)
dir.create(file.path(nc_dir, "LENGTH_AFTER"), showWarnings = FALSE)
all_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
file_LENGTH_AFTER <- all_files[grepl("__LENGTH_AFTER\\.nc$", basename(all_files))]
for (file in file_LENGTH_AFTER) {
  file.rename(file, file.path(nc_dir, "LENGTH_AFTER", basename(file)))
}
LENGTH_AFTER_dir<- file.path(nc_dir, "INPUT/LENGTH_AFTER")

#create LENGTH BEFORE folder
setwd(nc_dir)
dir.create(file.path(nc_dir, "LENGTH_BEFORE"), showWarnings = FALSE)
all_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
file_LENGTH_BEFORE <- all_files[grepl("__LENGTH_BEFORE\\.nc$", basename(all_files))]
for (file in file_LENGTH_BEFORE) {
  file.rename(file, file.path(nc_dir, "LENGTH_BEFORE", basename(file)))
}
LENGTH_BEFORE_dir<- file.path(nc_dir, "INPUT/LENGTH_BEFORE")

#create NOBS folder
setwd(nc_dir)
dir.create(file.path(nc_dir, "NOBS"), showWarnings = FALSE)
all_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
file_NOBS <- all_files[grepl("__NOBS\\.nc$", basename(all_files))]
for (file in file_NOBS) {
  file.rename(file, file.path(nc_dir, "NOBS", basename(file)))
}
NOBS_dir<- file.path(nc_dir, "INPUT/NOBS")

#create QFLAG folder
setwd(nc_dir)
dir.create(file.path(nc_dir, "QFLAG"), showWarnings = FALSE)
all_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
file_QFLAG <- all_files[grepl("__QFLAG\\.nc$", basename(all_files))]
for (file in file_QFLAG) {
  file.rename(file, file.path(nc_dir, "QFLAG", basename(file)))
}
QFLAG_dir<- file.path(nc_dir, "INPUT/QFLAG")

#create RMSE folder
setwd(nc_dir)
dir.create(file.path(nc_dir, "RMSE"), showWarnings = FALSE)
all_files <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)
file_RMSE <- all_files[grepl("__RMSE\\.nc$", basename(all_files))]
for (file in file_RMSE) {
  file.rename(file, file.path(nc_dir, "RMSE", basename(file)))
}
RMSE_dir<- file.path(nc_dir, "INPUT/RMSE")

#--------------------------------------------------------------------

#Variables folders
var_dirs <- list(
  LAI          = file.path(base_dir, "INPUT/nc_files/LAI"),  # the foldere where you have .nc files
  LENGTH_AFTER = file.path(base_dir, "INPUT/nc_files/LENGTH_AFTER"),
  LENGTH_BEFORE= file.path(base_dir, "INPUT/nc_files/LENGTH_BEFORE"),
  NOBS         = file.path(base_dir, "INPUT/nc_files/NOBS"),
  QFLAG        = file.path(base_dir, "INPUT/nc_files/QFLAG"),
  RMSE         = file.path(base_dir, "INPUT/nc_files/RMSE")
)
#if you don't have it, create output folder
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#move the directory to OUTPUT_LAI where I'll save the results
setwd(out_dir)

library(data.table)

# Read AOI and take to EPSG:4326
aoi_sf <- st_read(aoi_path, quiet = TRUE)
aoi_sf <- st_transform(aoi_sf, 4326)
aoi_vect <- terra::vect(aoi_sf)

##################################################################################################
# From .nc files, now I have to crate a .RData file with latitude, longitude, dates and LAI values
# It'll produce a very large table (n pixel x n file).
###########################################################################################

#WHAT THE SCRIP DOES
#Search all .nc in nc_dir and for each file:
#Reads the LAI variable with "terra", applies "scale_factor" and "_FillValue" if presents,
#trasforms the raster into a long table with columns x(lon), y(lat) and LAI, adds the "date"column
#extracts from the file name, it saves eache partial table as .rds (optional) to keep RAM low,
#finally it merge all tables (if possible) and saves a single .RData file or a compressed .rds

###############################################################################
# Script: nc_to_RData.R
# Objective: extract latitude, longitude, date e LAI from local .nc and save in .RData
# Dipendenze: terra, ncdf4, data.table
###############################################################################

# -------------------------- PARAMETERS to change -------------------------
#the directory of nc and out has already been set
out_dir <- file.path(base_dir, "OUTPUT")  # where save GeoTIFF, metadata and .rds / .RData
lai_dir<- file.path(base_dir, "INPUT/nc_files/LAI")  # the folder where there are .nc files

save_intermediate_rds <- TRUE                  # save every files like .rds (reduces RAM)
combine_and_save_RData <- TRUE                 # finally merge everything and save .RData
rds_chunk_prefix <- "lai_chunk_"               # prefx for intermediate .rds
final_RData_name <- file.path(out_dir, "LAI_pixels_YYYY_10y.RData") # change name
# If you want to restrict to AOI (if you have a shapefile), set "aoi_path"; if not leave NA
aoi_path <- NA # e.g. "C:/Progetti_LAI/input/aoi/aoi.shp" oppure NA
# ---------------------------------------------------------------------------

# function to extract the date from file name estrarre la data dal nome file (supports 12,8,6 digits)
extract_date_from_filename <- function(fname){
  m12 <- regmatches(fname, regexpr("\\d{12}", fname))
  if(length(m12) == 1 && nchar(m12) == 12){
    return(as.Date(substr(m12, 1, 8), format = "%Y%m%d"))
  }
  m8 <- regmatches(fname, regexpr("\\d{8}", fname))
  if(length(m8) == 1 && nchar(m8) == 8){
    return(as.Date(m8, format = "%Y%m%d"))
  }
  m6 <- regmatches(fname, regexpr("\\d{6}", fname))
  if(length(m6) == 1 && nchar(m6) == 6){
    return(as.Date(paste0(m6,"01"), format = "%Y%m%d"))
  }
  return(NA)
}

# list file .nc
lai_files <- list.files(path = lai_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
sapply(lai_files, extract_date_from_filename)
str(lai_dates)
range(lai_dates)
table(format(lai_dates, "%Y"))
cat("Trovati", length(lai_files), "file .nc in", lai_dir, "\n")
if(length(lai_files) == 0) stop("Nessun .nc trovato: controlla nc_dir.")

# if defined, load AOI and convert it to SpatVector of "terra"
aoi_vect <- NULL #verify that the file exist
if(!is.na(aoi_path)){
  if(!file.exists(aoi_path)) stop("aoi_path non trovato!")
  aoi_sf <- sf::st_read(aoi_path, quiet = TRUE)
  aoi_sf <- sf::st_transform(aoi_sf, 4326)
  aoi_vect <- terra::vect(aoi_sf)
  cat("AOI caricata e trasformata in EPSG:4326\n")
}

# function that process a single .nc file and returns data.table with x,y,LAI,date
process_single_nc_to_dt <- function(nc_path, aoi_vect = NULL){
  cat("Processing:", nc_path, "...\n")
  # extracts the base name and the date
  fname <- basename(nc_path)
  date_val <- extract_date_from_filename(fname)
  # reads attributes with "ncdf4" to find name of variable, scale, fill
  nc <- try(nc_open(nc_path), silent = TRUE)
  if(inherits(nc, "try-error")){
    warning("Impossibile aprire:", nc_path); return(NULL)
  }
  vars <- names(nc$var)
  # looks for variable LAI (case-insensitive)
  lai_varname <- vars[grepl("^LAI$|LAI", vars, ignore.case = TRUE)][1]
  if(is.na(lai_varname) || is.null(lai_varname)){
    warning("LAI non trovata in:", nc_path)
    nc_close(nc); return(NULL)
  }
  # reads scale_factor and _FillValue if present
  scale_att <- try(ncatt_get(nc, lai_varname, "scale_factor"), silent = TRUE)
  scale_factor <- if(is.list(scale_att) && !is.null(scale_att$value)) as.numeric(scale_att$value) else NA
  fill_att <- try(ncatt_get(nc, lai_varname, "_FillValue"), silent = TRUE)
  fillv <- if(is.list(fill_att) && !is.null(fill_att$value)) as.numeric(fill_att$value) else NA
  nc_close(nc)
  
  # Reads the LAI variable as raster with "terra"
  # use "subds = lai_varname" when it is possible; if faile, read first layer
  r <- try(rast(nc_path, subds = lai_varname), silent = TRUE)
  if(inherits(r, "try-error")){
    r <- try(rast(nc_path), silent = TRUE)
    if(inherits(r, "try-error")){
      warning("terra non riesce a leggere:", nc_path); return(NULL)
    }
    if(nlyr(r) > 1) r <- subset(r, 1)
  }
  
  # applies scale and fill
  if(!is.na(scale_factor)) r <- r * scale_factor
  if(!is.na(fillv)){
    r[r == fillv] <- NA
  }
  
  # if it is provided, run crop to the AOI (save memory)
  if(!is.null(aoi_vect)){
    r <- crop(r, aoi_vect, snap = "out")
    r <- mask(r, aoi_vect)
  }
  
  # converts the raster to data.frame with xy and value; "terra::as.data.frame" is efficient
  # use "na.rm = FALSE" before to can decide if remove NA
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE) # colnames: x, y, layer
  # if the column is called "layer" or other, take the first value col
  val_col <- setdiff(names(df), c("x","y"))[1]
  setDT(df) # trasforms to data.table
  dt <- df[, .(x = get("x"), y = get("y"), LAI = get(val_col))]
  # remove NA LAI to reduce size (optional: if you want to keep grid complete, comment)
  dt <- dt[!is.na(LAI)]
  # adds the date
  dt[, date := date_val]
  # free up memory
  rm(df); gc()
  return(dt)
}

library(data.table)
# Loop on .nc files: save a .rds for each file (hundle a chunk per RAM)
chunk_files <- character(0)
i <- 1
for(nc in lai_files){
  dt <- process_single_nc_to_dt(nc, aoi_vect = aoi_vect)
  if(is.null(dt)) next
  if(save_intermediate_rds){
    rds_name <- file.path(out_dir, paste0(rds_chunk_prefix, sprintf("%04d", i), ".rds"))
    saveRDS(dt, file = rds_name, compress = "xz") # good compress for the space
    chunk_files <- c(chunk_files, rds_name)
    cat("Saved chunk:", rds_name, " (rows:", nrow(dt), ")\n")
    # free memory
    rm(dt); gc()
  } else {
    # accumulate in list (be careful RAM)
    if(i == 1) all_list <- list(dt) else all_list[[length(all_list)+1]] <- dt
  }
  i <- i + 1
}

# if we saved chunk .rds, merge and save as .RData (or save as unique .rds)
if(save_intermediate_rds && combine_and_save_RData){
  cat("Unione dei chunk .rds (potrebbe richiedere tempo e spazio)...\n")
  # read every rds and rbind in data.table (use "rbindlist" with "fill=FALSE")
  dt_list <- vector("list", length(chunk_files))
  for(j in seq_along(chunk_files)){
    dt_list[[j]] <- readRDS(chunk_files[j])
    cat("Caricato chunk", j, "rows:", nrow(dt_list[[j]]), "\n")
  }
  big_dt <- data.table::rbindlist(dt_list, use.names = TRUE, fill = FALSE)
  cat("Dimensione tabella finale:", format(object.size(big_dt), units = "auto"), " - righe:", nrow(big_dt), "\n")
  # save as .RData or .rds
  final_RData_name <- file.path(out_dir, paste0("LAI_pixels_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData"))
  save(big_dt, file = final_RData_name, compress = "xz")
  cat("Salvato .RData:", final_RData_name, "\n")
  # optional: save also as rds
  saveRDS(big_dt, file = file.path(out_dir, "LAI_pixels_combined.rds"), compress = "xz")
  cat("Salvato .rds:", file.path(out_dir, "LAI_pixels_combined.rds"), "\n")
  # cleaning
  rm(dt_list); gc()
} else if(!save_intermediate_rds && combine_and_save_RData){
  # if no chunk have been saved, the "all_list" contains the dt (be careful RAM)
  big_dt <- data.table::rbindlist(all_list, use.names = TRUE, fill = FALSE)
  save(big_dt, file = file.path(out_dir, "LAI_pixels_combined.RData"), compress = "xz")
  cat("Salvato .RData (no chunks):", file.path(out_dir, "LAI_pixels_combined.RData"), "\n")
}

cat("FINE. Controlla la cartella:", out_dir, "\n")

#Generic function to read one NetCDF variable
process_variable_nc_to_dt <- function(nc_path, varname, aoi_vect = NULL){
  cat("Processing:", nc_path, "...\n")
  fname <- basename(nc_path)
  date_val <- extract_date_from_filename(fname)
  
  nc <- try(nc_open(nc_path), silent = TRUE)
  if(inherits(nc, "try-error")) return(NULL)
  if(!(varname %in% names(nc$var))) {
    warning("Variabile", varname, "non trovata in:", nc_path)
    nc_close(nc); return(NULL)
  }
  
  scale_att <- try(ncatt_get(nc, varname, "scale_factor"), silent = TRUE)
  scale_factor <- if(is.list(scale_att) && !is.null(scale_att$value)) as.numeric(scale_att$value) else NA
  fill_att <- try(ncatt_get(nc, varname, "_FillValue"), silent = TRUE)
  fillv <- if(is.list(fill_att) && !is.null(fill_att$value)) as.numeric(fill_att$value) else NA
  nc_close(nc)
  
  r <- try(rast(nc_path, subds = varname), silent = TRUE)
  if(inherits(r, "try-error")) return(NULL)
  if(!is.na(scale_factor)) r <- r * scale_factor
  if(!is.na(fillv)) r[r == fillv] <- NA
  if(!is.null(aoi_vect)) {
    r <- crop(r, aoi_vect, snap = "out")
    r <- mask(r, aoi_vect)
  }
  
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  val_col <- setdiff(names(df), c("x","y"))[1]
  setDT(df)
  dt <- df[, .(x = get("x"), y = get("y"), value = get(val_col))]
  dt <- dt[!is.na(value)]
  dt[, date := date_val]
  return(dt)
}

#Apply the function to QFLAG, NOBS, RMSE
qflag_dir <- file.path(base_dir, "INPUT/nc_files/QFLAG")
nobs_dir  <- file.path(base_dir, "INPUT/nc_files/NOBS")
rmse_dir  <- file.path(base_dir, "INPUT/nc_files/RMSE")

read_variable_folder <- function(folder, varname){
  files <- list.files(folder, pattern = "\\.nc$", full.names = TRUE)
  dt_list <- lapply(files, function(f) process_variable_nc_to_dt(f, varname, aoi_vect))
  dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  setnames(dt, "value", varname)
  return(dt)
}

qflag_dt <- read_variable_folder(qflag_dir, "QFLAG")
nobs_dt  <- read_variable_folder(nobs_dir, "NOBS")
rmse_dt  <- read_variable_folder(rmse_dir, "RMSE")

#merge all with LAI
merged_dt <- readRDS(file.path(out_dir, "LAI_pixels_combined.rds"))
range(merged_dt$date)               
table(format(merged_dt$date, "%Y")) 

#save
saveRDS(merged_dt, file = file.path(out_dir, "LAI_pixels_combined.rds"), compress = "xz")


big_dt <- readRDS("C:/Users/iacom/Desktop/COPERNICUS_LAI/OUTPUT/LAI_pixels_combined.rds")
big_dt[, LAI:=30*LAI]
print(big_dt)
str(big_dt)
head(big_dt)
names(big_dt)

# Create a unique ID for each unique coordinate pair
big_dt$ID <- match(paste(big_dt$x, big_dt$y), unique(paste(big_dt$x, big_dt$y)))

library(ggplot2)
# plot
lai_plot <- ggplot(big_dt, aes(x = date, y = LAI, group = ID, color = factor(ID))) +
  geom_line(size = 0.8)

lai_plot

#filter data
big_dt_filtered <- big_dt[LAI<6]

lai_plot_filtered <- ggplot(big_dt_filtered, aes(x = date, y = LAI, group = ID, color = factor(ID))) +
  geom_line(size = 0.8)

lai_plot_filtered

#remove duplicates
big_dt_filtered <- big_dt_filtered[order(ID, date)]
big_dt_filtered <- big_dt_filtered[, .SD[!duplicated(date)], by = ID]

# Create a daily date range for each ID
daily_dates <- big_dt_filtered[, .(date = seq(min(date), max(date), by = "1 day")), by = ID]

# Interpolate LAI for each ID
daily_lai <- big_dt_filtered[daily_dates, on = .(ID, date), roll = "nearest"]
daily_lai[, LAI := approx(x = big_dt_filtered[ID == .BY$ID]$date,
                          y = big_dt_filtered[ID == .BY$ID]$LAI,
                          xout = date)$y, by = ID]

big_dt_filtered[, .N, by = .(ID, date)][N>1]


# View result
lai_plot_daily <- ggplot(daily_lai, aes(x = date, y = LAI, group = ID, color = factor(ID))) +
  geom_line(size = 0.8)

lai_plot_daily

#TIME GRAPH (trend time of LAI)
lai_time <- daily_lai %>%
  group_by(date) %>%
  summarise(mean_LAI=mean(LAI, na.rm = TRUE))
ggplot(lai_time, aes (x = date, y = mean_LAI)) +
  geom_line(color="forestgreen") +
  geom_point(color="darkgreen") +
  labs(title = "Time Trend of LAI",
       x = "Year",
       y = "Mean LAI") +
  theme_minimal()
#it shows how LAI changes over time over area

#seasonal or annual Boxplot
lai <- daily_lai %>%
  mutate(year = format(date, "%Y"),
         month = format(date, "%m"))

ggplot(lai, aes(x = month, y = LAI, group = month)) +
  geom_boxplot(fill = "lightgreen", color = "darkgreen") +
  labs(title = "Monthly mean LAI values",
       x = "Month",
       y = "LAI") +
  theme_minimal()


#VERIFY
# Find the closest pixel(ID):
# Target coordinates
target_lat <- 39.94033
target_lon <- -5.77465

# Compute distance and find closest ID
closest <- big_dt_filtered[, .(ID, x, y,
                               dist = sqrt((x - target_lon)^2 + (y - target_lat)^2))][
                                 order(dist)][1]

# View result
print(closest)