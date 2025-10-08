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

aoi_path<- file.path(base_dir, "INPUT/AOI/Bounding_Box.shp")
out_dir <- file.path(base_dir, "OUTPUT_LAI")  # the folder where you save GeoTIFF and metadata
nc_dir<- file.path(base_dir, "INPUT/nc_files")  # the foldere where you have .nc files

#setting with resolution and LAI_scale factors
my_resolution <- 300
LAI_scaling_factor <- 30

file.exists(aoi_path)
file.exists(nc_dir)

setwd(nc_dir)
nc1 <- nc_open("c_gls_LAI300_201501100000_GLOBE_PROBAV_V1.0.1__LAI.nc")

print(nc1)

#if you don't have it, create output folder
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#move the directory to OUTPUT_LAI where I'll save the results
setwd(out_dir)

cat("Directory di lavoro impostato su:", getwd(), "\n")
cat("AOI", aoi_path, "\n")
cat("File .nc in:", nc_dir, "\n")
cat("OUTPUT_LAI in:", out_dir, "\n")

#Verify that AOI and file.nc exist
if(!file.exists(aoi_path)) stop("AOI non trovato! Controlla il percorso:", aoi_path)
if(!file.exists(nc_dir)) stop("Cartella nc_files non trovata! Controlla il percorso:", nc_dir)

# List file .nc 
#nc_dir is the folder where I put all the files .nc (YEARS)
#list.files() searches for all files ending in .nc
#full.names=TRUE returns the full path
#recursive=TRUE lokks in the sub folders
nc_files <- list.files(path = nc_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
cat("Trovati", length(nc_files)>0, "File .nc in", nc_dir, "\n")
head(nc_files, 3)

# Read AOI and take to EPSG:4326
aoi_sf <- st_read(aoi_path, quiet = TRUE)
aoi_sf <- st_transform(aoi_sf, 4326)
aoi_vect <- terra::vect(aoi_sf)

# utility: extracts the date from the filename (handles 12, 8, 6 digits)
extract_date_from_filename <- function(fname){
  # try pattern 12 (es: 201601100000), then 8 (20160110), and then 6 (201601)
  m12 <- regmatches(fname, regexpr("\\d{12}", fname))
  if(length(m12)==1 && nchar(m12)==12) {
    date_str <- substr(m12, 1, 8)
    return(as.Date(date_str, "%Y%m%d"))
  }
  m8 <- regmatches(fname, regexpr("\\d{8}", fname))
  if(length(m8)==1 && nchar(m8)==8){
    return(as.Date(m8, "%Y%m%d"))
  }
  m6 <- regmatches(fname, regexpr("\\d{6}", fname))
  if(length(m6)==1 && nchar(m6)==6){
    # Take YYYYMM -> before the month
    return(as.Date(paste0(m6, "01"), "%Y%m%d"))
  }
  return(NA)
}

# this function process a single netCDF
# it's a custom function that:
#1. receives the path to a .nc (nc_path) file as input
#2. reads data inside the file
#3. processes them
#4. the name of the generated GeoTIFF file returns
process_one_nc <- function(nc_path, aoi_vect, resolution, extract_qc = FALSE) {
  cat("Processing:", nc_path, "...\n")
  # 1) opens the NetCDF file with "nc_open()" and handles errors with "try()"
  nc <- try(nc_open(nc_path), silent = TRUE)
  if(inherits(nc, "try-error")){
    warning("Impossibile aprire:", nc_path) 
    return(NULL)
  }
  #identify variable LAI
  vars <- names(nc$var)
  # look for the variable containing LAI (case-insensitive) among the NetCDF variables
  lai_var <- vars[grepl("^LAI$|LAI", vars, ignore.case = TRUE)][1]
  if(is.na(lai_var)){
    warning("Variabile LAI non trovata in ", nc_path)
    nc_close(nc)
    return(NULL)
  }

#above it open the .nc file,looks at what variables there are inside (es. LAI, NOBS, RMSE,..)
  #look for the variable that contains the LAI dataand if it doesn't find it, it will warn you
  
  # Reads useful attributes (scale/_FillValue) using "ncatt_get()"
  scale_att <- try(ncatt_get(nc, lai_var, "scale_factor"), silent = TRUE)
  scale_factor <- if(is.list(scale_att) && !is.null(scale_att$value)) as.numeric(scale_att$value) else NA
  fill_att <- try(ncatt_get(nc, lai_var, "_FillValue"), silent = TRUE)
  fillv <- if(is.list(fill_att) && !is.null(fill_att$value)) as.numeric(fill_att$value) else NA
  nc_close(nc) #close the .nc file
  
  # 2) Reads the variable as raster with "terra" (read the LAI variables as raster,that is a numerical map)
  #each pixel of the raster rappresents a value of LAI
  r <- try(rast(nc_path, subds = lai_var), silent = TRUE)
  if(inherits(r, "try-error")){
    r <- try(rast(nc_path), silent = TRUE)
    if(inherits(r, "try-error")){
      warning("terra non riesce a leggere il file:", nc_path)
      return(NULL)
    }
    # se il file ha più layer, prova a scegliere il primo
    if(nlyr(r) > 1) r <- raster::subset(r, 1)
  }
  
  # 3) it applies "scale" and "fill" to correct because some .nc files save data multiplied by a factor (e.g. 0.001)
  #it applies the scale_factor and replaces "empty" (FillValue) values with NA (missing)
  if(!is.na(scale_factor)) r <- r * scale_factor
  if(!is.na(fillv)) r[r == fillv] <- NA
  
  # 4) if 300m (PROBA-V) it applies shift mezza-pixel (dx/dy= +or- 1/672°)
  if(resolution == my_resolution){
    # shift = 1/(336*2) gradi ≈ 1/672
    xshift <- 1/672
    yshift <- 1/672
    r <- shift(r, dx = -xshift, dy = yshift)
  }
  
  # 5) "crop" (cut the raster to AOI limits (reduce size) 
  #& "mask" (puts everything outside the AOI perimeter to NA)
  #Result: only data within my area
  r_crop <- crop(r, aoi_vect, snap = "out")
  r_mask <- mask(r_crop, aoi_vect)
  
  # 6) extracts the date from filename to nome the output
  fname <- basename(nc_path)
  date_val <- extract_date_from_filename(fname)
  if(is.na(date_val)){
    outname <- file.path(getwd(), paste0("LAI_", tools::file_path_sans_ext(fname), ".tif"))
  } else {
    outname <- file.path(getwd(), sprintf("LAI_%s.tif", format(date_val, "%Y%m%d")))
  }
 writeRaster(r, filename = outname, overwrite = TRUE)
 return(outname)
    }

cat("FINE pipeline.\n")

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
out_dir <- file.path(base_dir, "OUTPUT_LAI")  # where save GeoTIFF, metadata and .rds / .RData
nc_dir<- file.path(base_dir, "INPUT/nc_files")  # the folder where there are .nc files

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
  if(length(m12)==1 && nchar(m12)==12){
    return(as.Date(substr(m12,1,8), "%Y%m%d"))
  }
  m8 <- regmatches(fname, regexpr("\\d{8}", fname))
  if(length(m8)==1 && nchar(m8)==8){
    return(as.Date(m8, "%Y%m%d"))
  }
  m6 <- regmatches(fname, regexpr("\\d{6}", fname))
  if(length(m6)==1 && nchar(m6)==6){
    return(as.Date(paste0(m6,"01"), "%Y%m%d"))
  }
  return(NA)
}

# list file .nc
nc_files <- list.files(path = nc_dir, pattern = "\\.nc$|\\.NC$|\\.nc4$", full.names = TRUE, recursive = TRUE)
cat("Trovati", length(nc_files), "file .nc in", nc_dir, "\n")
if(length(nc_files) == 0) stop("Nessun .nc trovato: controlla nc_dir.")

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

# Loop on .nc files: save a .rds for each file (hundle a chunk per RAM)
chunk_files <- character(0)
i <- 1
for(nc in nc_files){
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

library(data.table)
big_dt <- readRDS("C:/Users/iacom/Desktop/COPERNICUS_LAI/OUTPUT_LAI/LAI_pixels_combined.rds")
big_dt[, LAI:=30*LAI]
print(big_dt)
str(big_dt)
head(big_dt)
names(big_dt)
#convertiamolo in un oggetto spaziale in R con "sf"
sf_lai <- st_as_sf(big_dt,
                   coords = c("x", "y"),
                   crs = 4326)
#ora sf_lai è un oggetto spaziale che posso visualizzare in R
plot(sf_lai["LAI"])

st_write(sf_lai, "C:/Users/iacom/Desktop/COPERNICUS_LAI/OUTPUT_LAI/LAI_points.gpkg", delete_dsn = TRUE)

#Create an unique ID for each coordinate pair
big_dt$ID <- match(paste(big_dt$x, big_dt$y), unique(paste(big_dt$x, big_dt$y)))

#Plot
lai_plot <- ggplot(big_dt, aes(x = date, y = LAI, group = ID, colour = factor(ID))) +
  geom_line(size = 0.8)
lai_plot

#GRAFICO TEMPORALE (andamento nel tempo del LAI)
library(ggplot2)
library(dplyr)

lai_time <- big_dt %>%
  group_by(date) %>%
  summarise(mean_LAI=mean(LAI, na.rm = TRUE))
ggplot(lai_time, aes (x = date, y = mean_LAI)) +
  geom_line(color="forestgreen") +
  geom_point(color="darkgreen") +
  labs(title = "Time Trend of LAI",
       x = "Year",
       y = "Mean LAI") +
  theme_minimal()
#mostra come cambia LAI (foglie/vegetazione) nel tempo sull'area

#Boxplot stagionali o annuali
lai <- big_dt %>%
  mutate(year = format(date, "%Y"),
         month = format(date, "%m"))

ggplot(lai, aes(x = month, y = LAI, group = month)) +
  geom_boxplot(fill = "lightgreen", color = "darkgreen") +
  labs(title = "Monthly mean LAI values",
       x = "Month",
       y = "LAI") +
  theme_minimal()

#if I want to filter years
big_dt <- big_dt %>%
  mutate(year = year(date))
big_dt_sub <- big_dt %>%
  filter(year >= 2020, year <= 2022)

lai_time_sub <- big_dt_sub %>%
  group_by(date) %>%
  summarise(mean_LAI=mean(LAI, na.rm = TRUE))

ggplot(lai_time_sub, aes (x = date, y = mean_LAI)) +
  geom_line(color="forestgreen") +
  geom_point(color="darkgreen") +
  labs(title = "Time Trend of LAI 2020-2022",
       x = "Year",
       y = "Mean LAI") +
  theme_minimal()

lai_sub <- big_dt_sub %>%
  mutate(year = format(date, "%Y"),
         month = format(date, "%m"))

ggplot(lai_sub, aes(x = month, y = LAI, group = month)) +
  geom_boxplot(fill = "lightgreen", color = "darkgreen") +
  labs(title = "Monthly mean LAI values 2020-2022",
       x = "Month",
       y = "LAI") +
  theme_minimal()
