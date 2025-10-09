library(ggplot2)
library(data.table)

#load data
load("C:/Users/checc/Downloads/big_st.rdata")

# Create a unique ID for each unique coordinate pair
big_dt$ID <- match(paste(big_dt$x, big_dt$y), unique(paste(big_dt$x, big_dt$y)))

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

# View result
lai_plot_daily <- ggplot(daily_lai, aes(x = date, y = LAI, group = ID, color = factor(ID))) +
  geom_line(size = 0.8)

lai_plot_daily




# Find the closest pixel(ID):
# Target coordinates
target_lat <- 39.94033
target_lon <- -5.77465

# Compute distance and find closest ID
closest <- big_dt[, .(ID, x, y,
                      dist = sqrt((x - target_lon)^2 + (y - target_lat)^2))][
                        order(dist)][1]

# View result
print(closest)
