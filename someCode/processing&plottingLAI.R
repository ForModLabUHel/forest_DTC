library(ggplot2)
library(data.table)

#load data
load("C:/Users/checc/Downloads/big_st.rdata")

# Create a unique ID for each unique coordinate pair
big_dt$ID <- match(paste(big_dt$x, big_dt$y), unique(paste(big_dt$x, big_dt$y)))

# plot
lai_plot <- ggplot(big_dt, aes(x = date, y = LAI, group = ID, color = factor(ID))) +
  geom_line(size = 0.8)

