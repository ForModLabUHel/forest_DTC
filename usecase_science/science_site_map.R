library(leaflet)

# Define site coordinates
sites <- data.frame(
  name = c("Guyaflux", "Hainich", "Hohes Holz", "Loobos", "Las Majadas",, "Sodankylä"),
  lon = c(-52.0, 10.4, 9.9, 5.77, -5.77, 26.60),
  lat = c(5.0, 51.2, 51.5, 52.20, 39.94, 67.37)
)

# Create the map
map <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldTopoMap) %>%  # Terrain background
  addCircleMarkers(data = sites, ~lon, ~lat,
                   label = ~name,
                   color = "red",
                   radius = nrow(sites),
                   fillOpacity = 0.9) %>%
  addLegend(position = "bottomright", colors = "red", labels = "Research Sites")

# Display the map
map

