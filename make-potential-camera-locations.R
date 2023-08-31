library(sf)
library(leaflet)
library(secr)

x <- readRDS("clt_surveyregion_anita.Rds")
class(x)

# plot regions to see what's happening
plot(st_geometry(x), axes= TRUE)
plot(st_geometry(x[1,]), axes= TRUE)
plot(st_geometry(x[2,]), add= TRUE) # area to remove
plot(st_geometry(x[1,]), axes= TRUE)
plot(st_geometry(x[3,]), add= TRUE) # area to remove
plot(st_geometry(x[1,]), axes= TRUE)
plot(st_geometry(x[4,]), add= TRUE) # area to remove
plot(st_geometry(x[1,]), axes= TRUE)
plot(st_geometry(x[5,]), add= TRUE) # area to add

# convert to UTM so can take intersections, unions etc
my_crs <- "+proj=utm +zone=34 +south +datum=WGS84 +units=m +no_defs"
x1_utm <- x[1,] %>% st_transform(my_crs)
x2_utm <- x[2,] %>% st_transform(my_crs)
x3_utm <- x[3,] %>% st_transform(my_crs)
x4_utm <- x[4,] %>% st_transform(my_crs)
x5_utm <- x[5,] %>% st_transform(my_crs)

# remove areas x2, x3, x4 and add area x5
potential_cams <- x1_utm %>% 
  st_difference(x2_utm) %>% 
  st_difference(x3_utm) %>% 
  st_difference(x4_utm) %>%
  st_union(x5_utm)
plot(potential_cams)

# plot on map
leaflet()%>%
  addProviderTiles(provider= "CartoDB.Positron") %>%
  addPolygons(data = potential_cams %>% st_transform(4326), color = "green", fillOpacity = 0.7, weight = 1, group = "editable")

# mesh of possible camera locations
# different spacing (1km, 2km, and 4km where last two are 2/3 sigma for sigmas of 3/6K 
potential_cams_mesh_1km <- make.mask(type = "polygon", poly = potential_cams, spacing = 1000)
potential_cams_mesh_SigmaF <- make.mask(type = "polygon", poly = potential_cams, spacing = 2000)
potential_cams_mesh_SigmaM <- make.mask(type = "polygon", poly = potential_cams, spacing = 4000)

#buffered versions
potential_cams_mesh_SigmaF.buff <- make.mask(type = "polybuffer", poly = potential_cams, spacing = 2000, buffer = 2000)

plot(potential_cams_mesh_1km)

save(potential_cams_mesh_1km, potential_cams_mesh_SigmaF, potential_cams_mesh_SigmaM, file = "potential_cams_mesh.Rdata")

# plot on map
potential_cams_mesh_sf <- st_as_sf(potential_cams_mesh_SigmaM, coords = c("x", "y"), crs = my_crs)
leaflet()%>%
  addProviderTiles(provider= "CartoDB.Positron") %>%
  addPolygons(data = potential_cams %>% st_transform(4326), color = "red", fill = NA, weight = 2) %>%
  addCircles(data = potential_cams_mesh_sf %>% st_transform(4326), color = "green", fillOpacity = 0.7, weight = 1, group = "editable")

