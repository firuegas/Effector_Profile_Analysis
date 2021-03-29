
### CREATE MAP WITH EFFECTOR PROFILE DISTRIBUTION, BY COUNTRY
### For usa, location can be broken down by state
### Fernand Iruegas Bocardo
### January 2020


### INSTALL AND LOAD PACKAGES

install.packages("ggmap")
install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
install.packages('rgeos')
install.packages("memisc")
install.packages('oz')
install.packages("scatterpie")
install.packages("maptools")

library(memisc)
library(assertthat)
library(sqldf)
library(magrittr)
library(dplyr)
library(reshape2)
library(ggplot2)
theme_set(theme_bw())
library(oz)
library(scatterpie)
library(rgdal)
library(maptools)
#library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggmap)
library(tidyverse)


### Upload input files

# Locations
location <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/All_analysis_062719/NMDS_Jan2020/Xp_locations_country.txt", 
                       sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = TRUE, blank=F)
str(location)
class(location)
#colnames(location) <- c("address")
locs <- unique(location$address)


#Metadata table
metadata <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/Xp_effProfile_2020/haplotypes_with_gev1063/metadata_coordinates_effclusters2.txt", 
                       sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = TRUE)
str(metadata)


#NMDS coodinates, to define effector profile clusters
nmds_coord <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/All_analysis_062719/NMDS_Jan2020/nmds_coordenates_annotated.txt", 
                           sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = TRUE, blank=F)
str(nmds_coord)
class(nmds_coord)

# OR
# Table with all the information already put together, Jan 2020
metadata <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/All_analysis_062719/NMDS_Jan2020/metadata_coordinates_effclusters.txt", 
                       sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = TRUE)
str(metadata)

# Wide datasset already formated
grouploc_wide <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/All_analysis_062719/NMDS_Jan2020/wide_table_effclusters.txt", 
                       sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = TRUE)
str(grouploc_wide)




### GET LON/LAT COORDENATES
# Get lon/lat coordenates from Google Maps

# Google Maps require registering
?register_google
# Personal Key API
# AIzaSyChuwqOXJczTqFjvCopTkrmBkHfYjV9qIk
#register_google(key = 'AIzaSyChuwqOXJczTqFjvCopTkrmBkHfYjV9qIk', write = TRUE)

# Find lon/lat coordinates for each location
# mutate_geocode() uses Google Maps to find longitud and latitude or each location
# requires two arguments: a dataframe and the name of the column with location data 
locations_df <- mutate_geocode(metadata, address)
head(locations_df)
#locations_df



### JOIN METADATA, LAT/LON AND EFFECTOR PROFILE CLUSTER GROUP IN A SINGLE TABLE

# Join metadata and lat/lon information in a single table metadata_loc
head(metadata)
head(location)

#location_order <- cbind(Sphap_nmds,locations_df$lon,locations_df$lat)
location_order <- locations_df[match(Sphap_nmds$strain, locations_df$strain),]
head(location_order)

metadata_loc <- merge(metadata, location_order, all=TRUE) # This table looses order from phylogenetic tree
head(metadata_loc)

# Join  metadata and effector profile clusters from nmds analysis in the table metadata_loc
head(nmds_coord)

nmds_order <- nmds_coord[match(metadata_loc$strain, nmds_coord$strain),]
head(nmds_order)

metadata_loc <- merge(metadata_loc, nmds_order, all=TRUE)
head(metadata_loc)

write.table(locations_df, file = "metadata_coordinates_effclusters2.txt", col.names = TRUE, row.names = TRUE, sep= "\t", quote = FALSE)


### FORMAT DATA  for map with piecharts
# columns by Effector Profile Cluster

# If the final metadata table is available, skip formating and input here
metadata_loc <- locations_df

grp <- metadata_loc %>% group_by( address, eff_cluster) %>%
  dplyr::summarise(count=n())
grp

#Spread long table into wide table
group_wide <- data.frame(t(spread(grp,address,count,fill = 0)))

#Get lon/lat coordinates of unique locations
uniq_loc <- metadata_loc[match(rownames(group_wide), metadata_loc$address),]
grouploc_wide <- cbind(group_wide, lon=uniq_loc$lon, lat=uniq_loc$lat, address=uniq_loc$address)
class(grouploc_wide)

colnames(grouploc_wide) <- c(LETTERS[1:9], "lon", "lat","address" )
grouploc_wide

grouploc_wide <- grouploc_wide[-1,]
grouploc_wide

# Make sure columns are numeric
u <- c(1:10)        # Select columns
grouploc_wide[ , u] <- apply(grouploc_wide[ , u], 2,            # Specify own function within apply
                    function(x) as.numeric(as.character(x)))
sapply(grouploc_wide, class)                           # Get classes of all columns


# Count no. strains per location, to use as scalar
strain <- c()
for(i in 1:nrow(grouploc_wide)){
  strain[i] <-  sum(grouploc_wide[i,1:9])
}
strain
grouploc_wide$strain <- strain

grouploc_wide       # Table ready to input in ggplot!!!

# Save to file
write.table(grouploc_wide, file = "Table_maps_Xp2.txt", col.names = TRUE, row.names = TRUE, sep= "\t", quote = FALSE)



##################################

### CREATE BASE MAPS

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

#USA <- ne_countries(country = "united states of america")
#plot(USA)



### PLOT MAP WIT PIECHART DISTRIBUTION OF EFFECTOR PROFILES

# Colors for effector clusters
mycoleff <- c("#FF7F50", "#0000FF", "#5CACEE", "#FFC125", "#B8DC8C", "#CD919E", "#FFFF00", "#FF0000", "#D2B48C")%>%
  setNames(ltr_ordr)

pdf("Map_effProfileClusters_Feb21.pdf")
ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  #coord_sf( ylim = c(70.77312286, -61.6092206), expand = FALSE) +
  geom_scatterpie(aes(x=lon, y = lat, group = address, r = (log10(strain+1))*4), data = grouploc_wide, cols=colnames(grouploc_wide)[1:9], legend_name = "Clusters", size = 0.05) +
  
  #geom_polygon(data = grouploc_wide, aes(x = lon, y = lat, group = address, size = 0.05)) +
  
  #geom_point(data = metadata_loc, aes(x = lon, y = lat, color = as.factor(group)), alpha= 0.7) +  # To plot a dot for each strain 
  scale_fill_manual(values= c(unique(mycoleff), "#0085C7","tomato", "#F4C300", "#0085C7","tomato","#F4C300")) +
  ggtitle("Xp Effector profile distribution") +
  guides( color =  guide_legend("phylogroup"),  shape =  guide_legend( "phylogroup")) +
  theme(panel.background = element_rect(fill = "white", 
                                        colour = "black", 
                                        size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, 
                                        linetype = 'solid',
                                        colour = "gray94"), 
        axis.text=element_text(size=14),
        axis.title=element_text(size=12), 
        plot.title = element_text(hjust = 0.5)) 

dev.off()



### TESTING 















