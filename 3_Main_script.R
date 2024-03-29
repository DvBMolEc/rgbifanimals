################################### Main Function ############################################
# Set working directory to directory where the R-script is saved
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # requires installation of package "rstudioapi"
#load in the functions
source("2_Functions.R")
# Create a rastered world
tr <- create_rastered_world("Inputs/new_tr.rdata")
# Read a species list
df <- read.csv("Output/Species_Location.csv")
# Read coordinates file
Coordinates <- read.csv("Inputs/Coordinates.csv")
# Find the shortest path for each sampling location to the species observation (as the crow flies)
sapply(df$Specieslist, function(species){
  print(species)
  tryCatch(find_closest_registered_place(species, Coordinates, tr, "Output/ShortestPath.csv"), error = function(e)return())
})

long <- pivot_longer(df, !Specieslist)
long <- long[long$value > 0, ]

apply(long, 1, function(row){
  samplelocation <- Coordinates[Coordinates$Observatory.ID == row[2], c("Longitude", "Latitude")]
  print(row[1])
  tryCatch({occurrence_data <- check_occurrence_data(row[1])
            find_shortest_route_in_sea(samplelocation, occurrence_data, tr, row, "Output/DistanceOverSea.csv")},
           error = function(errormessage){
             write.table(paste(c(row, 0, 0, NA), collapse = ","), file = "Output/DistanceOverSea.csv", 
                         append = TRUE, quote = FALSE, col.names = FALSE, row.names = FALSE)
           })
})

#### Clean R environment ####
rm(list = ls())

#### Additional, optional, cleaning step ####
results <- read.csv("Output/DistanceOverSea.csv")
results <- results[!is.na(results$distance),]
write.csv(results, "Output/DistanceOverSea.csv", quote = F, row.names = F)
