#Load all necessary packages
library("dplyr")
library("tidyverse")

#Set wd to wherever the script files are stored
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir <- getwd()

#Load BOLDigger output
data <- read.csv(file = "Inputs/BOLDigger_output.csv", sep = ";",
                 row.names = NULL)

#Extract locations of interest from main Data frame
locations_bold <- colnames(data)
locations_bold <- as.data.frame(stringr::str_split_fixed(string =
                                                           locations_bold,
                                                         pattern = "\\.", 2))
loi <- c("Toralla", "Getxo", "Vigo", "Roscoff", "Plymouth", "Galway",
         "BelgianCoast", "Bjorko", "Gbg", "Helsingborg", "Hjuvik", "Koster",
         "Laesoe1", "Laesoe2", "Laesoe3", "Limfjord", "Marstrand", "Preemraff",
         "Varberg", "Gdynia", "TZS")
locations_bold <- locations_bold %>%
  filter(V1 %in% loi)
#Paste colnames back together for later filtering steps
locations_bold$V3 <- paste0(locations_bold$V1, ".", locations_bold$V2)

#### Metadata ####
metadata <- read.csv("Inputs/MetaData.csv")
#Select applicable ARMS deployments
metadata <- left_join(metadata, locations_bold, by = c("Filename" = "V3"))
#Clean-up of data
metadata <- metadata %>% dplyr::select(-starts_with("X"), -V1, -V2)
metadata <- metadata %>% na_if("") %>% na.omit
metadata <- metadata[complete.cases(metadata), ]
metadata$Year <- format(as.Date(metadata$Deployment_date), "%Y")
metadata$Location_Year <- paste0(metadata$Observatory.ID, "_", metadata$Year)
#Add column in which ARMS fraction is removed
metadata$No_Fraction <- metadata$Filename
metadata$No_Fraction <- sub("\\.MF.00", "", metadata$No_Fraction)
metadata$No_Fraction <- sub("\\.MT.00", "", metadata$No_Fraction)
metadata$No_Fraction <- sub("\\.SF40", "", metadata$No_Fraction)

#Determine taxonomy order
taxonomy <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

#### Solve species name issues main df ####
#Merge Genus and species into Specieslist
data$Specieslist <- paste(data$Genus, data$Species)
#Replace , with . in Similarity, to be able to assign numeric variable
data$Similarity <- gsub(",", ".", x = data$Similarity, fixed = TRUE)
#Substitute species names with synonyms
synonyms <- read.csv("Inputs/Synonyms.csv")
old_name <- synonyms$Old_name
new_name <- synonyms$New_name
data$Specieslist[data$Specieslist %in% old_name] <-
  new_name[match(data$Specieslist, old_name, nomatch = 0)]
rm(old_name, new_name, synonyms)

#Filter dataset based on locations of interest
data_loi <- data[, c("sequence", taxonomy, "Specieslist", "Similarity",
                     metadata$Filename)]
rm(data)
# Remove contaminants ...
remove_species <- c("sapiens", "lupus", "scrofa")
data_loi <- data_loi[!(data_loi$Species %in% remove_species), ]
# ... and groups which are not of interest
remove_phyla <- c("Amoebozoa", "Ascomycota", "Bacillariophyta", "Chlorophyta",
                  "Heterokontophyta", "Ochrophyta", "Rhodophyta", "Zygomycota")
data_loi <- data_loi[!(data_loi$Phylum %in% remove_phyla), ]
# Remove sequences with no observations in our chosen subset
data_loi <- data_loi[rowSums(data_loi[, 10:ncol(data_loi)]) > 0, ]
#Reset rownumbers for easier subsetting
rownames(data_loi) <- c(seq_len(nrow(data_loi)))
rm(remove_phyla, remove_species)

#### Set of additional dataframes for additional insight ####
#Extract unique species with their corresponding lowest Similarity
uniques <- data_loi %>%
  dplyr::select(Specieslist, Similarity) %>%
  arrange(Specieslist) %>%
  distinct() %>%
  group_by(Specieslist) %>%
  slice_min(order_by = Similarity, with_ties = TRUE)

#ASV count for each species
asv_count <- data_loi[!duplicated(data_loi$sequence), ] %>%
  dplyr::count(Specieslist)

# Read count per ASV, per ARMS, per year
# (sum three fractions MF100, MF500, SF40)
# Create duplicate column names
df <- data_loi
colnames(df) <- c("sequence", taxonomy, "Specieslist", "Similarity",
                  metadata$No_Fraction)
column_order <- unique(colnames(df))
#Merge columns with duplicate names and sum contents
df2 <- sapply(unique(colnames(df)[duplicated(colnames(df))]), function(x) {
  rowSums(df[, grepl(paste(x, "$", sep = ""), colnames(df))])
})
read_count_asvs_arms <- cbind(df2, df[, !duplicated(colnames(df)) &
                                        !duplicated(colnames(df),
                                                    fromLast = TRUE)])
rm(df, df2)

#Read count per Species, per ARMS, per year
read_count_species_arms <- read_count_asvs_arms %>%
  dplyr::select(-sequence, -Phylum, -Order, -Class, -Family, -Genus, -Species,
                -Similarity) %>%
  dplyr::select(Specieslist, sort(colnames(.))) %>%
  group_by(Specieslist) %>%
  summarise_all(sum)
#Readable format
read_count_asvs_arms <- read_count_asvs_arms[, column_order]

#Read count per Location per Year
rc_species_loc_year <- data_loi[, c(taxonomy[1:4], "Specieslist",
                                    metadata$Filename)]
#Create duplicate column names
colnames(rc_species_loc_year) <- c(taxonomy[1:4], "Specieslist",
                                   metadata$Location_Year)
df <- rc_species_loc_year
#Merge columns with duplicate names and sum contents
df2 <- sapply(unique(colnames(df)[duplicated(colnames(df))]), function(x) {
  rowSums(df[, grepl(paste(x, "$", sep = ""), colnames(df))])
})
rc_species_loc_year <- as.data.frame(cbind(df2, df[, unique(colnames(df))]))
rm(df, df2)
#group by Specieslist and sum read counts locations
rc_species_loc_year <- rc_species_loc_year %>%
  relocate(Phylum:Specieslist, .before = BelgianCoast_2018) %>%
  mutate_at(c(6:ncol(rc_species_loc_year)), as.numeric) %>%
  group_by(Phylum, Class, Order, Family, Specieslist, ) %>%
  summarise(across(BelgianCoast_2018:Vigo_2019, sum))

#Presence/Absence matrix for Species per ARMS
pres_abs <- read_count_species_arms
pres_abs[2:ncol(pres_abs)][pres_abs[2:ncol(pres_abs)] > 0] <- 1

#Read count per species per fraction
read_count_species_fraction <- data_loi  %>%
  dplyr::select(-sequence, -Phylum, -Order, -Class, -Family, -Genus, -Species,
                -Similarity) %>%
  dplyr::select(Specieslist, sort(colnames(.))) %>%
  group_by(Specieslist) %>%
  summarise_all(sum)

#Export
dir.create(file.path(dir, "Output"))
setwd(file.path(dir, "Output"))
write.csv(uniques, "UniqueSpecies.csv", row.names = FALSE)
write.csv(asv_count, "ASV_Count.csv", row.names = FALSE)
write.csv(read_count_asvs_arms, "Read_Count_ASVs_ARMS.csv", row.names = FALSE)
write.csv(read_count_species_arms, "Read_Count_Species_ARMS.csv",
          row.names = FALSE)
write.csv(rc_species_loc_year, "Read_Count_Species_Loc_Year.csv",
          row.names = FALSE)
write.csv(read_count_species_fraction, "Read_Count_Species_Fraction.csv",
          row.names = FALSE)
write.csv(data_loi, "Data_LOI.csv", row.names = FALSE)

#### Build dataframe to be used in AlienIdentification.R ####
species_location <- data_loi  %>%
  dplyr::select(-sequence, -Phylum, -Order, -Class, -Family, -Genus, -Species,
                -Similarity) %>%
  dplyr::select(Specieslist, sort(colnames(.))) %>%
  group_by(Specieslist) %>%
  summarise_all(sum)
colnames(species_location) <- c("Specieslist", metadata$Observatory.ID)
df <- species_location
#Merge columns with duplicate names and sum contents
df2 <- sapply(unique(colnames(df)[duplicated(colnames(df))]), function(x) {
  rowSums(df[, grepl(paste(x, "$", sep = ""), colnames(df))])
})
species_location <- as.data.frame(cbind(df2, df[, unique(colnames(df))]))
species_location <- species_location %>%
  dplyr::select(Specieslist, everything())
rm(df, df2)
write.csv(species_location, "Species_Location.csv", row.names = FALSE)

# Export Coordinates, which are used in 3_Main_script and MetaData,
# which is used in 4_Visualisation
setwd(dir)
coordinates <- metadata %>%
  group_by(Observatory.ID) %>%
  dplyr::select(Observatory.ID, Longitude, Latitude) %>%
  slice_min(Longitude) %>%
  unique
write.csv(coordinates, "Inputs/Coordinates.csv", row.names = FALSE)
write.csv(metadata, "Inputs/MetaData_Adjusted.csv", row.names = FALSE)

#### Clean R environment ####
rm(list = ls())
