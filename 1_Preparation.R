#Load all necessary packages
library("dplyr")
library("tidyverse")

#Set wd to wherever the script files are stored
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
dir <- getwd()

#Load BOLDigger output
Data <- read.csv(file = "Inputs/BOLDigger_output.csv", sep = ";", row.names = NULL)

#Extract locations of interest from main Data frame
Locations_BOLD <- colnames(Data)
Locations_BOLD <- as.data.frame(stringr::str_split_fixed(string = Locations_BOLD, pattern = "\\.", 2))
LOI <- c("Toralla", "Getxo", "Vigo", "Roscoff", "Plymouth", "Galway", "BelgianCoast", "Bjorko", "Gbg", "Helsingborg", "Hjuvik", "Koster", "Laesoe1", "Laesoe2", "Laesoe3", "Limfjord", "Marstrand", "Preemraff", "Varberg", "Gdynia", "TZS")
Locations_BOLD <- Locations_BOLD %>% 
  filter(V1 %in% LOI)
#Paste colnames back together for later filtering steps
Locations_BOLD$V3 <- paste0(Locations_BOLD$V1, ".", Locations_BOLD$V2)

#### Metadata ####
MetaData <- read.csv("Inputs/MetaData.csv")
#Select applicable ARMS deployments
MetaData <- left_join(MetaData, Locations_BOLD, by = c("Filename" = "V3"))
#Clean-up of data
MetaData <- MetaData %>% dplyr::select(-starts_with("X"), -V1, -V2)
MetaData <- MetaData %>% na_if("") %>% na.omit
MetaData <- MetaData[complete.cases(MetaData), ]
MetaData$Year <- format(as.Date(MetaData$Deployment_date), "%Y")
MetaData$Location_Year <- paste0(MetaData$Observatory.ID, "_", MetaData$Year)
#Add column in which ARMS fraction is removed
MetaData$No_Fraction <- MetaData$Filename
MetaData$No_Fraction <- sub("\\.MF.00", "", MetaData$No_Fraction)
MetaData$No_Fraction <- sub("\\.MT.00", "", MetaData$No_Fraction)
MetaData$No_Fraction <- sub("\\.SF40", "", MetaData$No_Fraction)

#Determine taxonomy order
Taxonomy <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")

#### Solve species name issues main df ####
#Merge Genus and species into Specieslist
Data$Specieslist <- paste(Data$Genus, Data$Species)
#Replace , with . in Similarity, to be able to assign numeric variable
Data$Similarity <- gsub(",", ".", x = Data$Similarity, fixed = T)
#Substitute species names with synonyms
Synonyms <- read.csv("Inputs/Synonyms.csv")
Old_name <- Synonyms$Old_name
New_name <- Synonyms$New_name
Data$Specieslist[Data$Specieslist %in% Old_name] <- New_name[match(Data$Specieslist, Old_name, nomatch = 0)]
rm(Old_name, New_name, Synonyms)

#Filter dataset based on locations of interest
Data_LOI <- Data[, c("sequence", Taxonomy, "Specieslist", "Similarity", MetaData$Filename)]
rm(Data)
# Remove contaminants ...
RemoveSpecies <- c("sapiens", "lupus", "scrofa")
Data_LOI <- Data_LOI[!(Data_LOI$Species %in% RemoveSpecies), ]
# ... and groups which are not of interest
RemovePhyla <- c("Amoebozoa", "Ascomycota", "Bacillariophyta", "Chlorophyta", "Heterokontophyta", "Ochrophyta", "Rhodophyta", "Zygomycota")
Data_LOI <- Data_LOI[!(Data_LOI$Phylum %in% RemovePhyla), ]
# Remove sequences with no observations in our chosen subset
Data_LOI <- Data_LOI[rowSums(Data_LOI[, 10:ncol(Data_LOI)])>0, ]
#Reset rownumbers for easier subsetting
rownames(Data_LOI) <- c(1:nrow(Data_LOI))
rm(RemovePhyla, RemoveSpecies)

#### Set of additional dataframes for additional insight ####
#Extract unique species with their corresponding lowest Similarity
Uniques <- Data_LOI %>% 
  dplyr::select(Specieslist, Similarity) %>% 
  arrange(Specieslist) %>%
  distinct() %>%
  group_by(Specieslist) %>%
  slice_min(order_by = Similarity, with_ties = T)

#ASV count for each species
ASV_Count <- Data_LOI[!duplicated(Data_LOI$sequence), ] %>% dplyr::count(Specieslist)

#Read count per ASV, per ARMS, per year (sum three fractions MF100, MF500, SF40) 
#Create duplicate column names
df <- Data_LOI
colnames(df) <- c("sequence", Taxonomy, "Specieslist", "Similarity", MetaData$No_Fraction)
Column_Order <- unique(colnames(df))
#Merge columns with duplicate names and sum contents
df2 <- sapply(unique(colnames(df)[duplicated(colnames(df))]), function(x) rowSums(df[,grepl(paste(x, "$", sep=""), colnames(df))]))
Read_Count_ASVs_ARMS <- cbind(df2, df[,!duplicated(colnames(df)) & !duplicated(colnames(df), fromLast = TRUE)])
rm(df, df2)

#Read count per Species, per ARMS, per year
Read_Count_Species_ARMS <- Read_Count_ASVs_ARMS %>%
  dplyr::select(-sequence, -Phylum, -Order, -Class, -Family, -Genus, -Species, -Similarity) %>%
  dplyr::select(Specieslist, sort(colnames(.))) %>%
  group_by(Specieslist) %>%
  summarise_all(sum)
#Readable format
Read_Count_ASVs_ARMS <- Read_Count_ASVs_ARMS[, Column_Order]

#Read count per Location per Year
RC_Species_Loc_Year <- Data_LOI[, c(Taxonomy[1:4], "Specieslist", MetaData$Filename)]
#Create duplicate column names
colnames(RC_Species_Loc_Year) <- c(Taxonomy[1:4], "Specieslist", MetaData$Location_Year)
df <- RC_Species_Loc_Year
#Merge columns with duplicate names and sum contents
df2 <- sapply(unique(colnames(df)[duplicated(colnames(df))]), function(x) rowSums(df[,grepl(paste(x, "$", sep=""), colnames(df))]))
RC_Species_Loc_Year <- as.data.frame(cbind(df2, df[,!duplicated(colnames(df)) & !duplicated(colnames(df), fromLast = TRUE)]))
rm(df, df2)
#group by Specieslist and sum read counts locations
RC_Species_Loc_Year <- RC_Species_Loc_Year %>% 
  relocate(Phylum:Specieslist, .before = BelgianCoast_2018) %>% 
  mutate_at(c(6:ncol(RC_Species_Loc_Year)), as.numeric) %>%
  group_by(Phylum, Class, Order, Family, Specieslist,) %>% 
  summarise(across(BelgianCoast_2018:Vigo_2019, sum))

#Presence/Absence matrix for Species per ARMS
Pres_Abs <- Read_Count_Species_ARMS
Pres_Abs[2:ncol(Pres_Abs)][Pres_Abs[2:ncol(Pres_Abs)] > 0] <- 1

#Read count per species per fraction
Read_Count_Species_Fraction <- Data_LOI  %>%
  dplyr::select(-sequence, -Phylum, -Order, -Class, -Family, -Genus, -Species, -Similarity) %>%
  dplyr::select(Specieslist, sort(colnames(.))) %>%
  group_by(Specieslist) %>%
  summarise_all(sum)

#Export
dir.create(file.path(dir, "Output"))
setwd(file.path(dir,"Output"))
write.csv(Uniques, "UniqueSpecies.csv", row.names = FALSE)
write.csv(ASV_Count, "ASV_Count.csv", row.names = FALSE)
write.csv(Read_Count_ASVs_ARMS, "Read_Count_ASVs_ARMS.csv", row.names = FALSE)
write.csv(Read_Count_Species_ARMS, "Read_Count_Species_ARMS.csv", row.names = FALSE)
write.csv(RC_Species_Loc_Year, "Read_Count_Species_Loc_Year.csv", row.names = FALSE)
write.csv(Read_Count_Species_Fraction, "Read_Count_Species_Fraction.csv", row.names = FALSE)
write.csv(Data_LOI, "Data_LOI.csv", row.names = FALSE)

#### Build dataframe to be used in AlienIdentification.R ####
Species_Location <- Data_LOI  %>%
  dplyr::select(-sequence, -Phylum, -Order, -Class, -Family, -Genus, -Species, -Similarity) %>%
  dplyr::select(Specieslist, sort(colnames(.))) %>%
  group_by(Specieslist) %>%
  summarise_all(sum)
colnames(Species_Location) <- c("Specieslist", MetaData$Observatory.ID)
df <- Species_Location
#Merge columns with duplicate names and sum contents
df2 <- sapply(unique(colnames(df)[duplicated(colnames(df))]), function(x) rowSums(df[,grepl(paste(x, "$", sep=""), colnames(df))]))
Species_Location <- as.data.frame(cbind(df2, df[,!duplicated(colnames(df)) & !duplicated(colnames(df), fromLast = TRUE)]))
Species_Location <- Species_Location %>% dplyr::select(Specieslist, everything())
rm(df, df2)
write.csv(Species_Location, "Species_Location.csv", row.names = F)

#Export Coordinates, which are used in 3_Main_script and MetaData, which is used in 4_Visualisation
setwd(dir)
Coordinates <- MetaData %>% group_by(Observatory.ID) %>% dplyr::select(Observatory.ID, Longitude, Latitude) %>% slice_min(Longitude) %>% unique
write.csv(Coordinates, "Inputs/Coordinates.csv", row.names = F)
write.csv(MetaData, "Inputs/MetaData_Adjusted.csv", row.names = F)

#### Clean R environment ####
rm(list = ls())
