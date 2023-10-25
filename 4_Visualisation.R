#Load necessary packages
pacman::p_load(dplyr, tidyverse, phyloseq, rstatix, ggpubr, vegan, ggplot2)

#Set wd to wherever all files are stored
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Rarefaction curves ####
#Plot rarefaction results
read_count_species_arms <- read.csv("Output/Read_Count_Species_ARMS.csv")
otu_table_bar <- read_count_species_arms
otu_table_bar <- column_to_rownames(otu_table_bar, "specieslist")
otu_table_bar <- as.matrix(t(otu_table_bar))

rarefaction_df <- rarecurve(otu_table_bar, step = 20,
                            col = "blue",
                            cex = 0.6, tidy = TRUE)
ggplot(rarefaction_df, aes(x = Sample, y = species, group = Site)) +
  geom_line(color = "blue", alpha = 0.5) +
  geom_vline(xintercept = 100, linetype = "dashed", color = "black",
             alpha = 0.6, size = 0.6) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 80, 20)) +
  scale_x_continuous(expand = c(0, 500)) +
  ylab("Number of species") +
  xlab("Number of sequences")

#### alien read fraction and alien species fraction ####
gbif_distance <- read.csv("Output/DistanceOverSea.csv")
# Remove duplicates
gbif_distance <- unique(gbif_distance)
# Convert distance to km
gbif_distance$distance <- gbif_distance$distance / 1000
# Assign alien status
alien_distance_threshold <- 250
# Most NAs are assigned due to a bug in the ShortestPath function, where
# the sampling location and GBIF occurrence are too close to be measured.
# Therefore, change al NAs to the maximum potential error from the ShortestPath
# and assume these are native species
gbif_distance$distance[is.na(gbif_distance$distance)] <- 16
gbif_distance$alien_status <- ifelse(gbif_distance$distance >
                                       alien_distance_threshold, "YES", "NO")
rm(alien_distance_threshold)
#Remove lines where distance is not calculated
gbif_distance <- gbif_distance[complete.cases(gbif_distance), ]
#Create df with only aliens
only_aliens <- gbif_distance %>%
  filter(alien_status == "YES") %>%
  select(specieslist, name, value, distance)
only_natives <- gbif_distance %>%
  filter(alien_status == "NO") %>%
  select(specieslist, name, value, distance)
rm(gbif_distance)
#Create factor to be able to assign read count and metadata to alien status
alien_species_obs <- paste0(only_aliens$specieslist, "_", only_aliens$name)
#Spread df
only_aliens <- spread(data = only_aliens, key = name, value = value)
#Prepare main dataframe with read counts to assign alien status
rc_species_obs_year <- read_count_species_arms %>%
  ungroup() %>%
  select(specieslist:ncol(read_count_species_arms))
rc_species_obs_year <- rc_species_obs_year %>%
  gather(key = "No_fraction", value = "count", 2:ncol(rc_species_obs_year))
metadata <- read.csv("Inputs/MetaData_Adjusted.csv")
meta_rc <- metadata %>% select(-fraction, -Filename)
meta_rc <- left_join(rc_species_obs_year, meta_rc, by = "No_fraction")
meta_rc$species_Obs <- paste0(meta_rc$specieslist, "_", meta_rc$Observatory.ID)
#Filter main dataframe to one with only alien species
meta_rc_alien <- meta_rc %>% filter(species_Obs %in% alien_species_obs)
#Remove dupes in order to spread
meta_rc <- meta_rc[!duplicated(meta_rc), ]
meta_rc_alien <- meta_rc_alien[!duplicated(meta_rc_alien), ]
#Calculate alien species fraction
speciescount <- meta_rc %>%
  filter(count > 0) %>%
  group_by(No_fraction) %>%
  dplyr::summarise(Non_alien = n_distinct(specieslist))
alienspeciescount <- meta_rc_alien %>%
  filter(count > 0) %>%
  group_by(No_fraction) %>%
  dplyr::summarise(alien = n_distinct(specieslist))
alienspeciesfraction <- left_join(speciescount, alienspeciescount,
                                  by = "No_fraction")
alienspeciesfraction$alienspeciesfraction <- alienspeciesfraction$alien /
  (alienspeciesfraction$Non_alien + alienspeciesfraction$alien)
rm(speciescount, alienspeciescount)
#Spread dfs
meta_rc <- meta_rc %>%
  select(-species_Obs) %>%
  spread(specieslist, count)
meta_rc_alien <- meta_rc_alien %>%
  select(-species_Obs) %>%
  spread(specieslist, count)

#Fill empty cells with 0
meta_rc[is.na(meta_rc)] <- 0
meta_rc_alien[is.na(meta_rc_alien)] <- 0
#Calculate alien read fraction
totalreadcount <- rowSums(meta_rc[, 15:ncol(meta_rc)])
alienreadcount <- rowSums(meta_rc_alien[, 15:ncol(meta_rc_alien)])
alienreadfraction <- alienreadcount / totalreadcount
meta_rc$alienreadfraction <- alienreadfraction
meta_rc <- left_join(meta_rc, alienspeciesfraction, by = "No_fraction")
meta_rc <- meta_rc %>%
  select(alienreadfraction, alienspeciesfraction, Non_alien, alien,
         everything())
rm(totalreadcount, alienreadcount)

#Make dataframe which filters on at least 100 summed reads per ARMS
meta_rc_100 <- meta_rc
meta_rc_100 <- meta_rc_100[complete.cases(meta_rc_100), ]
meta_rc_100 <- filter(meta_rc_100,
                      rowSums(meta_rc_100[, 19:ncol(meta_rc_100)]) > 100)

##### NMDS plot BOLDigger ARMS #####
meta_rc[is.na(meta_rc)] <- 0
meta_rc_alien[is.na(meta_rc_alien)] <- 0
#Only keep complete cases
meta_rc_1000 <- meta_rc
meta_rc_1000 <- meta_rc_1000[complete.cases(meta_rc_1000), ]
#Filter on RC over 1000
minimum_rc <- 1000
meta_rc_1000 <- filter(meta_rc_1000,
                       rowSums(meta_rc_1000[, 19:ncol(meta_rc_1000)]) >
                         minimum_rc)
#Select species composition and environmental variables
com <- as.data.frame(meta_rc_100[, 19:ncol(meta_rc_100)])
env <- meta_rc_100 %>%
  select(alienreadfraction, alienspeciesfraction, country, Observatory.ID,
         Latitude, Longitude, Depth_m, Monitoring_area, Year,
         Sample_region_country)
com <- sapply(com, as.numeric)

#Calculate nmds
nmds <- metaMDS(com, distance = "bray")
en <- envfit(nmds, env, permutations = 999, na.rm = TRUE)
#Adjust dataframe for nicer plot
data_scores <- scores(nmds)
data_scores <- as.data.frame(data.scores$sites)
#Select one environmental variable to focus on
data_scores$Sample_region_country <- env$Sample_region_country
en_coord_cont <- as.data.frame(scores(en, "vectors")) * ordiArrowMul(en) * 3
en_coord_cat <- as.data.frame(scores(en, "factors"))
palette <- pals::cols25(length(unique(data.scores$Sample_region_country)))
#Plot NMDS focussing on the chosen variable
ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = data.scores, aes(colour = Sample_region_country),
             size = 5, alpha = 0.7) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey30"),
        axis.ticks = element_blank(), axis.text = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(size = 10, face = "bold",
                                    colour = "grey30"),
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Sample region (country)") +
  scale_color_manual(values = palette)


##### STATISTICAL ANALYSES AND VISUALISATION ALIEN READ/SPECIES FRACTION #####
comparisons_mon <- list(c("Industrial port", "Marina"),
                        c("Industrial port", "LHI"),
                        c("Industrial port", "MPA"),
                        c("Marina", "LHI"), c("Marina", "MPA"),
                        c("LHI", "MPA"))
ggline(data = meta_rc_1000, x = "Monitoring_area", y = "alienspeciesfraction",
       add = c("mean_se", "jitter")) +
  stat_compare_means(ref.group = "MPA")

ggplot(meta_rc_1000, aes(Monitoring_area, alienspeciesfraction,
                         color = country)) +
  geom_boxplot() +
  geom_jitter() +
  stat_compare_means(ref.group = "Marina") +
  theme_pubr() +
  scale_color_manual(values = palette)

##### alien species versus native #####
#Measure difference of alienspeciesfraction between Monitoring_area
stat_test <- wilcox_test(meta_rc, alienspeciesfraction ~ Monitoring_area) %>%
  add_xy_position()
#plot
ggboxplot(meta_rc, x = "Monitoring_area", y = "alienspeciesfraction",
          palette = palette, order = c(unique(stat_test$group1), "MPA"), ) +
  stat_pvalue_manual(stat_test, hide.ns = FALSE) +
  scale_y_continuous(expand = c(0, 0.01))

#Stretch df to fit barplot
barplot_df <- meta_rc_100 %>%
  select(No_fraction, Non_alien, alien, Monitoring_area, ) %>%
  arrange(Monitoring_area, Non_alien)
#Save PlotOrder
barplot_df <- barplot_df %>%
  arrange(factor(Monitoring_area,
                 levels = c("Industrial port", "Marina", "LHI", "MPA")))
plot_order_ma <- barplot_df$Monitoring_area
barplot_df <- pivot_longer(data = barplot_df, cols = c(Non_alien, alien),
                           names_to = "species_status")
#Create the plot
ggplot(barplot_df, aes(x = No_fraction, y = value, fill = species_status)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_x_discrete(labels = plot_order_ma) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.2),
    text = element_text(size = 30)
  ) +
  aes(x = forcats::fct_inorder(No_fraction)) +
  xlab("Monitoring area") +
  ylab("Number of species") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 80, 20)) +
  labs(fill = "species status") +
  scale_fill_manual(labels = c("alien", "Native"),
                    values = c("grey60", "grey40"))


##### Phyloseq #####
#OTU Table
library("tidyverse")
read_count_species_fraction <-
  read.csv("Output/Read_Count_Species_Fraction.csv")
otu_table <- read_count_species_fraction
otu_table <- column_to_rownames(otu_table, "specieslist")
otu_table <- as.matrix(t(otu_table))

meta_phylo <- metadata
rownames(metadata) <- metadata$Filename

#Taxonomic info
data_loi <- read.csv("Output/Data_LOI.csv")
tax_info <- data_loi %>%
  select(Phylum:specieslist, -species) %>%
  distinct()
tax_info <- tax_info[!duplicated(tax_info$specieslist), ]
rownames(tax_info) <- tax_info$specieslist
tax_info <- as.matrix(tax_info)

#Combine into Phyloseq object
physeq <- phyloseq(otu_table(otu_table, taxa_are_rows = FALSE),
                   tax_table(tax_info),
                   sample_data(metadata))

#Add Plot Order
plot_order_phy <- metadata %>%
  arrange(Sample_region_country, Observatory.ID, Year) %>%
  select(Filename) %>%
  distinct()
plot_order_phy <- as.list(plot_order_phy)
sample_data(physeq)$NewID <- factor(sample_names(physeq))
sample_data(physeq)$NewID <- factor(sample_data(physeq)$NewID,
                                    levels = plot_order_phy$Filename)

#Remove samples with less than 100 reads
physeq <- prune_samples(sample_sums(physeq) > 100, physeq)

#Merge samples based on variable by summing read counts
palette_src <- pals::cols25(length(unique(meta_rc$Sample_region_country)))
#Merge by Sample_region_country...
mergedphyseq_src <- merge_samples(physeq, "Sample_region_country")
plot_heatmap(mergedphyseq_src, taxa.order = "Phylum", max.label = 500,
             taxa.label = "Phylum")
#... or Monitoring_area
mergedphyseq_ma <- merge_samples(physeq, "Monitoring_area")
plot_heatmap(mergedphyseq_ma, taxa.order = "Phylum", max.label = 500)

#### Diversity of native species ####
#Plot diversity
richness <- estimate_richness(physeq,
                              measures = c("Observed", "Chao1", "Shannon"))
richness <- left_join(rownames_to_column(richness),
                      rownames_to_column(metadata), by = "rowname")
richness <- richness %>%
  arrange(factor(Monitoring_area,
                 levels = c("Industrial port", "Marina", "LHI", "MPA")))
#Measure difference of Observed between Monitoring_area
stat_test <- wilcox_test(richness, Observed ~ Monitoring_area) %>%
  add_xy_position()
#plot diversity of entire species pool
ggboxplot(richness, x = "Monitoring_area", y = "Observed", palette = palette) +
  stat_pvalue_manual(stat_test, hide.ns = TRUE, label = "p.adj", size = 10) +
  theme(text = element_text(size = 30))

#plot diversity of just native species
stat_test <- wilcox_test(meta_rc_100, Non_alien ~ Monitoring_area) %>%
  add_xy_position()
ggboxplot(meta_rc_100, x = "Monitoring_area", y = "Non_alien",
          palette = palette) +
  stat_pvalue_manual(stat_test, hide.ns = TRUE, label = "p.adj", size = 10)

#### alien phyloseq prep ####
aliens_otu <- read_count_species_fraction
aliens_otu <- aliens_otu %>%
  filter(specieslist %in% unique(only_aliens$specieslist))
aliens_otu <- column_to_rownames(aliens_otu, "specieslist")
aliens_otu[aliens_otu > 0] <- 1
aliens_otu <- as.matrix(aliens_otu)

#Add alien_Type to Taxonomic info
as_the_crow_flies <- read.csv("Output/ShortestPath.csv")
#Extract observation count, to re-add after gather function
observation_count <- select(as_the_crow_flies,
                            c("Speciesname", "total.observations",
                              "Unique.locations"))
as_the_crow_flies <- as_the_crow_flies %>%
  select(-total.observations, -Unique.locations, -Errors) %>%
  gather(key = "Observatory.ID", value = "distance_atcf", 2:13)
as_the_crow_flies <- left_join(as_the_crow_flies, observation_count,
                               by = "Speciesname")
#Convert meters to kilometres for easier reading
as_the_crow_flies$distance_atcf <- as.numeric(as_the_crow_flies$distance_atcf)
as_the_crow_flies$distance_atcf <- ceiling(as_the_crow_flies$distance_atcf /
                                             1000)
#Find species with no observation on GBIF
cryptogenic <- as_the_crow_flies[is.na(as_the_crow_flies$total.observations), ]
#Calculate mean distance to closes observation in order to determine alien_Type
mean <- as_the_crow_flies %>%
  na.omit %>%
  group_by(speciesname) %>%
  summarise_at(vars(distance_atcf), list(name = mean))
as_the_crow_flies <- left_join(as_the_crow_flies, mean, by = "Speciesname")
as_the_crow_flies <- as_the_crow_flies[as_the_crow_flies$speciesname %in%
                                         only_aliens$specieslist, ]
as_the_crow_flies$alien_Type <- ifelse(as_the_crow_flies$total.observations
  < 50 |  is.na(as_the_crow_flies$total.observations),
  as_the_crow_flies$alien_Type <- "Cryptogenic",
  ifelse(as_the_crow_flies$name > 2500,
         as_the_crow_flies$alien_Type <-
           "Hitchhiker", "Range expander")
)
#DF with alien types
alien_type <- as_the_crow_flies %>%
  select(speciesname, alien_type) %>%
  unique() %>%
  arrange(speciesname)
#Add species with no gbif occurences
cryptogenic <- as.data.frame(unique(cryptogenic$speciesname))
cryptogenic$alien_Type <- "Cryptogenic"
colnames(cryptogenic) <- c("Speciesname", "Alien_Type")

# Plot native biodiversity and alien species fraction
library("SciViews")
palette_4 <- c("#440154", "#12c5ed", "#fde725", "#0fb859")
ggplot(meta_rc_100, aes(Non_alien, alienspeciesfraction,
                        color = Monitoring_area)) +
  geom_point(size = 8, alpha = 0.7) +
  stat_smooth(method = lm, formula = y ~ ln(x), inherit.aes = FALSE,
              data = meta_rc_100,
              mapping = aes(Non_alien, alienspeciesfraction)) +
  stat_regline_equation(formula = y ~ ln(x), inherit.aes = FALSE,
                        data = meta_rc_100, label.x = 60, label.y = 0.33,
                        mapping = aes(Non_alien, alienspeciesfraction,
                                      label = ..adj.rr.label..),
                        size = 8) +
  stat_regline_equation(formula = y ~ ln(x), inherit.aes = FALSE,
                        data = meta_rc_100, label.x = 60, label.y = 0.37,
                        mapping = aes(Non_alien, alienspeciesfraction,
                                      label = ..eq.label..), size = 8) +
  theme_pubclean(base_size = 30) +
  scale_color_manual(values = palette_4) +
  xlab("species with confirmed presence") +
  ylab("New alien species fraction")

#Add alien_type to Taxonomix info, prepare for phyloseq
aliens_tax <- as.data.frame(tax_info)
aliens_tax <- aliens_tax %>% filter(specieslist %in% alien_Type$speciesname)
aliens_tax <- left_join(aliens_tax, alien_Type,
                        by = c("specieslist" = "Speciesname"))
rownames(aliens_tax) <- aliens_tax$specieslist
alien_order <- aliens_tax %>% arrange(alien_Type) %>% select(specieslist)
aliens_tax <- as.matrix(aliens_tax)

#### Phyloseq of aliens ####
alien_physeq <- phyloseq(otu_table(aliens_otu, taxa_are_rows = TRUE),
                         tax_table(aliens_tax),
                         sample_data(metadata))

#Add plotorder of samples
plot_order_al <- metadata %>%
  arrange(factor(Monitoring_area, levels =
                   c("Industrial Port", "Marina", "LHI", "MPA"))) %>%
  select(Filename)
sample_data(alien_physeq)$NewID <- factor(sample_names(alien_physeq))
sample_data(alien_physeq)$NewID <- factor(sample_data(alien_physeq)$NewID,
                                          levels = plot_order_al$Filename)

#Remove samples with less than 100 alien reads
alien_physeq <- prune_samples(sample_sums(alien_physeq) > 100, alien_physeq)

#Plot heatmap of SRC
alienphyseq_src <- merge_samples(alien_physeq, "Sample_region_country")
plot_heatmap(alienphyseq_src, taxa.order = alien_order$specieslist,
             taxa.label = "Alien_Type", )

#Plot heatmap of MA
alienphyseq_ma <- merge_samples(alien_physeq, "Monitoring_area")
plot_heatmap(alienphyseq_ma, taxa.order = alien_order$specieslist,
             taxa.label = "Alien_Type")

#Merge alien types
range_expander <- alien_Type[alien_Type$alien_Type == "Range expander", ]
hitchhiker <- alien_Type[alien_Type$alien_Type == "Hitchhiker", ]
cryptogenic <- alien_Type[alien_Type$alien_Type == "Cryptogenic", ]

#Merge alien types and Monitoring areas
x1 <- merge_taxa(alienphyseq_ma, range_expander$speciesname, 2)
x2 <- merge_taxa(x1, hitchhiker$speciesname, 2)
x3 <- merge_taxa(x2, cryptogenic$speciesname, 2)
#plot heatmap of merged monitoring type and merged alien type
plot_heatmap(x3, taxa.order = c("Acartia hudsonica", "Cephalothrix spiralis",
                                "Amblyosyllis madeirensis"),
             sample.order = c("Industrial port", "Marina", "LHI", "MPA")) +
  scale_y_discrete(labels = c("Range expanders", "Hithchikers", "Cryptogenic"))

#Statistics on alien type
alien_type_stats <- as.data.frame(otu_table(x3))
colnames(alien_type_stats) <- c("Range expanders", "Cryptogenic", "Hitchhikers")
chisq.test(alien_type_stats)
alien_type_stats <- rownames_to_column(alien_type_stats, "Monitoring_area")
alien_type_stats <- pivot_longer(alien_type_stats, cols = c("Range expanders",
                                                            "Cryptogenic",
                                                            "Hitchhikers"))
