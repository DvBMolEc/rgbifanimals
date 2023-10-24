#Input data
straight_line <- read.csv("Output/ShortestPath.csv")
through_sea <- read.csv("Output/DistanceOverSea.csv")

test <- merge(through_sea, straight_line, by.y = "Speciesname",
              by.x = "Specieslist")
for (item in levels(test$name)) {
  test$straight[test$name == item] <- test[test$name == item, item]
}

#remove unused sampling locations
test <- test[, !names(test) %in% c(levels(test$name), "Galway", "Gdynia")]

#plot
ggplot(test) +
  geom_point(aes(x = straight, y = distance), col = "green") +
  geom_abline()

ggplot(test) +
  geom_density(aes(x = log10(straight), fill = is.na(distance)), alpha = 0.5)
