library(ggplot2)
library(dplyr)

################ OMEGA = 0.2 ##############


################ OMEGA = 0.5 ##############


################ OMEGA = 0.7 ##############

km <- read.csv("kmeans_index.txt", sep= ",", header = F)
km["Metodo"] <- "K-Means" 
mclust <- read.csv("mclust_index.txt", sep= ",", header = F)
mclust["Metodo"] <- "Model Based"
spec <- read.csv("spec_index.txt", sep= ",", header = F)
spec["Metodo"] <- "Espectral"
dbscan <- read.csv("dbscan_index.txt", sep= ",", header = F)
dbscan["Metodo"] <- "DBSCAN"
fcm <- read.csv("fcm_index.txt", sep= ",", header = F)
fcm["Metodo"] <- "Fuzzy C-Means"

df <- as.data.frame(rbind(km, mclust, spec, fcm))

colnames(df) <- c("n", "k", "p", "Omega", "RI", "Prop", "VI", "Metodo")




df %>% ggplot(aes(x = Omega, y = RI, color = Metodo)) +
  geom_line(size=1) +
  theme_classic()
#stat_smooth(method = "lm",
#            formula = y ~ poly(x, 10), se = FALSE)


df %>% ggplot(aes(x = Omega, y = Prop, color = Metodo)) +
  geom_line(size=1) +  theme_classic()
#stat_smooth(method = "lm",
#            formula = y ~ poly(x, 10), se = FALSE)

df %>% ggplot(aes(x = Omega, y = VI, color = Metodo)) +
  #geom_line(size=1) +  theme_classic() +
  stat_smooth(method = "lm",
              formula = y ~ poly(x, 10), se = FALSE)+  
  theme_classic()

