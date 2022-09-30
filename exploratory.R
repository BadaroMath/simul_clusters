library(ggplot2)
library(dplyr)

################  VARIACAO DO OMEGA ##############


baromg <- read.csv("clusters_index_ BarOmega .txt", 
                   sep= ",", 
                   header = F)


head(baromg)
colnames(baromg) <- c("Obs", "Clusters", "Comp", 
                      "BarOmega", "RI", "ClassProp", 
                      "VI","TimeElapsed", "Method")




baromg %>% ggplot(aes(x = BarOmega, y = RI, color = Method)) +
  #geom_line(size=1) +
  theme_classic()+
  stat_smooth(method = "lm",
              formula = y ~ poly(x, 10), se = FALSE)


baromg %>% ggplot(aes(x = BarOmega, y = ClassProp, color = Method)) +
  #geom_line(size=1) +  theme_classic()+
  stat_smooth(method = "lm",
              formula = y ~ poly(x, 10), se = FALSE)

baromg %>% ggplot(aes(x = BarOmega, y = VI, color = Method)) +
  #geom_line(size=1) +  theme_classic() +
  stat_smooth(method = "lm",
              formula = y ~ poly(x, 10), se = FALSE)+  
  theme_classic()



################  VARIACAO DE CLUSTERS ##############


clusters <- read.csv("clusters_index_ Clusters .txt", 
                   sep= ",", 
                   header = F)


head(clusters)
colnames(clusters) <- c("Obs", "Clusters", "Comp", 
                      "BarOmega", "RI", "ClassProp", 
                      "VI","TimeElapsed", "Method")




clusters %>% ggplot(aes(x = Clusters, y = RI, color = Method)) +
  geom_line(size=1) +
  theme_classic()
  #stat_smooth(method = "lm",
  #            formula = y ~ poly(x, 8), se = FALSE)
#

clusters %>% ggplot(aes(x = Clusters, y = ClassProp, color = Method)) +
  geom_line(size=1) +  
  theme_classic()
  #stat_smooth(method = "lm",
  #            formula = y ~ poly(x, 8), se = FALSE)

clusters %>% ggplot(aes(x = Clusters, y = VI, color = Method)) +
  geom_line(size=1) +  theme_classic() +
  #stat_smooth(method = "lm",
  #            formula = y ~ poly(x, 10), se = FALSE)+  
  theme_classic()


################  VARIACAO DE CLUSTERS ##############


comps <- read.csv("clusters_index_ Components .txt", 
                     sep= ",", 
                     header = F)


head(comps)
colnames(comps) <- c("Obs", "Clusters", "Comp", 
                        "BarOmega", "RI", "ClassProp", 
                        "VI","TimeElapsed", "Method")




comps %>% ggplot(aes(x = Comp, y = RI, color = Method)) +
  geom_line(size=1) +
  theme_classic()
#stat_smooth(method = "lm",
#            formula = y ~ poly(x, 8), se = FALSE)
#

comps %>% ggplot(aes(x = Comp, y = ClassProp, color = Method)) +
  geom_line(size=1) +  
  theme_classic()
#stat_smooth(method = "lm",
#            formula = y ~ poly(x, 8), se = FALSE)

comps %>% ggplot(aes(x = Comp, y = VI, color = Method)) +
  geom_line(size=1) +  theme_classic() +
  #stat_smooth(method = "lm",
  #            formula = y ~ poly(x, 10), se = FALSE)+  
  theme_classic()

################  VARIACAO DE Obs ##############


obs <- read.csv("clusters_index_ N .txt", 
                  sep= ",", 
                  header = F)


head(obs)
colnames(obs) <- c("Obs", "Clusters", "Comp", 
                     "BarOmega", "RI", "ClassProp", 
                     "VI","TimeElapsed", "Method")




obs %>% ggplot(aes(x = Obs, y = RI, color = Method)) +
  #geom_line(size=1) +
  theme_classic()+
  stat_smooth(method = "lm",
              formula = y ~ poly(x, 8), se = FALSE)
  

obs %>% ggplot(aes(x = Obs, y = ClassProp, color = Method)) +
  geom_line(size=1) +  
  theme_classic()
#stat_smooth(method = "lm",
#            formula = y ~ poly(x, 8), se = FALSE)

obs %>% ggplot(aes(x = Obs, y = VI, color = Method)) +
  geom_line(size=1) +  theme_classic() +
  #stat_smooth(method = "lm",
  #            formula = y ~ poly(x, 10), se = FALSE)+  
  theme_classic()


