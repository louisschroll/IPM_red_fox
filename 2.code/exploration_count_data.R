
library(readxl)
library(tidyverse)

GIC_Domagne_fox <- read_excel(
  path = "0.raw_data/comptages/GICDOMAGNETOT.xls",
  sheet = "distren"
)

# How many transect?
GIC_Domagne_fox %>% pull(transect) %>% unique()

# How many obs / transect?
GIC_Domagne_fox %>% 
  group_by(transect) %>% 
  count()

mean_count <- GIC_Domagne_fox %>% 
  group_by(transect) %>% 
  count() %>% 
  pull(n) %>% 
  mean()

GIC_Domagne_fox %>% 
  ggplot(aes(x = transect)) +
  geom_bar() +
  geom_hline(yintercept = mean_count, color = "red") +
  theme_bw()

GIC_Domagne_fox %>% 
  group_by(transect, annee) %>% 
  count()  %>% 
  mutate(annee = as.factor(annee)) %>% 
  ggplot(aes(x = transect, y = n, fill = annee, label = n)) +
  geom_bar(stat = "identity") +
  theme_bw()

# Is the effort always the same?
GIC_Domagne_fox %>% pull(effort) %>% summary()
