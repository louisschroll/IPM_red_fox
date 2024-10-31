
library(readxl)
library(tidyverse)

### Look at the line transect data from ille-et-Vilaine (35)

GIC_Domagne_fox <- read_excel(
  path = "0.raw_data/comptages/GICDOMAGNETOT.xls",
  sheet = "distren"
)

GIC_Fougeres_fox <- read_excel(
  path = "0.raw_data/comptages/GICFOUGERESTOT.xls",
  sheet = "distren"
) %>% 
  rename(effort = eff)

GIC_Vendelais_fox <- read_excel(
  path = "0.raw_data/comptages/GICVENDELAISTOT.xls",
  sheet = "distren"
) %>% 
  rename(effort = eff) 

# Indicate which GIC you want to look at
GIC_fox <- GIC_Vendelais_fox

# How many transect?
GIC_fox %>% pull(transect) %>% unique()

# How many obs / transect?
GIC_fox %>% 
  group_by(transect) %>% 
  count()

mean_count <- GIC_fox %>% 
  group_by(transect) %>% 
  count() %>% 
  pull(n) %>% 
  mean()

GIC_fox %>% 
  ggplot(aes(x = transect)) +
  geom_bar() +
  geom_hline(yintercept = mean_count, color = "red") +
  theme_bw()

GIC_fox %>% 
  group_by(transect, annee) %>% 
  count()  %>% 
  mutate(annee = as.factor(annee)) %>% 
  ggplot(aes(x = transect, y = n, fill = annee, label = n)) +
  geom_bar(stat = "identity") +
  theme_bw()

# Is the effort always the same?
GIC_fox %>% pull(effort) %>% summary()

# How many obs / year ?
GIC_fox %>% 
  group_by(annee) %>% 
  count() %>% 
  left_join(GIC_fox %>% 
              group_by(annee) %>% 
              summarise(efftot = sum(effort)),
            by = join_by(annee)) %>% 
  mutate(n_adjusted = n/efftot*600) %>% 
  pivot_longer(-c(annee, efftot)) %>% 
  ggplot(aes(x = annee, y = value, color = name)) +
  geom_point(pch = 4) +
  geom_line(linewidth=0.8) +
  theme_classic()
  
### look at the point transect data from Aube (10)

GIC_Barrois_fox <- read_excel(
  path = "0.raw_data/comptages/pointtransectrenardsAUBE_aout13.xls",
  sheet = "Barrois"
)

GIC_Sarce_fox <- read_excel(
  path = "0.raw_data/comptages/pointtransectrenardsAUBE_aout13.xls",
  sheet = "Sarce"
)

GIC_fox <- GIC_Sarce_fox

# What are the columns and what do they contain
GIC_Barrois_fox %>% colnames()

for (coln in colnames(GIC_fox)){
  print(coln)
  print(GIC_fox %>% pull(coln) %>% unique())
  print(" ")
}


# How many obs / year ?
GIC_fox %>% 
  group_by(annee) %>% 
  count() %>% 
  left_join(GIC_fox %>% 
              group_by(annee) %>% 
              summarise(efftot = sum(eff150m)),
            by = join_by(annee)) %>% 
  mutate(n_adjusted = n/efftot*120) %>% 
  pivot_longer(-c(annee, efftot)) %>% 
  ggplot(aes(x = annee, y = value, color = name)) +
  geom_point(pch = 4) +
  geom_line(linewidth=0.8) +
  theme_classic()
