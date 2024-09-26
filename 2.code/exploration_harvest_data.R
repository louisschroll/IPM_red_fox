
library(readxl)
library(tidyverse)


harvest_fox_10 <- read_excel("0.raw_data/prelevements/prvtsCONNUSFDC10_copie.xlsx") %>% 
  janitor::clean_names() %>% 
  mutate(total = as.numeric(total))

colnames(harvest_fox_10)

# How many foxes were killed / season?

harvest_fox_10 %>% pull(saison) %>% 
  unique()

harvest_fox_10 %>% 
  filter(comm != "Plus") %>% 
  group_by(saison) %>% 
  summarise(sum_kill = sum(total),
            mean_kill = mean(total),
            med_kill = median(total)) %>% 
  ggplot(aes(x = saison, y = sum_kill)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "none")

harvest_fox_10 %>% 
  filter(comm != "Plus") %>% 
  ggplot(aes(x = saison, group = saison, y = total)) +
  geom_violin(trim = TRUE, fill='#A4A4A4', color="darkred") +
  geom_boxplot(width = 0.07) +
  geom_jitter(position = position_jitter(0.2), size = 1) +
  theme_minimal()


# How many foxes / places?
harvest_fox_10 %>% pull(comm) %>% unique()
harvest_fox_10 %>% pull(com) %>% unique()
harvest_fox_10 %>% pull(gic) %>% unique()

all(harvest_fox_10$com == harvest_fox_10$comm)

harvest_fox_10 %>% 
  ggplot(aes(x = comm, y = total, group = comm)) +
  geom_boxplot() +
  theme_minimal()

harvest_fox_10 %>% 
  ggplot(aes(x = com, y = total, group = com)) +
  geom_boxplot() +
  theme_minimal()

harvest_fox_10 %>% 
  ggplot(aes(x = gic, y = total, group = gic)) +
  geom_boxplot() +
  theme_minimal()

harvest_fox_10 %>% 
  filter(comm != "Plus") %>% 
  ggplot(aes(x = saison, y = total, group = comm, color = comm)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "none")


# Foxes killed by each methods

harvest_fox_10 %>% 
  select(piegeur:route_ni) %>% 
  pivot_longer(cols = everything(), names_to = "method", values_to = "nb_kill") %>% 
  mutate(nb_kill = as.numeric(nb_kill)) %>% 
  filter(nb_kill > 0) %>% 
  ggplot(aes(x = method, y = nb_kill, fill = method)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(0.2), size = 1) +
  theme_minimal()

harvest_fox_10 %>% 
  select(piegeur:route_ni) %>% 
  pivot_longer(cols = everything(), names_to = "method", values_to = "nb_kill") %>% 
  mutate(nb_kill = as.numeric(nb_kill)) %>% 
  group_by(method) %>% 
  summarise(sum_kill = sum(nb_kill),
            mean_kill = mean(nb_kill),
            mean_kill2 = mean(nb_kill[nb_kill > 0]))
