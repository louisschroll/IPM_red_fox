
library(readxl)
reprod_age_data_10 <- read_excel("0.raw_data/reprod-age/reprodFDC10_copie.xlsx") %>% 
  janitor::clean_names() %>% 
  mutate(mode_prelt = case_when(
    mode_prelt == "A" ~ "A",
    mode_prelt == "C" ~ "Chasse",
    mode_prelt == "D" ~ "Deterrage",
    mode_prelt == "G" ~ "GP",
    mode_prelt == "P" ~ "Piegeage",
    mode_prelt == "R" ~ "R",
    mode_prelt == "S" ~ "Sauteux",
    mode_prelt == "T" ~ "Tir",
    mode_prelt == "V" ~ "Voiture"
  ),
  agedents = as.numeric(agedents))

reprod_age_data_35 <- read_excel("0.raw_data/reprod-age/reprodFDC35_copie.xlsx") %>% 
  janitor::clean_names() %>% 
  mutate_at(across(poids,age), as.numeric) %>% 
  rename(mode_prelt = mode_de_prel)

colnames(reprod_data_10)
colnames(reprod_data_35)

nrow(reprod_data_10)
nrow(reprod_data_35)

# What are the sex-ratio
reprod_data_10$sexe %>% unique()
reprod_data_35$sexe %>% unique() # female only !

reprod_age_data_10 %>% 
  mutate(sexe = as.factor(sexe)) %>% 
  ggplot(aes(x = sexe)) +
  geom_bar() +
  theme_bw()


reprod_age_data_10 %>% 
  mutate(mode_prelt = as.factor(mode_prelt)) %>% 
  ggplot(aes(x = mode_prelt, fill = sexe)) +
  geom_bar() +
  theme_bw()

# What are the age
reprod_age_data_10 %>% pull(age) %>% unique()
reprod_age_data_10 %>% pull(agedents) %>% unique()

reprod_age_data_10 %>% select(age, agedents) %>% filter(!is.na(age)) %>% print(n=20)
all(reprod_age_data_10$age == reprod_age_data_10$agedents)
reprod_age_data_10 %>% select(age, agedents) %>% filter(age != agedents)

reprod_age_data_10 %>% 
  ggplot(aes(x = age)) +
  geom_bar()

reprod_age_data_35 %>% pull(age) %>% unique()

# How are age class distributed for each culling method ?
reprod_age_data_10 %>% 
  mutate(mode_prelt = as.factor(mode_prelt)) %>% 
  ggplot(aes(x = mode_prelt, fill = as.factor(age))) +
  geom_bar() +
  theme_bw()

reprod_age_data_10 %>% 
  mutate(mode_prelt = as.factor(mode_prelt)) %>% 
  filter(!is.na(age)) %>% 
  ggplot(aes(x = mode_prelt, fill = as.factor(age))) +
  geom_bar() +
  theme_bw()

reprod_age_data_10 %>% 
  select(age, mode_prelt) %>% 
  mutate(mode_prelt = as.factor(mode_prelt)) %>% 
  filter(!is.na(age)) %>% 
  group_by(age, mode_prelt) %>% 
  count() %>% 
  group_by(mode_prelt) %>% 
  mutate(total_n = sum(n),
         prop = n / total_n) %>% 
  ungroup() %>%
  ggplot(aes(x = mode_prelt, y = prop, fill = as.factor(age))) +
  geom_bar(stat = "identity") +  #, position = "dodge"
  labs(x = "Mode de prélèvement", y = "Proportion", fill = "Age Class") +
  theme_minimal() +                                   
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75))


reprod_age_data_35 %>% 
  mutate(mode_prelt = as.factor(mode_prelt)) %>% 
  filter(!is.na(age)) %>% 
  ggplot(aes(x = mode_prelt, fill = as.factor(age))) +
  geom_bar() +
  theme_bw()

reprod_age_data_35 %>% 
  select(age, mode_prelt) %>% 
  mutate(mode_prelt = as.factor(mode_prelt)) %>% 
  filter(!is.na(age)) %>% 
  group_by(age, mode_prelt) %>% 
  count() %>% 
  group_by(mode_prelt) %>% 
  mutate(total_n = sum(n),
         prop = n / total_n) %>% 
  ungroup() %>%
  ggplot(aes(x = mode_prelt, y = prop, fill = as.factor(age))) +
  geom_bar(stat = "identity") +  #, position = "dodge"
  labs(x = "Mode de prélèvement", y = "Proportion", fill = "Age Class") +
  theme_minimal() +                                   
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75))


